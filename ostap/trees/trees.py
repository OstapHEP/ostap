#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/trees.py
#  Module with decoration of Tree/Chain objects for efficient use in python
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (  
    'Chain'              , ## helper class , needed for multiprocessing 
    'Tree'               , ## helper class , needed for multiprocessing
    'ActiveBranches'     , ## context manager to activate certain branches 
    'active_branches'    , ## context manager to activate certain branches
    'UseAliases'         , ## context manager to redefine aliases 
    'use_aliases'        , ## context manager to redefine aliases
    'numpy_buffer_types' , ## allowed buffer types from numpy  
    'array_buffer_types' , ## allowed array buffer types 
  ) 
# =============================================================================
from   ostap.core.ostap_types    import ( integer_types      , long_type      ,
                                          string_types       , sequence_types ,
                                          sized_types        , num_types      ,
                                          dictlike_types     , list_types     )
from   ostap.core.core           import ( std , Ostap , VE   , WSE ,
                                          hID , fID   ,  
                                          rootException      ,  
                                          ROOTCWD , strings  , 
                                          valid_pointer      , rootError  ) 
from   ostap.logger.utils        import print_args  
from   ostap.math.reduce         import root_factory
from   ostap.histos.histos       import histo_book2, histo_keys
from   ostap.stats.statvars      import data_decorate , data_range 
from   ostap.trees.cuts          import vars_and_cuts , order_warning
from   ostap.utils.basic         import ( isatty , terminal_size ,
                                          NoContext, loop_items  , typename )
from   ostap.utils.strings       import ( split_string           ,
                                          split_string_respect   ,
                                          var_separators         )
from   ostap.utils.cidict        import cidict, cidict_fun
from   ostap.utils.progress_bar  import progress_bar
from   ostap.utils.progress_conf import progress_conf
from   ostap.utils.scp_copy      import scp_copy
from   ostap.utils.utils         import chunked
from   ostap.math.base           import np2raw  
from   ostap.math.base           import FIRST_ENTRY, LAST_ENTRY, evt_range, all_entries
from   ostap.logger.symbols      import tree           as tree_symbol
from   ostap.logger.symbols      import branch         as branch_symbol
from   ostap.logger.symbols      import leaves         as leaves_symbol
from   ostap.logger.symbols      import tape_cartridge as files_symbol
# 
import ostap.trees.treereduce 
import ostap.trees.param
import ostap.trees.funcs 
import ostap.io.root_file 
import ROOT, os, math, array, sys, numpy    
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.trees' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for Tree/Chain objects')
# =============================================================================
## check validity/emptiness  of TTree/TChain
#  require non-zero poniter and non-empty Tree/Chain
def _tt_nonzero_ ( tree ) :
    """ Check validity/emptiness  of TTree/TChain
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
#    >>> for entry, weight in tree.withCuts ( 'pt>5' ) : print ( entry.y , weight ) 
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
#  for entry, weight  in tree.withCuts ( 'pt>10' , active = ( 'pt' , 'y') ) : sum_y += entry.y 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-06
def _iter_cuts_ ( tree , cuts = '' , first = 0 , last = LAST_ENTRY , progress = False , active = () ) :
    """ Iterator over ``good events'' in TTree/TChain:
    
    >>> tree = ... # get the tree
    >>> for entry, weight  in tree.withCuts ( 'pt>5' ) : print ( i.y , weight ) 
    
    Attention: TTree::GetEntry is already invoked for accepted events,
    no need in second call

    - If only (small) fraction of branches is used in 'cuts' and/or
    only small fraction of branches will be used in the loop,
    the processing can be speed up siginificantly
    by specification of `active' branches:

    >>> tree = ...
    >>> sum_y = 0 
    >>> for entry, weight in tree.withCuts ( 'pt>10' , active = ( 'pt' , 'y') ) : sum_y += entry.y 
    """

    ## show progress only for tty 
    progress = progress and isatty()

    ## redefine first/last 
    first, last = evt_range ( len ( tree ) , first , last ) 
    ##
    #
    if not cuts : cuts = '1'
    #
    if active :
        
        cvar    = tree.the_variables ( cuts , *active )
        abrs    = tuple ( set ( [ c for c in cvar ] ) ) 
        context = ActiveBranches  ( tree , *abrs )
        
    else : context = NoContext () 

    from ostap.utils.progress_conf import progress_conf
    with context :
        
        if progress : pit = Ostap.PyIterator ( tree , progress_conf () , cuts , first , last )
        else        : pit = Ostap.PyIterator ( tree ,                    cuts , first , last )
        ##
        if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
        
        ## the tree is advanced 
        mytree = pit.tree   ()
        weight = pit.weight () 
        while valid_pointer ( mytree ) :
            yield mytree, weight 
            mytree = pit.next   ()
            weight = pit.weight ()
            
        del pit

    ## set the tree at the initial position 
    ievt = tree.GetEntryNumber ( 0 )
    if 0 <= ievt : tree.GetEntry ( ievt )

    
ROOT.TTree .withCuts  = _iter_cuts_ 
ROOT.TChain.withCuts  = _iter_cuts_ 

ROOT.TTree. __len__   = lambda s : s.GetEntries()


# =============================================================================
## Iterator over `good events' in TTree/TChain:
#  @code 
#    >>> tree = ... # get the tree
#    >>> for entry, wweight  in tree( 0, 100, 'pt>5' ) : print ( entry.y, weight) 
#  @endcode
#  @see Ostap::PyIterator
#  @see Ostap::Formula
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-06
def _tc_call_ ( tree , first = 0 , last = LAST_ENTRY  , cuts = None , progress = False , active = () ) :
    """ Iterator over ``good events'' in TTree/TChain:
    
    >>> tree = ... # get the tree
    >>> for entry, weight in tree(0, 100 , 'pt>5' ) : print ( entry.y, weight) 
    
    """

    ## show progress obly for tty 
    progress = progress and isatty()

    first, last = evt_range ( len ( tree ) , first , last )
    
    if active :
        
        cvar    = tree.the_variables ( cuts , *active )
        abrs    = tuple ( set ( [ c for c in cvar ] ) ) 
        context = ActiveBranches  ( tree , *abrs )
        
    else :
        
        from ostap.utils.basic import NoContext 
        context = NoContext () 

    from ostap.utils.progress_conf import progress_conf
    with context : 
        
        if cuts : ## use Ostap.PyIterator 
            
            if progress : pit = Ostap.PyIterator ( tree , progress_conf () , cuts , first , last )
            else        : pit = Ostap.PyIterator ( tree ,                    cuts , first , last )
            
            if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
            
            mytree = pit.tree   ()
            weight = pit.weight () 
            while valid_pointer ( mytree ) :                
                yield mytree, weight            ## YIELD HERE                  
                mytree = pit.next   ()          ## advance to the next entry  
                weight = pit.weight () 
            del pit
            
        else : ## trivial loop 

            ## explicit loop over entries 
            for entry in progress_bar ( range ( firts , last ) , silent = not progress , description = 'Entries:' )  :
                
                ievt = tree.GetEntryNumber ( entry  )
                
                if ievt < 0 :
                    logger.error ( "Invalid entry/1 %s< skip it!" % entry )
                    continue
                
                ievt = tree.GetEntry ( ievt )            
                if ievt < 0 :
                    logger.error ( "Invalid entry/2 %s< skip it!" % entry )
                    continue
                
                yield tree, 1.0         ## YIELD HERE 
            
    ## set the tree at the initial position 
    ievt = tree.GetEntryNumber ( 0 )
    if 0 <= ievt : tree.GetEntry ( ievt )
    
ROOT.TTree .__call__  = _tc_call_ 
ROOT.TChain.__call__  = _tc_call_

# =============================================================================
if numpy : 
    # =========================================================================
    def get_result ( vct ) : return numpyfrombuffer ( vct.data() , count = len ( vct ) , dtype = float )
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    def get_result ( vct ) : return array.array ( 'd' , vct )
    # =========================================================================
# =============================================================================
##  Iterate over tree entries and get a row/array of values for each good entry
#   @code
#   tree = ...
#   for row, weight  in tree.rows ( 'a a+b/c sin(d)' , 'd>0' ) :
#      print ( row , weight ) 
#   @code 
def _tt_rows_ ( tree     , 
               variables , 
               cuts      = ''          , 
               first     = FIRST_ENTRY , 
               last      = LAST_ENTRY  , 
               progress  = False       , 
               active    = ()          ) :
    """ Iterate over tree entries and get a row/array of values for each good entry
    >>> tree = ...
    >>> for row, weight  in tree.rows ( 'a a+b/c sin(d)' , 'd>0' ) :
    >>>    print ( row , weight ) 
    """
    
    ## show progress obly for tty 
    progress = progress and isatty()

    ## redefine first/last 
    first, last = evt_range ( tree , first , last ) 
    
    vars , cuts, _ = vars_and_cuts ( variables , cuts )
    
    if active : context = ActiveBranches  ( tree , *active  )
    else      : context = NoContext () 
        
        
    vars = strings ( vars ) 

    from ostap.utils.progress_conf import progress_conf
    with context :
        
        getter = Ostap.Trees.Getter ( vars , tree ) 
        result = std.vector('double')()
        
        if cuts : ## more efficient loop using Ostap.PyIterator 
            
            if progress : pit = Ostap.PyIterator ( self , progress_conf () , cuts , first , last )
            else        : pit = Ostap.PyIterator ( self ,                    cuts , first , last )
            
            assert pit and pit.ok() , 'ROWS: Invalid formula %s' % cuts  
            
            mytree = pit.tree   ()
            weight = pit.weight () 
            while valid_pointer ( mytree ) :
                
                sc = getter.eval ( result , mytree )
                if sc.isFailure () :
                    logger.error('ROWS: Error status from getter %s' % sc  ) 
                    break
                
                yield get_result ( result ) , weight                  
                mytree = pit.next   ()
                weight = pit.weight () 
            del pit
        
        else : ## trivial loop
          
            for event in progress_bar ( range ( first, last ) , silent = not progress , description = 'Entries:' ) : 
                
                ## tt      = getter.tree()
                tt      = tree
            
                ievent  = tt.GetEntryNumber ( event  )
                if ievent < 0 :
                    logger.error('ROWS: Cannot read entry %s' % event ) 
                    break
                
                ientry  = tt.GetEntry ( ievent )
                if ientry < 0 :
                    logger.error('ROWS: Cannot get entry  %s' % event ) 
                    break
                
                sc = getter.eval ( result , tt )
                if sc.isFailure () :
                    logger.error('ROWS: Error status from getter %s' % sc  ) 
                    break 
                
                yield get_result ( result ) , 1.0                  

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
#    >>> tree.Project ( h1.GetName() , 'm', cuts = 'chi2<10' ) ## standart ROOT 
#    
#    >>> h1   = ROOT.TH1D(... )
#    >>> tree.project ( h1.GetName() , 'm', cuts = 'chi2<10' ) ## ditto 
#    
#    >>> h1   = ROOT.TH1D(... )
#    >>> tree.project ( h1           , 'm', cuts = 'chi2<10' ) ## use histo
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
def tree_project ( tree                     ,
                   histo                    ,
                   what                     ,
                   cuts       = ''          ,
                   first      = FIRST_ENTRY ,          
                   last       = LAST_ENTRY  , 
                   native     = False       ,
                   use_frame  = False       , ## use DataFrame ? 
                   parallel   = False       , ## use parallel stuff?
                   progress   = False       ) :
    """ Helper project method
    
    >>> tree = ...
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.Project ( h1.GetName() , 'm', cuts2<10' ) ## standart ROOT 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1.GetName() , 'm', 'chi2<10' ) ## ditto 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1           , 'm', 'chi2<10' ) ## use histo

    - histo : the histogram (or histogram name)
    - what  : variable/expression to project. It can be expression or list/tuple of expression or comma (or semicolumn) separated expression
    - cuts  : selection criteria/weights
    
    - for 2D&3D cases the natural order of varibales is used.

    """

    ## 1) adjust the first/last
    first , last = evt_range ( tree , first , last  )

    ## 2) if the histogram is specified by the name, try to locate it in ROOT memory  
    if isinstance ( histo , string_types ) :
        groot = ROOT.ROOT.GetROOT()
        h     = groot.FindObject ( histo )
        assert h , "Cannot get locate histogram by name:`%s'" % histo       
        assert isinstance ( h , ROOT.TH1 ) , "Object `%s' exists, but not ROOT.TH1" % typename ( h ) 
        histo = h
        
    ## 3) redirect to the appropriate method
    
    ## chain helper object?
    if   isinstance ( tree , Chain           ) : tree = tree.chain
    elif isinstance ( tree , ROOT.RooAbsData ) :
        from ostap.fitting.dataset import ds_project as _ds_project_
        return _ds_project_ ( tree , histo , what , cuts = cuts , first = first , last = last , progress = progress )
    
    assert isinstance ( tree , ROOT.TTree ) , "Invalid type of 'tree': %s" % type ( tree ) 
    
    ## 4) target        
    target = histo    
    
    ## input histogram?
    input_histo = isinstance ( target , ROOT.TH1 )
    from ostap.trees.param import param_types_nD        
    assert input_histo or  isinstance ( target , param_types_nD ) , \
        'Invalid target/histo type %s' % typename ( target ) 
        
    ## 3) parse input expressions
    varlst, cuts, input_string = vars_and_cuts  ( what , cuts )

    ## copy/clone the target 
    def target_copy   ( t ) : return t.Clone() if isinstance ( t , ROOT.TH1 ) else type ( t ) ( t )
    ## reset the target 
    def target_reset  ( t ) :
        if isinstance ( t , ROOT.TH1 ) : 
            t.Reset ()
            if not t.GetSumw2() : t.Sumw2 () 
        else : t.reset ()
        return t
        
    ## reset the target 
    target = target_reset ( target ) 

    ## dimension of the target 
    dim = target.dim ()
    assert 1 <= dim <= 4 , 'Invalid dimension of target: %s' % dim  

    nvars = len ( varlst )
    assert ( 1 == dim and dim < nvars and not native ) or dim == nvars , \
        'Mismatch between the target/histo dimension %d and input variables %s' % ( dim , varlst )
 
    if input_string and 2 <= dim and order_warning :
        vv = ' ; '.join  ( varlst  ) 
        logger.attention ("project: from v1.10.1.9 variables are in natural order [x;y;..]=[ %s ]" % vv  )
            
    ## avoid looping 
    if first == last : return target

    ## Native ROOT processing
    if native and isinstance ( target , ROOT.TH1 ) and isinstance ( tree , ROOT.TTree ) and \
       ROOT.ROOT.GetROOT().FindObject ( target.GetName() ) is target : 
        
        ## ATTENTION: 
        ## here the inverse/contrintuitive ROOT convention
        ## is used for the ordering of variables
        the_vars = ' : '.join ( '( %s )' % v for v in reversed ( varlst ) )
        ## Use the native ROOT machinery
        rr = tree.Project ( target.GetName () , the_vars , cuts , '' , last - first , first )
        if rr < 0  : logger.error ("Error from TTree::Project %+d" % rr )
        else :
            nn = target.GetEntries ()
            if nn != rr : logger.error ("Mismath for #entries: %+g vs %+d" % ( nn , rr) )
                
        return target
    
    ## use frame if requested and if/when possible 
    if use_frame and all_entries ( tree , first , last ) : 
        import ostap.frames.frames as F 
        frame  = F.DataFrame ( tree )
        if progress : frame , _ = F.frame_progress ( frame , len ( tree ) )            
        return F.frame_project ( frame , target , expressions = varlst , cuts = cuts , lazy = False  )

    ## use parallel processing if requested and if/when possible  
    if parallel and isinstance ( tree , ROOT.TChain ) and all_entries ( tree , first , last ) and 1 < tree.nFiles () :
        if input_histo :
            from ostap.parallel.parallel_project import parallel_project
            return parallel_project ( tree , target , what , cuts , use_frame = use_frame , progress = progress ) 

    tail = cuts , first , last
    
    ## Use our own loop/fill/project machinery
    hp = Ostap.Project ( progress_conf ( progress ) ) 
    
    ## get the list of active branches 
    with ActiveBranches  ( tree , cuts , *varlst ) :
        ## very special case of projection of several expressions into the same 1D-target 
        if 1 == dim and dim < nvars : 
            ## very special case of projections of several expressions into the same 1D-target 
            htmp  = target_copy ( target )  ## prepare the temporary object 
            for var in varlst :
                htmp = target_reset ( htmp ) ## rest the temporary object 
                sc   = hp.project1  ( tree , htmp , var , cuts , *tail )
                if not sc.isSuccess() : logger.error ( "Error from Ostap.Project.project1(%s) %s" %  ( sc , var ) )
                ## update results 
                target += htmp
            del htmp 
        elif 1 == dim :
            sc =  hp.project1 ( tree , target , *varlst , *tail )
            if not sc.isSuccess() : logger.error ( "Error from Ostap.Project.project1 %s" % sc )
        elif 2 == dim : 
            sc =  hp.project2 ( tree , target , *varlst , *tail )
            if not sc.isSuccess() : logger.error ( "Error from Ostap.Project.project2 %s" % sc )
        elif 3 == dim : 
            sc =  hp.project3 ( tree , target , *varlst , *tail )
            if not sc.isSuccess() : logger.error ( "Error from Ostap.Project.project3 %s" % sc )
        elif 4 == dim : 
            sc =  hp.project4 ( tree , target , *varlst , *tail )
            if not sc.isSuccess() : logger.error ( "Error from Ostap.Project.project4 %s" % sc )

        ## return None on error 
        return target if sc.isSuccess() else None 

ROOT.TTree .project = tree_project
ROOT.TChain.project = tree_project

# ======================================================================
## Draw the variables/expressions from TTree object
def tree_draw ( tree                     , 
                what                     ,
                cuts       = ''          ,
                opts       = ''          , 
                first      = FIRST_ENTRY , 
                last       = LAST_ENTRY  , 
                use_frame  = False       , ## use DataFrame ? 
                parallel   = False       , ## use parallel processing?
                native     = False       , ## use native ROOT processinng
                delta      = 0.01        ,
                progress   = False       , **kwargs ) :  
 
    ## check type of opts 
    assert isinstance ( opts , string_types ) , "Invalid type of `opts' : %s" % type ( opts )
    
    ## adjust first/last indices 
    first , last = evt_range ( tree , first , last )
    
    ## decode variables/cuts 
    varlst, cuts, input_string = vars_and_cuts  ( what , cuts )
    if input_string and 2 <= len ( varlst ) and order_warning :
        vv = ' ; '.join  ( varlst ) 
        logger.attention ("draw: from v1.10.1.9 variables are in natural order [x;y;..]=[ %s ]" % vv  )
    
    nvars = len ( varlst ) 
    assert 1 <= nvars <= 3 , "Invalid number of variables: %s" % str ( varlst )
    
    kw = cidict ( transform = cidict_fun , **kwargs )
        
    if native and 1 == nvars :
        with ROOTCWD () , ActiveBranches ( tree , cuts , *varlst ) :
            groot  = ROOT.ROOT.GetROOT()
            groot.cd()            
            hname  = hID ()
            varexp = '(%s) >> %s' % ( varlst [ 0 ] , hname )
            tree.Draw ( varexp , cuts , opts , last - first , first )
            cdir   = ROOT.gDirectory()
            histo  = cdir.Get ( hname )
            assert histo and isinstance ( histo , ROOT.TH1 ) and histo.GetName() == hname , \
                "Cannot retrive the histogram %s from %s" %  ( hname , cdir.GetName () ) 
            ## remove keys related to the booking 
            for k in histo_keys : kw.pop( k , None )
            ## draw the histogram 
            histo.draw ( opts , **kw ) 
            return histo
            
    ## get the suitable ranges for the variables
    args  = first, last 
    ranges = data_range ( tree   ,
                          varlst ,
                          cuts   , *args , 
                          use_frame = use_frame , 
                          parallel  = parallel  , 
                          delta     = delta     , 
                          progress  = progress  )
    if not ranges :
        logger.warning ( 'tree_draw: nothing to draw, return None' )
        return None
    
    histos = []
    for var in varlst :
        mn, mx = ranges [ var ]
        item   = var, ( mn, mx) 
        histos.append ( item ) 

    ## book the histogram 
    histo = histo_book2  ( histos , kw )

    ## fill the histogram 
    histo = tree_project ( tree                  , 
                           histo                 , 
                           varlst                , 
                           cuts      = cuts      , 
                           first     = first     , 
                           last      = last      ,
                           native    = native    ,  
                           use_frame = use_frame , 
                           parallel  = parallel  , 
                           progress  = progress  )
    
    ## draw the histogram 
    histo.draw ( opts , **kw )
    return histo

ROOT.TTree .draw = tree_draw 
ROOT.TChain.draw = tree_draw 

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
    """ Check if object is in tree/chain  :
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
    """ Get min/max for the certain variable in chain/tree
    >>> chain = ...
    >>> mn,mx = chain.vminmax('pt')
    >>> mn,mx = chain.vminmax('pt','y>3')
    """
    if hasattr ( tree , 'pstatVar' ) : 
        if cuts : s = tree.pstatVar ( var , cuts = cuts )
        else    : s = tree.pstatVar ( var )
    else :
        if cuts : s = tree.statVar  ( var , cuts = cuts )
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
    vlst =  tuple ( v.GetName()  for v in t.GetListOfLeaves() )
    if not vlst : return tuple()

    if not pattern :

        lst  = [ v for v in vlst  ]
        lst.sort ()
        return tuple ( lst ) 

    if isinstance ( pattern , string_types ) : pattern  = [ pattern ]

    lst = set()
    for p in pattern : 
        try : 
            import re
            c    =  re.compile ( p , *args )
            vars = [ v for v in vlst if c.match ( v ) ]
            lst  = lst | set ( vars  ) 
        except :
            logger.error ('leaves("%s"): exception is caught, use all ' % p  , exc_info = True ) 
            lst  = lst | set ( [ v for v in vlst  ] )
            
    lst = list ( lst )
    lst.sort () 
    return tuple ( lst ) 

ROOT.TTree.leaves   = _rt_leaves_

# ==============================================================================
## Get the leaf with the certain name 
def _rt_leaf_ ( tree , leaf ) :
    """ Get the leaf with certain name:
    >>> tree = ...
    >>> l = tree.leaf('pt') 
    """
    lst = tuple ( v for v in tree.GetListOfLeaves() )
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
def tree_branches ( t , pattern = '' , *args ) :
    """ Get the list of branch names
    
    >>> tree = ...
    >>> lst = tree.branches()
    >>> for b in lst : print b
    >>> lst = tree.branches( '.*(Muon).*' , re.I )
    >>> for b in lst : print b
    
    """
    vlst =  [ b.GetName() for b in t.GetListOfBranches() ]
    if not vlst : return tuple()

    if pattern :
        # ======================================================================
        try : # ================================================================
            # ==================================================================
            import re
            c  =  re.compile ( pattern , *args )
            lst  = [ v for v in vlst if c.match ( v  ) ]
            lst.sort()
            return tuple ( lst )
            # ================================================================
        except : # ===========================================================
            # ================================================================
            logger.error ('branches: exception is caught, skip it' , exc_info = True ) 
            
    lst  = [ v for v in vlst  ]
    lst.sort()
    return tuple ( lst ) 

ROOT.TTree.branches = tree_branches

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
    """ Simplified print out for tree/chain

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

    result = res.replace ('\n','\n# ')
    return result 
 

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
def _rt_table_0_ ( tree , 
                   pattern = None ,
                   cuts    = ''   ,
                   prefix  = ''   ,
                   title   = ''   ,
                   style   = ''   , *args ) :
    """ Print tree as table 
    """
    ## get list of branches/leaves  
    brs = tree.leaves ( pattern )
    if 'TObject' in brs :
        brs = list  ( brs )
        brs.remove  ( 'TObject' ) 
        brs = tuple ( brs )

    ## collect information
    _vars = []

    ## 
    s0 = tree. statVar ( '1' , *args , cuts = cuts , use_frame = True , parallel = True )    
    n0 = s0.nEntries   ()

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
        
        tn        = l.GetTypeName ()

        selvars += 1 
        if not _in_types ( tn ) : continue
        
        bbs.append ( b ) 


    ## get the statistic 
    bbstats = tree. statVar ( bbs , *args , cuts = cuts , use_frame = True , parallel= True )

    from ostap.stats.counters import WSE, SE  
    if   isinstance ( bbstats , ( WSE , SE ) )  : bbstats = { bbs [ 0 ] : bbstats } 
    
    for b in brs :
        
        l = tree.leaf ( b )

        if not l :
            logger.warning ("table: can't get the leaf  \"%s\"" % b )
            continue
        
        type_name = l.get_type()
        
        rr = [ b , l.get_type () ]
        
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
            title  = '%s("%s","%s"): %d#,' % ( typename ( tree ) , tree.path , tt , len ( tree ) )
        else :
            title  = '%s("%s"): %d#,'      % ( typename ( tree ) , tree.path ,      len ( tree ) )
        if tree_symbol : title = '%s %s' %  ( tree_symbol , title )
        
        nb = len ( tree.branches () )
        title += '%d%s' % ( nb , branch_symbol if branch_symbol else 'branches' ) 
        nl = len ( tree.leaves   () )
        if nl != nb :            
            title += '%d%s'  % ( nl , leaves_symbol if leaves_symbole else  'leaves' )  
        
        if isinstance ( tree , ROOT.TChain ) :
            nfiles = len ( tree.files() )
            if 1 < nfiles : title += ',%d%s' %  ( nfiles , files_symbol if files_symbol else 'files' )
             

    import ostap.logger.table as T
    if table_data : table_data = T.remove_empty_columns ( table_data ) 
    t  = T.table ( table_data , title , prefix = prefix , style = style )
    w  = T.table_width ( t )
    return t , w 
    
# ==============================================================================
## print tree as table 
def _rt_table_1_ ( tree , 
                   variables      ,
                   cuts    = ''   ,
                   prefix  = ''   ,
                   title   = ''   ,
                   style   = ''   , *args ) :
    """ Print tree as table 
    """
    if isinstance ( variables , string_types ) :
        variables = split_string ( variables , strip = True , respect_groups = True )

    bbs = tuple ( sorted ( variables ) ) 

    bbstats = tree. statVar ( bbs , *args , cuts = cuts , use_frame = True , parallel = True ) 

    from ostap.stats.counters import WSE 
    if isinstance ( bbstats , WSE )  : bbstats = { bbs[0] : bbstats } 

    _vars = []
    
    for v in bbstats :

        rr = [ v ]        
        stat = bbstats [ v ] 
        n    = stat.nEntries() 
        mnmx = stat.minmax ()
        mean = stat.mean   () 
        rms  = stat.rms    ()
        rr += [ ( '%+.5g' % mean.value() ).strip() , ## 1
                ( '%.5g'  % rms          ).strip() , ## 2 
                ( '%+.5g' % mnmx[0]      ).strip() , ## 3
                ( '%+.5g' % mnmx[1]      ).strip() ] ## 4
            
        _vars.append ( tuple  ( rr ) )
        
    _vars.sort()
    
    name_l  = len ( 'Variable' )  
    mean_l  = len ( 'mean' ) 
    rms_l   = len ( 'rms'  ) 
    min_l   = len ( 'min'  )  
    max_l   = len ( 'max'  )  
    for v in _vars :
        name_l = max ( name_l , len ( v [ 0 ] ) )
        mean_l = max ( mean_l , len ( v [ 1 ] ) )
        rms_l  = max ( rms_l  , len ( v [ 2 ] ) )
        min_l  = max ( min_l  , len ( v [ 3 ] ) )
        max_l  = max ( max_l  , len ( v [ 4 ] ) )
    

    index_l =   int ( math.ceil ( math.log10( len ( _vars ) + 1 ) ) )
        
    fmt_name  = '%%%ds. %%-%ds' % ( index_l , name_l )
    fmt_mean  = '%%%ds'         % mean_l
    fmt_rms   = '%%-%ds'        % rms_l
    fmt_min   = '%%%ds'         % mean_l
    fmt_max   = '%%-%ds'        % rms_l

    title_l = index_l + 2 + name_l 
    header = (
        ( '{:^%d}' % title_l ).format ( 'Variable' ) ,
        ( '{:^%d}' % mean_l  ).format ( 'mean'     ) ,
        ( '{:^%d}' % rms_l   ).format ( 'rms'      ) ,
        ( '{:^%d}' % min_l   ).format ( 'min'      ) ,
        ( '{:^%d}' % max_l   ).format ( 'max'      ) )
               
    table_data = [ header ] 
    for i , v in enumerate ( _vars ) :
        table_data.append ( ( fmt_name  % ( i + 1 , v [ 0 ] ) ,
                              fmt_mean  %           v [ 1 ] ,
                              fmt_rms   %           v [ 2 ] ,
                              fmt_min   %           v [ 3 ] ,
                              fmt_max   %           v [ 4 ] ) )

    if not title :
        
        tt = tree.GetTitle()
        if tt and tt != tree.GetName() : 
            title  = '%s("%s","%s") %d entries,' % ( typename ( tree ) , tree.path , tt , len ( tree ) )
        else :
            title  = '%s("%s") %d entries,'      % ( typename ( tree ) , tree.path ,      len ( tree ) )
        
        if isinstance ( tree , ROOT.TChain ) :
            nfiles = len ( tree.files() )
            if 1 < nfiles : title += '/%d files ' % nfiles 
        
    import ostap.logger.table as T
    if table_data : table_data = T.remove_empty_columns ( table_data ) 
    t  = T.table ( table_data , title , prefix = prefix , style = style )
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
    """ Get a type for TLeaf object
    >>> tree = ...
    >>> leaf = t.leaf ( 'QQQ' )
    >>> print leaf.get_type ( )
    """
    
    if not leaf : return 'NULL'
    
    branch    = leaf.GetBranch   ()
    type_name = leaf.GetTypeName () 
    
    name     = branch.GetTitle() 
    p1       = name. find ( '[' ) 
    p2       = name.rfind ( ']' )
    if   0 < p1 < p2 :
        type_name = '%s [%s]' % ( type_name , name [ p1 + 1 : p2 ] )
    elif 0 < p1 :
        type_name = '%s [%s]' % ( type_name , name [ p1 : ] )

    type_name = type_name.replace ( 'Float_t'  , 'float'  ) 
    type_name = type_name.replace ( 'Double_t' , 'double' ) 
    type_name = type_name.replace ( 'Bool_t'   , 'bool'   )
    
    return type_name 


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
    """ Get a type for TLeaf object
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
## print root-tree in a form of the table
#  @code
#  data = ...
#  print dat.table() 
#  @endcode
def _rt_table_ (  dataset        ,
                  variables = [] ,
                  cuts      = '' ,
                  prefix    = '' ,
                  title     = '' ,
                  style     = '' , *args ) :
    """ Print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _rt_table_0_ ( dataset          ,
                          variables        ,
                          cuts   = cuts    ,
                          prefix = prefix  ,
                          title  = title   ,
                          style  = style   , *args )[0]

# ==============================================================================
## print root-tree in a form of the table
#  @code
#  data = ...
#  print dat.table2 () 
#  @endcode
def _rt_table2_ (  dataset      ,
                   variables    ,
                   cuts    = '' ,
                   prefix  = '' ,
                   title   = '' ,
                   style   = '' , *args ) :
    """ Print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _rt_table_1_ ( dataset          ,
                          variables        ,
                          cuts   = cuts    ,
                          prefix = prefix  ,
                          title  = title   ,
                          style  = style   , *args )[0]

# =============================================================================
##  print DataSet
def _rt_print2_ ( data  , prefix = '' ) :
    """ Print TTree/TChain"""
    
    br = len ( data.branches () ) + len ( data.leaves() )  
    l  = len ( data             )
    if 10000000 < br * l : return _rt_print_ ( data )
    
    from   ostap.utils.basic import isatty
    if not isatty() : return _rt_table_ ( data )
    tw  , th   = terminal_size()
    rep , wid  = _rt_table_0_ ( data , prefix = prefix ) 
    if wid < tw  : return rep
    ##
    return _rt_print_ ( data )


ROOT.TTree.__repr__ = _rt_print2_
ROOT.TTree.__str__  = _rt_print2_
ROOT.TTree.table    = _rt_table_ 
ROOT.TTree.table2   = _rt_table2_ 

# =============================================================================
## get list of files used for the given chain
#  @code
#  >>> chain = ... ## get the files 
#  >>> files = chain.files() 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rc_files_ ( chain ) :
    """ Get the list of files used for the chain
    
    >>> chain = ... ## get the files 
    >>> files = chain.files()
    """
    return tuple ( i.GetTitle() for i in chain.GetListOfFiles() )


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
    return len ( _rc_files_ ( chain ) )

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
## iterator over individual trees in the chain
#  @code
#  chain = ...
#  for tree in  chain.trees ()  :
#       print len(tree)
#  @endcode 
def _rc_itrees_   ( self ) :
    """ Iterator over individual trees in the echain
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
#  vars , weight = tree.slice ( ['Pt', 'mass'] , 'eta>3' )
#  @endcode 
#  @see numpy.array 
#  @author Albert BURSCHE
#  @date 2015-07-08
def tree_slice ( tree                     ,
                 expressions              ,
                 cuts       = ''          , 
                 structured = True        ,
                 transpose  = True        , 
                 progress   = False       , 
                 use_frame  = False       ,
                 parallel   = False       , 
                 first      = FIRST_ENTRY ,
                 last       = LAST_ENTRY  ) :
    
    """ Get `slice' from TTree in a form of numpy.array
    >>> tree = ...
    >>> varr , _  = tree.slice ( ['Pt','mass'] , 'eta>3' )
    >>> print ( varr )  
    """
    
    ## adjust first/last indices 
    first , last = evt_range ( tree , first , last ) 
    if last <= first :
        return () , None 

    if parallel or use_frame : 
        from ostap.stats.data_statvars import data_slice
        return data_slice ( tree         ,
                            expressions  ,  
                            cuts         , first , last , 
                            structured   = structured   ,
                            progress     = progress     ,
                            use_frame    = use_frame    ,
                            parallel     = parallel     ) 
    
    ## decode cuts & the expressions 
    varlst, cuts , _ = vars_and_cuts  ( expressions , cuts )
    
    ## preliminary result: list of arrays 
    result = [] 

    nvars = len ( varlst ) 
    vars  = ' : '.join ( '(%s)' % v for v in varlst )

    import numpy
    
    with rootException () :
        
        ge   = tree.GetEstimate()
        
        ## redefine the expected estimate
        tree.SetEstimate ( last - first )
        
        ## run internal ROOT machinery 
        num   = tree.Draw ( vars , cuts , "goff" , last - first , first )
        num  = tree.GetSelectedRows ()    
        if tree.GetEstimate() < num :
            tree.SetEstimate ( num )
            logger.debug ( 'Re-run Draw(goff) machinery' ) 
            num = tree.Draw ( vars , cut , "goff" , last - first , first )
            num = tree.GetSelectedRows ()
            
        assert num <= tree.GetEstimate  () , "Something wrong with SetEstimate/GetEstimate/GetSelectedRows"
            
        for k, v  in enumerate ( varlst ) :
            result.append ( numpy.array ( numpy.frombuffer ( tree.GetVal ( k ) , count = num , dtype = numpy.float64 ) , copy = True ) )
        if cuts   :
            result.append ( numpy.array ( numpy.frombuffer ( tree.GetW   (   ) , count = num , dtype = numpy.float64 ) , copy = True ) )
        
        ## reset estimate to the previosus value
        tree.SetEstimate ( ge ) 
            
    if not result :
        return () , None 

    if not cuts :
        weights = None
    else : 
        weights = result [ -1]
        result  = result [:-1]
        if numpy.all ( weights == 1 ) : weights = None 
        
    if structured :
        
        dt   = numpy.dtype ( [ ( v , numpy.float64 ) for v in varlst ] )
        part = numpy.zeros ( num  , dtype = dt )
        for v , a in zip ( varlst , result ) : part [ v ] = a 
        result = part
        
    else :
        
        result = numpy.stack ( result )
        if transpose : result = numpy.transpose ( result )
        
    return result, weights

ROOT.TTree .slice  = tree_slice
ROOT.TTree .slices = tree_slice

# =============================================================================
## extending the existing chain 
def _tc_iadd_ ( self ,  other ) :
    """ Add elements (files,  chains) to existing chain
    >>>  chain  = ...
    >>>  chain += 'myfile.root'
    >>>  chain += ( 'myfile1.root' , 'myfile2.root' )    
    """
    if   self == other                             : return self
    
    elif isinstance ( other , ROOT.TChain ) :
        return _tc_iadd_ ( self , other.files() ) 
    
    elif isinstance ( other , strint_types ) :        
        if not other in self.files () : self.Add ( other )
        return self

    elif isinstance ( o , ( list , tuple , set ) ) :        
        for f in other : _tc_iadd_ ( self , f )
        return  self
        
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
## get all variables needed to evaluate the expressions for the given tree
#  @code
#  tree = 
#  vars = tree.the_variables ( [ 'x>0&& y<13' , 'zzz*15' ] )   
#  vars = the_variables ( tree , [ 'x>0&& y<13' , 'zzz*15' ] ) ## ditto
#  @endcode 
def the_variables ( tree , expression , *args ) :
    """ Get all variables needed to evaluate the expressions for the given tree
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
        assert tf.ok () , "the_variables: Invalid expression : `%s'" % e 
        if not tf.ok()  :
            logger.error ("the_variables: Invalid expression : `%s'" % e )
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
## valid formula expression?
#  @code
#  tree =
#  if not tree.valid_expression ( 'QQ>1' ) : ...
#  @endcode 
def _rt_valid_formula_ ( tree , expression ) :
    """ Valid formula expression?
    >>> tree =
    >>> if not tree.valid_expression ( 'QQ>1' ) : ...
    """
    with rootError () :
        ff   = Ostap.Formula ( fID() , expression , tree ) 
        result = ff.ok ()
        del ff 
        return result
    
ROOT.TTree.valid_formula    = _rt_valid_formula_ 
ROOT.TTree.valid_expression = _rt_valid_formula_ 

# ===============================================================================
## Get all ``size''-variables
#  @code
#  tree  = ... 
#  sizes = tree.size_vars() 
#  @endcode 
def _rt_size_vars_ ( tree ) :
    """ Get all ``size''-variables
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
    """ Get all ``array-like''-variables
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
    """ Context manager to activate only certain branches in the tree.
    - It drastically speeds up the iteration over the tree.
    >>> tree = ...
    >>> with ActiveBraches( tree , '*_Lb', 'eta_Lc') :
    ...    for i in range(1000000) :
    ...    tree.GetEntry ( i )
    ...    print tree.pt_Lb, tree.eta_Lc
    """
    def __init__ ( self , tree , *vars ) :
        
        ## assert tree and vars , 'ActiveBranches: both tree and vars must be valid!'
        
        self.__tree = tree if  tree and isinstance ( tree , ROOT.TTree )          else None 
        self.__vars = self.__tree.the_variables ( *vars ) if vars and self.__tree else () 
        
    ## context manager: ENTER 
    def __enter__ ( self ) :
        
        if self.__vars : ## and self.__tree : 
            ## deactivate the all branches
            self.__tree.SetBranchStatus ( '*' , 0 )     ## deactivate *ALL* branches
            for var in self.__vars :
                ##  activate certain branches 
                self.__tree.SetBranchStatus ( var , 1 )  ## activate only certain branches
        return self.__tree 
    
    ## context manager: EXIT 
    def __exit__ ( self , *_ ) :
        ## reactivate all branches again
        if self.__vars : ## and self.__tree :
            self.__tree.SetBranchStatus ( '*' , 1 )     ## reactivate *ALL* branches
            
    @property
    def vars ( self ) :
        """`vars' : list of activated Tree variabes, or () """
        return self.__vars
    
    @property
    def tree ( self ) :
        """`tree` : the actual TTree or None """
        return self.__tree 

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
    """ Context manager to activate only certain branches in the tree.
    - It drastically speeds up the iteration over the tree.
    >>> tree = ...
    >>> with active__braches ( tree , '*_Lb', 'eta_Lc') :
    ...    for i in range(1000000) :
    ...    tree.GetEntry ( i )
    ...    print tree.pt_Lb, tree.eta_Lc
    """
    return ActiveBranches ( tree , *vars ) 
    
# =============================================================================
## get the list of actiev branches
#  @code
#  tree = ...
#  active = tree.active_branches() 
#  @endcode 
def tree_active_branches ( tree ) :
    """ Get the list of actiev branches
    >>> tree = ...
    >>> active = tree.active_branches() 
    """
    return tuple ( sorted ( b.GetName() for b in tree.GetListOfBranches() if tree.GetBranchStatus ( b.GetName () ) ) ) 

ROOT.TTree.active_branches = tree_active_branches

# =============================================================================
## files and utilisties for TTree/TChain "serialization"
# =============================================================================
## get some file info for the given path 
def file_info ( fname ) :
    """ Get some file info for the given path
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
    """ Simple class to keep definition of tree/chain ``pickable''
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
        self.__nevents = nevents if 0 <= nevents < LAST_ENTRY else -1 
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
        if 0 < self.__first or 0 < nevents < LAST_ENTRY :
            
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
        """ Get/calculate the lengths of the individual trees
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

        if chunk_size <= 0 : chunk_size = LAST_ENTRY 
        
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
        """ Split the tree for several trees with chunk_size entries
        >>> tree = ....
        >>> trees = tree.split ( chunk_size = 1000000 ) 
        """
        if chunk_size <= 0 : chunk_size = LAST_ENTRY 
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
    """ Simple class to keep definition of tree/chain
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
        last = min ( ll , self.first + self.nevents if 0 <= self.nevents else LAST_ENTRY ) 

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

# ==============================================================================
## get a full path of the TTree object 
#  @code
#  rdir = ...
#  path = rdirt.full_path 
#  @endcode
def tree_path ( tree ) :
    """Get a full path of the directory
    >>> tree = ...
    >>> path = tree.full_path 
    """
    rdir = tree.GetDirectory()
    if not rdir : return tree.GetName()
    return os.path.join ( rdir.full_path , tree.GetName() ) 

ROOT.TTree.full_path = property ( tree_path , None , None )


# =============================================================================
## Some decoration for Ostap.Trees.Branches
# =============================================================================
def _brs_iter_  ( brs ) :
    """ Simple iterator over known valid names of branches """
    for name in brs.names():
        branch = brs.branch ( name ) 
        if name and branch : yield name  
# =============================================================================
def _brs_items_  ( brs ) :
    """ simple iterator over valid (name,function) items"""
    for name in brs.names() :
        branch = brs.branch ( name ) 
        if name and branch  : yield name , branch
# =============================================================================
def _brs_getitem_  ( brs , name ) :
    branch = brs.branch ( name )
    if branch : return branch
    raise KeyError ( "Unknown branch name:%s" % name )
# =============================================================================
## Reconstruct Ostap::Trees::Branches objetc 
def _brs_factory_  ( *items ) :
    """ Reconstruct `Ostap::Trees::Branches` object
    """
    branches = Ostap.Trees.Branches () 
    for n, b in items : branches.add ( n , b ) 
    return branches ;
# =============================================================================
## Reduce Ostap.Trees.Brnaches object
def _brs_reduce_  ( brs ) :
    """ Reduce `Ostap.Trees.Branches` object
    """
    return _brs_factory_ , tuple ( ( n , b ) for ( n , b ) in brs.items () )
# =============================================================================
def _brs_str_      ( brs ) :
    if not brs : return "{}"
    return '{ %s }' % ( ', '.join ( "'%s': %s" % item for item in brs.items() ) )
def _brs_contains_ ( brs , name  ) :
    return brs.has_key (  name )

Ostap.Trees.Branches.__len__     = lambda s : s.size()
Ostap.Trees.Branches.__iter__    = _brs_iter_ 
Ostap.Trees.Branches.items       = _brs_items_ 
Ostap.Trees.Branches.iteritems   = _brs_items_ 
Ostap.Trees.Branches.__getitem__ = _brs_getitem_ 
Ostap.Trees.Branches.__contains__= _brs_contains_ 
Ostap.Trees.Branches.__reduce__  = _brs_reduce_ 
Ostap.Trees.Branches.__str__     = _brs_str_ 
Ostap.Trees.Branches.__repr__    = _brs_str_

Ostap.IFuncTree.__str__          = lambda s : typename ( s )


store = []

# ===============================================================================
func_keywords = 'function' , 'callable' , 'what' , 'how' , 'calculator'
# ===============================================================================
## parse&prepare branches
#  @see Ostap::Trees::add_branch 
def prepare_branches ( tree , branch , / , **kwargs ) :
    """ parse& prepare new branches  
    - see `Ostap.Trees.add_branch` 
    """
    assert isinstance ( tree , ROOT.TTree ) and valid_pointer ( tree ) , \
        "prepare_branches: TTree* is invalis!"

    logger.debug ( 'prepare_branches: start %s' % typename ( branch ) ) 
    
    ## names for new branches 
    new_branches = set ()

    ## keep all information 
    keeper = [ branch ]
    ## for k , v in loop_items ( kwargs ) : keeper.append ( v ) 
    
    keeper   = []

    ## the simplest case: branches are already prepared 
    if   isinstance ( branch , Ostap.Trees.Branches        ) and not 'name' in kwargs :        
        ## already prepared ....  not to much action
        for name  in branch : new_branches.add ( name )  
        args = branch ,
        logger.debug ( 'prepare_branches: case [1] %s' % typename ( branch ) ) 

    elif isinstance ( branch , Ostap.Utils.SPLOT ) and not 'name' in kwargs :        
        ## (very) special case
        
        prefix   = kwargs.pop ( 'prefix'  , ''    )
        suffix   = kwargs.pop ( 'suffix'  , '_sw' )
        mapping  = kwargs.pop ( 'mapping' , {}    )
        the_map  = Ostap.Trees.DCT()
        for key, value in loop_items ( mapping ) : the_map [ key ] = value
        for c in branch.coefficients() : new_branches.add ( prefix + c.name + suffix )        
        args     = branch , prefix , suffix , the_map 
        logger.debug ( 'prepare_branches: case [2] %s' % typename ( branch ) ) 

    elif isinstance ( branch , Ostap.Utils.COWs ) and ( not 'name' in kwargs ) : 
    
        ## (very) special case
        names    = tuple ( 'CMP%02d_cw' % i for i in range ( branch.size() ) ) 
        names    = kwargs.pop ( 'names' , names ) 
        names    = strings    ( names ) 
        mapping  = kwargs.pop ( 'mapping' , {} )
        the_map  = Ostap.Trees.DCT() 
        for key, value in loop_items ( mapping ) : the_map [ key ] = value
        for name in names : new_branches.add ( name )        
        args     = branch , names , the_map
        logger.debug ( 'prepare_branches: case [3] %s' % typename ( branch ) ) 

    elif isinstance ( branch , Ostap.IFuncTree ) and 'name' in kwargs :
        ## a simple function        
        args     = branch ,
        logger.debug ( 'prepare_branches: case (IFuncTree) %s' % typename ( branch ) ) 

    elif isinstance ( branch , string_types    ) and 'name' in kwargs :
        ## single branch/formula expression 
        args     = branch , 
        logger.debug ( 'prepare_branches: case [4] %s' % typename ( branch ) ) 

    elif isinstance ( branch , string_types ) and \
         1 == sum ( key in kwargs for key in func_keywords ) and not 'name' in kwargs :
        ## branch is actually the name of the branch
        for key in func_keywords :
            if key in kwargs : 
                func = kwargs.pop ( key )
                logger.debug ( 'prepare_branches: case [5->delegate] %s' % typename ( branch ) )
                args , nb, kw , keep = prepare_branches ( tree , func , name = branch , **kwargs )
                return args , nb , kw , keeper + keep 
            
    elif isinstance ( branch , string_types ) and \
         ( 'formula' in kwargs and isinstance ( kwargs['formula'] , string_types ) ) and not 'name' in kwargs :        
        ## branch is actually the name of the branch 
        formula = kwargs.pop ( 'formula' )
        logger.debug ( 'prepare_branches: case [6->delegate] %s' % typename ( branch ) )        
        args , nb , kw , keep = prepare_branches ( tree , formula , name = branch , **kwargs )
        return args , nb , kw , keeper + keep  
        
    elif isinstance ( branch , Ostap.Utils.RooFun )  and 'name' in kwargs :
        ## RooFit construction 
        mapping  = kwargs.pop ( 'mapping' , {}  )
        the_map  = Ostap.Trees.DCT()
        for key, value in loop_items ( mapping ) : the_map[ key ] = value
        args     = branch , the_map
        logger.debug ( 'prepare_branches: case [7] %s' % typename ( branch ) ) 

    elif isinstance ( branch , ROOT.RooAbsReal ) and 'name' in kwargs and 'observables' in kwargs :        
        ## add informaton from RooFit object 
        ## need to have a little bit of RooFit here 
        import ostap.fitting.roocoelctions 
        ## RooFit construction 
        observables   = kwargs.pop ( 'obsevables'    )
        normalization = kwargs.pop ( 'normalization' , ROOT.nullptr )
        mapping       = kwargs.pop ( 'mapping' , {}  )
        the_map  = Ostap.Trees.DCT()
        for key, value in loop_items ( mapping ) : the_map[ key ] = value
        keeper.append ( observables   )
        keeper.append ( normalization )
        args     = branch , observables , normalization , the_map 
        logger.debug ( 'prepare_branches: case [8] %s' % typename ( branch ) ) 

    elif isinstance ( branch , Ostap.Math.Histo3D ) and all ( ( k in kwargs ) for k in ( 'xname' , 'yname' , 'zname' ) ) :
        
        names = kwargs.pop ( 'xname' ) , kwargs.pop ( 'yname' ) , kwargs.pop ( 'zname' )
        for n in names : new_branches.add ( n )
        args   = names +  ( branch , )
        logger.debug ( 'prepare_branches: case [9] %s' % typename ( branch ) ) 

    elif isinstance ( branch , Ostap.Math.Histo2D ) and all ( ( k in kwargs) for k in ( 'xname' , 'yname' ) ) :
        
        names = kwargs.pop ( 'xname' ) , kwargs.pop ( 'yname' )
        for n in names : new_branches.add ( n )
        args  = names + ( branch , )
        logger.debug ( 'prepare_branches: case [10] %s' % typename ( branch ) ) 

    elif isinstance ( branch , Ostap.Math.Histo1D ) and 'xname' in kwargs :
        
        names = kwargs.pop ( 'xname' ) ,
        for n in names : new_branches.add ( n )
        args  = names + ( branch , )
        logger.debug ( 'prepare_branches: case [11] %s' % typename ( branch ) ) 
        
    elif isinstance ( branch , ROOT.TH1 ) and 1 <= branch.GetDimension() <=3 : 
        ## sample from 1D,2D,3D historgam
        
        names = () 
        if   3 == branch.GetDimension () : names = kwargs.pop ( 'xname' ) , kwargs.pop ( 'yname' ) , kwargs.pop ( 'zname' )
        elif 2 == branch.GetDimension () : names = kwargs.pop ( 'xname' ) , kwargs.pop ( 'yname' ) 
        elif 1 == branch.GetDimension () : names = kwargs.pop ( 'xname' ) ,
        
        for n in names : new_branches.add ( n )
        keeper.append ( branch ) 
        args     = names +  ( branch , )
        logger.debug ( 'prepare_branches: case [12] %s' % typename ( branch ) ) 
        
    elif isinstance ( branch , dictlike_types ) and not 'name' in kwargs :        
        ## general dict-like stuff, converted to Ostap.Trees.Branches 

        branches = Ostap.Trees.Branches()

        for key , value in loop_items ( branch ) :
            keeper.append ( value ) 
            
            assert isinstance ( key , string_types ) and Ostap.Trees.valid_name_for_branch ( key ) ,\
                "Invalid branch name type/value" % ( typename ( key ) , key )
            
            if   isinstance ( value , string_types    ) : branches .add ( key , value , tree )
            elif isinstance ( value , Ostap.IFuncTree ) : branches .add ( key , value , tree )
            elif callable   ( value ) :
                from ostap.trees.funcs import PyTreeFunction as PTF
                fun = PTF     ( value , tree )
                branches.add  ( key , fun , tree ) 
            else :
                raise TypeError ( "Invalid branch type:%s for key=%s" % ( typename ( value ) , key ) )
            
            new_branches.add ( key )
            
        keeper.append ( branches )  
        args     = branches ,
        logger.debug ( 'prepare_branches: case [13] %s' % typename ( branch )  )

    elif callable ( branch  ) and 'name' in kwargs and 'arguments' in kwargs : 
        ## generic callable that takes 1-3 arguments from the tree 
        
        vars = kwargs.pop( 'arguments' , () )
        if isinstance ( vars , string_types ) :
            vars = split_string ( vars , var_separators , strip = True , respect_groups = True )
            
        assert 1 <= len ( vars ) <= 3 and all ( ( isinstance ( a , string_types ) and a.strip() ) for a in vars ) , \
            "1-to-3 non-empty strings must be specified as `arguments' for `branch'=%s" % typename ( branch )

        fvars = ( branch , ) + vars + ( tree , )        
        keeper.append ( fvars[:-1]  ) 
        
        args = vars + ( branch , )
        ## if   ( 6 , 26 ) <= root_info : args = vars + ( branch , )
        ## else :
        ## 
        ##    fvars = ( branch , ) + vars + ( tree , )
        ##    if   1 == len ( vars ) : args = Ostap.Functions.Func1D ( *fvars ) , 
        ##    elif 2 == len ( vars ) : args = Ostap.Functions.Func2D ( *fvars ) , 
        ##    elif 3 == len ( vars ) : args = Ostap.Functions.Func3D ( *fvars ) ,

        logger.debug ( 'prepare_branches: case [14] %s' % typename ( branch ) )
        
    elif callable ( branch  ) and 'name' in kwargs  and not any ( key in  kwargs for key in func_keywords ) :
        
        ## generic callable  ( function takes TTree as argument ) 
        name      = kwargs.pop ( 'name' )
        newbranch = { name : branch }
        keeper.append (    branch )
        keeper.append ( newbranch ) 
        logger.debug ( 'prepare_branches: case [15->delegate] %s' % typename ( branch ) )
        args , nb , kw , keep = prepare_branches ( tree , newbranch , **kwargs )
        return args , nb , kw , keeper + keep 
    
    else :
        
        table = print_args ( branch , **kwargs ) 
        logger.error    ( 'Invalid/inconsistent branch/kwargs structure:\n%s' % table )  
        raise TypeError ( 'Invalid/inconsistent branch?kwargs structure') 

    if not new_branches :
        name = kwargs.pop ( 'name' )
        assert name , "Branch name `name=...' must be provided!'"
        new_branches.add ( name )
        args = ( name , ) + args

    ## check names for new brnaches 
    for name in new_branches :
        assert name and isinstance ( name , string_types ) , "Invalid new brach name!"
        assert Ostap.Trees.valid_name_for_branch ( name  ) , "Invalid name for bew brnach:'%s'" % name
        assert not name in tree , "Branch/leave `%s' is already in the Treee!" % name   

    for a in args : keeper.append ( a )
    
    logger.debug ( 'prepare_bramnches, end...' ) 
    return args , new_branches , kwargs , keeper 


# ===============================================================================
## Add new branch to the tree
#  @see Ostap::Trees::adD_branch 
def add_new_branch ( tree , branch , / , **kwargs ) :
    """ Add new branch to the tree
    - see `Ostap,Trees.add_branch`
    """
    assert valid_pointer ( tree ) and isinstance ( tree , ROOT.TTree ) , \
        "Tree* is invalid!"

    verbose = kwargs.pop ( 'verbose ' , False ) ## ATTENTION! 

    ## parse argumensts
    args , expected , kw , keeper = prepare_branches ( tree , branch , **kwargs ) 

    progress = kw.pop ( 'progress' , True  )
    report   = kw.pop ( 'report'   , True  ) 

    if verbose :
        logger.info ( 'Initial   arguments:\n%s' % print_args ( branch , **kwargs ) )
        logger.info ( 'Processed arguments:\n%s' % print_args ( *args  , **kw     ) ) 
        
    assert args     , "No arguments for Ostap.Trees.add_branch are collected!"
    assert expected , "No expected branche detected!!"

    kw.pop ( 'prefix' , '' ) 
    kw.pop ( 'title'  , '' ) 
    if kw : 
        title1 = 'add_new_branch: Unknown/unprocessed arguments'
        title2 = 'Unknown/uprocessed arguments'
        logger.warning  ( '%s:\n%s' % ( title1 , print_args ( prefix = '# ' , title = title2 , **kw ) ) ) 

    ## start the actual processing 
    chain    = push_2chain ( tree , *args ,  progress = progress , report = report )

    missing  = sorted ( br for br in expected if not br in chain  )
    if missing : logger.warning ( 'Extected but missing branches: %s' % ( ', '.join ( m for m in missing ) ) ) 

    del keeper 
    return chain
    
# =============================================================================
## Add new branch(s) to (single) TTree
#  @param tree input tree
#  @param params parameters
#  @see Ostap::Trees::add_branch
def push_2tree ( tree , *config , progress = True , report = True ) :
    """ Add new brnaches to (snigole) TTree according to specificaions
    - @see `Ostap.Trees.add_branch`
    """
    assert isinstance ( tree , ROOT.TTree ) and valid_pointer ( tree  ) , \
        "push_2tree: `tree' is invalid"
    
    if isinstance ( tree , ROOT.TChain ) and 1 < len ( tree.files() ) :
        return push_2chain ( tree , *config , progress = progress , report = report ) 

    treepath = tree.path
    the_file = tree.topdir
    groot = ROOT.ROOT.GetROOT() 
    assert treepath and the_file and ( not the_file is groot ) and isinstance ( the_file , ROOT.TFile ) , \
        'This is not the file-resident TTree* object! addition of new branch is not posisble!'
    the_file = the_file.GetName() 
    
    tname = tree.GetName      ()
    tdir  = tree.GetDirectory ()
    tpath = tree.path

    ## list of existing brranches/leaves 
    branches = ( set ( tree.branches() ) | set ( tree.leaves() ) ) if report else set() 

    from ostap.utils.progress_conf import progress_conf
    
    args = tuple ( a for a in config  ) 
    if progress : args += ( progress_conf ( progress ) , ) 

    from ostap.io.root_file  import REOPEN     
    with ROOTCWD() , REOPEN ( tdir ) as tfile :
        
        tfile.cd()
        ttree    = tfile.Get ( tpath )
        assert valid_pointer ( ttree ) , 'Invalid TTree:%s in file:%s' % ( tpath , the_file )

        assert tfile.IsWritable() , 'The file:%s is not writeable!'
        
        ## Add the branch!
        ## table = print_args ( *args ) 
        ## logger.always ( 'ARGUMETS:\n%s' % table )
        sc    = Ostap.Trees.add_branch ( ttree , *args  )
        assert sc.isSuccess () , "Error from Ostap.Trees.add_branch %s" % sc

        from ostap.utils.root_utils import implicitMT 
        with implicitMT ( False ) : tfile.Write ( "" , ROOT.TObject.kOverwrite )
         
        ttree = ROOT.nullptr 

    chain = ROOT.TChain ( tpath )
    chain.Add  ( the_file )

    if report :
        new_branches = sorted ( ( set ( chain.branches () ) | set ( chain.leaves () ) ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 == n : title = 'Added %s branch to TChain'   % n 
            else      : title = 'Added %s branches to TChain' % n 
            table = chain.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            chain = ROOT.TChain ( tpath )
            chain.Add  ( the_file )     
            
    return chain

# =============================================================================
## Add new branch(s) to chain 
#  @param tree input tree
#  @param params parameters
#  @see Ostap::Trees::add_branch
def push_2chain ( chain , *config , progress = True , report = True ) :
    """ Add new brnaches to (snigole) TTree accoring to specificaions
    - @see `Ostap.Trees.add_branch` 
    """
    assert isinstance  ( chain , ROOT.TTree ) and valid_pointer ( chain ) , \
        "push_2chain: `chain' is invalid"

    ## delegate to TTree: 
    if not isinstance ( chain , ROOT.TChain ) or 2 > len ( chain.files() ) :
        return push_2tree ( chain , *config , report = report , progress = progress ) 

    ## name & path 
    cname, cpath = chain.GetName () , chain.path

    ## list of existing brranches/leaves 
    branches = ( set ( chain.branches() ) | set ( chain.leaves() ) ) if report else set() 

    ## files to be processed 
    files = chain.files()
    
    chain_progress = progress and ( 5 <= len ( files ) ) 
    tree_progress  = progress and ( 5 >  len ( files ) ) 
    
    for fname in progress_bar ( files , silent = not chain_progress , description = 'Files:' ) :
        tree = ROOT.TChain ( cname )
        tree.Add ( fname ) 
        ## treat the tree 
        push_2tree ( tree , *config , report = False , progress = tree_progress ) 
        
    ## reconstruct the resulting chain 
    chain = ROOT.TChain ( cname )
    for fname in files : chain.Add ( fname )
            
    if report :        
        new_branches = sorted ( ( set ( chain.branches () ) | set ( chain.leaves () ) ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 >= n : title = 'Added %s branch to TChain'   % n 
            else      : title = 'Added %s branches to TChain' % n 
            table = chain.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            chain = ROOT.TChain ( cname )
            for fname in files : chain.Add ( fname )
    
    return chain

# ==============================================================================
try : # ========================================================================
    # ==========================================================================
    import numpy
    numpy_buffer_types = ( numpy.float16 ,
                           numpy.float32 ,
                           numpy.float64 ,
                           numpy.byte    ,   
                           numpy.int8    , 
                           numpy.int_    ,
                           numpy.int16   ,
                           numpy.int32   ,
                           numpy.int64   ,
                           numpy.ubyte   ,
                           numpy.uint8   ,
                           numpy.uint    ,
                           numpy.uint16  ,
                           numpy.uint32  ,
                           numpy.uint64  )    
    # ==========================================================================
except ImportError : # =========================================================
    # ==========================================================================
    numpy = None
    numpy_buffer_types = () 

# =========================================================================
## valid array-types for adding to TTree
# =========================================================================
array_buffer_types = ( 'f' , 'd' , 'h' , 'i' , 'l' , 'H' , 'I' , 'L' , 'b' , 'B' , 'q' , 'Q' )

# ============================================================================
## Add new buyffer to TTree
#  @see Ostap::Trees::add_buffer 
def add_new_buffer ( tree , name , buffer , **kwargs ) :
    """ Add new buyffer to TTree
    - see `Ostap.Trees.add_buffer` 
    """
    
    assert valid_pointer ( tree ) and isinstance ( tree , ROOT.TTree ) , \
        "add_new_buffer: Tree* is invalis!"

    assert name and isinstance ( name , string_types ) , \
        "add_new_buffer: Invalid name for new branch: :%s/%s" % ( typoname ( name ) , name )

    ## assert not name in tree , "Name `'%s' is already in TTree!" % name
    assert Ostap.Trees.valid_name_for_branch ( name ) , \
        "add_new_buffer: name `%s' is not valid!" % name

    keep  = [ buffer ] 
    for k , v in loop_items ( kwargs ) : keep.append ( ( k , v ) )

    extra_floats = () if not numpy else ( numpy.float16 , )
    extra_ints   = () if not numpy else ( numpy.int8    , )
    
    if numpy and np2raw and isinstance ( buffer , numpy.ndarray ) and buffer.dtype in numpy_buffer_types :
        ## numpy array of valid types

        keep.append ( buffer )
        
        ## recast 
        if   buffer.dtype == numpy.float16 :
            buffer = numpy.asarray ( buffer , dtype = numpy.float32 )
            return add_new_buffer ( tree , name , buffer , **kwargs ) 

        buflen = len ( buffer ) 
        ## convert ndarray into raw C++ buffer 
        raw_buffer, _ = np2raw ( buffer ) 

        if   buffer.dtype in ( numpy.int8  , numpy.byte  ) :
            the_buffer = Ostap.Utils.schar_buffer ( raw_buffer , buflen ) ## SCHAR_MAKE 
        elif buffer.dtype in ( numpy.uint8 , numpy.ubyte ) :
            the_buffer = Ostap.Utils.uchar_buffer ( raw_buffer , buflen ) ## CCHAR_MAKE 
        else :
            the_buffer = Ostap.Utils.make_buffer  ( raw_buffer , buflen  )
                        
        keep.append ( raw_buffer ) 
        value      = kwargs.pop ( 'value' , 0 )
        if value : the_buffer.setValue ( value )
            
    elif isinstance ( buffer , array.array ) and buffer.typecode in array_buffer_types :
        ## arra.array f valied types

        if   'b' == buffer.typecode : the_buffer = Ostap.Utils.schar_buffer ( buffer , len ( buffer ) )
        elif 'B' == buffer.typecode : the_buffer = Ostap.Utils.uchar_buffer ( buffer , len ( buffer ) )
        else                        : the_buffer = Ostap.Utils.make_buffer  ( buffer , len ( buffer ) )

        value      = kwargs.pop ( 'value' , 0  )
        if value : the_buffer.setValue ( value )

    elif numpy and isinstance ( buffer ,  ( bytes , bytearray , memoryview ) ) :

        ## construct numpyarray from raw buffer 
        the_array = numpy.frombuffer ( buffer , dtype = numpy.byte )
        
        keep.append ( the_array ) 
        return add_new_buffer ( tree , name , the_array , **kwargs )

    elif isinstance ( buffer ,  ( bytes , bytearray , memoryview ) ) :

        ## construct arra from raw  buffer 
        the_array = array.array ( 'b' , buffer )
        
        keep.append ( the_array )                
        return add_new_buffer ( tree , name , the_array , **kwargs )
    
    elif numpy and isinstance ( buffer , sequence_types ) :
        ## seqence convertivel to numpy

        ## construct arrat from the sequence 
        the_array = numpy.asarray ( buffer )
        
        keep.append ( the_array )
        return add_new_buffer ( tree , name , the_array , **kwargs )

    elif isinstance ( buffer , sequence_types ) :
        ## seqence convertivel to array-array 

        typecode  = kwargs.pop ( 'typecode' , 'd'    ) 
        the_array = array.array ( typecode  , buffer )

        keep.append ( the_array ) 
        return add_new_buffer ( tree , name , the_array , **kwargs )

    else :
        
        table = print_args ( name , buffer , **kwargs ) 
        logger.error    ( 'Invalid/inconsistent name/buffer/kwargs structure:\n%s' % table )  
        raise TypeError ( 'Invalid/inconsistent name/buffer/kwargs structure') 

    keep.append ( the_buffer )
    
    progress = kwargs.pop ( 'progress' , True  )
    report   = kwargs.pop ( 'report'   , True  ) 

    if kwargs :
        title1 = 'add_new_branch: Unknown arguments'
        title2 = 'Unknown arguments'
        logger.warning  ( '%s:\n%s' % ( title1 , print_args ( prefix = '# ' , title = title2 , **kwargs ) ) ) 

    ## start the actual processing 
    return buffer_2chain ( tree , name , the_buffer , progress = progress , report = report )

# =============================================================================
## Add new branch(s) to (single) TTree
#  @param tree input tree
#  @param params parameters
#  @see Ostap::Trees::add_branch
def buffer_2tree ( tree , name , buffer , progress = True , report = True ) :
    """ Add new brnaches to (snigole) TTree accoring to specificaions
    - @see Ostap.Trees.add_branch 
    """
    assert isinstance ( tree , ROOT.TTree ) and valid_pointer ( tree  ) , \
        "buffer_2tree: `tree' is invalid"
    
    if isinstance ( tree , ROOT.TChain ) and 1 < len ( tree.files() ) :
        return  buffer_2chain ( tree , config , progress = progress , report = report ) 

    treepath = tree.path
    the_file = tree.topdir
    groot    = ROOT.ROOT.GetROOT()
    assert treepath and the_file and ( not the_file is groot ) and isinstance ( the_file , ROOT.TFile ) , \
        'This is not the file-resident TTree* object! addition of new branch is not posisble!'
    the_file = the_file.GetName() 
    
    tname = tree.GetName      ()
    tdir  = tree.GetDirectory ()
    tpath = tree.path

    ## list of existing brranches/leaves 
    branches = ( set ( tree.branches() ) | set ( tree.leaves() ) ) if report else set() 
    from ostap.utils.progress_conf import progress_conf

    if progress : args = name , buffer , progress_conf ( progress )
    else        : args = name , buffer 

    from ostap.io.root_file  import REOPEN     
    with ROOTCWD() , REOPEN ( tdir ) as tfile :
        
        tfile.cd()
        ttree    = tfile.Get ( tpath )
        assert valid_pointer ( ttree ) , 'Invalid TTree:%s in file:%s' % ( tpath , the_file )

        assert tfile.IsWritable() , 'The file:%s is not writeable!'
        
        ## Add the branch!
        sc = Ostap.Trees.add_buffer ( ttree , *args  )
        assert sc.isSuccess () , "Error from Ostap.Trees.add_branch %s" % sc

        from ostap.utils.root_utils import implicitMT 
        with implicitMT ( False ) : tfile.Write ( "" , ROOT.TObject.kOverwrite )

    ## recostrut the tree/chain 
    chain = ROOT.TChain ( tpath )
    chain.Add  ( the_file )

    if report :
        new_branches = sorted ( ( set ( chain.branches () ) | set ( chain.leaves () ) ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 == n : title = 'Added %s branch to TChain'   % n 
            else      : title = 'Added %s branches to TChain' % n 
            table = chain.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
          
    return chain

# =============================================================================
## Add new branch(s) to chain 
#  @param tree input tree
#  @param params parameters
#  @see Ostap::Trees::add_branch
def buffer_2chain ( chain , name , buffer , progress = True , report = True ) :
    """ Add new brnaches to (snigole) TTree accoring to specificaions
    - @see Ostap.Trees.add_branch 
    """
    assert isinstance  ( chain , ROOT.TTree ) and valid_pointer ( chain ) , \
        "buffer_2chain: `chain' is invalid"

    ## delegate to TTree: 
    if not isinstance ( chain , ROOT.TChain ) or 2 > len ( chain.files() ) :
        return buffer_2tree ( chain , name , buffer , report = report , progress = progress ) 

    ## name & path 
    cname, cpath = chain.GetName () , chain.path

    ## list of existing brranches/leaves 
    branches = ( set ( chain.branches() ) | set ( chain.leaves() ) ) if report else set() 

    ## files to be processed 
    files = chain.files()
    
    chain_progress = progress and ( 5 <= len ( files ) ) 
    tree_progress  = progress and ( 5 >  len ( files ) ) 
    
    for fname in progress_bar ( files , silent = not chain_progress , description = 'Files:' ) :
        tree = ROOT.TChain ( cname )
        tree.Add ( fname ) 
        ## treat the tree
        
        tree   = buffer_2tree ( tree , name , buffer , report = False , progress = tree_progress )
        offset = len ( tree )         
        buffer = buffer.offset ( offset )

    ## reconstruct the chain 
    chain = ROOT.TChain ( cname )
    for fname in files : chain.Add ( fname )

    if report :        
        new_branches = sorted ( ( set ( chain.branches () ) | set ( chain.leaves () ) ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 >= n : title = 'Added %s branch to TChain'   % n 
            else      : title = 'Added %s branches to TChain' % n 
            table = chain.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
          
    return chain

ROOT.TTree.add_new_branch  = add_new_branch 
ROOT.TTree.add_new_buffer  = add_new_buffer 

# =============================================================================
## Context manager to temporariliy redefine aliases
#  @code
#  tree = ...
#  with UseAliases ( tree , px = 'PX/1000' , py = 'PY/1000' ) :
#  ... 
#  @endcode 
class UseAliases(object) :
    """ Context manager to temporariliy redefine aliases
    >>> tree = ...
    >>> with UseAliases ( tree , px = 'PX/1000' , py = 'PY/1000' ) :
    """
    def __init__ ( self , tree , **aliases ) :

        self.__aliases  = {}
        self.__aliases.update ( aliases )
        slef.__tree     = tree 
        self.__previous = {}
        
    def __enter__ ( self ) :

        if self.__tree and self.__aliases :
            for a in self.__aliases :
                ## 1) save previosu alias
                ea = tree.GetAlias ( a )
                if ea : self.__previous [ a ] = ea
                ## 1) define new alias 
                tree.SetAlias ( a , self.__aliases [ a ] )
                
        return self
    
    def __exit__ ( self , *_ ) :
        
        if self.__tree and self.__previous :
            # reset all aliases 
            for a in self.__previous  :
                tree.SetAlias ( a , self.__previous [ a ] )
            
# =============================================================================
## Context manager to temporariliy redefine aliases
#  @code
#  tree = ...
#  with use_aliases ( tree , px = 'PX/1000' , py = 'PY/1000' ) :
#  ... 
#  @endcode 
def use_aliases ( tree , **aliases ) :
    """ Context manager to temporariliy redefine aliases
    >>> tree = ...
    >>> with use_aliases ( tree , px = 'PX/1000' , py = 'PY/1000' ) :
    """
    return UseAliases ( tree , **aliases ) 

# =============================================================================
_new_methods_   = data_decorate ( ROOT.TTree )
del data_decorate 
_new_methods_  += (
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
    ROOT.TTree .statCovs  ,
    ROOT.TChain.statCovs  ,
    #
    ROOT.TTree .vminmax   ,
    ROOT.TChain.vminmax   ,
    #  
    ROOT.TTree.branches        , 
    ROOT.TTree.active_branches , 
    ## 
    ROOT.TTree .__repr__     , 
    ROOT.TTree .__str__      ,
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
    ROOT.TTree.valid_formula    ,
    ROOT.TTree.valid_expression ,
    #
    ROOT.TTree.the_variables    ,
    ROOT.TTree.add_new_branch   ,
    ##
    ROOT.TLeaf.get_type         ,
    ROOT.TLeaf.get_type_short   ,
    ROOT.TLeaf.get_short_type   ,
)


_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,   
    ROOT.TLeaf      
)
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
