#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/frames/tree_reduce.py
#  Helper module to "Reduce" tree using frames
#  @see Ostap::DataFrame
#  @see ROOT::RDataFrame
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-06-16
# =============================================================================
"""Helper module to ``reduce'' tree using frames
- see Ostap.DataFrame
- see ROOT.ROOT.RDataFrame
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'ReduceTree' ,
    'reduce'     ,    
    ) 
# =============================================================================
import ROOT, os 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.frames.tree_reduce' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( "``Reduce'' TTree using ROOT::RDataFrame object")
# =============================================================================
import ostap.trees.trees
from   ostap.core.core     import cpp, Ostap 
from   ostap.utils.cleanup import CleanUp 
# =============================================================================
## @class ReduceTree
#  Reduce TTree object using intermediate (temporary
#  @code
#  tree    = ...
#  r       = ReduceTree ( tree , cuts , [ 'px', 'py', 'pz' ] , 'new_file.root' )
#  reduced = t.tree 
#  @endcode
class ReduceTree(CleanUp):
    """Reduce ROOT.TTree object 
    >>> tree    = ...
    >>> r       = ReduceTree ( tree , cuts , [ 'px', 'py', 'pz' ]
    >>> reduced = r.tree    
    """
    def __init__ ( self               ,
                   chain              ,   ## input TChain/TTree 
                   selection  = {}    ,   ## selection/cuts  
                   save_vars  = ()    ,   ## list of variables to save 
                   new_vars   = {}    ,   ## new variables
                   no_vars    = ()    ,   ## exclude these variables
                   ##
                   output     = ''    ,   ## output file name
                   addselvars = False ,   ## add varibles from selections?
                   tmp_keep   = False ,   ## keep the temporary file 
                   silent     = False ):  ## silent processing 
        
        from   ostap.frames.frames import DataFrame
        frame  = DataFrame ( chain )
        report = None

        nvars = [] 
        ## new variables 
        for nv in new_vars :
            frame = frame.Define ( nv , new_vars [ nv] )
            nvars.append ( nv )

        from ostap.core.ostap_types import  ( string_types   ,
                                              listlike_types ,
                                              dictlike_types )
        
        cut_types = string_types + ( ROOT.TCut , )
        
        Lmax = 30 

        selections = [] 
        if    selection and isinstance ( selection , cut_types ) :
            ss = str ( selection ).strip() 
            if len ( ss )  < Lmax : filter_name = ss          
            else                  : filter_name = 'SELECTION'
            logger.info  ( 'Add filter/1 %s :  %s' % ( ss , filter_name ) )
            frame = frame.Filter ( ss  , filter_name )
            selections.append ( ss ) 
        elif selection and isinstance ( selection , dictlike_types  ) :
            for filter_name in selection :
                s = selection [ filter_name ]
                assert isinstance ( s , cut_types ),\
                       'Invalid selection type %s/%s' % ( s , type ( s ) )
                ss = str ( s ).strip()
                logger.info ( 'Add filter/2 %s :  %s' % ( ss , filter_name ) )
                frame = frame.Filter ( ss , str ( filter_name ) ) 
                selections.append ( ss ) 
        elif selection and isinstance ( selection , listlike_types ) :
            for i , s in enumerate  ( selection ) :
                assert isinstance ( s , cut_types ),\
                       'Invalid selection type %s/%s' % ( s , type ( s ) )
                ss = str( s ).strip()
                ##
                if len ( ss ) < Lmax          : filter_name = ss
                else                          : filter_name = 'SELECTION%d' % i
                #
                logger.info ( 'Add filter/3 %s :  %s' % ( ss , filter_name ) )
                frame = frame.Filter ( ss , filter_name )
                selections.append ( ss ) 
        elif selection :
             raise TypeError('Invalid  selection type %s/%s' %  ( selection , type ( selection ) ) )

        if not output : 
            output = self.tempfile ( prefix = 'frame-' , suffix = '.root' )
            logger.debug ( 'ReduceTree: output file is %s' % output )  
            if not tmp_keep : self.trash.add ( output  )

        if  selections : report = frame.Report()
            
        if selections and addselvars :
            bvars     = chain.the_variables ( selections )
            save_vars = list  ( bvars ) + [ v for v in save_vars if not v in bvars ]
            save_vars = tuple ( save_vars ) 

        ## exclude some variables 
        if no_vars and not save_vars :
            bvars     = list ( chain.branches () )
            all_vars  = list ( bvars ) + [ v for v in nvars if not v in bvars ]
            save_vars = tuple ( [ v for v in all_vars if not v in no_vars ] )
        elif no_vars :
            bvars    = chain.the_variables ( *save_vars )
            all_vars  = list ( bvars ) + [ v for v in nvars if not v in bvars ]
            save_vars = tuple ( [ v for v in all_vars if not v in no_vars ] ) 
            
        nb_ = len ( chain.branches () )
        ne_ = len ( chain             )
        
        ## chain name:
        ## FIXME!
        # cname = chain.GetName()  ## produces ROOT error 
        _ , _ , cname = chain.GetName().rpartition ( '/' ) 
        
        if not save_vars : 
            snapshot = frame.Snapshot ( cname , output )
        else :
            bvars    = chain.the_variables ( *save_vars )
            all_vars = list ( bvars ) + [ v for v in nvars if not v in bvars ]
            from ostap.core.core import strings as _strings
            all_vars = _strings ( all_vars  ) 
            snapshot = frame.Snapshot ( cname , output , all_vars )

        assert os.path.exists ( output ) and\
               os.path.isfile ( output ) , 'Invalid file %s' % fname

        self.__chain = ROOT.TChain ( cname )
        self.__chain.Add ( output ) 
        self.__output = output 

        self.__report = 'Tree -> Frame -> Tree filter/transformation'
        self.__table  = [] 
        if report :
            from ostap.frames.frames import report_print, report_as_table 
            title = self.__report 
            self.__report += '\n%s' % report_print ( report , title , '# ')
            self.__table   = report_as_table ( report )
                
        fs = os.path.getsize ( self.__output )        
        gb , r = divmod ( fs ,  1024 * 1024 * 1024 )
        mb , r = divmod ( r  ,  1024 * 1024 )
        kb , r = divmod ( r  ,  1024 )
        

        if   gb : fs = '%.1fGB' % ( float ( fs ) / 1024 / 1024 / 1024 )
        elif mb : fs = '%.1fMB' % ( float ( fs ) / 1024 / 1024 )
        elif kb : fs = '%.1fkB' % ( float ( fs ) / 1024 ) 
        else    : fs = '%sB'    %           fs
        
        nb = len ( self.__chain.branches () )
        ne = len ( self.__chain             )

        self.__report += '\n# Reduce %d -> %d branches, %d -> %d entries' % ( nb_ , nb , ne_ , ne ) 
        self.__report += '\n# Output:%s size:%s'                          % ( self.__output , fs  )
        self.__report += '\n# %s' % str ( self.__chain ) 
        
    def __str__  ( self ) : return self.__report 
    def __repr__ ( self ) : return self.__report 
        
    @property
    def output   ( self ) :
        """``output'' : the output file name"""
        return self.__output
    
    @property
    def chain  ( self ) :
        """``chain'': the reduced chain/tree (same as tree)"""
        return self.__chain
    
    @property
    def tree   ( self ) :
        """``tree'': the reduced chain/tree (same as chain)"""
        return self.__chain

    @property
    def table ( self ) :
        """``table'' : get the statitics as table"""
        return self.__table
    
# ===============================================================================
## Reduce the tree/chain
#  @code
#  tree = ....
#  reduced1 = tree.reduce  ( 'pt>1' )
#  reduced2 = tree.reduce  ( 'pt>1' , vars = [ 'p', 'pt' ,'q' ] )
#  reduced3 = tree.reduce  ( 'pt>1' , no_vars = [ 'Q', 'z' ,'x' ] )
#  reduced4 = tree.reduce  ( 'pt>1' , new_vars = { 'pt2' : 'pt*pt' } )
#  reduced5 = tree.reduce  ( 'pt>1' , new_vars = { 'pt2' : 'pt*pt' } , output = 'OUTPUT.root' )
#  @endcode
def reduce  ( tree               ,
              selection          ,
              save_vars  = ()    , 
              new_vars   = {}    ,
              no_vars    = ()    ,
              output     = ''    ,
              addselvars = False ,
              silent     = False ) :
    
    """ Reduce the tree/chain
    >>> tree = ....
    >>> reduced1 = tree.reduce  ( 'pt>1' )
    >>> reduced2 = tree.reduce  ( 'pt>1' , vars = [ 'p', 'pt' ,'q' ] )
    >>> reduced3 = tree.reduce  ( 'pt>1' , no_vars = [ 'Q', 'z' ,'x' ] )
    >>> reduced4 = tree.reduce  ( 'pt>1' , new_vars = { 'pt2' : 'pt*pt' } )
    >>> reduced5 = tree.reduce  ( 'pt>1' , new_vars = { 'pt2' : 'pt*pt' } , output = 'OUTPUT.root' )
    """

    reduced = ReduceTree ( tree                    ,
                           selection  = selection  ,
                           save_vars  = save_vars  ,
                           new_vars   = new_vars   ,
                           no_vars    = no_vars    ,
                           output     = output     ,
                           addselvars = addselvars ,
                           tmp_keep   = True       , 
                           silent     = silent     )

    from ostap.trees.trees import Chain
    
    chain = Chain ( reduced.chain )
    if not output : chain.trash.add ( reduced.output )  

    return chain


ROOT.TTree. reduce = reduce

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree ,    
    )

_new_methods_       = (
    ROOT.TTree.reduce ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
