#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Baisc ROOT-bases utils 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Basi cROOT-based utils 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    ##
    'TakeIt'             , ## take and later delete ...
    'takeIt'             , ## take and later delete ...
    ##
    'Batch'              , ## context manager to keep/force certain ROOT ``batch''-mode
    'batch'              , ## context manager to keep/force certain ROOT ``batch''-mode
    'batch_env'          , ## chek&set the bacth from environment 
    ##
    'KeepCanvas'         , ## context manager to keep the current ROOT canvas
    'keepCanvas'         , ## context manager to keep the current ROOT canvas
    'InvisibleCanvas'    , ## context manager to use the invisible current ROOT canvas
    'invisibleCanvas'    , ## context manager to use the invisible current ROOT canvas
    ##
    'implicitMT'         , ## context manager to enable/disable implicit MT in ROOT 
    ##
    'ImplicitMT'         , ## context manager to enable/disable implicit MT in ROOT    
    ##
    'hadd'               , ## merge ROOT files using command `hadd`
    'hadd2'              , ## merge ROOT files using command `hadd`
    )

# =============================================================================
from   ostap.utils.timing import Wait 
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.root_utils' )
else                       : logger = getLogger( __name__                 )
del getLogger
# =============================================================================
## @class TakeIt
#  Take some object, keep it and delete at the exit
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2014-08-03    
class TakeIt(object):
    """ Take some object, keep it and delete at the exit
    
    >>> ds = dataset.reduce('pt>1')
    >>> with takeIt ( ds ) :
    ...
    
    """
    def __init__  ( self , other ) :
        self.other = other
        
    def __enter__ ( self ) :
        ROOT.SetOwnership ( self.other , True )
        return self.other
    
    def __exit__  ( self , *_ ) :

        o = self.other

        ## delete it! 
        del self.other
        
        if o and hasattr ( o , 'reset'  ) : o.reset  ()
        if o and hasattr ( o , 'Reset'  ) : o.Reset  ()
        if o and hasattr ( o , 'Delete' ) : o.Delete ()
        
        if o : del o
                            
# =============================================================================
## Take some object, keep it and delete at the exit
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2014-08-03    
def takeIt (  other ):
    """Take some object, keep it and delete at the exit
    >>> ds = dataset.reduce('pt>1')
    >>> with takeIt ( ds ) :
    ...    
    """
    return TakeIt ( other ) 


# =============================================================================
## context manager to keep ROOT ``batch'' state
#  @code
#  with Batch() :
#  ... do something here 
#  @endcode 
class Batch(object) :
    """ Context manager to keep ROOT ``batch'' state
    >>> with Batch() :
    ... do something here 
    """
    def __init__  ( self , batch = True ) :
        self.__batch = batch 
    ## contex manahger: ENTER
    def __enter__ ( self ) :
        groot = ROOT.ROOT.GetROOT()
        self.old_state = groot.IsBatch()
        if self.old_state != self.__batch : groot.SetBatch ( self.__batch ) 
        return self
    ## contex manager: EXIT
    def __exit__  ( self , *_ ) :
        groot = ROOT.ROOT.GetROOT()
        if self.old_state != groot.IsBatch() : groot.SetBatch( self.old_state ) 

# =============================================================================
## check batch from environmen variables, set it ans issue the message
def batch_env  ( logger = logger ) :
    """ chek&set the bacth from environment 
    - Check batch environmen variable
    - set ROOT.TROOT.SetBatch(True) 
    - issue the message 
    """
    groot = ROOT.ROOT.GetROOT()
    if not groot.IsBatch () :
        from   ostap.utils.env        import get_env , OSTAP_BATCH
        if get_env ( OSTAP_BATCH , '' ).lower() not in ( '' , '0' , 'not' , 'off' , 'false' ) :        
            groot.SetBatch ( True )
            logger.attention ( "BATCH processing is activated (environment)  " )                    
    return groot.IsBatch() 
        
# =============================================================================
## context manager to keep ROOT ``batch'' state
#  @code
#  with batch() :
#  ... do something here 
#  @endcode 
def batch( batch = True ) :
    """Context manager to keep ROOT ``batch'' state
    >>> with batch() :
    ... do something here 
    """
    return Batch ( batch )


# =============================================================================
## @class KeepCanvas
#  helper class to keep the current canvas
#  @code
#  with KeepCanvas() :
#  ... do something here 
#  @endcode 
class KeepCanvas(Wait) :
    """ Helper class to keep the current canvas
    >>> with KeepCanvas() :
    ... do something here 
    """
    def __init__ ( self , wait = 0 ) :
        Wait.__init__ ( self , after = wait ) 
        self.__old_canvas  = None 
    def __enter__ ( self ) :
        Wait.__enter__ ( self ) 
        import ROOT
        ## pad  = ROOT.TVirtualPad.Pad()
        pad     = ROOT.Ostap.Utils.get_pad()  
        cnv     = pad.GetCanvas()  if pad  else None 
        self.__old_canvas = cnv if cnv else None 
    def __exit__  ( self , *_ ) :
        Wait.__exit__ ( self , *_ )         
        if self.__old_canvas:
            self.__old_canvas.cd()
        self.__old_canvas = None             
    @property
    def old_canvas ( self ) :
        """``old_canvas'': canvas to be preserved"""
        return self.__old_canvas
    
# =============================================================================
#  Keep the current canvas
#  @code
#  with keepCanvas() :
#  ... do something here 
#  @endcode
def keepCanvas() :
    """Keep the current canvas
    >>> with keepCanvas() :
    ... do something here
    """
    return KeepCanvas()



# =============================================================================
## @class InvisibleCanvas
#  Use context ``invisible canvas''
#  @code
#  with InvisibleCanvas() :
#  ... do somehing here 
#  @endcode
class InvisibleCanvas(KeepCanvas) :
    """Use context ``invisible canvas''
    >>> with InvisibleCanvas() :
    ... do something here 
    """
    ## context manager: ENTER 
    def __enter__ ( self ) :
        ## start from keeping the current canvas 
        KeepCanvas.__enter__ ( self )
        ## create new canvas in batch mode 
        with Batch( True ) : 
            import ROOT 
            self.batch_canvas = ROOT.TCanvas()
            self.batch_canvas.cd ()
            return self.canvas

    ## context manager: EXIT
    def __exit__ ( self , *_ ) :
        if self.batch_canvas :
            self.batch_canvas.Close() 
            del self.batch_canvas             
        KeepCanvas.__exit__ ( self , *_ )

# =============================================================================
## Use context ``invisible canvas''
#  @code
#  with invisibleCanvas() :
#  ... do something here 
#  @endcode
def invisibleCanvas() :
    """ Use context ``invisible canvas''
    >>> with invisibleCanvas() :
    ... do something here 
    """
    return InvisibleCanvas() 

# =============================================================================
## EnableImplicitMT
#  Context manager to enable/disable implicit MT in ROOT 
#  @see ROOT::EnableImplicitMT 
#  @see ROOT::DisableImplicitMT 
#  @see ROOT::IsImplicitMTEnabled
#  @code
#  with ImplicitMT( True ) :
#  ...
#  @endcode
class ImplicitMT(object) :
    """ Context manager to enable/disable implicit MT in ROOT 
    >>> with ImplicitMT( True ) :
        ...
    - see ROOT::EnableImplicitMT 
    - see ROOT::DisableImplicitMT 
    - see ROOT::IsImplicitMTEnabled
    """
    def __init__  ( self , enable = True ) :

        if   isinstance ( enable , bool ) : 
            self.__enable   =        enable
            self.__nthreads =        0
        elif isinstance ( enable , int  ) and 0 <= enable : 
            self.__enable   = bool ( enable ) 
            self.__nthreads =        enable 
        else :
            raise  TypeError ( "ImplicitMT: invalid ``enable'' flag :%s/%s" % ( enable , type ( enable ) ) )

    @property
    def enable   ( self ) : return self.__enable
    @property
    def nthreads ( self ) : return self.__nthreads
    
    ## Context manager: ENTER 
    def __enter__ ( self ) :
            
        self.__initial = ROOT.ROOT. IsImplicitMTEnabled ()
        
        if bool ( self.__initial ) == bool ( self.enable ) : pass 
        elif self.enable : ROOT.ROOT.EnableImplicitMT  ( self.__nthreads )
        else             : ROOT.ROOT.DisableImplicitMT ()

        return self
    
    ## Context manager: EXIT
    def __exit__ ( self , *_ ) :

        _current = ROOT.ROOT.IsImplicitMTEnabled()

        if   _current == self.__initial : pass
        elif _current                   : ROOT.ROOT.DisableImplicitMT ()
        else                            : ROOT.ROOT.EnableImplicitMT  ()
            
# =============================================================================
## Context manager to enable/disable implicit MT in ROOT 
#  @see ROOT::EnableImplicitMT 
#  @see ROOT::DisableImplicitMT 
#  @see ROOT::IsImplicitMTEnabled
#  @code
#  with implicitMT( True ) :
#  ...
#  @endcode
def implicitMT ( enable = True ) :
    """ Context manager to enable/disable implicit MT in ROOT 
    >>> with implicitMT( True ) :
        ...
    - see ROOT::EnableImplicitMT 
    - see ROOT::DisableImplicitMT 
    - see ROOT::IsImplicitMTEnabled
    """
    return ImplicitMT ( enable ) 

# =============================================================================

# =========================================================================
## merge all files using <code>hadd</code> script from ROOT
#  @param output  name of the output merged file, if None,
#                 the temporary name will be generated,
#                 that will be deleted at the end of the session
#  @param opts   options for command <code>hadd</code>
#  @return the name of the merged file
# OPTIONS:
# -a                                   Append to the output
# -k                                   Skip corrupt or non-existent files, do not exit
# -T                                   Do not merge Trees
# -O                                   Re-optimize basket size when merging TTree
# -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
# -j                                   Parallelize the execution in multiple processes
# -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
# -d                                   Carry out the partial multiprocess execution in the specified directory
# -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
# -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
# -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
# -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
# -fk                                  Sets the target file to contain the baskets with the same compression
#                                      as the input files (unless -O is specified). Compresses the meta data
#                                      using the compression level specified in the first input or the
#                                      compression setting after fk (for example 206 when using -fk206)
# -ff                                  The compression level use is the one specified in the first input
# -f0                                  Do not compress the target file
# -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)  
def hadd ( files , output = None , dir = None , opts = "-ff -O" ) :
    """ Merge all files using <code>hadd</code> script from ROOT
    - `output`  name of the output merged file
    - `opts`   options for command <code>hadd</code>
    It returns the name of the merged file
    
    If no output file name is specified, the temporary name
    will be generate and the temporary file will be deleted
    at the end of the session
    
    OPTIONS:
    # -a                                   Append to the output
    # -k                                   Skip corrupt or non-existent files, do not exit
    # -T                                   Do not merge Trees
    # -O                                   Re-optimize basket size when merging TTree
    # -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
    # -j                                   Parallelize the execution in multiple processes
    # -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
    # -d                                   Carry out the partial multiprocess execution in the specified directory
    # -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
    # -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
    # -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
    # -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
    # -fk                                  Sets the target file to contain the baskets with the same compression
    #                                      as the input files (unless -O is specified). Compresses the meta data
    #                                      using the compression level specified in the first input or the
    #                                      compression setting after fk (for example 206 when using -fk206)
    # -ff                                  The compression level use is the one specified in the first input
    # -f0                                  Do not compress the target file
    # -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)                            
    """

    if isinstance ( files , string_types ) : files = [ files ]

    import glob
    all_files = []
    for p in files :
        all_files += [ f for f in glob.iglob ( p ) ]

    all_files = list  ( set ( all_files ) )
    all_files.sort()
    all_files = tuple ( all_files )
    
    if not output :
        import ostap.utils.cleanup as CU
        suffix = '.root'
        if files :
            base , suffix = os.path.splitext  ( all_files [0] )
            if base : base    = os.path.basename ( base   )
            if base : suffix  = '-%s%s' % ( base , suffix )            
        output = CU.CleanUp.tempfile ( prefix = 'ostap-hadd-merged-' , suffix = suffix , dir = dir )
            
                            
    cargs    = [ 'hadd' ] + opts.split() + [ output ] + [ f for f in all_files ]

    import subprocess
    subprocess.check_call ( cargs )
    
    if os.path.exists ( output ) and os.path.isfile ( output ) :
        return output 
    
    raise IOError ( "The output file %s does not exist!" % output )


# =============================================================================
def hadd2 ( args ) :

    if   isinstance ( args , string_types     ) : return hadd ( args )
    elif isinstance ( args , listlike_types   ) \
         and all ( ( isinstance ( i , string_types ) and i ) for i in args ) :
        return hadd ( args )    
    elif isinstance ( args , dictlike_types   ) :
        return hadd ( **args )

    return hadd  ( *args ) 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================
