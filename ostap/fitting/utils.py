#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/utils.py
#  Set of useful technical utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful technical utilities to build various fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2018-08-14"
__all__     = (
    'cov_qual'          , ## Quality of covariance matrix from Minuit
    'fit_status'        , ## Fit status from Minuit
    #
    'RangeVar'          , ## Helper class to temporary change a range for the variable 
    #
    'component_similar' , ## Should one use "similar" component?
    'component_clone'   , ## Should one use "cloned" component?
    'numcpu'            , ## number of CPUs
    'ncpu'              , ## fuction to build ROOT.RooFit.NumCPU
    #
    ## add/remove RooFit topic
    'remove_topic'      , ## remove topic from RooMsgService
    'add_topic'         , ## add    topic from RooMsgService
    'suppress_topics'   , ## suppress topics from RooMsgService 
    #
    'roo_gaussian'      , ## generate gaussian random number 
    'roo_poisson'       , ## generate Poisson random number
)
# =============================================================================
from   ostap.core.core         import ( Ostap   , rootID     , VE ,
                                        isequal , roo_silent )
from   ostap.core.ostap_types  import ( num_types      , list_types     ,
                                        integer_types  , string_types   ,
                                        is_good_number , sequence_types ,
                                        is_integer     )
from   ostap.math.random_ext   import ve_gauss, poisson
from   ostap.fitting.variables import SETVAR 
from   ostap.utils.basic       import numcpu , items_loop
import ostap.fitting.variables 
import ostap.fitting.roocollections
import ROOT, math, random 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.utils' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
ZERO = ROOT.RooFit.RooConst ( 0 ) 
ONE  = ROOT.RooFit.RooConst ( 1 ) 
# =============================================================================
## Generate Poisson random number coherentry  with RooFit
def roo_poisson ( mu ) :
    """ Generate Poisson random number coherentry  with RooFit
    """
    return ROOT.RooRandom.randomGenerator().Poisson ( mu )
# ==============================================================================
## Generate Gaussian random number coherentry  with RooFit
def roo_gaussian ( mu , sigma ) :
    """ Generate Gaussian random number coherentry  with RooFit
    """
    return ROOT.RooRandom.randomGenerator().Gaus   ( mu , sigma )
# =================================================================================
## make a name from prefix, name and suffix 
def make_name ( prefix , name = '' , suffix = '' ) :
    """ Make a name from prefix, name and suffix
    """
    
    prefix = prefix.strip ( '_ ' ) 
    suffix = suffix.strip ( '_ ' )
    name   = name  .strip ( '_ ' )

    if   prefix and name and suffix : return "%s_%s_%s" % ( prefix , name , suffix ) 
    elif prefix and suffix          : return "%s_%s"    % ( prefix ,        suffix ) 
    elif prefix and name            : return "%s_%s"    % ( prefix , name          ) 
    elif suffix and name            : return "%s_%s"    % (          name , suffix )

    return "%s" % ( name or prefix or suffix ) 
    
# =============================================================================
## MINUIT covariance matrix status:
# - status = -1 :  not available (inversion failed or Hesse failed)
# - status =  0 : available but not positive defined
# - status =  1 : covariance only approximate
# - status =  2 : full matrix but forced pos def
# - status =  3 : full accurate matrix
_cov_qual_ = {
    -1 :  '-1/not available (inversion failed or Hesse failed or externally provided)' ,
    0  :  ' 0/available but not positive defined',
    1  :  ' 1/covariance only approximate',
    2  :  ' 2/full matrix but forced pos def',
    3  :  ' 3/full accurate matrix',
    }
# =============================================================================
## MINUIT covariance matrix status:
# - status = -1 : not available (inversion failed or Hesse failed)
# - status =  0 : available but not positive defined
# - status =  1 : covariance only approximate
# - status =  2 : full matrix but forced pos def
# - status =  3 : full accurate matrix
def cov_qual ( status ) : return _cov_qual_.get( status , "%s" % status )
# =============================================================================
## Miniut::minimize status code
# - status = 1    : Covariance was made pos defined
# - status = 2    : Hesse is invalid
# - status = 3    : Edm is above max
# - status = 4    : Reached call limit
# - status = 5    : Any other failure
_fit_status_ = {
    0    : ' 0/success' ,
    1    : ' 1/Covariance was made pos defined',
    2    : ' 2/Hesse is invalid',
    3    : ' 3/Edm is above max',
    4    : ' 4/Reached call limit',
    5    : ' 5/Any other failure',
    }
# =============================================================================
## convert fit stats ointo the string 
def fit_status ( status ) :
    """ convert fit stats ointo the string """
    return _fit_status_.get( status ,"%s" % status )
                  
# =============================================================================
_nemax = 1000 ## number of events per CPU-core 
_ncmax =   16 ## maximal number of CPUs: there are some problems with >= 7
              ## @see https://sft.its.cern.ch/jira/browse/ROOT-4897
# ==============================================================================
_ncpus = []

# =============================================================================
## prepare "NumCPU" argument with reasonable choice of #cpu, depending on
#  number of events in dataset 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-03-31
def ncpu ( events ) :
    """Prepare 'NumCPU' argument with reasonable choice of #cpu, depending on
    the number of events in dataset 
    """
    #
    n_cores = numcpu() 
    if n_cores <= 1 : return ROOT.RooFit.NumCPU ( 1 ) ## fake!!! 
    #
    n  = events // _nemax
    if n       <= 1 : return ROOT.RooFit.NumCPU ( 1 ) ## fake!!! 
    #
    num = min ( n , n_cores , _ncmax )
    if not _ncpus : _ncpus.append ( num )   
    #
    return ROOT.RooFit.NumCPU ( num )

# =============================================================================
## helper class to temporary change a range for the variable 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class RangeVar(object) :
    """Helper class to temporary change a range for the variable 
    """
    def __init__( self , var , vmin , vmax ) :
        self.var  = var
        self.vmin = min ( vmin , vmax ) 
        self.vmax = max ( vmin , vmax )
        self.omin = self.var.getMin ()
        self.omax = self.var.getMax ()
        
    def __enter__ ( self ) :
        self.omin = self.var.getMin ()
        self.omax = self.var.getMax ()
        self.var.setMin ( self.vmin ) 
        self.var.setMax ( self.vmax )
        return self
    
    def __exit__  ( self , *_ ) :        
        self.var.setMin ( self.omin ) 
        self.var.setMax ( self.omax )



# ==============================================================================
## Should one use "similar" component?
def component_similar ( same ) :
    """Should one use 'similar' component?
    """
    if   same is Ellipsis           : return True
    elif same is NotImplemented     : return True
    elif isinstance ( same , str  ) \
         and same.strip().lower() in ( 'ditto' , 'similar' ) : return True
    return False

# =============================================================================
## Should one      "clone"  component?
def component_clone  ( same ) :
    """Should one use 'cloned' component?
    """    
    if isinstance ( same , str ) \
       and same.strip().lower() in ( 'clone' , 'cloned' , 'same' ) : return True        
    return False 
# =============================================================================

    
# =============================================================================
##  get <code>i</code>-th component from <code>what</code>
def get_i ( what , i , default = None ) :
    """
    """
    if   isinstance ( what , ROOT.RooArgList ) and i in what             : return what[i]
    elif isinstance ( what , ROOT.RooAbsReal ) and 0 == i                : return what 
    elif isinstance ( what , num_types       ) and 0 == i                : return what  
    elif isinstance ( what , list_types      ) and 0 <= i < len ( what ) : return what[i] 

    return default

# =============================================================================
levels = ( 'DEBUG' , 'INFO' , 'PROGRESS' , 'WARNING' , 'ERRROR' , 'FATAL' )
        
# =============================================================================
## consruct MsgTopic
#  @see RooFit::MsgTopic
#  @code
#  topic = msgTopic ( ROOT.RooFit.Fitting ) 
#  topic = msgTopic ( ROOT.RooFit.Fitting , ROOT.RooFit.Caching )
#  topic = msgTopic ( 'Fitting' , 'Caching' )
#  @endcode
def msg_topic ( *topics ) :
    """Consruct MsgTopic
    >>> topic = msgTopic ( ROOT.RooFit.Fitting ) 
    >>> topic = msgTopic ( ROOT.RooFit.Fitting , ROOT.RooFit.Caching )
    >>> topic = msgTopic ( 'Fitting' , 'Caching' )
    """
    topic = 0
    for i in  topics : 
        if   isinstance ( i , integer_types )  : topic |= i
        elif isinstance ( i , string_types  )  :

            tt = i
            
            a, sep, b = i.partition(':')
            if   sep and b and a.upper() in levels : tt = b
            elif sep and a and b.upper() in levels : tt = a

            ii = tt.lower() 
            if   ii == 'generation'            : topic |=  ROOT.RooFit.Generation 
            elif ii == 'minimization'          : topic |=  ROOT.RooFit.Minimization
            elif ii == 'minization'            : topic |=  ROOT.RooFit.Minimization
            elif ii == 'plotting'              : topic |=  ROOT.RooFit.Plotting
            elif ii == 'fitting'               : topic |=  ROOT.RooFit.Fitting 
            elif ii == 'integration'           : topic |=  ROOT.RooFit.Integration 
            elif ii == 'linkstatemgmt'         : topic |=  ROOT.RooFit.LinkStateMgmt
            elif ii == 'eval'                  : topic |=  ROOT.RooFit.Eval
            elif ii == 'caching'               : topic |=  ROOT.RooFit.Caching
            elif ii == 'optimization'          : topic |=  ROOT.RooFit.Optimization
            elif ii == 'optimisation'          : topic |=  ROOT.RooFit.Optimization
            elif ii == 'objecthandling'        : topic |=  ROOT.RooFit.ObjectHandling
            elif ii == 'inputarguments'        : topic |=  ROOT.RooFit.InputArguments
            elif ii == 'tracing'               : topic |=  ROOT.RooFit.Tracing
            elif ii == 'contents'              : topic |=  ROOT.RooFit.Contents
            elif ii == 'datahandling'          : topic |=  ROOT.RooFit.DataHandling
            elif ii == 'numintegration'        : topic |=  ROOT.RooFit.NumIntegration
            elif ii == 'numericintegration'    : topic |=  ROOT.RooFit.NumIntegration
            elif ii == 'numericalintegration'  : topic |=  ROOT.RooFit.NumIntegration
            elif ii == 'fastevaluations'       : topic |=  ROOT.RooFit.FastEvaluations
            else : logger.error ( 'MsgTopic/1: unknown topic %s, skip' % i )
        else : logger.error ( 'MsgTopic/2: unknown topic %s/%s, skip' % ( i , type ( i ) ) )
        
    return topic 
    
# =============================================================================
# upgraded constructor for class Ostap::Utils::RemoveTopicsd
# @code
# with RemoveTopic ( [ 'Fitting' , 'Plotting' ] ) :
#    ... do something ...
# with RemoveTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
#    ... do something ...
# @endcode
# @see Ostap::Utils::AddTopic
# @see Ostap::Utils::RemoveTopic
def _rt_new_init_ ( self , topics , level = ROOT.RooFit.INFO , streams = -1  ) :
    """ Upgraded constructor for class Ostap::Utils::RemoveTopics
    >>> with RemoveTopic ( [ 'Fitting' , 'Plotting' ] ) :
    ...    ... do something ...
    >>> with RemoveTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
    ...    ... do something ...
    - see Ostap::Utils::AddTopic
    - see Ostap::Utils::RemoveTopic
    """
    if isinstance ( topics , integer_types ) and 0 < topics and topics <= 2**16 :
        return self._old_init_ ( topics , level , streams )    
    if isinstance ( topics , string_types  ) : topics = topics.split()
    topic = msg_topic ( *topics )
    return self._old_init_ ( topic , level , streams )

if not hasattr ( Ostap.Utils.RemoveTopic , '_old_init_' ) :
    Ostap.Utils.RemoveTopic._old_init_ = Ostap.Utils.RemoveTopic.__init__
    Ostap.Utils.RemoveTopic.__init__   = _rt_new_init_
    Ostap.Utils.RemoveTopic.__enter__  = lambda s : s
    Ostap.Utils.RemoveTopic.__exit__   = lambda s,*_ : s.exit() 
    
# =============================================================================
# upgraded constructor for class Ostap::Utils::AddTopic
# @code
# with AddTopic ( [ 'Fitting' , 'Plotting' ] ) :
#    ... do something ...
# with AddTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
#    ... do something ...
# @endcode
# @see Ostap::Utils::AddTopic
# @see Ostap::Utils::RemoveTopic
def _at_new_init_ ( self , topics , streams = -1  ) :
    """ Upgraded constructor for class Ostap::Utils::AddTopics
    >>> with RemoveTopic ( [ 'Fitting' , 'Plotting' ] ) :
    ...    ... do something ...
    >>> with RemoveTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
    ...    ... do something ...
    - see Ostap::Utils::AddTopic
    - see Ostap::Utils::RemoveTopic
    """
    if isinstance ( topics , integer_types ) and 0 < topics and topics <= 2**16 :
        return self._old_init_ ( topics , streams )
    
    if isinstance ( topics , string_types  ) : topics = [ topics ]

    topic = msg_topic ( *topics )
    return self._old_init_ ( topic , streams )
                
if not hasattr (  Ostap.Utils.AddTopic , '_old_init_' ) :
    Ostap.Utils.AddTopic._old_init_ = Ostap.Utils.AddTopic.__init__
    Ostap.Utils.AddTopic.__init__   = _at_new_init_
    Ostap.Utils.AddTopic.__enter__  = lambda s : s
    Ostap.Utils.AddTopic.__exit__   = lambda s,*_ : s.exit() 
    
# ================================================================================
## remove topic from Roofit message streams
#  @see RooMsgService
#  @code
#  with remove_topic ( ROOT.RooFit.Fitting ) :
#    ...
#  with remove_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
#    ...
#  with remove_topic ( [ 'Fitting' , 'Plotting' ] ) :
#    ...
#  @endcode
#  @see Ostap::Utils::RemoveTopic
#  @see Ostap::Utils::AddTopic
def remove_topic ( topics , level = ROOT.RooFit.INFO , stream  = -1 ) :
    """Remove topic from Roofit message streams
    - see RooMsgService
    >>> with remove_topic ( ROOT.RooFit.Fitting ) :
    ...  ...
    >>> with remove_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
    ... ...
    >>> with remove_topic ( [ 'Fitting' , 'Plotting' ] ) :
    ... ...
    - see Ostap::Utils::RemoveTopic
    - see Ostap::Utils::AddTopic
    """
    return Ostap.Utils.RemoveTopic ( topics , level , stream ) 


# ================================================================================
## add topic from RooFit message streams
#  @see RooMsgService
#  @code
#  with add_topic ( ROOT.RooFit.Fitting ) :
#    ...
#  with add_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
#    ...
#  with add_topic ( [ 'Fitting' , 'Plotting' ] ) :
#    ...
#  @endcode
#  @see Ostap::Utils::RemoveTopic
#  @see Ostap::Utils::AddTopic
def add_topic ( topics , stream  = -1 ) :
    """Add topic to RooFit message streams
    - see RooMsgService
    >>> with add_topic ( ROOT.RooFit.Fitting ) :
    ...  ...
    >>> with add_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
    ... ...
    >>> with add_topic ( [ 'Fitting' , 'Plotting' ] ) :
    ... ...
    - see Ostap::Utils::RemoveTopic
    - see Ostap::Utils::AddTopic
    """
    return Ostap.Utils.AddTopic ( topics , level , stream ) 



# =============================================================================
## suppress certain message topics
#  @code
#  suppress_topics ( 'Fitting'  , 'Caching' ) 
#  @endcode 
def suppress_topics ( *topics ) :
    """suppress certain message topics
    >>> suppress_topics ( 'Fitting'  , 'Caching' ) 
    """
    if topics and 1 == len( topics ) :
        t = str ( topics [ 0 ] ).lower()
        if 'config' == t : return suppress_topics() 

    if not topics :
        newtopics = [] 
        import ostap.core.config as CONFIG
        if 'RooFit' in CONFIG.config :
            import string
            ws     = string.whitespace 
            node   = CONFIG.config [ 'RooFit' ]
            data   = node.get('RemoveTopics','' )
            topics = tuple ( i.strip ( ws ) for i in data.split ( ',' ) if i.strip ( ws ) ) 

    if topics :
        
        svc = ROOT.RooMsgService.instance()

        for topic in topics :

            level = ROOT.RooFit.INFO 
            tt    = topic 
            a, sep, b = topic.partition(':')
            if   sep and b and a.upper() in levels :
                tt    = b
                level = levels.index ( a.upper() )
            elif sep and a and b.upper() in levels :
                tt    = a
                level = levels.index ( b.upper() )

            svc.saveState ()
            
            num   = svc.numStreams()
            tt    = msg_topic ( tt  )
            
            for i in range ( num ) :
                ok = Ostap.Utils.remove_topic ( i , tt , level ) 

# =============================================================================
## and finally suppress exra RooFit topics! 
suppress_topics ()


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
#                                                                       The END 
# =============================================================================
