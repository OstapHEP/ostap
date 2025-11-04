#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/chopping.py
#  `TMVA Chopper' - helper utility to train/use  TMVA using `chopping'
#
# `Chopping' is a jargon name for k-fold cross-validation technique
# @see https://machinelearningmastery.com/k-fold-cross-validation
# 
# The most frequest case:
# - TMVA is trained using the simulated events as `Signal' and realatively
# limited sample of data events (e.g. sidebands)  as `Backrgound'.
# To avoid using the same backrgound events for training and the final evaluation,
# `background' sample is split into `N' independent categories and
# `N'-independent TMVA trainings are performed.
# - For each training the corresponding category of backroung events is
# not used for this training.
# - For the final TMVA evaluation, for the events from category `i'
# the corresponding trained TMVA is used [By construction,  these events have
# not been used for the training of corresponding TMVA].
# 
# @attention For the large number of categories `N' it could be rather slow,
#            since too many TMVA's need to be trained)
#
# The interface is very similaer to TMVATrainer/TMVAReader, but one needs to specify
#  -  For training:
#  * <code>N</code>: number of categories 
#  * <code>category</code>: the string/expression with TTree variables
#                           that used to construct the category,
#                            e.g.  <code>'event'</code>. The actual expression used
#                            to get the category number is constructed as 
#                             <code>'(category)%N'</code>
#  - For reading
#  * <code>N</code> must be the same as above
#  * <code>categoryfunc</code>: python callable that gets category number, e.g.
#    <code>categoryfunc = lambda s : int(s.event)%N </code>
#
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-09-10
# =============================================================================
"""`TMVAChopper' - helper utility to train/use  TMVA using `chopping'

`Chopping' is a jargon name for k-fold cross-validation technique.
 - see e.g. https://machinelearningmastery.com/k-fold-cross-validation

Most frequest case:

TMVA is trained using the simulated events as `Signal' and realatively
limited sample of data events (e.g. sidebands)  as `Backrgound'.
To avoid using the same backrgound events for training and the final evaluation,
`background' sample is split into `N' independent categories and
`N'-independent TMVA trainings are performed.
For each training the corresponding category of backroung events is not
used for this training.
For the final TMVA evaluation, for the events from category `i' the corresponding
trained TMVA is used [By construction,  these events have not been used for
the training of corresponding TMVA].

- For the large number of categories it could be rather slow, since too many
TMVA's need to be trained

The interface is very similar to TMVATrainer/TMVAReader, but one needs to specify:

-  For training:
*  N        : number of categories 
*  category : the string/expression with TTree variables 
that used to construct the category, e.g.  'event'.

The actual expression used to get the category number is constructed as '(category)%N'

- For reading
* N            :   number of  categories (must be the same as above)
* categoryfunc : python callable that gets category number from TTree, e.g.: 

>>> N = 10 
>>> categoryfunc = lambda s : int(s.event)%N 

"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2017-09-10"
__all__     = (
    ##
    "Trainer"              , ## the `chopper' trainer for TMVA 
    "Reader"               , ## the `chopper' reader  for TMVA
    'addChoppingResponse'  , ## add `chopping' response to TTree or RooDataSet
    )
# =============================================================================
from   ostap.core.meta_info      import python_info, root_info
from   ostap.core.ostap_types    import ( integer_types  , string_types   ,
                                          num_types      , 
                                          sequence_types , dictlike_types )  
from   ostap.core.core           import VE, WSE 
from   ostap.core.pyrouts        import hID, h1_axis, Ostap 
from   ostap.utils.cleanup       import CleanUp
from   ostap.tools.tmva          import Trainer as TMVATrainer
from   ostap.tools.tmva          import Reader  as TMVAReader
from   ostap.tools.tmva          import ( dir_name      , good_for_negative ,
                                          trivial_opts  , make_tarfile      ,
                                          decode_vars   , NO_PROCESSING     )
from   ostap.utils.basic         import items_loop, typename 
from   ostap.stats.statvars      import data_statistic
from   ostap.utils.progress_conf import progress_conf 
from   ostap.utils.progress_bar  import progress_bar
from   ostap.utils.root_utils    import ImplicitMT 
from   ostap.utils.timing        import timing 
import ostap.trees.trees 
import ostap.trees.cuts
import ostap.utils.utils         as     Utils 
import ROOT, os, math, shutil, tarfile  
# =============================================================================
from ostap.logger.logger      import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.chopping' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
## @class Trainer
#  The `chopping'  trainer. The interface is very similar to TMVA Trainer
#  with two   additional mandatory parameters:
#  -  <code>N</code>        : number of categories 
#  -  <code>category</code> : the string/expression with TTree variables 
#  that used to construct the category, e.g.
#  <code>'event'</code> or <code>'137*event+813*run'</code> or 
#  The actual expression used to get the category number is constructed as
#  <code>'(category)%N'</code>
#
# - Book the trainer
# @code 
# >>> N = 11 
# >>> trainer = Trainer (
# ... category = '137*evt+813*run' ,
# ... N        = N                 , 
# ... methods =  [ # type                   name   configuration
# ...      ( ROOT.TMVA.Types.kMLP        , 'MLP'        , 'H:!V:EstimatorType=CE:VarTransform=N:NCycles=200:HiddenLayers=N+3:TestRate=5:!UseRegulator' ) ,
# ...      ( ROOT.TMVA.Types.kBDT        , 'BDTG'       , 'H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2' ) , 
# ...      ( ROOT.TMVA.Types.kCuts       , 'Cuts'       , 'H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart' ) ,
# ...      ( ROOT.TMVA.Types.kFisher     , 'Fisher'     , 'H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10' ),
# ...      ( ROOT.TMVA.Types.kLikelihood , 'Likelihood' , 'H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50' ) ] ,
# ... variables  = [ 'var1' , 'var2' ,  'var3' ] ,  ## Variables to use in the training
# ... signal     = signal_tree      , ## TTree/TChain with `signal' sample   
# ... background = background_tree  , ## TTree/TChain with `background' sample   
# ... name       = 'TMVAChopper'    ,
# ... verbose    = False )
# @endcode
#
# - Use the trainer
# @code
# >>> trainer.train()
# @endcode
# 
# - Get  results from the  trainer
# @code
# >>> weights_files = trainer.weights_files ## weights files (XML) 
# >>> class_files   = trainer.  class_files ## class files (C++)
# >>> output_files  = trainer. output_files ## output ROOT files 
# >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
# @endcode
class Trainer(object) :
    """ The `chopping'  trainer. The interface is very similar to TMVA Trainer
    with two   additional mandatory parameters:
    1. `N'        : number of categories 
    2. `category' : the string/expression with TTree variables 
    that used to construct the category, e.g.  'event'.
    The actual expression used to get the category number is constructed as
    '(category)%N'

    - Book the trainer: 
    >>> N = 11 
    >>> trainer = Trainer (
    ... category = '137*evt+813*run' ,
    ... N        = N                 , 
    ... methods =  [ # type                   name   configuration
    ...      ( ROOT.TMVA.Types.kMLP        , 'MLP'        , 'H:!V:EstimatorType=CE:VarTransform=N:NCycles=200:HiddenLayers=N+3:TestRate=5:!UseRegulator' ) ,
    ...      ( ROOT.TMVA.Types.kBDT        , 'BDTG'       , 'H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2' ) , 
    ...      ( ROOT.TMVA.Types.kCuts       , 'Cuts'       , 'H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart' ) ,
    ...      ( ROOT.TMVA.Types.kFisher     , 'Fisher'     , 'H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10' ),
    ...      ( ROOT.TMVA.Types.kLikelihood , 'Likelihood' , 'H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50' ) ] ,
    ... variables  = [ 'var1' , 'var2' ,  'var3' ] ,  ## Variables to use in the training
    ... signal     = signal_tree      , ## TTree/TChain with `signal' sample   
    ... background = background_tree  , ## TTree/TChain with `background' sample   
    ... name       = 'TMVAChopper'    ,
    ... verbose    = False )

    - Use the trainer
    >>> trainer.train()

    - Get  results from the  trainer
    >>> weights_files = trainer.weights_files ## weights files (XML) 
    >>> class_files   = trainer.  class_files ## class files (C++)
    >>> output_files  = trainer. output_files ## output ROOT files 
    >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
    """
    def __init__ ( self                              ,
                   category                          ,  ## accessor to category 
                   N                                 ,  ## number of categories 
                   methods                           ,  ## list of TMVA methods
                   variables                         ,  ## list of variables
                   ##
                   signal                            ,  ## signal tree
                   background                        ,  ## background tree
                   ##
                   signal_vars       = {}            ,  ## dictionary with new variables for signal sample 
                   background_vars   = {}            ,  ## dictionary with new variables for background sample 
                   ##                   
                   signal_cuts       = ''            ,  ## signal cuts 
                   background_cuts   = ''            ,  ## background cuts 
                   spectators        = []            ,
                   bookingoptions    = "Transformations=I;D;P;G,D" , 
                   configuration     = "SplitMode=Random:NormMode=NumEvents:!V" ,
                   signal_weight     = None          ,                
                   background_weight = None          ,
                   ## 
                   prefilter            = ''      ,  ## prefilter cuts before TMVA data loader
                   prefilter_signal     = ''      ,  ## separate prefilter for signal data 
                   prefilter_background = ''      ,  ## separate prefilter for background data                    
                   ##
                   signal_train_fraction     = -1 , ## fraction of signal events used for training     : 0<=f<1 
                   background_train_fraction = -1 , ## fraction of background events used for training : 0<=f<1
                   ##
                   prescale_signal      = 1       , ## prescale factor for signal 
                   prescale_background  = 1       , ## prescale factor for background 
                   ##
                   signal_add_vars      = {}      , ## add signal     vars, presumably fake/misisng  vars 
                   background_add_vars  = {}      , ## add backgriund vars, presumably fake/misisng  vars 
                   ##                    
                   name              = 'TMVAChopper'         ,   # the name 
                   verbose           = False                 ,   # verbose ? 
                   chop_signal       = False                 ,   # chop the signal     ?
                   chop_background   = True                  ,   # chop the background ?
                   logging           = True                  ,   # create log-files    ?
                   make_plots        = False                 ,   # make standard plots ?
                   ## 
                   control_plots_signal     = ()  , ## control plots for signal 
                   control_plots_background = ()  , ## control plot for background 
                   ##                    
                   workdir           = ''                    ,   # working directory   
                   multithread       = True                  ,   # use multithreading  ?
                   parallel          = True                  ,   # parallel training   ? 
                   parallel_conf     = {}                    ,   # parallel configuration ?
                   logger            = None                  ) : # logger to be used 
        
        """Create TMVA `chopping' trainer
        
        >>> N = 11 
        >>> trainer = Trainer (
        ... category = '137*evt+813*run' ,
        ... N        = N                 , 
        ... methods =  [ # type                   name   configuration
        ...      ( ROOT.TMVA.Types.kMLP        , 'MLP'        , 'H:!V:EstimatorType=CE:VarTransform=N:NCycles=200:HiddenLayers=N+3:TestRate=5:!UseRegulator' ) ,
        ...      ( ROOT.TMVA.Types.kBDT        , 'BDTG'       , 'H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2' ) , 
        ...      ( ROOT.TMVA.Types.kCuts       , 'Cuts'       , 'H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart' ) ,
        ...      ( ROOT.TMVA.Types.kFisher     , 'Fisher'     , 'H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10' ),
        ...      ( ROOT.TMVA.Types.kLikelihood , 'Likelihood' , 'H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50' ) ] ,
        ... variables  = [ 'var1' , 'var2' ,  'var3' ] ,  ## Variables to use in the training
        ... signal     = signal_tree      , ## TTree/TChain with `signal' sample   
        ... background = background_tree  , ## TTree/TChain with `background' sample   
        ... name       = 'TMVAChopper'    ,
        ... verbose    = False )
        
        """

        
        assert isinstance ( N , integer_types ) and 1 < N , "Invalid number of categories"

        self.__name   = name
        self.__logger = logger if logger else getLogger ( self.name )

        self.__verbose         = True if verbose    else False 
        self.__make_plots      = True if make_plots else False

        self.__chop_signal     = True if chop_signal     else False 
        self.__chop_background = True if chop_background else False 
        self.__parallel        = True if parallel        else False 
        self.__logging         = True if logging         else False

        self.__signal_train_fraction     = signal_train_fraction     if 0 <= signal_train_fraction     < 1 else -1.0 
        self.__background_train_fraction = background_train_fraction if 0 <= background_train_fraction < 1 else -1.0 

        assert ( isinstance ( prescale_signal     , integer_types ) and 1 <= prescale_signal         ) or \
               ( isinstance ( prescale_signal     , float         ) and 0 <  prescale_signal     < 1 )  , \
               "Invalid 'prescale_signal'"
        assert ( isinstance ( prescale_background , integer_types ) and 1 <= prescale_background     ) or \
               ( isinstance ( prescale_background , float         ) and 0 <  prescale_background < 1 )  , \
               "Invalid 'prescale_background'"
        
        self.__prescale_signal     = prescale_signal
        self.__prescale_background = prescale_background
        
        if self.parallel :
            from ostap.parallel.parallel import has_dill 
            if not has_dill :
                self.logger.attention  ("Disable parallel chopping due to old `dill` version ")
                self.__parallel = False
                
        self.__parallel_conf   = {}
        self.__parallel_conf.update ( parallel_conf )
        
        assert  self.__chop_signal or self.__chop_background, "Neither signal nor background chopping"         
        self.__category  = category 
        self.__N         = N

        ##
        ## Signal/background sources 
        ##

        ## RooDataSet as signal source 
        if isinstance ( signal     , ROOT.RooDataSet ) :
            self.logger.info ( 'Signal     converted from RooDataSet to TTree')
            signal.convertToTreeStore ()
            if signal.isWeighted() : 
                ws = Ostap.Utils.getWeight ( signal ) 
                if ws :
                    sw = ws if not signal_weight else signal_weight * ROOT.TCut ( ws )
                    signal_weight = sw
                    self.logger.info ( 'Redefine Signal     weight to be %s' % signal_weight )
            signal     = Chain ( signal.tree () ) 

        ## RooDataSet as background source 
        if isinstance ( background , ROOT.RooDataSet ) :
            self.logger.info ( 'Background converted from RooDataSet to TTree')
            background.convertToTreeStore ()
            if background.isWeighted() : 
                from ostap.core.core import Ostap             
                ## try to get the weight from dataset 
                ws = Ostap.Utils.getWeight ( background ) 
                if ws :
                    bw = ws if not background_weight else background_weight * ROOT.TCut ( ws )
                    backround_weight = bw 
                    self.logger.info ( 'Redefine Background weight to be %s' % background_weight )
            background = Chain ( background.tree () )

        from ostap.trees.trees import Chain
        source_types = Chain, ROOT.TTree
        
        assert isinstance ( signal     , source_types   ) or  \
            (  isinstance ( signal     , sequence_types ) and \
               all ( isinstance ( s    , source_types ) for s in signal ) ) , \
               "Invalid `signal` type %s" % typename ( signal ) 
        
        assert isinstance ( background , source_types   ) or  \
            ( isinstance  ( background , sequence_types ) and \
              all ( isinstance ( s , source_types ) for s in background  ) ) , \
              "Invalid `signal` type %s" % typename ( backdound ) 

        if  isinstance ( signal     , source_types ) : signal     = signal     , 
        if  isinstance ( background , source_types ) : background = background , 
        
        
        self.__signals           = tuple ( Chain ( o ) for o in signal ) 
        self.__signal            = self.__signals [ 0  ] 
        self.__more_signals      = self.__signals [ 1: ]

        self.__backgrounds       = tuple ( Chain ( o ) for o in background ) 
        self.__background        = self.__backgrounds [ 0  ] 
        self.__more_backgrounds  = self.__backgrounds [ 1: ]

        self.__methods              = tuple ( methods)
        
        self.__signal_weight        = signal_weight 
        self.__signal_cuts          = ROOT.TCut ( signal_cuts )     
        self.__background_weight    = background_weight 
        self.__background_cuts      = ROOT.TCut ( background_cuts ) 
        
        self.__prefilter            = ROOT.TCut ( prefilter   )
        self.__prefilter_signal     = prefilter_signal     if prefilter_signal     else '' 
        self.__prefilter_background = prefilter_background if prefilter_background else ''    
        
        variables                   = list ( variables ) 
        the_vars, _ , _sig_vars, _bkg_vars = decode_vars ( variables )
        
        if signal_vars     : _sig_vars.update ( signal_vars     )
        if background_vars : _bkg_vars.update ( background_vars )

        signal_vars      = _sig_vars 
        background_vars  = _bkg_vars 

        self.__signal_vars          = {}
        if signal_vars     : self.__signal_vars.update     ( signal_vars )
        
        self.__background_vars      = {}
        if background_vars : self.__background_vars.update ( background_vars ) 

        ## addtional signal variables 
        self.__signal_add_vars      = {}
        if signal_add_vars :
            addvars = signal_add_vars 
            assert isinstance ( addvars , dictlike_types ) , \
                "Invalid `signal_add_vars` %s` " % typename ( addvars) 
            assert all ( isinstance ( k , string_types ) for k in addvars.keys() ) , \
                "Invalid keys  in `signal_add_vars`"
            assert all ( isinstance ( v , string_types + num_types ) for v in addvars.values() ) , \
                "Invalid value in `signal_add_vars`"
            
            for k , v in addvars.items () :
                if    isinstance ( v , string_types  ) : vv =          v 
                elif  isinstance ( v , integer_types ) : vv = '%d'   % v 
                elif  isinstance ( v , num_types     ) : vv = '%.9g' % v
                self.__signal_add_vars [ k ] = vv 
                
        ## addtional background variables 
        self.__background_add_vars  = {}
        if background_add_vars :
            addvars = background_add_vars 
            assert isinstance ( addvars , dictlike_types ) , \
                "Invalid `background_add_vars` %s` " % typename ( addvars) 
            assert all ( isinstance ( k , string_types ) for k in addvars.keys() ) , \
                "Invalid keys  in `background_add_vars`"
            assert all ( isinstance ( v , string_types + num_types ) for v in addvars.values() ) , \
                "Invalid value in `background_add_vars`"
            
            for k , v in addvars.items () :
                if    isinstance ( v , string_types  ) : vv =          v 
                elif  isinstance ( v , integer_types ) : vv = '%d'   % v 
                elif  isinstance ( v , num_types     ) : vv = '%.9g' % v
                self.__background_add_vars [ k ] = vv 

                
        self.__variables         = tuple ( variables  )
        
        self.__spectators        = tuple ( s for s in spectators )
        self.__bookingoptions    = bookingoptions
        self.__configuration     = configuration        
        
        self.__sig_histos        = ()
        self.__bkg_histos        = ()
        
        dirname                  = dir_name ( self.name ) 
        self.__dirname           = dirname 
        self.__trainer_dirs      = [] 

        ##
        ## Signal & Background control plots
        ##
        
        self.__control_plots_signal     = control_plots_signal 
        self.__control_plots_background = control_plots_background 
        
        self.__multithread   = True if isinstance ( multithread , bool ) or ( isinstance ( multithread , int ) and 0 <= multithread ) else False 

        if not workdir : workdir = os.getcwd()

        if not os.path.exists ( workdir ) :
            from ostap.utils.basic import make_dirs
            make_dirs ( workdir )
            
        assert os.path.exists ( workdir ) and os.path.isdir ( workdir ), \
               'No valid working directory!'
        
        self.__workdir = os.path.abspath ( workdir ) 
        self.logger.info ("Working directory is %s" % self.__workdir )

        # =====================================================================
        ## prefilter if required 
        # =====================================================================
        
        ##if self.verbose :
        all_vars = list ( the_vars ) 
        ## for v in the_vars :
        ##     all_vars.append ( varexp )            
        for v in self.spectators :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )
            all_vars.append ( vv [ 0 ] )
            
        if self.prefilter         : all_vars.append ( self.prefilter )
        
        ## if self.signal_weight     : all_vars.append ( self.signal_weight     )
        ## if self.background_weight : all_vars.append ( self.background_weight )
        
        ## if self.signal_cuts       : all_vars.append ( self.signal_cuts       )
        ## if self.background_cuts   : all_vars.append ( self.background_cuts   )
        
        ## do not forget to process the chopping category index!
        all_vars.append ( self.category ) 

        if self.signal_weight :
            assert self.signal    .good_variables ( self.signal_weight )  or self.signal_add_vars        , \
                "Do no know how to treat `signal_weight`: %s" %  self.signal_weight
            assert self.background.good_variables ( self.signal_weight )  or self.background_add_vars    , \
                "Do no know how to treat `signal_weight`: %s" %  self.signal_weight
            
        if self.background_weight :
            assert self.signal    .good_variables ( self.background_weight ) or self.signal_add_vars     , \
                "Do no know how to treat `background_weight`: %s" %  self.background_weight
            assert self.background.good_variables ( self.background_weight ) or self.background_add_vars , \
                "Do no know how to treat `background_weight`: %s" %  self.background_weight
                    
        ## prefilter signal if required 
        if self.prefilter_signal     or \
           self.prefilter            or \
           self.signal_cuts          or \
           self.signal_vars          or \
           self.signal_add_vars      or \
           self.__more_signals       or 1 != self.prescale_signal :
            
            the_vars = list ( all_vars ) 
            for weight in ( self.signal_weight , self.background_weight ) :
                if not weight             : continue
                if   self.signal_add_vars : continue ## assume it defines the missnig/nesessary variabes 
                elif self.signal.good_variables ( weight ) : the_vars.append ( weight  )

            from ostap.trees.cuts import vars_and_cuts 
            for _ , w in self.control_plots_signal :  
                vv , _ , _ = vars_and_cuts ( w  , '' )                
                the_vars += vv 
                
            import ostap.trees.trees
            
            for k, v in items_loop ( self.signal_vars ) : the_vars.append ( v ) 
            avars = self.signal.the_variables ( the_vars )
            
            keys  = set   ( self.background_vars.keys () ) - set ( self.signal_vars    .keys () )
            avars = tuple ( sorted ( set ( avars ).union ( keys ) ) ) 

            scuts = {}
            if self.prefilter_signal : scuts.update ( { 'PreFilter/Signal' : self.prefilter_signal } )                    
            if self.prefilter        : scuts.update ( { 'PreFilter/Common' : self.prefilter        } )
            if self.signal_cuts      : scuts.update ( { 'Signal'           : self.signal_cuts      } )
            
            import ostap.frames.frames 
            import ostap.frames.tree_reduce       as TR
                        
            ## reduced signals
            with timing ( 'Pre-filter Signal     before processing' , logger = self.logger ) , ImplicitMT ( True ) :                 
                inputs   = ( self.signal, )  + self.__more_signals
                silent   = not self.verbose or not self.category in ( 0, -1 )
                new_vars = { }
                new_vars.update ( self.signal_vars     )
                new_vars.update ( self.signal_add_vars )
                self.__RS = tuple ( TR.reduce ( i                    ,
                                                selection = scuts    ,
                                                save_vars = avars    ,
                                                name      = 'SIGNAL' ,
                                                output    = CleanUp.tempfile ( suffix = '.root' , prefix = 'ostap-chopping-SIGNAL-' ) , 
                                                new_vars  = new_vars , 
                                                prescale  = self.prescale_signal ,  
                                                silent    = False    ) for i in inputs )
            ## merge signals 
            files = ()
            for rs in self.__RS : files += rs.files
            self.__SigTR  = Chain ( name = 'SIGNAL' , files = files )                 
            
            self.__signal          = self.__SigTR
            self.__signal_cuts     = ROOT.TCut() 
            self.__more_signals    = () 
            self.__signal_add_vars = {} 
     
        missvars = [ v for v in self.signal_vars if not v in self.signal ]
        assert not missvars , "Variables %s are not in signal sample!" % missvars                          
        
        assert not self.signal_weight     or self.signal.good_variables ( self.signal_weight     ) , \
            "The weight %s is not accessible in signal sample" % self.signal_weight 
        assert not self.background_weight or self.signal.good_variables ( self.background_weight ) , \
            "The weight %s is not accessible in signal sample" % self.backgound_weight 
  
        # =================================================================
        ## prefilter background if required 
        if self.prefilter_background     or \
           self.prefilter                or \
           self.background_cuts          or \
           self.backgruund_vars          or \
           self.background_add_vars      or \
           self.__more_backgrounds       or 1 != self.prescale_background  :
            
            the_vars = list ( all_vars ) 
            for weight in ( self.signal_weight , self.background_weight ) :
                if not weight : continue 
                if   self.background_add_vars  : continue ## assume it defines missing/nesessary variables 
                elif self.background.good_variables ( weight ) : the_vars.append ( weight  )

            from ostap.trees.cuts import vars_and_cuts 
            for _ , w in self.control_plots_background :
                vv , _ , _ = vars_and_cuts ( w , '' )                
                the_vars += vv 
                                
            import ostap.trees.trees

            for k, v in items_loop ( self.background_vars ) : the_vars.append ( v )             
            bvars = self.background.the_variables ( the_vars )
                        
            keys  = set   ( self.signal_vars.keys () ) - set ( self.background_vars    .keys () )
            bvars = tuple ( sorted ( set ( bvars ).union ( keys ) ) ) 

            bcuts = {}
            if self.prefilter_background : bcuts.update ( { 'PreFilter/Background' : self.prefilter_background } )                    
            if self.prefilter            : bcuts.update ( { 'PreFilter/Common'     : self.prefilter            } )
            if self.background_cuts      : bcuts.update ( { 'Background'           : self.background_cuts      } )
            
            import ostap.frames.frames 
            import ostap.frames.tree_reduce       as TR
                    
            ## reduced backgrounds 
            with timing ( 'Pre-filter Background before processing' , logger = self.logger ) , ImplicitMT ( True ) :
                silent   = not self.verbose or not self.category in ( 0, -1 )
                inputs   = ( self.background , ) + self.__more_backgrounds 
                new_vars = {}
                new_vars.update ( self.background_vars     ) 
                new_vars.update ( self.background_add_vars )                
                self.__RB = tuple ( TR.reduce ( i                  ,
                                                selection = bcuts  ,
                                                save_vars = bvars  ,
                                                name      = 'BACKGROUND'             ,
                                                output    = CleanUp.tempfile ( suffix = '.root' , prefix = 'ostap-chopping-BACKGROUND-' ) ,
                                                new_vars  = new_vars                 , 
                                                prescale  = self.prescale_background ,  
                                                silent    = False  ) for i in inputs )
            ## "merge" backgrounds 
            files = ()
            for rb in self.__RB : files += rb.files
            self.__BkgTR  = Chain ( name = 'BACKGROUND' , files = files ) 
            
            self.__background          = self.__BkgTR
            self.__background_cuts     = ROOT.TCut() 
            self.__more_backgrounds    = () 
            self.__background_add_vars = {}

        missvars = [ v for v in self.background_vars if not v in self.background ]
        assert not missvars , "Variables %s are not in background sample!" % missvars                         
                              
        assert not self.signal_weight     or self.background.good_variables ( self.signal_weight     ) , \
            "The weight %s is not accessible in background sample!" % self.signal_weight 
        assert not self.background_weight or self.background.good_variables ( self.background_weight ) , \
            "The weight %s is not accessible in background sample!" % self.background_weight 
          
        # =====================================================================
        ## check for signal weigths
        # =====================================================================
        if self.signal_weight :
            from ostap.stats.statvars import data_statistic
            sw = data_statistic ( self.signal , 
                                  self.signal_weight ,
                                  cuts      = self.signal_cuts ,
                                  progress  = False ,  
                                  use_frame = True  , 
                                  parallel  = True  )
            mn , mx = sw.minmax() 
            if mn < 0 : 
                ## there are negative weights :
                for m in self.methods :
                    if not m[0] in good_for_negative :
                        self.logger.error ( "Method `%s' does not support negative (signal) weights" % m[1] )
            
        # =====================================================================
        ## check for background weights
        # =================================================================
        if self.background_weight :
            from ostap.stats.statvars import data_statistic
            bw = data_statistic ( self.background , 
                                  self.background_weight , 
                                  cuts      = self.background_cuts , 
                                  progress  = False , 
                                  use_frame = True  , 
                                  parallel  = True  )
            mn , mx = bw.minmax() 
            if mn < 0 : 
                ## there are negative weights :
                for m in self.methods :
                    if not m[0] in good_for_negative :
                        self.logger.error ( "Method `%s' does not support negative (background) weights" % m[1] )
                        
        # =====================================================================
        ## Category population:
        # =====================================================================
        cat = '(%s)%%%d' % ( self.category , self.N  )        
        if self.chop_signal      :
            hs1 = ROOT.TH1F( hID() , 'Signal categories' , self.N * 5 , -0.5 , self.N - 1 ) 
            hs2 = h1_axis ( [ -0.5+i for i in range(   self.N + 1 ) ] , title = hs1.GetTitle() ) 
            self.signal.project     ( hs1 , cat , self.signal_cuts )
            self.signal.project     ( hs2 , cat , self.signal_cuts )
            self.__sig_histos = hs1   , hs2
            st = hs2.stat()
            if 0 >=  st.min()  : self.logger.warning ("Some signal categories are empty!")                 
            self.logger.info('Signal     category population mean/rms: %.1f/%-.1f min/max: %.1f/%-.1f' % ( st.mean() , st.rms() , st.min() , st.max() ) )
                        
        if self.chop_background  :
            hb1 = ROOT.TH1F( hID() , 'Background categories' , self.N * 5 , -0.5 , self.N - 1 ) 
            hb2 = h1_axis ( [ -0.5+i for i in range(   self.N + 1 ) ] , title = hb1.GetTitle() ) 
            self.background.project ( hb1 , cat , self.background_cuts )
            self.background.project ( hb2 , cat , self.background_cuts )
            self.__bkg_histos = hb1 , hb2
            ##
            st = hb2.stat()
            if 0 >=  st.min()  : self.logger.warning ("Some background categories are empty!")                 
            self.logger.info('Background category population mean/rms: %.1f/%-.1f min/max: %.1f/%-.1f' % ( st.mean() , st.rms() , st.min() , st.max() ) )
        
        ##  trick to please Kisa 
        from ostap.trees.trees import Chain
        if isinstance ( self.__signal     , ROOT.TTree ) :  self.__signal     = Chain ( self.__signal     ) 
        if isinstance ( self.__background , ROOT.TTree ) :  self.__background = Chain ( self.__background ) 

        ##  trick to please Kisa for ROOT 6.31/01  (enum TMVA.Types.ETMVA is not pickable...
        if ( 6 , 31 ) <= root_info : ## < ( 6 , 36 ) :
            ## @see https://github.com/root-project/root/issues/15104
            ms = list ( self.__methods )
            for i, e in enumerate  ( ms ) :
                e = list ( e )
                e = tuple ( [ int ( e [ 0 ] ) ] + e [ 1 : ]) 
                ms [ i ] = e
            self.__methods = ms
            
        ## book the trainers 
        self.__trainers      = () 
        self.__weights_files = []
        self.__class_files   = []
        self.__output_files  = []
        self.__tar_file      = None 
        self.__log_file      = None 

        self.__AUCs          = {} 
        
        if self.verbose :
            
            rows = [ ( 'Item' , 'Value' ) ]
            
            row  = 'Name'      , self.name
            rows.append ( row )

            row = 'Category'  , self.category
            rows.append ( row )
            
            row = 'N'         , '%d' % self.N
            rows.append ( row )

            row = 'Variables'      , ' '.join ( self.variables )
            rows.append ( row )

            if self.signal_vars :
                row = 'Signal vars'  , str ( self.signal_vars  ) 
                rows.append ( row )
                
            if self.background_vars :
                row = 'Background vars'  , str ( self.background_vars  ) 
                rows.append ( row )

            if self.spectators :
                row = 'Spectators' , ' '.join ( self.spectators )
                rows.append ( row )
            
            for i , o in enumerate ( [ o for o in self.bookingoptions.split ( ':' )  if not o in trivial_opts ] ) :
                if 0 == i : row = 'Booking options' , o
                else      : row = ''                , o 
                rows.append ( row )
                
            for i , o in enumerate ( [ o for o in self.configuration.split ( ':' ) if not o in trivial_opts  ] ) :
                if 0 == i : row = 'TMVA configuraton'    , o
                else      : row = ''                     , o 
                rows.append ( row )
  
            if self.prefilter :
                row = 'Prefilter'  , str ( self.prefilter ) 
                rows.append ( row )                
                              
            if self.prefilter_signal :
                row = 'Prefilter Signal' , str ( self.prefilter_signal ) 
                rows.append ( row )
                
            if self.prefilter_background :
                row = 'Prefilter Background' , str ( self.prefilter_background ) 
                rows.append ( row )
     
            if self.signal_cuts : 
                row = 'Signal cuts' , str ( self.signal_cuts ) 
                rows.append ( row )

            if self.signal_weight : 
                row = 'Signal weight' , str ( self.signal_weight ) 

            if 0 < self.signal_train_fraction < 1 :
                row = 'Signal train fraction' , '%.1f%%' % ( 100 *  self.signal_train_fraction ) 
                rows.append ( row )
                
            if 1 != self.prescale_signal :
                row = 'Signal prescale'       , '%s' % self.prescale_signal 
                rows.append ( row )
  
            if self.background_cuts : 
                row = 'Background cuts' , str ( self.background_cuts )
                rows.append ( row )

            if self.background_weight : 
                row = 'Background weight', str ( self.background_weight ) 
                rows.append ( row )

            if 0 < self.background_train_fraction < 1 :
                row = 'Backgroundtrain fraction' , '%.1f%%' % ( 100 *  self.background_train_fraction ) 
                rows.append ( row )
                
            if 1 != self.prescale_background:
                row = 'Backgronud prescale'      , '%s' % self.prescale_background
                rows.append ( row )

            ms  = [ m[1] for m in self.methods ]
            row = 'Methods' , ' '.join ( ms ) 
            rows.append ( row )
            
            from ostap.logger.colorized import allright 
            for m in self.methods :
                row = 'Method' , '%s Id:%s' % ( allright ( m[1] ) , m[0] )
                rows.append ( row )
                for i , o in enumerate ( [ o for o in m[2].split(':' ) if not o in trivial_opts ] ) :
                    if 0 == i : row = 'Method configruration' , o
                    else      : row =   ''                    , o 
                    rows.append ( row ) 
                
            import ostap.logger.table as T
            title = "Chopping Trainer %s created" % self.name 
            table = T.table (  rows , title = title , prefix = "# " , alignment = "lw" )
            self.logger.info ( "%s\n%s" % ( title , table ) )

        ## one more table (a short one with input samples summary:
        if self.verbose or self.signal_weight or self.background_weight :

            ES, EB, SW, BW = None, None , None , None
            cnf = {'use_frame' : True, 'parallel' : True , 'progress' :  False }

            from ostap.stats.statvars import data_statistic
            sc = ROOT.TCut      ( self.    signal_cuts )
            bc = ROOT.TCut      ( self.background_cuts )
            ss = data_statistic ( self.signal     , '1' , cuts = sc , **cnf )
            sb = data_statistic ( self.background , '1' , cuts = bc , **cnf )
            NS = ss.nEntries()
            NB = sb.nEntries()
            
            if self.    signal_weight :
                scw  = sc * self.    signal_weight
                sw   = data_statistic ( self.signal , '1' , cuts = scw , **cnf )
                sw   = WSE ( sw ) 
                ES   = sw.nEff () 
                SW   = VE  ( sw.sumw () , sw.sumw2() )
                
            if self.background_weight :
                bcw  = bc * self.background_weight                    
                bw   = data_statistic ( self.background  , '1' , cuts = bcw , **cnf )
                bw   = WSE ( bw ) 
                EB   = bw.nEff () 
                BW   = VE  ( bw.sumw () , bw.sumw2() ) 
                
            from ostap.logger.symbols import times, sum_symbol        
            rows = [ ( 'Sample' , '#Events' , '%sw' % sum_symbol , '' , '#Eff' , '' ) ]
            
            if isinstance ( SW , VE ) or ES :
                SW = VE ( SW )
                ES = VE ( ES ) 
                s1 , e1  = SW.pretty_print ( precision = 4 , width = 6 , parentheses = False )
                s2 , e2  = ES.pretty_print ( precision = 4 , width = 6 , parentheses = False )
                row = 'Signal' , '%d' % NS , \
                    s1 , '%s10^{%+d}' % ( times , e1 ) if e1 else '' , \
                    s2 , '%s10^{%+d}' % ( times , e2 ) if e2 else '' 
            else :
                row = 'Signal' , '%d' % NS
                
            rows.append ( row )
            
            if isinstance ( BW , VE ) or EB :
                BW = VE ( BW )
                EB = VE ( EB )
                s1 , e1  = BW.pretty_print ( precision = 4 , width = 6 , parentheses = False )
                s2 , e2  = EB.pretty_print ( precision = 4 , width = 6 , parentheses = False )
                row = 'Background' , \
                    s1 , '%s10^{%+d}' % ( times , e1 ) if e1 else '' , \
                    s2 , '%s10^{%+d}' % ( times , e2 ) if e2 else '' 
            else :
                row = 'Background' , '%d' % NB
                
            rows.append ( row )

            import ostap.logger.table as T
            rows  = T.remove_empty_columns ( rows ) 
            title = "Input samples"
            table = T.table ( rows , prefix = '# ' , title = title , alignment = "lccc" )
            self.logger.info ( '%s:\n%s' % ( title , table ) ) 

    ## create all trainers 
    def __create_trainers ( self ) :
        if self.trainers : self.logger.debug ('Remove existing trainers ')
        self.__trainers = [] 
        for i in  range ( self.N ) : self.__trainers.append ( self.create_trainer ( i ) )
        self.__trainers     = tuple ( self.__trainers ) 
        self.__trainer_dirs = [ t.dirname for t in self.trainers ]
        
    ## create the trainer for category "i"
    def create_trainer ( self , i , verbose = True ) :
        """ Create the trainer for category `i'
        """
        cat       = '(%s)%%%d' % ( self.category , self.N  )
        nam       =  '%s_%03d' % ( self.name , i )
        scuts     = self.    signal_cuts 
        bcuts     = self.background_cuts 
        icategory = "(%s)!=%d" % ( cat , i ) 
        if self.chop_signal     :
            scuts = icategory * scuts if scuts else icategory
        if self.chop_background :
            bcuts = icategory * bcuts if bcuts else icategory

        if self.parallel :
            mp = self.make_plots and ( self.verbose or 0 == i )        
            vb = self.verbose    and (                 0 == i ) 
        else :
            mp = self.make_plots 
            vb = self.verbose   
        
        vv = set ( self.variables )
        for v in self.signal_vars     : vv.add ( v ) 
        for v in self.background_vars : vv.add ( v )
        vv = sorted ( vv )
        
        t  = TMVATrainer ( methods           = self.methods            ,
                           variables         = vv                      ,
                           ##
                           signal            = self.__signal           ,
                           background        = self.__background       ,
                           ##
                           spectators        = self.spectators         ,
                           bookingoptions    = self.bookingoptions     ,
                           configuration     = self.configuration      ,
                           signal_weight     = self.signal_weight      ,
                           background_weight = self.background_weight  ,
                           ##
                           prefilter         = ''                      , ## attention! 
                           ##
                           signal_train_fraction     = self.signal_train_fraction     ,  
                           background_train_fraction = self.background_train_fraction ,
                           ##
                           signal_add_vars     = self.signal_add_vars     ,
                           background_add_vars = self.background_add_vars ,                           
                           ##
                           output_file       = ''                      , 
                           ##                           
                           signal_cuts       = scuts                   , 
                           background_cuts   = bcuts                   ,
                           ##
                           control_plots_signal     = self.control_plots_signal     , ## control plot for signal 
                           control_plots_background = self.control_plots_background , ## control plot for background 
                           ##                            
                           name              = nam                     ,
                           verbose           = vb                      , ## verbose 
                           logging           = self.logging            ,
                           make_plots        = mp                      , ## make plots 
                           multithread       = self.multithread        ,
                           category          = i                       ,
                           workdir           = self.workdir            )
        
        return t
    
    @property
    def name    ( self ) :
        """`name'    : the name of TMVA chopper"""
        return self.__name

    @property
    def logger ( self ) :
        """`logger' : logger instance for the trainer"""
        return self.__logger
    
    @property
    def dirname  ( self ) :
        """`dirname'  : the output directory name"""
        return str(self.__dirname) 

    @property
    def trainers ( self ) :
        """`trainers' - the  actual list of N-TMVA  trainers for N-categories"""
        return self.__trainers 

    @property
    def trainer_dirs ( self ) :
        """`trainer_dirs' - list of direcorries with trainer's output"""
        return tuple(self.__trainer_dirs)

    @property
    def category ( self ) :
        """`category' -  the accessor(string) to the category"""
        return self.__category 

    @property
    def parallel ( self ) :
        """`parallel' : use parallelisation for training"""
        return self.__parallel

    @property
    def parallel_conf ( self ) :
        """`parallel_conf' : configuration for parallel processing"""
        return self.__parallel_conf
    
    @property
    def logging  ( self ) :
        """`logging' : create the log-files"""
        return self.__logging
    
    @property
    def make_plots ( self ) :
        """`make_plots' : make standard TMVA plots?"""
        return self.__make_plots
    
    @property
    def N        ( self ) :
        """`N' - number of categories for chopping"""
        return self.__N

    @property
    def chop_signal ( self ) :
        """`chop_signal' : use chopping for signal?"""
        return self.__chop_signal
    
    @property
    def chop_background ( self ) :
        """`chop_background' : use chopping for background?"""
        return self.__chop_background
    
    @property
    def methods ( self ) :
        """`methods' : the list of TMVA methods to be used"""
        return tuple(self.__methods)
    
    @property
    def variables ( self ) :
        """`variables' : the list of variables  to be used for training"""
        return tuple(self.__variables)

    @property
    def spectators ( self ) :
        """`spectators' : the list of spectators to be used"""
        return tuple(self.__spectators)

    @property
    def signal ( self ) :
        """`signal' :  TTree for signal events"""
        return self.__signal.chain
    
    @property
    def signal_cuts ( self ) :
        """`signal_cuts' :  cuts to be applied for `signal' sample"""
        return ROOT.TCut(self.__signal_cuts)

    @property
    def signal_weight ( self ) :
        """`signal_weight' : weight to be applied for `signal' sample"""
        return self.__signal_weight
 
    @property
    def background ( self ) :
        """`background' :  TTree for background events"""
        return self.__background.chain 
    
    @property
    def background_cuts ( self ) :
        """`background_cuts' :  cuts to be applied for `backgroud' sample """
        return ROOT.TCut(self.__background_cuts)
    
    @property
    def background_weight ( self ) :
        """`background_weight' : weight to be applied for `background' sample"""
        return self.__background_weight

    @property
    def signal_vars ( self ) :
        """'signal_vars' :  variables for 'signal' sample"""
        result = {}
        if self.__signal_vars : result.update ( self.__signal_vars ) 
        return result
    
    @property
    def background_vars ( self ) :
        """'background_vars' :  variables for `backgroun' sample"""
        result = {}
        if self.__background_vars : result.update ( self.__background_vars ) 
        return result 

    @property
    def signal_add_vars ( self ) :
        """'signal_add_vars' :  variables to be added into signal sampoe"""
        return self.__signal_add_vars
    
    @property
    def background_add_vars ( self ) :
        """'background_add_vars' :  variables to be added into background sample"""
        return self.__background_add_vars

    @property
    def prefilter ( self ) :
        """`prefilter' : cuts ot be applied/prefilter before processing"""
        return self.__prefilter
    
    @property
    def prefilter_signal ( self ) :
        """'prefilter_signal' : cuts to be applied/prefilter for signal before processing"""
        return self.__prefilter_signal
    
    @property
    def prefilter_background ( self ) :
        """'prefilter_background' : cuts to be applied/prefilter for background before processing"""
        return self.__prefilter_background 

    @property
    def bookingoptions ( self ) :
        """`bookingoptions' : options used to book TMVA::Factory"""
        return str(self.__bookingoptions)
    
    @property
    def configuration ( self ) :
        """`configuration' : options used to book TMVA"""
        return str(self.__configuration)
    
    @property
    def verbose ( self ) :
        """`verbose' : verbosity  flag"""
        return self.__verbose

    @property
    def silent ( self ) :
        """`silent' : verbosity  flag"""
        return not self.verbose

    @property
    def signal_categories ( self ) :
        """`signal_categories' - two histograms(different binning) with signal category population"""
        return self.__sig_histos
    
    @property
    def signal_train_fraction ( self ) :
        """`signal_train_fraction': if non-negative, fraction of signal events used for training"""
        return self.__signal_train_fraction

    @property
    def background_train_fraction ( self ) :
        """`background_train_fraction': if non-negative, fraction of background events used for training"""
        return self.__background_train_fraction
    
    @property
    def prescale_signal ( self ) :
        """`prescale_signal': prescale the signal sample"""
        return self.__prescale_signal
    
    @property
    def prescale_background ( self ) :
        """`prescale_background': prescale the background sample"""
        return self.__prescale_background

    @property
    def background_categories ( self ) :
        """`background_categories' - two histograms(different binning) with background category population"""
        return self.__bkg_histos
    
    @property
    def weights_files ( self ) :
        """`weights_files': the list/tuple of final files with TMVA weights"""
        return tuple(self.__weights_files)
    @property
    def class_files ( self ) :
        """`class_files' : the list/tuple of final files with TMVA classes"""
        return tuple(self.__class_files)
    @property
    def output_files ( self ) :
        """`output_files': the output files """
        return tuple(self.__output_files)

    @property
    def tar_file ( self ) :
        """`tar_file': the compressed tar file"""
        return str(self.__tar_file) if self.__tar_file else None
    
    @property
    def log_file ( self ) :
        """`log_file': the compressed tar file with tarinng log-files"""
        return str(self.__log_file) if self.__log_file else None 

    @property 
    def multithread ( self ) :
        """`multithread' : make try to use multithreading in TMVA"""
        return self.__multithread 

    @property
    def workdir ( self ) :
        """`workdir' : working directory"""
        return self.__workdir

    @property
    def control_plots_signal     ( self ) :
        """`control_plots_signal` : control plots for signal sample 
        """
        return self.__control_plots_signal
    
    @property
    def control_plots_background     ( self ) :
        """`control_plots_background` : control plots for background sample 
        """
        return self.__control_plots_background

    @property
    def AUCs ( self ) :
        """`AUCs` : dictionary of dictionaries { category : { method : AUC } } , where
        - AUC is "Area Under ROC curce 
        """
        return self.__AUCs


    # =========================================================================
    ## The main method: training of all subsamples 
    #  - Use the trainer
    # @code 
    # >>> trainer.train()
    # @endcode
    #
    # - Get  results from the  trainer
    # @code
    # >>> weights_files = trainer.weights_files ## weights files (XML) 
    # >>> class_files   = trainer.  class_files ## class files (C++)
    # >>> output_files  = trainer. output_files ## output ROOT files 
    # >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
    # @endcode
    def train ( self ) :
        """ The main method: training of all subsamples 
        - Use the trainer
        >>> trainer.train()
        
        - Get  results from the  trainer
        >>> weights_files = trainer.weights_files ## weights files (XML) 
        >>> class_files   = trainer.  class_files ## class files (C++)
        >>> output_files  = trainer. output_files ## output ROOT files 
        >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
        """
        from ostap.utils.utils import keepCWD
        with keepCWD ( self.workdir ) :
            result = self._train()
            
        rows = [ ( 'Item' , 'Value' ) ]

        vv = set ( self.variables )
        for k in self.signal_vars     : vv.add ( k )
        for k in self.background_vars : vv.add ( k )
        vv = sorted ( vv ) 
        
        row  = 'Variables', ', '.join ( vv )
        rows.append ( row )
        
        if self.signal_vars :
            row = 'Signal vars'  , str ( self.signal_vars  ) 
            rows.append ( row )
            
        if self.background_vars :
            row = 'Background vars'  , str ( self.background_vars  ) 
            rows.append ( row )
            
        if self.spectators : 
            row  = 'Spectators', ', '.join( self.spectators )
            rows.append ( row )
            
        row = 'Methods' ,  ', '.join( [ i[1] for i in  self.methods ] )
        rows.append ( row )

        row = 'Category' , self.category 
        rows.append ( row )

        row = 'N' , "%d" % self.N
        rows.append ( row )
        
        if self.tar_file and os.path.exists ( self.tar_file  ) and tarfile.is_tarfile ( self.tar_file ) : 
            row = "Tar file" , self.tar_file 
            rows.append ( row )
            
        if self.log_file and os.path.exists ( self.log_file  ) and tarfile.is_tarfile ( self.log_file ) : 
            row = "Log file" , self.log_file 
            rows.append ( row )
            
        def flatten(t):
            return [ item for sublist in t for item in sublist]
        
        if self.weights_files : 
            wfiles = flatten ( self.weights_files )
            from ostap.utils.basic import commonpath
            wd = commonpath ( wfiles ) if 1 < len ( wfiles  ) else ''
            if wd :
                row = "Weights directory" , wd
                rows.append ( row )
            lw = len ( wd )
            for w in sorted ( wfiles ) :
                row = "Weights file" , w if not lw else w[lw+1:]
                rows.append ( row )

        if self.class_files : 
            cfiles = flatten ( self.class_files )
            from ostap.utils.basic import commonpath
            wd = commonpath ( cfiles ) if 1 < len ( cfiles  ) else ''
            if wd :
                row = "Class directory" , wd
                rows.append ( row )
            lw = len ( wd )
            for w in sorted ( cfiles ) :
                row = "Class  file" , w if not lw else w[lw+1:]
                rows.append ( row )

        import ostap.logger.table as T
        title = "Chopping %s outputs" % self.name 
        table = T.table (  rows , title = title , prefix = "# " , alignment = "lw" )
        self.logger.info ( "%s\n%s" % ( title , table ) ) 

        
        if self.verbose                          and \
               self.tar_file                     and \
               os.path.exists ( self.tar_file  ) and \
               tarfile.is_tarfile ( self.tar_file ) :
            self.logger.info ( 'Content of the tar file %s' % self.tar_file )
            with tarfile.open ( self.tar_file , 'r' ) as tar : tar.list ()

        if self.verbose                          and \
               self.log_file                     and \
               os.path.exists ( self.log_file  ) and \
               tarfile.is_tarfile ( self.log_file ) :
            self.logger.info ( 'Content of the tar/log file %s' % self.log_file )
            with tarfile.open ( self.log_file , 'r' ) as tar : tar.list ()
            
        return result 
        
            
    # =========================================================================
    ## The main method: training of all subsamples 
    #  - Use the trainer
    # @code 
    # >>> trainer.train()
    # @endcode
    #
    # - Get  results from the  trainer
    # @code
    # >>> weights_files = trainer.weights_files ## weights files (XML) 
    # >>> class_files   = trainer.  class_files ## class files (C++)
    # >>> output_files  = trainer. output_files ## output ROOT files 
    # >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
    # @endcode
    def _train ( self ) :
        """ The main method: training of all subsamples 
        - Use the trainer
        >>> trainer.train()
        
        - Get  results from the  trainer
        >>> weights_files = trainer.weights_files ## weights files (XML) 
        >>> class_files   = trainer.  class_files ## class files (C++)
        >>> output_files  = trainer. output_files ## output ROOT files 
        >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
        """
        result = self.p_train () if self.parallel else self.s_train ()

        if os.path.exists ( self.dirname ) and os.path.isdir ( self.dirname ) :
            try :
                shutil.rmtree ( self.dirname )
            except :
                pass
            
        try :
            os.mkdir ( self.dirname )
        except :
            pass

        if os.path.exists ( self.dirname ) and os.path.isdir ( self.dirname ) :
            
            if self.tar_file and tarfile.is_tarfile ( self.tar_file ) :            
                try :
                    shutil.move ( self.tar_file , self.dirname )
                    ntf = os.path.join ( self.dirname , self.tar_file )
                    if os.path.exists ( ntf ) and tarfile.is_tarfile ( ntf ) :
                        self.__tar_file = os.path.abspath ( ntf )
                except :
                    pass
                
            if self.log_file and tarfile.is_tarfile ( self.log_file ) :            
                try :
                    shutil.move ( self.log_file , self.dirname )
                    ntf = os.path.join ( self.dirname , self.log_file )
                    if os.path.exists ( ntf ) and tarfile.is_tarfile ( ntf ) :
                        self.__log_file = os.path.abspath ( ntf )
                except :
                    pass

            for tdir in self.trainer_dirs :
                if os.path.exists ( tdir ) and os.path.isdir ( tdir ) :
                    try :
                        shutil.move ( tdir , self.dirname )
                    except :
                        pass


        if self.verbose and self.AUCs :
            title = 'ROC/AUC summary'
            table = self.AUC_table ( prefix = '# ' , title = title )
            self.logger.info ( '%s:\n%s' % ( title , table ) )
            
        return self.tar_file


    # =========================================================================
    ## The main method: training of all subsamples sequentially 
    #  - Use the trainer for sequential training:
    # @code 
    # >>> trainer.s_train()
    # @endcode
    #
    # - Get  results from the  trainer
    # @code
    # >>> weights_files = trainer.weights_files ## weights files (XML) 
    # >>> class_files   = trainer.  class_files ## class files (C++)
    # >>> output_files  = trainer. output_files ## output ROOT files 
    # >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
    # @endcode
    def s_train ( self ) :
        """The main method: training of all subsamples sequentially 
        - Use the trainer for   sequential training 
        >>> trainer.train()
        
        - Get  results from the  trainer
        >>> weights_files = trainer.weights_files ## weights files (XML) 
        >>> class_files   = trainer.  class_files ## class files (C++)
        >>> output_files  = trainer. output_files ## output ROOT files 
        >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
        """
        ## create the trainers 
        self.__create_trainers()

        assert 1<= self.N and self.N == len ( self.trainers ), 'Invalid trainers!'
        
        weights  = []
        classes  = []
        outputs  = [] 
        tarfiles = [] 
        logfiles = []
        
        from ostap.utils.progress_bar import progress_bar
        for t in progress_bar ( self.trainers , silent = self.verbose ) :
            if self.verbose : self.logger.info  ( "Train the trainer `%s'" % ( t.name ) ) 
            t.train() 
            weights  += [ t.weights_files ] 
            classes  += [ t.  class_files ] 
            outputs  += [ t. output_file  ] 
            tarfiles += [ t.    tar_file  ] 
            logfiles += [ t.    log_file  ] if t.log_file and os.path.exists ( t.log_file ) else [] 
            self.__AUCs [ t.category ] = t.AUC
            
        self.__weights_files = tuple ( weights ) 
        self.__class_files   = tuple ( classes )
        self.__output_files  = tuple ( outputs )

        ## create the final tar-file 
        self.__tar_file      = make_tarfile ( output  = '.'.join ( [ self.name , 'tgz' ] ) ,
                                              files   = tarfiles     ,
                                              verbose = self.verbose , 
                                              tmp     = True         )
        
        if logfiles : 
            self.__log_file  = make_tarfile ( output  = '.'.join ( [ '%s_logs' % self.name , 'tgz' ] ) ,
                                              files   = logfiles     ,
                                              verbose = self.verbose , 
                                              tmp     = True         )  
        else : 
            self.__log_file  = ''

        return self.tar_file 

    # =========================================================================
    ## The main method: training of all subsamples in parallel  
    #  - Use the trainer for parallel training 
    # @code 
    # >>> trainer.p_train()
    # @endcode
    #
    # - Get  results from the  trainer
    # @code
    # >>> weights_files = trainer.weights_files ## weights files (XML) 
    # >>> class_files   = trainer.  class_files ## class files (C++)
    # >>> output_files  = trainer. output_files ## output ROOT files 
    # >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
    # @endcode

    ## use the parallel training 
    def p_train ( self ) :
        """ The main method: training of all subsamples in parallel  
        - Use the trainer for parallel training 
        >>> trainer.p_train()
        
        - Get  results from the  trainer
        >>> weights_files = trainer.weights_files ## weights files (XML) 
        >>> class_files   = trainer.  class_files ## class files (C++)
        >>> output_files  = trainer. output_files ## output ROOT files 
        >>> tar_file      = trainer.    tar_file  ## tar-file (XML&C++)
        """
        from ostap.parallel.parallel_chopping import chopping_training as _training_
        
        ##  train it!
        results = _training_ ( self , **self.parallel_conf )
        
        assert self.N == len ( results [ 0 ] ) , 'Invalid number of weights files'
        assert self.N == len ( results [ 1 ] ) , 'Invalid number of   class files'
        assert self.N == len ( results [ 2 ] ) , 'Invalid number of  output files'
        assert self.N == len ( results [ 3 ] ) , 'Invalid number of     tar files'
        assert self.N == len ( results [ 4 ] ) , 'Invalid number of     dir files'
        assert self.N == len ( results [ 5 ] ) , 'Invalid number of     log files'
        
        ## assert self.N == len ( results [ 6 ] ) , 'Invalid AUC results'
        
        weights  = [ i [ 1 ] for i in results [ 0 ]            ]
        classes  = [ i [ 1 ] for i in results [ 1 ]            ]
        outputs  = [ i [ 1 ] for i in results [ 2 ]            ]
        tarfiles = [ i [ 1 ] for i in results [ 3 ]            ]
        dirnames = [ i [ 1 ] for i in results [ 4 ]            ]
        logfiles = [ i [ 1 ] for i in results [ 5 ] if i [ 1 ] ]
        
        self.__weights_files = tuple ( weights ) 
        self.__class_files   = tuple ( classes )
        self.__output_files  = tuple ( outputs )

        ## get the AUCs from all trainers 
        for c, auc in results [ 6 ] : self.__AUCs [ c ] = auc
        
        ## create the final tar-file 
        self.__tar_file      = make_tarfile ( output  = '.'.join ( [ self.name , 'tgz' ] ) ,
                                              files   = tarfiles     ,
                                              verbose = self.verbose , 
                                              tmp     = True         )
        
        if logfiles and False : 
            self.__log_file  = make_tarfile ( output  = '.'.join ( [ '%s_logs' % self.name , 'tgz' ] ) ,
                                              files   = logfiles     ,
                                              verbose = self.verbose , 
                                              tmp     = True         )  
        else : 
            self.__log_file  = ''

        self.__trainer_dirs = tuple ( dirnames ) 
        
        return self.tar_file

    # =========================================================================
    ## build a summary AUC table
    #  @code
    #  trainer = ...
    #  print ( trainer.AUC_table() )
    #  @endcode 
    def AUC_table ( self , title = '' , prefix = '' , style = '' ) :
        """ Build a summary AUC table  
        >>> trainer = ...
        >>> print ( trainer.AUC_table() )
        """
        from ostap.logger.pretty  import pretty_float
        from ostap.logger.symbols import times
        header = 'Method' , '#' , 'AUC' , '1-AUC' , '' , '-log10(1-AUC)' 
        rows   = [] 
        for category , aucs in self.AUCs.items()  :
            for key, auc in aucs.items() : 
                da = 1 - auc 
                a1 , e1             = pretty_float ( da                 , precision = 4 , width = 6 )
                if 0 < da : a2 , e2 = pretty_float ( -math.log10 ( da ) , precision = 3 , width = 5 )
                else      : a2 , e2 = '' , 0             
                row = key , '%d' % category , '%.6f' % auc                              , \
                    a1 ,  '%s10^{%+d}' % ( times , e1 ) if e1 else '' , \
                    a2 ,  '%s10^{%+d}' % ( times , e2 ) if e2 else ''             
                rows.append ( row )
        rows = [ header ] + sorted ( rows ) 
        import ostap.logger.table as T
        title = title if title else "ROC/AUC compare"
        rows  = T.remove_empty_columns ( rows ) 
        return T.table ( rows , prefix = prefix , title = title , alignment = "lrcc" , style = '' )

# =============================================================================
## @class WeightFiles
#  helper structure  to deal with weights files
class WeightsFiles(CleanUp) :
    """ Helper structure  to deal with weights files
    """
    def __init__ ( self , weights_files ) :
        
        if isinstance ( weights_files , str  ) :
            
            wf  = weights_files
            import tarfile, os
            assert os.path.exists  ( wf ) and tarfile.is_tarfile ( wf ) , "Non-existing or invalid tarfile %s "  % wf
            
            with tarfile.open ( wf , 'r' ) as tar :
                logger.debug ( "Open tarfile %s" % wf )
                ## tar.list()
                tmpdir = self.tempdir ( prefix = 'ostap-chopping-weights-' )  
                self.trash.add ( tmpdir )
                args = { 'path' : tmpdir }
                if ( 3 , 12 ) <= python_info : args [ 'filter' ] = 'data'
                tar.extractall ( **args )                
                logger.debug ('Un-tar into temporary directory %s' % tmpdir ) 
                weights_files  = [ os.path.join ( tmpdir , i ) for i in tar.getnames() ]
                self.tmpfiles += weights_files
                
        self.__weights_files = weights_files

    @property
    def files   ( self ) :
        """`files': the weights file
        """
        import copy
        return copy.deepcopy ( self.__weights_files ) 
# =============================================================================
## @class Reader
#  The `chopping' TMVA reader.
#  The interface is very similar to TMVA Trainer
#  with two   additional mandatory parameters:
#  - <code>N</code>           : number of categories (must be the  same as for training) 
#  - <code>categoryfunc</code>: python callable that gets category number from TTree, e.g.:
#  @code 
#  >>> N = 11 
#  >>> categoryfun = lambda s : int(137*s.evt+813*s.run)%N
#  @endcode
#
# - Book the reader
# @code 
# >>> reader = Reader ( 
# ...    name          = 'CHOPPER' ,
# ...    categoryfunc  = categoryfun 
# ...    N             = N         ,
# ...    variables     = [ ('var1' , lambda s : s.var1 )   ,
# ...                      ('var2' , lambda s : s.var2 )   ,
# ...                      ('var3' , lambda s : s.var3 ) ] ,
# ...     weights_files = weights_files  )
# @endcode
# 
# Here <code>weights_files</code> can be :
# - single tar/tgz/tar.gz-file with weights files (output from <code>Trainer.tar_file</code>)
# - the structure of xml-files with weights       (output from <code>Trainer.weights_files</code>)
#
# - Use the reader
# @code 
# >>> tree =  ....  ## TTree/TChain/RooDataSet with data
# >>> for entry in tree :
# ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
# ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
# ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   )
# @endcode 
#
# - A bit more efficient form is :
# @code 
# >>> tree =  ....  ## TTree/TChain/RooDataSet with data
# >>> mlp_fun  =  reader.MLP
# >>> bdgt_fun =  reader.BDTG
# >>> for entry in tree :
# ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
# ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
# ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg))
# @endcode 
# - It it natually merges with Ostap's <code>SelectorWithVars</code> utility
class Reader(object) :
    """ The `chopping' TMVA reader.
    The interface is very similar to TMVA Trainer
    with two   additional mandatory parameters:
    1. `N'           : number of categories (must be the  same as for training) 
    2. `categoryfun' : python callable that gets category number from TTree, e.g.:

    >>> N = 11 
    >>> categoryfun = lambda s : int(137*s.evt+813*s.run)%N

    - Book the reader:
    >>> reader = Reader ( 
    ...    name          = 'CHOPPER' ,
    ...    categoryfunc  = categoryfun 
    ...    N             = N         ,
    ...    variables     = [ ('var1' , lambda s : s.var1 )   ,
    ...                      ('var2' , lambda s : s.var2 )   ,
    ...                      ('var3' , lambda s : s.var3 ) ] ,
    ..     weights_files = weights_files  )
    
    `weights_files' can be :
    - single tar/tgz/tar.gz-file with weights files (output from `Trainer.tar_file')
    - the structure of xml-files with weights       (output from `Trainer.weights_files')
    
    - Use the reader
    >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    >>> for entry in tree :
    ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
    ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
    ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   )

    - A bit more efficient form is :
    >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    >>> mlp_fun  =  reader.MLP
    >>> bdgt_fun =  reader.BDTG
    >>> for entry in tree :
    ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
    ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
    ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg))
    
    - It it natually merges with Ostap's `SelectorWithVars' utility     
    """
    def __init__ ( self                           , * , 
                   categoryfunc                   ,
                   N                              , 
                   variables                      ,
                   weights_files                  ,
                   spectators   = ()              , 
                   name         = 'ChopperReader' ,
                   options      = ''              ,
                   logger       = None            ,
                   verbose      = False           ) :
        """ Book the reader:
        >>> reader = Reader ( 
        ...    name          = 'CHOPPER' ,
        ...    categoryfunc  = categoryfun 
        ...    N             = N         ,
        ...    variables     = [ ('var1' , lambda s : s.var1 )   ,
        ...                      ('var2' , lambda s : s.var2 )   ,
        ...                      ('var3' , lambda s : s.var3 ) ] ,
        ..     weights_files = weights_files  )
        
        `weights_files' can be :
        - single tar/tgz/tar.gz-file with weights files (output from `Trainer.tar_file')
        - the structure of xml-files with weights       (output from `Trainer.weights_files')
        """

        assert isinstance ( N , integer_types ) and 1 <= N , "`N' is illegal %s/%s"  % ( N , type(N) )

        self.__name          = str(name)
        self.__logger        = logger if logger else getLogger ( self.name )
  
        self.__categoryfunc  = categoryfunc, 
        self.__N             = N
        
        self.__variables     = tuple ( sorted ( variables  ) ) 
        self.__spectators    = tuple ( sorted ( spectators ) ) 
        self.__methods       = []

        import copy
        self.__weights        = WeightsFiles ( weights_files )
        files = self.weights.files

        self.__weights_files = copy.deepcopy ( files )
        assert len ( self.weights_files ) == N , "Invalid length of `weights_files'"

        self.__readers   = []
        for i in range ( self.N ) :
            
            inam = '%s_%03d'   % ( self.name , i )            
            self.__readers.append ( TMVAReader ( name          = inam                  ,
                                                 variables     = self.variables        ,
                                                 spectators    = self.spectators       , 
                                                 weights_files = self.weights_files[i] ,
                                                 options       = options               ,
                                                 verbose       = verbose and i == 0    ) )
                
        self.__readers  = tuple   ( self.__readers )
        self.__histo    = h1_axis ( [ -0.5 + i for i in range ( self.N + 1 ) ] ,
                                    title ="Category population" )
        for r in self.__readers :
            mr = list ( r.methods )
            if not self.__methods : self.__methods = mr
            m1 = set ( mr )
            m2 = set ( self.__methods ) 
            assert m1 == m2 , "Inconsistent configuration of readers is detected"

        self.__methods = tuple (  self.__methods ) 
        
    @property
    def name ( self ) :
        """`name' - the name of Chopper Reader"""
        return self.__name

    @property
    def logger ( self ) :
        """`logger' : the logger instace for this Trainer"""
        return self.__logger
    
    @property
    def N ( self ) :
        """`N' - number of categories"""
        return self.__N
     
    @property
    def methods ( self ) :
        """`methods' - the list/tuple of booked TMVA methods"""
        return tuple (self.__methods ) 

    @property
    def variables ( self ) :
        """`variables' - the list variables with accessor functions, e.g.
        >>> variables = [ ## name      accessor  
        ...              ( 'pt'   , lambda s : s.pt ) ,
        ...              ( 'ip'   , lambda s : s.ip ) ,
        ...                'var1'                     ,   ## s.var1 will be used 
        ...                'var2'                     ] , ## s.var2 will be used 
        """
        return self.__variables
    
    @property
    def spectators ( self ) :
        """'spectators' : the list of spectators to be used"""
        return self.__spectators

    @property
    def weights_files( self ) :
        """`weight_files' - TMVA weight files"""
        return tuple(self.__weights_files)

    @property
    def weights ( self ) :
        """`weights' : input structure of weigth  files """
        return self.__weights
    
    @property
    def readers ( self ) :
        """`readers' -  the actual list/tuple of readers"""
        return self.__readers

    @property
    def categoryfunc ( self ) :
        """`categoryfunc' - the actual callable for the category classification, 
        e.g. for 11 categories:
        >>> categoryfun = lambda s : int(137*s.evt+813*s.run)%11          
        """
        return self.__categoryfunc[0]

    @property
    def histo ( self ) :
        """`histo': histogram with the category populations statistic"""
        return self.__histo 

    # =========================================================================
    ## Helper class to get the decision of `chopper'
    #  @code 
    #  reader = ...
    #  var = reader[ method ]
    #  val = var ( entry )
    #  @endcode
    class Method(TMVAReader.Var) :
        """ Helper class to get the decision of `chopper'
        >>> reader = ...
        >>> var = reader[ method ]
        >>> val = var ( entry )
        """
        def __init__ ( self , reader , method ) :
            TMVAReader.Var.__init__ (  self ,   reader , method )
            self.__N = reader.N 
        # =====================================================================
        ## Evaluate the chopper 
        #  @code
        #  tree = ...
        #  method = reader.MLP
        #  print('Response %s' % method ( tree ))
        #  @endcode
        #  or using categroy and parameters
        #  @code
        #  category  = 2
        #  pt , y, phi = ...
        #  print('Response %s' % method ( category , pt , y , phi ))
        #  @endcode 
        def __call__ ( self , arg ,  *args ) :
            """ Evaluate the chopper from TTree/TChain/RooAbsData:
            >>>tree = ...
            >>> method = reader.MLP
            >>> print('Response %s' % method ( tree ))
            Using category and parameters
            >>> category  = 2
            >>> pt , y, phi = ...
            >>> print('Response %s' % method ( category , pt , y , phi ))
            """
            if isinstance ( arg , int ) and args :
                category = arg
                return self.evaluate ( category , *args ) 
            return self.eval ( arg , *args )
        # =====================================================================
        ## Evaluate the method from parameters 
        # @code 
        # method       = ...
        # pt, eta, phi = 5 ,  3.0 , 0  ## variables
        # category     = 2 
        # print('Response is %s'    % method.evaluate ( category , pt ,  eta , phi ) )
        # @endcode 
        def evaluate ( self , category , *args ) :
            """ Evaluate the method from parameters 
            >>> method       = ...
            >>> pt, eta, phi = 5 ,  3.0 , 0  ## variables
            >>> category     = ...
            >>> print('Response is %s'    % method.evaluate ( category ,  pt ,  eta , phi ) )
            """
            return self.reader.evaluate ( category , self.method , *args )

        # ====================================================================
        ## Get the  mean over all categories
        #  @code
        #  method = ...
        #  pt, eta, phi = 5 ,  3.0 , 0  ## variables
        #  print ('Mean  response is %s' % method.mean (  pt , eta , phi ) )
        #  @endcode
        def mean ( self , *args ) :
            """ Get the  mean over all categories
            >>> method = ...
            >>> pt, eta, phi = 5 ,  3.0 , 0  ## variables
            >>> print('Mean  response is %s' % method (  pt , eta , phi ) )
            """
            sum = 0.0
            for i in range(self.__N) : sum += self.evaluate ( i , *args )
            return sum / float ( self.__N ) 
        
        # ====================================================================
        ## Get the full statistic over all categories
        #  @code
        #  method = ...
        #  pt, eta, phi = 5 ,  3.0 , 0  ## variables
        #  print('Response statistic is %s' % method.stat (  pt , eta , phi ) )
        #  @endcode        
        def stat ( self , *args ) :
            """ Get the full statistic over all categories
            >>> method = ...
            >>> pt, eta, phi = 5 ,  3.0 , 0  ## variables
            >>> print('Response statistic is %s' % method.stat (  pt , eta , phi ) )
            """
            from ostap.stats.counters import SE
            se = SE()
            for i in range(self.__N) : se += self.evaluate ( i , *args )
            return se 
                        
    ## =======================================================================
    ## helper utility to  get the correspondig function from the  reader:
    #  @code 
    #  >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    #  >>> mlp_fun  =  reader['MLP']  ## <-- here!
    #  >>> bdgt_fun =  reader['BDTG'] ## <-- here!
    #  >>> for entry in tree :
    #  ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
    #  ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
    #  ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg))
    # @endcode        
    def __getitem__ ( self , method ) :
        """ Helper utility to  get the correspondig function from the  reader:
        - Use the reader
        >>> tree =  ....  ## TTree/TChain/RooDataSet with data
        >>> mlp_fun  =  reader['MLP']  ## <-- here! 
        >>> bdgt_fun =  reader['BDTG'] ## <-- here!
        >>> for entry in tree :
        ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
        ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
        ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   )
        """
        if not method in self.__methods :
            return KeyError( 'No method %s is booked!' %  method )
        return Reader.Method  ( self , method )
    
    ## =======================================================================
    ## helper utility to  get the correspondig function from the  reader:
    #  @code 
    #  >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    #  >>> mlp_fun  =  reader.MLP  ## <-- here! 
    #  >>> bdgt_fun =  reader.BDTG ## <-- here!
    #  >>> for entry in tree :
    #  ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
    #  ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
    #  ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg))
    # @endcode        
    def __getattr__ ( self , method ) :
        """ Helper utility to  get the correspondig function from the  reader:
        - Use the reader
        >>> tree =  ....  ## TTree/TChain/RooDataSet with data
        >>> mlp_fun  =  reader.MLP  ## <-- here! 
        >>> bdgt_fun =  reader.BDTG ## <-- here!
        >>> for entry in tree :
        ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
        ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
        ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   )
        """                
        if not method in self.__methods :
            raise AttributeError( 'No method %s is booked!' %  method )
        return Reader.Method  ( self , method ) 

    # =========================================================================
    ## the main method - evaluate of TMVA from the certain category reader 
    #  - Use the reader
    #  @code 
    #  >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    #  >>> for entry in tree :
    #  ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
    #  ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
    #  ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   )
    #  @endcode 
    #  @attention it is *not* CPU efficient
    #  Ugly trick with arrays is needed due to some technical problems
    #  (actually TMVA reader needs the address of `float'(in C++ sense) variable
    def __call__ ( self , method , entry , cut_efficiency = 0.90 ) :
        """ The main method - evaluate of TMVA from the certain category reader 
        
        - Use the reader
        >>> tree =  ....  ## TTree/TChain/RooDataSet with data
        >>> for entry in tree :
        ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
        ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
        ...     print ( 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   ) 
        """
        icatfunc = self.__categoryfunc[0]
        ic       = icatfunc ( entry )
                 
        assert isinstance ( ic , integer_types ) and 0 <= ic < self.__N, \
               "Invalid `category' %s/%s" % ( ic ,  type ( ic ) )
        return self.__readers[ ic ] ( method ,  entry , cut_efficiency ) 
        
    # ========================================================================
    ## evaluate TMVA
    #  @code
    #  reader   = ...
    #  pt, y    = ...  ##
    #  category = ... 
    #  print( 'MLP response is: ', reader.evaluate ( category , 'MLP' , pt , y ))
    #  @endcode
    def evaluate ( self , category , method , *args ) :
        """ Evaluate TMVA
        >>> reader = ...
        >>> pt, y  = ...  ##
        >>> print('MLP response is: ', reader.evaluate ( category , 'MLP' , pt , y ))
        """
        assert isinstance ( category , integer_types ) and 0 <= category < self.__N, \
               "Invalid `category' %s/%s" % ( category ,  type ( category ) )
        return self.__readers[ category ].evaluate ( method , *args ) 
                                


# =============================================================================
## Specific action to ROOT.TTree
def _add_response_tree_ ( tree , * , 
                          chopper                  , 
                          category                 ,
                          N                        ,
                          inputs                   ,
                          weights                  ,
                          options    = ""          ,
                          prefix     = 'tmva_'     ,
                          suffix     = '_response' ,
                          spectators = ()   ,                          
                          aux        = 0.9  ,
                          progress   = True ,
                          report     = True ) :
    """ Specific action to ROOT.TTree
    """
    
    import ostap.trees.trees
    from   ostap.core.core   import   valid_pointer
    
    assert isinstance ( tree , ROOT.TTree ) and valid_pointer ( tree ) , \
        "Invalid TTree/TChain object"
    
    ## delegate to chain 
    if isinstance ( tree , ROOT.TChain ) and 2 <= tree.nFiles :
        return _add_response_chain_ ( tree     , 
                                      chopper    = chopper    , 
                                      category   = category   ,
                                      N          = N          ,
                                      inputs     = inputs     ,
                                      weights    = weights    ,
                                      spectators = spectators , 
                                      options    = options    ,
                                      prefix     = prefix     ,
                                      suffix     = suffix     ,
                                      aux        = 0.9        , 
                                      report     = report     ,
                                      progress   = progress   ) 
    
    
    from   ostap.core.core           import Ostap, ROOTCWD
    from   ostap.io.root_file        import REOPEN
    from   ostap.utils.root_utils    import implicitMT
    
    
    treepath = tree.fullpath
    the_file = tree.topdir    
    groot    = ROOT.ROOT.GetROOT() 
    assert treepath and the_file and ( not the_file is groot ) and isinstance ( the_file , ROOT.TFile ) , \
        'This is not the file-resident TTree* object! addition of new branch is not posisble!'

    the_file = the_file.GetFile ()    
    assert the_file and isinstance ( the_file, ROOT.TFile )  , \
        'This is not the file-resident TTree* object! addition of new branch is not posisble!'
    
    filename = the_file.GetName()
    
    branches = set ( tree.branches() ) | set (  tree.leaves () ) 

    matched = []    
    if prefix or suffix :
        matched = sorted ( v for v in branches if v.startswith ( prefix ) and v.endswith ( suffix ) ) 
        
    if matched or category in branches :
        matched = ','.join ( matched )         
        logger.warning ( "add_response_tree: Variables/Category '%s/%s' already in TTree, skip" % ( matched , category ) )
        return tree
    
    ## (5) display progress ? 
    progress = progress_conf ( progress )
    adder    = Ostap.AddTMVA ( progress ) 

    from ostap.math.base import strings
    if isinstance ( spectators , string_types ) : spectators = spectators,    
    _spectators = strings ( *spectators ) 
    
    with ROOTCWD () , REOPEN ( the_file  ) as tfile : 
        
        tfile.cd()
        assert tfile.IsWritable() , "The file `%s` is not writable!" % flename

        the_tree = tfile.Get ( treepath )
        assert valid_pointer ( the_tree ) and isinstance ( the_tree , ROOT.TTree ) , \
            'Invalid TTree:%s in file:%s' % ( treepath , filename  )
        
        ## run it! 
        sc = adder.addChoppingResponse ( the_tree    ,
                                         chopper     ,
                                         category    ,
                                         N           ,
                                         inputs      ,
                                         weights     ,
                                         _spectators , 
                                         options     , 
                                         prefix      ,
                                         suffix      ,
                                         aux         ) 

        assert sc.isSuccess () , 'Error from Ostap::AddTMVA::addChoppingResponse %s' % sc
        
        ## tfile.Write() ##  "" ) ## , ROOT.TFile.kOverwrite )
        with implicitMT ( False  ) : tfile.Write( "" , ROOT.TFile.kOverwrite )
        
        the_tree = ROOT.nullptr 

    chain = ROOT.TChain ( treepath )
    chain.Add  ( filename )

    if report :
        new_branches = sorted ( ( set ( chain.branches () ) | set ( chain.leaves () ) ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 >= n : title = "Added %s branch to TTree(%s)"   % ( n , treepath ) 
            else      : title = "Added %s branches to TTree(%s)" % ( n , treepath )  
            table = chain.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            chain = ROOT.TChain ( treepath )
            chain.Add  ( filename )     
            
    return chain

        
# =============================================================================d
## Specific action to ROOT.TChain
def _add_response_chain_ ( chain      , * ,
                           chopper    , 
                           category   ,
                           N          ,
                           inputs     ,
                           weights    ,
                           options    ,
                           prefix     = 'tmva_'      ,
                           suffix     = '_response'  ,
                           spectators = ()   ,                          
                           aux        = 0.9  , 
                           progress   = True ,
                           report     = True ) :
    """ Specific action to ROOT.TChain
    """
    
    import ostap.trees.trees
    from   ostap.core.core   import   valid_pointer

    assert isinstance ( chain , ROOT.TTree ) and valid_pointer ( chain ) , \
        "Invalid TChain/TTree!"

    if not isinstance ( chain , ROOT.TChain ) or chain.nFiles <= 1 :
        return _add_response_tree_ ( chain      ,
                                     chopper    = chopper    ,   
                                     category   = category   ,
                                     N          = N          ,
                                     inputs     = inputs     ,
                                     weights    = weights    ,
                                     spectators = spectators , 
                                     options    = options    ,
                                     prefix     = prefix     ,
                                     suffix     = suffix     ,
                                     aux        = 0.9        , 
                                     report     = report     ,
                                     progress   = progress   ) 
        
    files    = chain.files   
    treepath = chain.fullpath
    
    branches = set ( chain.branches () ) | set ( chain.leaves() ) if report else set()
    
    tree_progress  = progress and      len ( files ) < 5
    chain_progress = progress and 5 <= len ( files )
 
    for f in progress_bar ( files , len ( files ) , silent = not chain_progress ) :

        ch = ROOT.TChain ( treepath  )
        ch.Add ( f )
        ## treat the single tree 
        _add_response_tree_ ( ch  ,
                              chopper    = chopper       ,     
                              category   = category      ,
                              N          = N             ,
                              inputs     = inputs        ,
                              weights    = weights       ,
                              spectators = spectators    , 
                              options    = options       ,
                              prefix     = prefix        ,
                              suffix     = suffix        ,
                              aux        = aux           , 
                              report     = False         ,
                              progress   = tree_progress ) 
        
        
    chain = ROOT.TChain ( treepath )
    chain.Add  ( filename )

    if report :
        new_branches = sorted ( ( set ( chain.branches () ) | set ( chain.leaves () ) ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 >= n : title = "Added %s branch to TChain(%s)"   % ( n , treepath ) 
            else      : title = "Added %s branches to TChain(%s)" % ( n , treepath )  
            table = chain.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            chain = ROOT.TChain ( treepath )
            chain.Add  ( filename )     
            
    return chain
  

# =============================================================================
## Helper function to add TMVA/chopping response into dataset
#  @code
#  tar_file = trainer.tar_file
#  dataset  = ...
#  inputs   = [ 'var1' , 'var2' , 'var2' ] ## input variables to TMVA
#  dataset.addTMVAResponse ( dataset , chopper , inputs , tar_file , prefix = 'tmva_' )
#  @endcode
#  @param dataset input dataset to be updated
#  @param chopper       chopping category/formula
#  @param N             number of categories
#  @param inputs        input variables
#  @param weights_files files with TMVA weigths (tar/gz or xml)
#  @param category_name the category
#  @param prefix        prefix for TMVA-variable
#  @param suffix        suffix for TMVA-variable
#  @param options       options to be used in TMVA Reader
#  @param verbose       verbose operation?
#  @param aux           obligatory for the cuts method, where it represents the efficiency cutoff 
def addChoppingResponse ( dataset                     , ## input dataset to be updated
                          * , 
                          chopper                     , ## chopping category/formula 
                          N                           , ## number of categrories
                          inputs                      , ## input variables 
                          weights_files               , ## files with TMVA weigths (tar/gz or xml)
                          spectators    =  ()         , 
                          category_name = 'chopping'  , ## category name 
                          prefix        = 'tmva_'     , ## prefix for TMVA-variable         
                          suffix        = '_response' , ## suffix for TMVA-variable 
                          options       =  ''         , ## TMVA-reader options
                          verbose       = True        , ## verbosity flag
                          progress      = True        ,   ## verbosity flag
                          aux           = 0.9         ,   ## for Cuts method : efficiency cut-off                      
                          report        = True        ) : ## final report?
    """ Helper function to add TMVA/chopping  response into dataset
    >>> tar_file = trainer.tar_file
    >>> dataset  = ...
    >>> inputs   = [ 'var1' , 'var2' , 'var2' ] ## input variables to TMVA 
    >>> dataset.addChoppingResponse ( dataset , chopper ,  inputs , tar_file , prefix = 'tmva_' )
    """
    assert isinstance ( N , int ) and 1 < N < 10000 , 'Invalid "N" %s' % N

    assert category_name and ( prefix or suffix ) , \
        'addChoppingResponse: invalid category/prefix/suffix %s/%s' % ( category_name , prefix , suffix ) 
    
    if isinstance ( dataset , ROOT.TTree ) :
        import ostap.trees.trees        
        vars    = set ( dataset.branches() ) | set( dataset.leaves () ) 
        matched = sorted ( v for v in vars if v.startswith ( prefix ) and v.endswith (  suffix ) ) 
        if matched or category_name in vars :
            matched = ','.join ( matched )
            logger.warning ( "addChoppingResponse: Variables/Category '%s/%s' already in TTree, skip" % ( matched , category_name ) )
            return dataset
     
    ## decode inputs&weights
    
    from ostap.tools.tmva import _inputs2map_ , _weights2map_ , opts_replace
    
    _inputs = _inputs2map_  ( inputs )
    
    weights_files = WeightsFiles  ( weights_files )
    files         = weights_files.files
    files__       = [ _weights2map_ ( f ) for f in files ]
    files_        = [ f[0] for f in files__ ]
    weights_      = [ f[1] for f in files__ ]
    del files__
    
    from ostap.core.core import cpp, std, Ostap
    MAP   = std.map    ( 'std::string', 'std::string' )
    MAPS  = std.vector ( MAP ) 
    _maps = MAPS()
    for m in files_ : _maps.push_back( m ) 

    options = opts_replace ( options , 'V:'      ,     verbose )
    options = opts_replace ( options , 'Silent:' , not verbose )
    
    from ostap.utils.basic import isatty
    options = opts_replace ( options , 'Color:'  , verbose and isatty() )

    # =========================================================================
    ## TTree or TChain here
    # =========================================================================
    if   isinstance ( dataset , ROOT.TTree ) :        
        return _add_response_chain_ ( dataset    ,
                                      chopper    = chopper       ,
                                      category   = category_name ,
                                      N          = N             ,
                                      inputs     = _inputs       ,
                                      weights    = _maps         ,
                                      spectators = spectators    , 
                                      options    = options       ,
                                      prefix     = prefix        ,
                                      suffix     = suffix        ,
                                      aux        = aux           ,
                                      report     = report        ,
                                      progress   = progress      )

    # =========================================================================
    ## RooFit ?
    # =========================================================================
    assert isinstance ( dataset , ROOT.RooDataSet ) , "Invalid type of `dataset` %s" %  typename ( dataset ) 
   
    ## here we deal with RooAbsData
    
    if isinstance ( chopper , string_types ) :
        
        from   ostap.fitting.variables   import make_formula
        varset  = dataset.get()
        chopper = make_formula ( 'chopper' , chopper , varset )
            
    assert isinstance ( chopper , ROOT.RooAbsReal ), \
        'Invalid chopping type %s' % typename ( chopper ) 

    category = ROOT.RooCategory ( category_name , 'Chopping category: (%s)%%%d' %  ( chopper.GetTitle() , N ) )

    
    if   N <    10 : fmt = category_name + '_%d'    
    if   N <   100 : fmt = category_name + '_%02d'  
    elif N <  1000 : fmt = category_name + '_%03d'  
    elif N < 10000 : fmt = category_name + '_%04d'  
    else           : fmt = category_name + '_%d'   
    ## 
    for i in range ( N ) : category.defineType ( fmt % i  , i )

    
    progress = progress_conf ( progress )
    adder    = Ostap.AddTMVA ( progress )
    
    from ostap.math.base import strings
    if isinstance ( spectators , string_types ) : spectators = spectators,    
    _spectators = strings ( *spectators ) 
    
    sc = adder.addChoppingResponse ( dataset     ,
                                     chopper     ,
                                     category    ,
                                     N           ,
                                     _inputs     ,
                                     _maps       ,
                                     _spectators , 
                                     options     ,
                                     prefix      ,
                                     suffix      ,
                                     aux         )  
    assert sc.isSuccess () , 'Error from Ostap::AddTMVA::addChoppingResponse %s' % sc 

    return dataset 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
