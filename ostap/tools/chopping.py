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
from   ostap.core.meta_info    import root_info
from   ostap.tools.tmva        import Trainer as TMVATrainer
from   ostap.tools.tmva        import Reader  as TMVAReader
from   ostap.tools.tmva        import ( dir_name     , good_for_negative,
                                        trivial_opts ) 
from   ostap.core.core         import WSE 
from   ostap.core.pyrouts      import hID, h1_axis, Ostap 
from   ostap.core.ostap_types  import integer_types 
from   ostap.utils.cleanup     import CleanUp
import ostap.trees.trees 
import ostap.trees.cuts
import ostap.utils.utils       as     Utils 
from   ostap.fitting.variables import make_formula 
import ROOT, os, shutil, tarfile  
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
    """The `chopping'  trainer. The interface is very similar to TMVA Trainer
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
                   category                          ,   # accessor to category 
                   N                                 ,   # number of categories 
                   methods                           ,   # list of TMVA methods
                   variables                         ,   # list of variables
                   ##
                   signal                            ,   # signal tree
                   background                        ,   # background tree
                   ##
                   signal_vars       = {}       ,  ## dictionary with new variables for signal sample 
                   background_vars   = {}       ,  ## dictionary with new variables for background sample 
                   ##                   
                   signal_cuts       = ''            ,   # signal cuts 
                   background_cuts   = ''            ,   # background cuts 
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
                   name              = 'TMVAChopper'         ,   # the name 
                   verbose           = False                 ,   # verbose ? 
                   chop_signal       = False                 ,   # chop the signal     ?
                   chop_background   = True                  ,   # chop the background ?
                   logging           = True                  ,   # create log-files    ?
                   make_plots        = False                 ,   # make standard plots ?
                   workdir           = ''                    ,   # working directory   
                   multithread       = False                 ,   # use multithreading  ?
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
            from ostap.parallel.parallel import DILL_PY3_issue
            if DILL_PY3_issue :
                self.logger.warning ("Disable parallel chopping due to DILL_PY3_issue")
                self.__parallel = False
                
        if self.parallel and root_info < (6,15) :
            self.logger.warning ( "Parallel chopping is activated only for ROOT>6.15")
            self.__parallel = False
            
        self.__parallel_conf   = {}
        self.__parallel_conf.update ( parallel_conf )
        
        assert  self.__chop_signal or self.__chop_background, "Neither signal nor background chopping"         
        self.__category  = category 
        self.__N         = N

        assert signal     , 'Invalid Signal     is specified!'
        assert background , 'Invalid Background is specified!'

        from ostap.trees.trees import Chain      
        if   isinstance ( signal     , Chain           ) : pass 
        elif isinstance ( signal     , ROOT.TTree      ) : signal = Chain ( signal ) 
        elif isinstance ( signal     , ROOT.RooAbsData ) :
            signal.convertToTreeStore ()
            if signal.isWeighted() : 
                from ostap.core.core import Ostap
                ## try to get the weight from dataset 
                ws = Ostap.Utils.getWeight ( signal ) 
                if ws :
                    sw = ws if not signal_weight else signal_weight * ROOT.TCut ( ws )
                    signal_weight = sw
                    self.logger.info ( 'Redefine Signal     weight to be %s' % signal_weight )
            signal     = Chain ( signal.tree () )
                        
        if   isinstance ( background , Chain           ) : pass 
        if   isinstance ( background , ROOT.TTree      ) : background = Chain ( background ) 
        elif isinstance ( background , ROOT.RooAbsData ) :
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

        self.__signal            = signal     
        self.__background        = background 

        self.__methods           = tuple ( methods) 
        self.__signal_weight     = signal_weight 
        self.__signal_cuts       = ROOT.TCut ( signal_cuts )     
        
        self.__prefilter            = ROOT.TCut ( prefilter   )
        self.__prefilter_signal     = prefilter_signal     if prefilter_signal     else '' 
        self.__prefilter_background = prefilter_background if prefilter_background else ''    
        
        self.__background_weight = background_weight 
        self.__background_cuts   = ROOT.TCut ( background_cuts ) 


        variables                = list ( variables ) 
        
        _sig_vars = {}
        _bkg_vars = {}
        
        vars = []
        for v in variables  :
            a , s1 , b  = v.partition ( ':' )    ## a : b 
            a = a.strip ()
            b = b.strip ()
            if a and s1 and b :
                c , s2 , d = b.partition ( '?' ) ## a : c ? d 
                c = c.strip()
                d = d.strip()
                if c and s2 and d :
                    vars.append ( a )
                    signal_vars    .update ( { a : c } ) ## ATTENTION HERE 
                    background_vars.update ( { a : d } ) ## ATTENTION HERE 
                else :
                    vars.append ( a )
                    signal_vars    .update ( { a : b } )  ## ATTENTION HERE 
                    background_vars.update ( { a : b } )  ## ATTENTION HERE 
            else :
                vars.append ( v )
                
        variables = vars 
        variables.sort ()
        
        if signal_vars     : _sig_vars.update ( signal_vars     )
        if background_vars : _bkg_vars.update ( background_vars )

        signal_vars      = _sig_vars 
        background_vars  = _bkg_vars 

        self.__signal_vars          = {}
        if signal_vars     : self.__signal_vars.update     ( signal_vars )        
        self.__background_vars      = {}
        if background_vars : self.__background_vars.update ( background_vars ) 
        

        self.__variables         = tuple ( variables  )
        
        self.__spectators        = tuple ( spectators )

        self.__bookingoptions    = bookingoptions
        self.__configuration     = configuration
        
        self.__verbose           = True if verbose    else False 
        self.__make_plots        = True if make_plots else False
        
        self.__sig_histos        = ()
        self.__bkg_histos        = ()
        
        dirname                  = dir_name ( self.name ) 
        self.__dirname           = dirname 
        self.__trainer_dirs      = [] 

        self.__multithread       = multithread and (6,18)<= root_info

        if not workdir : workdir = os.getcwd()

        if not os.path.exists ( workdir ) :
            from ostap.utils.basis import make_dirs
            make_dirs ( workdir )
            
        assert os.path.exists ( workdir ) and os.path.isdir ( workdir ), \
               'No valid working directory!'
        
        self.__workdir = os.path.abspath ( workdir ) 
        self.logger.info ("Working directory is %s" % self.__workdir )

        
        # =====================================================================
        ## prefilter if required 
        # =====================================================================
        
        
        ##if self.verbose :
        all_vars = []           
        for v in self.variables  :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )
            
            varexp = vv [ 0 ].strip()  
            nick , sep , expr = varexp.partition ( ":=" )
            nick = nick.strip()
            expr = expr.strip() 
            if nick and sep and expr : varexp = expr
            else                     : nick   = varexp 
            
            if   varexp in self.signal_vars     or ( nick and nick in self.signal_vars     ) : continue
            elif varexp in self.background_vars or ( nick and nick in self.background_vars ) : continue 
            else  : all_vars.append ( varexp )
            
        for v in self.spectators :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )
            all_vars.append ( vv[0] )
            
        if self.prefilter         : all_vars.append ( self.prefilter         )                
        ## if self.signal_cuts       : all_vars.append ( self.signal_cuts       )
        ## if self.signal_weight     : all_vars.append ( self.signal_weight     )
        ## if self.background_cuts   : all_vars.append ( self.background_cuts   )
        ## if self.background_weight : all_vars.append ( self.background_weight )
        
        ## do not forget to process the chopping category index!
        all_vars.append ( self.category ) 
        
        ## prefilter signal if required 
        if self.prefilter_signal or self.prefilter or 1 != self.prescale_signal or self.signal_vars :

            if self.signal_weight : all_vars.append ( self.signal_weight )

            import ostap.trees.trees
            avars = self.signal.the_variables ( all_vars )
            
            keys  = set   ( self.background_vars.keys () ) - set ( self.signal_vars    .keys () )
            avars = tuple ( sorted ( set ( avars ).union ( keys ) ) ) 

            scuts = {}
            if self.prefilter_signal : scuts.update ( { 'PreFilter/Signal' : self.prefilter_signal } )                    
            if self.prefilter        : scuts.update ( { 'PreFilter/Common' : self.prefilter        } )
            if self.signal_cuts      : scuts.update ( { 'Signal'           : self.signal_cuts      } )
            
            if ( 6 , 24 ) <= root_info :
                import ostap.frames.frames 
                import ostap.frames.tree_reduce       as TR
            else :
                import ostap.parallel.parallel_reduce as TR
                
            silent = not self.verbose or not self.category in ( 0, -1 )
            self.logger.info ( 'Pre-filter Signal     before processing' )
            self.__SigTR = TR.reduce ( self.signal        ,
                                       selection = scuts  ,
                                       save_vars = avars  ,
                                       new_vars  = self.signal_vars     , 
                                       prescale  = self.prescale_signal ,  
                                       silent    = False  )

            self.__signal      = self.__SigTR
            self.__signal_cuts = ROOT.TCut() 
            
        # =================================================================
        ## prefilter background if required 
        if self.prefilter_background or self.prefilter or 1 != self.prescale_background or self.background_vars :
            
            if self.background_weight : all_vars.append ( self.background_weight )

            import ostap.trees.trees
            bvars = self.background.the_variables ( all_vars )
            
            keys  = set   ( self.signal_vars.keys () ) - set ( self.background_vars    .keys () )
            bvars = tuple ( sorted ( set ( bvars ).union ( keys ) ) ) 

            bcuts = {}
            if self.prefilter_background : bcuts.update ( { 'PreFilter/Background' : self.prefilter_background } )                    
            if self.prefilter            : bcuts.update ( { 'PreFilter/Common'      : self.prefilter            } )
            if self.background_cuts      : bcuts.update ( { 'Background'           : self.background_cuts      } )
            
            if ( 6 , 24 ) <= root_info :
                import ostap.frames.frames 
                import ostap.frames.tree_reduce       as TR
            else :
                import ostap.parallel.parallel_reduce as TR
                    
            silent = not self.verbose or not self.category in ( 0, -1 )
            self.logger.info ( 'Pre-filter Background before processing' )
            self.__BkgTR = TR.reduce ( self.background    ,
                                       selection = bcuts  ,
                                       save_vars = bvars  ,
                                       new_vars  = self.background_vars     , 
                                       prescale  = self.prescale_background ,  
                                       silent    = False  )
            
            self.__background      = self.__BkgTR
            self.__background_cuts = ROOT.TCut() 
            
        # =====================================================================
        ## check for signal weigths
        # =====================================================================
        if self.signal_weight :
            if ( 6 , 20 ) <= root_info : 
                from ostap.frames.frames import frame_statVar 
                sw = frame_statVar ( self.signal , self.signal_weight , self.signal_cuts )
            else :
                sw = self.signal    .statVar ( self.signal_weight     , self.signal_cuts )
            if isinstance ( sw , WSE ) : sw = sw.values()
            mn , mx = sw.minmax() 
            if mn < 0 : 
                ## there are negative weights :
                for m in self.methods :
                    if not m[0] in good_for_negative :
                        self.logger.error ( "Method `%s' does not support negative (signal) weights" % m[1] )
            
        # =====================================================================
        ## check for background weigths
        # =================================================================
        if self.background_weight :
            if ( 6 , 20 ) <= root_info :
                from ostap.frames.frames import frame_statVar 
                bw = frame_statVar ( self.background , self.background_weight , self.background_cuts )
            else :
                bw = self.background.statVar ( self.background_weight , self.background_cuts )
            if isinstance ( bw , WSE ) : bw = bw.values()
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
        
        ## book the trainers 
        self.__trainers      = () 
        self.__weights_files = []
        self.__class_files   = []
        self.__output_files  = []
        self.__tar_file      = None 
        self.__log_file      = None 

        
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

    ## create all trainers 
    def __create_trainers ( self ) :
        if self.trainers : self.logger.debug ('Remove existing trainers ')
        self.__trainers = [] 
        for i in  range ( self.N ) : self.__trainers.append ( self.create_trainer ( i ) )
        self.__trainers     = tuple ( self.__trainers ) 
        self.__trainer_dirs = [ t.dirname for t in self.trainers ]
        
    ## create the trainer for category "i"
    def create_trainer ( self , i , verbose = True ) :
        """Create the trainer for category `i'
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
            
        t  = TMVATrainer ( methods           = self.methods          ,
                           variables         = self.variables         ,
                           signal            = self.__signal          ,
                           background        = self.__background      ,
                           spectators        = self.spectators        ,
                           bookingoptions    = self.bookingoptions    ,
                           configuration     = self.configuration     ,
                           signal_weight     = self.signal_weight     ,
                           background_weight = self.background_weight ,
                           ##
                           prefilter         = ''                     , ## attention! 
                           ##
                           signal_train_fraction     = self.signal_train_fraction     ,  
                           background_train_fraction = self.background_train_fraction ,
                           ##
                           output_file       = ''                     , 
                           ##                           
                           signal_cuts       = scuts                  , 
                           background_cuts   = bcuts                  ,
                           ##
                           name              = nam                    ,
                           verbose           = vb                     , ## verbose 
                           logging           = self.logging           ,
                           make_plots        = mp                     , ## make plots 
                           multithread       = self.multithread       ,
                           category          = i                      ,
                           workdir           = self.workdir           )
        
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
        
        row  = 'Variables', ', '.join( self.variables )
        rows.append ( row )
        
        if self.spectators : 
            row  = 'Spectators', ', '.join( self.spectators )
            rows.append ( row )
            
        row = 'Methdos' ,  ', '.join( [ i[1] for i in  self.methods ] )
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


        self.__weights_files = tuple ( weights ) 
        self.__class_files   = tuple ( classes )
        self.__output_files  = tuple ( outputs )

        self.make_tarfile ( tarfiles , logfiles )

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
        """The main method: training of all subsamples in parallel  
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
        
        assert self.N == len ( results [0] ) , 'Invalid number of weights files '
        assert self.N == len ( results [1] ) , 'Invalid number of   class files '
        assert self.N == len ( results [2] ) , 'Invalid number of  output files '
        assert self.N == len ( results [3] ) , 'Invalid number of     tar files '
        assert self.N == len ( results [4] ) , 'Invalid number of     dir files '
        assert self.N == len ( results [5] ) , 'Invalid number of     log files '
        
        weights  = [ i[1] for i in results [0]         ]
        classes  = [ i[1] for i in results [1]         ]
        outputs  = [ i[1] for i in results [2]         ]
        tarfiles = [ i[1] for i in results [3]         ]
        dirnames = [ i[1] for i in results [4]         ]
        logfiles = [ i[1] for i in results [5] if i[1] ]
        
        self.__weights_files = tuple ( weights ) 
        self.__class_files   = tuple ( classes )
        self.__output_files  = tuple ( outputs )

        self.make_tarfile ( tarfiles , logfiles )

        self.__trainer_dirs = tuple ( dirnames ) 
        
        return self.tar_file
        
    ## create the tarfile from the list of tarfiles 
    def make_tarfile ( self , tarfiles , logfiles = [] ) : 
    
        import tarfile, os
        
        tfile = '%s.tgz' % self.name 
        if os.path.exists ( tfile ) :
            self.logger.verbose ( "Remove existing tar-file %s" % tfile ) 
            try :
                os.remove ( tfile )
            except :
                pass

        # create temporary tar-file
        tmptar = CleanUp.tempfile ( prefix = 'ostap-tmp-tarfile-' , suffix = '.tgz' )
            
        with tarfile.open ( tmptar , 'w:gz' ) as tar :
            for x in  tarfiles: tar.add ( x )
            self.logger.debug ( "Tar/gz    file  : %s" % tmptar ) 
            if self.verbose : tar.list ()

        assert os.path.exists ( tmptar ) and os.path.isfile ( tmptar ) and tarfile.is_tarfile ( tmptar ) , \
               'Non-existing or invalid temporary tar-file!'

        ## copy it 
        shutil.copy ( tmptar , tfile )
        
        assert os.path.exists ( tfile ) and os.path.isfile ( tfile ) and tarfile.is_tarfile ( tfile ) , \
               'Non-existing or invalid temporary tar-file!'
            
        ## finally set the tar-file 
        if os.path.exists ( tfile ) and tarfile.is_tarfile( tfile ) :
            self.__tar_file = tfile ; ## os.path.abspath ( tfile ) 
            self.logger.debug ( "Tar/gz    file  : %s" % tfile ) 
        
        if not logfiles : return tfile
        
        lfile = '%s_logs.tgz' % self.name 
        if os.path.exists ( lfile ) :
            self.logger.verbose ( "Remove existing tar-logfile %s" % lfile )
            try :
                os.remove ( lfile )
            except :
                pass
            
        # create temporary tar-file
        tmptar = CleanUp.tempfile ( prefix = 'ostap-tmp-logfile-' , suffix = '.tgz' )
        with tarfile.open ( tmptar , 'w:gz' ) as tar :
            for x in  logfiles:
                if os.path.exists ( x ) : tar.add ( x )
            self.logger.debug ( "Tar/gz logfile  : %s" % tmptar ) 
            if self.verbose : tar.list ()
                
        if os.path.exists ( tmptar ) and os.path.isfile ( tmptar ) and tarfile.is_tarfile ( tmptar ) :
            shutil.copy ( tmptar , lfile )                   

        ## finally set the tar/log-file 
        if os.path.exists ( lfile ) and tarfile.is_tarfile( lfile ) :
            self.__log_file = lfile ## os.path.abspath ( lfile )
            self.logger.debug ( "Tar/gz logfile  : %s" % lfile  ) 
            

        return tfile


# =============================================================================
## @class WeightFiles
#  helper structure  to deal with weights files
class WeightsFiles(CleanUp) :
    """Helper structure  to deal with weights files
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
                tar.extractall ( path = tmpdir )                
                logger.debug ('Un-tar into temporary directory %s' % tmpdir ) 
                weights_files  = [ os.path.join ( tmpdir , i ) for i in tar.getnames() ]
                self.tmpfiles += weights_files
                
        self.__weights_files = weights_files

    @property
    def files   ( self ) :
        "`files': the weights file"
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
    """The `chopping' TMVA reader.
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
    def __init__ ( self                           ,
                   categoryfunc                   ,
                   N                              , 
                   variables                      ,
                   weights_files                  , 
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
        
        variables            = list  ( variables ) ; variables.sort ()        
        self.__variables     = tuple(variables)
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
        """Helper class to get the decision of `chopper'
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
            """Evaluate the chopper from TTree/TChain/RooAbsData:
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
            """Evaluate the method from parameters 
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
            """Get the  mean over all categories
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
        """Helper utility to  get the correspondig function from the  reader:
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
        """Helper utility to  get the correspondig function from the  reader:
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
        """The main method - evaluate of TMVA from the certain category reader 
        
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
        """Evaluate TMVA
        >>> reader = ...
        >>> pt, y  = ...  ##
        >>> print('MLP response is: ', reader.evaluate ( category , 'MLP' , pt , y ))
        """
        assert isinstance ( category , integer_types ) and 0 <= category < self.__N, \
               "Invalid `category' %s/%s" % ( category ,  type ( category ) )
        return self.__readers[ category ].evaluate ( method , *args ) 
                                


# =============================================================================
def _add_response_tree ( tree , verbose , *args ) :
    """Specific action to ROOT.TTree
    """
    
    import ostap.trees.trees
    from   ostap.core.core           import Ostap, ROOTCWD
    from   ostap.io.root_file        import REOPEN
    from   ostap.utils.progress_conf import progress_conf

    tdir     = tree.GetDirectory()
    branches = tree.branches () if verbose else set() 

    with ROOTCWD () , REOPEN ( tdir )  as tfile  : 
        
        tdir.cd()

        ## add progress bar 
        if verbose : sc = Ostap.TMVA.addChoppingResponse ( tree , progress_conf () , *args  )
        else       : sc = Ostap.TMVA.addChoppingResponse ( tree ,                    *args  )
        
        if sc.isFailure() :
            logger.error ( 'Error from Ostap::TMVA::addChoppingResponse %s' % sc )
        
        if tfile.IsWritable() :
            tfile.Write( "" , ROOT.TFile.kOverwrite )
        else :
            logger.error ( "Can't write TTree back to the file" )

            
    return sc , tdir.Get ( tree.GetName() )    ## RETURN

## status = sc
## newt   = tdir.Get ( tree.GetName() )
## if verbose :
##     new_branches = set ( newt.branches () ) | set ( newt.leaves() )
##     new_branches = sorted ( new_branches - set ( branches ) )
##     if new_branches : 
##         n = len ( new_branches )  
##         if 1 == n  : title = 'Added %s branch to TChain'   % n
##         else       : title = 'Added %s branches to TChain' % n
##         table = newt.table ( new_branches , title = title , prefix = '# ' )
##         logger.info ( '%s:\n%s' % ( title , table ) ) 
        
##         return sc , tree                               ## RETURN
    
    


# =============================================================================
def _add_response_chain ( chain , verbose , *args ) :
    """Specific action to ROOT.TChain
    """
    
    import ostap.trees.trees
    
    files    = chain.files    ()
    cname    = chain.GetName  () 
    branches = chain.branches () if verbose else set() 
    
    if not files :
        logger.warning ( 'addChoppingResponse: empty chain (no files)' )
        return Ostap.StatusCode ( 900 ) , chain 

    status  = None 

    tree_verbose  = verbose and      len ( files ) < 5
    chain_verbose = verbose and 5 <= len ( files )
 
    from ostap.utils.progress_bar import progress_bar
    for f in progress_bar ( files , len ( files ) , silent = not chain_verbose ) :

        with  ROOT.TFile.Open ( f , 'UPDATE' ) as ff  :
            ## get the tree 
            tt      = ff.Get(cname)
            ## treat the tree 
            sc , nt = _add_response_tree ( tt , tree_verbose , *args )
            if status is None or sc.isFailure() : status = sc 
            
    newc = ROOT.TChain ( cname )
    for f in  files : newc.Add ( f  )

    if verbose :
        new_branches = set ( newc.branches () ) | set ( newc.leaves() )
        new_branches = sorted ( new_branches - set ( branches ) ) 
        if new_branches : 
            n = len ( new_branches )  
            if 1 == n  : title = 'Added %s branch to TChain'   % n
            else       : title = 'Added %s branches to TChain' % n
            table = newc.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
   
    return status, newc
  

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
                          chopper                     , ## chopping category/formula 
                          N                           , ## number of categrories
                          inputs                      , ## input variables 
                          weights_files               , ## files with TMVA weigths (tar/gz or xml)
                          category_name = 'chopping'  , ## category name 
                          prefix        = 'tmva_'     , ## prefix for TMVA-variable         
                          suffix        = '_response' , ## suffix for TMVA-variable 
                          options       =  ''         , ## TMVA-reader options
                          verbose       = True        , ## verbosity flag 
                          aux           = 0.9         ) :
    """
    Helper function to add TMVA/chopping  response into dataset
    >>> tar_file = trainer.tar_file
    >>> dataset  = ...
    >>> inputs   = [ 'var1' , 'var2' , 'var2' ] ## input variables to TMVA 
    >>> dataset.addChoppingResponse ( dataset , chopper ,  inputs , tar_file , prefix = 'tmva_' )
    """
    assert isinstance ( N , int ) and 1 < N < 10000 , 'Invalid "N" %s' % N

    
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

    args = chopper , category_name , N , _inputs , _maps , options , prefix , suffix , aux
    
    if   isinstance ( dataset , ROOT.TChain  ) :
        
        sc , newdata = _add_response_chain ( dataset , verbose , *args ) 
        if sc.isFailure() : logger.error ( 'Error from Ostap::TMVA::addChoppingResponse %s' % sc )
        return newdata
    
    elif isinstance ( dataset , ROOT.TTree   ) :
        
        sc , newdata = _add_response_tree  ( dataset , verbose , *args )
        if sc.isFailure() : logger.error ( 'Error from Ostap::TMVA::addChoppingResponse %s' % sc )
        
        return newdata 

    ## here we deal with RooAbsData
    
    if isinstance ( chopper , str ) :
        
        if chopper in dataset : chopper = getattr ( dataset , chopper )            
        else :
            
            varset  = dataset.get()
            chopper = make_formula ( 'chopping' , chopper , varset )
            
    assert isinstance ( chopper , ROOT.RooAbsReal ), 'Invalid chopper type %s' % chopper 


    category = ROOT.RooCategory ( category_name ,
                                  'Chopping category: (%s)%%%d' %  ( chopper.GetTitle() , N ) )

    
    if   N <    10 : fmt = category_name + '_%d'    
    if   N <   100 : fmt = category_name + '_%02d'  
    elif N <  1000 : fmt = category_name + '_%03d'  
    elif N < 10000 : fmt = category_name + '_%04d'  
    else           : fmt = category_name + '_%d'   
    ## 
    for i in range ( N ) : category.defineType ( fmt % i  , i )

    args = chopper , category , N , _inputs , _maps , options , prefix , suffix , aux
    
    from   ostap.utils.progress_conf import progress_conf
    ## add progress bar 
    if verbose : sc = Ostap.TMVA.addChoppingResponse ( dataset , progress_conf () , *args )
    else       : sc = Ostap.TMVA.addChoppingResponse ( dataset ,                    *args )

    if sc.isFailure() : logger.error ( 'Error from Ostap::TMVA::addChoppingResponse %s' % sc )

    return dataset 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
