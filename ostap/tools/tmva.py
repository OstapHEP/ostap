#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/tools/tmva.py
#  Python interface to basic TMVA functionality: Trainer and Reader 
#
#  Actually for the Trainer, it is a bit simplified version of Albert's code 
#   - thanks to Albert PUIG
#  Inspired from
#  @see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py
#
#  @date   2013-10-02
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
#  - thanks to Albert PUIG
#
# =============================================================================
""" Python interface to two major TMVA classes
- Trainer
- Reader

Actually for the Trainer, it is a bit simplified version of Albert's code [thanks Albert Puig],
inspired from
http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py

"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2013-10-02"
__version__ = '$Revision$'
__all__     = (
    "Trainer"         , ## the basic TMVA trainer 
    "Reader"          , ## the basic TMVA reader
    "addTMVAResponse" , ## add TMVA response to RooDataSet/TTree/TChain
    "tmvaGUI"
    )
# =============================================================================
from   ostap.core.core         import items_loop, WSE, Ostap, rootWarning 
from   ostap.core.ostap_types  import num_types, string_types, integer_types 
from   ostap.core.meta_info    import root_version_int, root_info  
from   ostap.utils.cleanup     import CleanUp
import ostap.io.root_file
import ROOT, os, glob, math, tarfile, shutil, itertools 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger     import getLogger, attention
# =============================================================================
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.tmva' )
else                       : logger = getLogger ( __name__           )
# =============================================================================
pattern_XML   = "%s/weights/%s*.weights.xml"
pattern_CLASS = "%s/weights/%s*.class.C" 
pattern_PLOTS = "%s/plots/*.*"
# =============================================================================
trivial_opts = (
    'V', '!V' , 'H' , '!H' , 'S' , '!S' ,
    'Verbose' , '!Verbose' ,
    'Silent'  , '!Silent'  , ''
    )
# =============================================================================
good_for_negative = (
    ROOT.TMVA.Types.kLikelihood ,
    ROOT.TMVA.Types.kPDERS      ,
    ROOT.TMVA.Types.kPDEFoam    ,
    ROOT.TMVA.Types.kKNN        ,
    ROOT.TMVA.Types.kSVM        ,
    ROOT.TMVA.Types.kBDT
    )
# =============================================================================
def dir_name ( name ) :
    name = str( name )
    for s in ' %!><\n?(){}[]+:.,;-^&|$#@="\'/' :
        while s in name : name = name.replace ( ' ' , '_' )
    while '__' in name : name = name.replace ('__','_')
    return  name 
# =============================================================================
## @class WeightFiles
#  helper structure to deal with weights files
class WeightsFiles(CleanUp) :
    """Helper structure  to deal with weights files
    """
    def __init__ ( self , weights_files ) :

        self.__trash = None 
        ## string ? treat it a a single tar-file 
        if isinstance ( weights_files , str  ) :
            
            assert os.path.exists  ( weights_files ) , \
                   "Non-existing `weights_file'  %s"   %  weights_files 
            
            if tarfile.is_tarfile ( weights_files ) :
                def xml_files ( archive ) :
                    for tarinfo in archive:
                        if os.path.splitext(tarinfo.name)[1] == ".xml":
                            yield tarinfo
                            
                with tarfile.open ( weights_files , 'r' ) as tar :
                    ## tar.list() 
                    xmls   = [ f for f in xml_files ( tar ) ] 
                    tmpdir = self.tempdir ( prefix = 'ostap-tmva-weights-' )
                    tar.extractall ( path = tmpdir , members = xml_files ( tar ) )
                    logger.debug ('Un-tar into temporary directory %s' % tmpdir ) 
                    weights_files  = [ os.path.join ( tmpdir, x.name ) for x  in xmls ]
                    self.trash.add ( tmpdir ) 
            else :
                weights_files = [ weights_files ]

        ## list/tuple/etc 
        if not isinstance ( weights_files , dict ) : 

            wfs = {} 
            for wf in weights_files :
                
                if isinstance ( wf , str ) : method , xml = None, wf 
                else                       : method , xml =       wf
                
                assert os.path.exists ( xml ) and os.path.isfile ( xml ), \
                       "No weights file '%s'" %  xml 
                
                if not method :
                    f      = os.path.split ( xml ) [-1]
                    p,s,m  = f.rpartition('_')
                    method = m[:m.find('.weights.xml')]
                    
                if not method :
                    if not s : raise AttributeError("Can't extract method name from %s" % xml )
                    
                wfs[method] = xml
                
            weights_files   = wfs

        ## dictionary
        assert isinstance ( weights_files , dict  ), \
               "Invalid type of 'weight_files'  %s "    % weights_files

        for method , xml in items_loop ( weights_files ) :            
            assert os.path.exists ( xml ) and os.path.isfile ( xml ), \
                   "No weights file '%s'" %  xml 
            
        self.__methods       = tuple ( [ i for  i in  weights_files.keys() ] )
        import copy
        self.__weights_files = copy.deepcopy ( weights_files ) 

    @property
    def methods ( self ) :
        "'methods': the known methods from weights-file"
        return self.__methods
    @property
    def files   ( self ) :
        "'files': the weights file"
        import copy
        return copy.deepcopy ( self.__weights_files )

# =============================================================================
## some manipulations with TMVA options 
def opts_replace ( opts , expr , direct = True ) :
    """some manipulations with TMVA options"""
    if  direct : 
        if   0 <= opts.find ( '!' + expr ) : opts  = opts.replace ( '!' + expr , expr )
        elif 0 <= opts.find (       expr ) : pass
        else                               : opts += ':'  + expr
    else :
        if   0 <= opts.find ( '!' + expr ) : pass 
        elif 0 <= opts.find (       expr ) : opts  = opts.replace ( expr , '!' + expr )
        else                               : opts += ':!' + expr
        #
    while 0<= opts.find ( '::' ) : opts = opts.replace ( '::' , ':' )
    return opts 

# =============================================================================
## Create the tar file from components, optionally create it as tmp,
#  and later copy to the final destination.
#  (Sometime for unreliable file systems (like EOS via fsmount)
#  normal creation of tar-gile raises OSError
#  @code
#  files =...
#  t = make_tarfile ( 'outptu.tgz' , files , varbose = True , tmp = True )  
#  @endcode
def make_tarfile ( output , files , verbose = False , tmp = False ) :
    """Create the tar file from components, optionally create it as tmp,
    and later copy to the final destination.
    (Sometime for unreliable file systems (like EOS via fsmount)
    normal creation of tar-gile raises OSError
    >>> files =...
    >>> t = make_tarfile ( 'outptu.tgz' , files , varbose = True , tmp = True )  
    """
    
    assert not os.path.exists ( output ) or os.path.isfile ( outptu ) , \
           "Invalid destination for output tar-file: `%s'" % output
    
    if os.path.exists ( output ) and os.path.isfile ( output ) :
        try :
            os.remove ( output )
        except :
            pass

    if tmp :
        
         # 1) create & fill the temporary tar-file
         tmptar = CleanUp.tempfile ( prefix = 'ostap-tmp-tarfile-' , suffix = '.tgz' )         
         with tarfile.open ( tmptar , 'w:gz' ) as tar :
             for x in files :
                 if os.path.exists ( x ) and os.path.isfile ( x ) : tar.add ( x )
             if verbose : tar.list()
         # 2) check it 
         assert os.path.exists ( tmptar ) and os.path.isfile ( tmptar ) and tarfile.is_tarfile ( tmptar ) , \
                'Non-existing or invalid temporary tar-file!'
         # 3) move to the final destination 
         shutil.move ( tmptar , output )
         
    else :

        ## create & fill the tar-file
        with tarfile.open ( output, 'w:gz' ) as tar :
            for x in files : 
                 if os.path.exists ( x ) and os.path.isfile ( x ) : tar.add ( x )
            if verbose : tar.list()

    ## check the result 
    assert os.path.exists ( output ) and os.path.isfile ( output ) and tarfile.is_tarfile ( output ) , \
           "Non-existing or invalid tar-file:`%s'" % output 

    return output 
    
# =============================================================================
## @class Trainer
#  Helper class to train TMVA
# 
#  - Book the trainer: 
#  @code
#  >>> from ostap.tools.tmva import Trainer
#  >>> t = Trainer(
#  ...   ## TMVA methods  
#  ...   methods =  [ ( ROOT.TMVA.Types.kMLP ,
#  ...                'MLP',
#  ...                'H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator' ) ] , 
#  ...   variables = [ 'dtfchi2' , 'ctau', 'ptb' , 'vchi2' ] , ## list  of variables 
#  ...   signal          = treeSignal            ,  ## TTree for 'signal' sample  
#  ...   background      = treeBackgrund         ,  ## TTree for 'background'  sample 
#  ...   signal_cuts     = cuts_for_signal       ,
#  ...   background_cuts = cuts_for_background   )
#  @endcode 
#
#  - Use the trainer: 
#  @code 
#  >>> t.train()     ## Use the trainer
#  @endcode 
#
#  - Get the results from the trainer:
#  @code
#  >>> xml      = t.weights_files ## get the weights XML   files 
#  >>> classes  = t.class_files   ## get the classes (C++) files
#  >>> output   = t.output_file   ## get the output ROOT file  
#  >>> tar_file = t.tar_file      ## get the tar-file with XML&C++
#  @endcode
#
#  For more detailes
#  @see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py
# 
#  @date   2013-10-02
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
#  Thanks to Albert PUIG
class Trainer(object):
    """Helper class to train TMVA:  
    
    >>> from ostap.tools.tmva import Trainer
    >>> t = Trainer(
    ...   ## TMVA methods  
    ...   methods =  [ ( ROOT.TMVA.Types.kMLP ,
    ...                'MLP',
    ...                'H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator' ) ] , 
    ...   variables = [ 'dtfchi2' , 'ctau', 'ptb' , 'vchi2' ] , ## list  of variables 
    ...   signal          = treeSignal            ,  ## TTree for 'signal' sample  
    ...   background      = treeBackgrund         ,  ## TTree for 'background'  sample 
    ...   signal_cuts     = cuts_for_signal       ,
    ...   background_cuts = cuts_for_background   )

    Use the trainer: 
    >>> t.train()     ## Use the trainer

    Get the results:
    >>> xml      = t.weights_files ## get the weights XML   files 
    >>> classes  = t.class_files   ## get the classes (C++) files
    >>> output   = t.output_file   ## get the output ROOT file  
    >>> tar_file = t.tar_file      ## get the tar-file with XML&C++
    
    Actually it is a bit simplified version of the original code by Albert PUIG,
    inspired from
    http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py
    """
    # =========================================================================
    ## constructor
    #
    #  - Book the trainer: 
    #  @code
    #  >>> from ostap.tools.tmva import Trainer
    #  >>> t = Trainer(
    #  ...   ## TMVA methods  
    #  ...   methods =  [ ( ROOT.TMVA.Types.kMLP ,
    #  ...                'MLP',
    #  ...                'H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator' ) ] , 
    #  ...   variables = [ 'dtfchi2' , 'ctau', 'ptb' , 'vchi2' ] , ## list  of variables 
    #  ...   signal          = treeSignal            ,  ## TTree for 'signal' sample  
    #  ...   background      = treeBackgrund         ,  ## TTree for 'background'  sample 
    #  ...   signal_cuts     = cuts_for_signal       ,
    #  ...   background_cuts = cuts_for_background   )
    #  @endcode 
    #  For more detailes
    #  @see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
    def __init__(  self                           ,
                   methods                        ,
                   variables                      ,  ## list of variables
                   signal                         ,  ## signal sample/tree
                   background                     ,  ## background sample/tree
                   ## 
                   signal_vars          = {}      ,  ## dictionary with new variables for signal sample 
                   background_vars      = {}      ,  ## dictionary with new variables for background sample 
                   ## 
                   signal_cuts          = ''      ,  # signal cuts 
                   background_cuts      = ''      ,  # background cuts                   
                   spectators           = []      ,
                   bookingoptions       = "Transformations=I;D;P;G,D" , 
                   configuration        = "SplitMode=Random:NormMode=NumEvents" ,
                   ##
                   signal_weight        = None    ,                
                   background_weight    = None    ,
                   ##
                   prefilter            = ''      ,  ## prefilter cuts before TMVA data loader
                   prefilter_signal     = ''      ,  ## separate prefilter for signal data 
                   prefilter_background = ''      ,  ## separate prefilter for background data 
                   ##
                   signal_train_fraction     = -1 , ## fraction of signal events used for training     : 0<=f<1 
                   background_train_fraction = -1 , ## fraction of background events used for training : 0<=f<1
                   ##
                   prescale_signal      = 1       , ## prescale factor for the signal 
                   prescale_background  = 1       , ## prescale factor for the background
                   ## 
                   output_file          = ''      , ## the name of output file
                   verbose              = True    ,
                   logging              = True    ,
                   name                 = 'TMVA'  ,
                   make_plots           = True    ,
                   workdir              = ''      , 
                   category             = -1      ,
                   multithread          = False   ,
                   logger               = None    ) :
        
        """Constructor with list of methods
        
        >>> from ostap.tools.tmva import Trainer
        >>> t = Trainer(
        ...   ## TMVA methods  
        ...   methods =  [ ( ROOT.TMVA.Types.kMLP ,
        ...                'MLP',
        ...                'H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator' ) ] , 
        ...   variables = [ 'dtfchi2' , 'ctau', 'ptb' , 'vchi2' ] , ## list  of variables 
        ...   signal          = treeSignal            ,  ## TTree for 'signal' sample  
        ...   background      = treeBackgrund         ,  ## TTree for 'background'  sample 
        ...   signal_cuts     = cuts_for_signal       ,
        ...   background_cuts = cuts_for_background   )        
        - For more detailes
        see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
        """

        self.__name   = name 
        self.__logger = logger if logger else getLogger ( self.name )

        dirname         = dir_name        ( self.name )
        self.__dirname  = dirname

        self.__logging           = False
        if  logging :
            if isinstance ( logging , string_types ) : self.__logging = logging 
            else                                     : self.__logging = "%s.log" % self.name


        
        for i , e in enumerate ( methods ) :
            t = e [ 0 ] 
            assert ROOT.TMVA.Types.kVariable <= t < ROOT.TMVA.Types.kMaxMethod, \
                   'Invalid TMVA.Types.EMVA %s' % t 
                  
        self.__methods           = tuple ( methods    )
        
        variables                = list  ( variables  ) 

        _sig_vars = {}
        _bkg_vars = {}
        
        vars = []
        for v in variables  :
            a , s1 , b  = v.partition ( ':' )
            a = a.strip ()
            b = b.strip ()
            if a and s1 and b :
                c , s2 , d = b.partition ( '?' )
                c = c.strip()
                d = d.strip()
                if c and s2 and d :
                    vars.append ( a )
                    _sig_vars.update ( { a : c } ) ## ATTENTION HERE 
                    _bkg_vars.update ( { a : d } ) ## ATTENTION HERE 
                else :
                    vars.append ( a )
                    _sig_vars.update ( { a : b } )  ## ATTENTION HERE 
                    _bkg_vars.update ( { a : b } )  ## ATTENTION HERE 
            else :
                vars.append ( v )

        if signal_vars     : _sig_vars.update ( signal_vars     )
        if background_vars : _bkg_vars.update ( background_vars )

        signal_vars      = _sig_vars 
        background_vars  = _bkg_vars 

        variables = vars 
        variables.sort ()
        self.__variables         = tuple ( variables )

        self.__configuration     = configuration
        
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
        
        from ostap.trees.trees import Chain
        
        if   isinstance ( signal     , Chain           ) : pass 
        elif isinstance ( background , ROOT.TTree      ) : signal = Chain ( signal ) 
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
        elif isinstance ( background , ROOT.TTree      ) : background = Chain ( background ) 
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
            

        self.__signal               = signal
        self.__signal_cuts          = signal_cuts  
        self.__signal_weight        = signal_weight

        self.__background           = background
        self.__background_cuts      = background_cuts 
        self.__background_weight    = background_weight

        self.__signal_vars          = {}
        if signal_vars     : self.__signal_vars.update     ( signal_vars )        
        self.__background_vars      = {}
        if background_vars : self.__background_vars.update ( background_vars ) 
        

        self.__prefilter            = ROOT.TCut ( prefilter ) 
        self.__prefilter_signal     = prefilter_signal     if prefilter_signal     else '' 
        self.__prefilter_background = prefilter_background if prefilter_background else '' 
        
        self.__category             = int ( category )
        self.__make_plots           = True if make_plots else False
        
        self.__spectators           = [] 
        
        ## self.__verbose = True if verbose else False 
        self.__verbose = True if ( verbose and  self.category <= 0 ) else False 


        self.__bookingoptions   = bookingoptions
                
        import os 
        if not workdir : workdir = os.getcwd()
        if not os.path.exists ( workdir ) :
            from ostap.utils.basis import make_dirs
            make_dirs ( workdir )
        assert workdir and os.path.exists ( workdir ) and os.path.isdir ( workdir ), \
               'No valid working directory!'
        
        self.__workdir = os.path.abspath ( workdir ) 
        if self.verbose : self.logger.debug ("Working directory is %s" % self.workdir )
   
        ##
        ## outputs :
        ## 
        self.__weights_files = []
        self.__class_files   = []
        self.__output_file   = output_file if output_file else os.path.join ( self.workdir , '%s.root' % self.name )
        self.__tar_file      = None 
        self.__log_file      = '' 
        self.__plots         = []
        
        self.__output_file   = os.path.abspath ( self.output_file ) if self.output_file else None  

        #
        ## minor adjustment
        #
        
        opts = self.__bookingoptions
        
        if 0 > opts.find ( "AnalysisType=" ) :
            opts += ":AnalysisType=Classification"
            self.logger.debug('Booking options are appended with ":AnalysisType=Classification"' )

        opts = opts_replace ( opts , 'V:'               ,     self.verbose )
        opts = opts_replace ( opts , 'Silent:'          , not self.verbose )
        
        from ostap.utils.basic import isatty
        OK1  = self.verbose and self.category in ( 0 , -1 )
        OK2  = OK1 and isatty () 
        
        opts = opts_replace ( opts , 'DrawProgressBar:' , OK1 ) 
        opts = opts_replace ( opts , 'Color:'           , OK2 ) 

        _methods  = [] 
        for _m in self.methods :
            _m    = list( _m ) 
            _m[2] = opts_replace ( _m[2] , 'H:' , OK1 )
            _m[2] = opts_replace ( _m[2] , 'V:' , OK1 )
            _methods.append (  tuple ( _m ) )
        self.__methods = tuple (  _methods ) 

        self.__bookingoptions = str ( opts )

        
        pattern_xml = pattern_XML   % ( self.dirname ,  self.dirname )
        pattern_C   = pattern_CLASS % ( self.dirname ,  self.dirname )

        rf = []
        import glob
        for f in glob.glob ( pattern_xml ) :
            rf.append ( f ) 
            os.remove ( f ) 
        for f in glob.glob ( pattern_C   ) :
            rf.append ( f ) 
            os.remove ( f )
        if rf : self.logger.debug ( "Remove existing xml/class-files %s" % rf ) 
        
        self.__pattern_xml   = pattern_xml 
        self.__pattern_C     = pattern_C 
        self.__pattern_plots = pattern_PLOTS % self.name

        self.__multithread = multithread and 61800 <= root_version_int 

        if self.verbose :
            
            rows = [ ( 'Item' , 'Value' ) ]
            
            row  = 'Name'      , self.name
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
            
            for i , o in enumerate ( [ o for o in self.bookingoptions.split ( ':' ) if not o in trivial_opts ] ) :
                if 0 == i : row = 'Booking options' , o
                else      : row = ''                , o 
                rows.append ( row )
                
            for i , o in enumerate ( [ o for o in self.configuration.split ( ':' ) if not o in trivial_opts ] ) :
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
                row = 'Background weight'    , str ( self.background_weight ) 
                rows.append ( row )

            if 0 < self.background_train_fraction < 1 :
                row = 'Background train fraction' , '%.1f%%' % ( 100 *  self.background_train_fraction ) 
                rows.append ( row )

            if 1 != self.prescale_background:
                row = 'Background prescale'      , '%s' % self.prescale_background
                rows.append ( row )

            ms  = [ m[1] for m in self.methods ]
            row = 'Methods' , ' '.join ( ms ) 
            rows.append ( row )
            
            from ostap.logger.colorized import allright 
            for m in self.methods :
                row = 'Method' , '%s Id:%s' % ( allright ( m[1] ) , m[0] )
                rows.append ( row )
                for i , o in enumerate ( [ o for o in m[2].split(':' ) if not o in trivial_opts ] ) :
                    if 0 == i : row = 'Method configuration' , o
                    else      : row =   ''                    , o 
                    rows.append ( row ) 
                
            import ostap.logger.table as T
            title = "TMVA Trainer %s created" % self.name 
            table = T.table (  rows , title = title , prefix = "# " , alignment = "lw" )
            self.logger.info ( "%s\n%s" % ( title , table ) ) 

    @property
    def name    ( self ) :
        """'name'    : the name of TMVA trainer"""
        return self.__name
    
    @property
    def methods ( self ) :
        """'methods' : the list of TMVA methods to be used"""
        return tuple(self.__methods)

    @property
    def method_names ( self ) :
        """'method_names' : tuple of method names"""
        return tuple ( m[1] for m in self.__methods ) 
        
    @property
    def logger ( self ) :
        """'logger' : the logger instace for this Trainer"""
        return self.__logger
    
    @property
    def variables ( self ) :
        """'variables' : the list of variables  to be used for training"""
        return tuple(self.__variables)

    @property
    def spectators ( self ) :
        """'spectators' : the list of spectators to be used"""
        return tuple(self.__spectators)

    @property
    def signal ( self ) :
        """'signal' :  TTree for signal events"""
        return self.__signal.chain 
    
    @property
    def signal_cuts ( self ) :
        """'signal_cuts' :  cuts to be applied for 'signal' sample"""
        return str(self.__signal_cuts)

    @property
    def signal_weight ( self ) :
        """'signal_weight' : weight to be applied for 'signal' sample"""
        return self.__signal_weight
    
    @property
    def background ( self ) :
        """'background' :  TTree for background events"""
        return self.__background.chain
    
    @property
    def background_cuts ( self ) :
        """'background_cuts' :  cuts to be applied for 'backgroud' sample """
        return str(self.__background_cuts)

    @property
    def background_weight ( self ) :
        """'background_weight' : weight to be applied for 'background' sample"""
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
        """'prefilter' : cuts to be applied/prefilter before processing"""
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
        """'bookingoptions' : options used to book TMVA::Factory"""
        return str(self.__bookingoptions)

    @property
    def configuration ( self ) :
        """'configuration' : options used to book TMVA"""
        return str(self.__configuration)

    @property
    def signal_train_fraction ( self ) :
        """'signal_train_fraction': if non-negative, fraction of signal events used for training"""
        return self.__signal_train_fraction

    @property
    def background_train_fraction ( self ) :
        """'background_train_fraction': if non-negative, fraction of background events used for training"""
        return self.__background_train_fraction

    @property
    def prescale_signal ( self ) :
        """'prescale_signal': prescale the signal sample"""
        return self.__prescale_signal

    @property
    def prescale_background ( self ) :
        """'prescale_background': prescale the background sample"""
        return self.__prescale_background

    @property
    def verbose ( self ) :
        """'verbose' : verbosity  flag"""
        return self.__verbose

    @property
    def silent ( self ) :
        """'silent' : verbosity flag"""
        return not self.verbose
    
    @property
    def logging ( self ) :
        """'logging' : logging flag : produce log-file?"""
        return self.__logging

    @property
    def make_plots ( self ) :
        """'make_plots' : make standard TMVA plots?"""
        return self.__make_plots
    
    @property
    def dirname  ( self ) :
        """'dirname'  : the output directiory name"""
        return str(self.__dirname) 

    @property
    def weights_files ( self ) :
        """'weights_files' : the list/tuple of final files with TMVA weights"""
        return tuple(self.__weights_files)
    @property
    def class_files ( self ) :
        """'class_files' : the list/tuple of final files with TMVA classes"""
        return tuple(self.__class_files)
    @property
    def output_file ( self ) :
        """'output_file'  : the output file (if any)"""
        return str(self.__output_file) if  self.__output_file else None
    
    @property
    def tar_file ( self ) :
        """'tar_file'  : the compressed (gz) tar file with results"""
        return str(self.__tar_file) if self.__tar_file else None 
    @property
    def log_file ( self ) :
        """'log_file'  : the name of log-file """
        return str(self.__log_file) if self.__log_file else '' 

    @property
    def category    ( self ) :
        """'category'  : chopping category"""
        return self.__category

    @property 
    def multithread ( self ) :
        """'multithread' : make try to use multithreading in TMVA"""
        return self.__multithread
    
    @property
    def workdir ( self ) :
        """'workdir' : working directory"""
        return self.__workdir

    @property
    def plots ( self ) :
        """'plots': list of produced plots"""
        return tuple ( self.__plots ) 

    @property
    def show_plots ( self ) :
        """'show_plots': show plots?"""
        return self.verbose and ( self.category in ( 0 , -1 ) )
        
    # =========================================================================
    ## train TMVA 
    #  @code
    #  trainer.train ()
    #  @endcode  
    #  @return the name of output XML file with the weights 
    def train ( self )  :
        """ train TMVA 
        >>> trainer.train ()
        return the name of output XML files with the weights 
        """
        from ostap.utils.utils import keepCWD
        with keepCWD ( self.workdir ) :
            return self._train() 

    # =========================================================================
    ## train TMVA 
    #  @code
    #  trainer.train ()
    #  @endcode  
    #  @return the name of output XML file with the weights 
    def _train ( self )  :
        """ train TMVA 
        >>> trainer.train ()
        return the name of output XML files with the weights 
        """

        log = self.logging 

        self.__log_file = None

        if  log :
            
            try :
                if os.path.exists ( log ) and os.path.isfile ( log ) : os.remove ( log )
            except :
                pass

            from ostap.logger.utils import TeeCpp , OutputC  
            context  = TeeCpp ( log ) if self.verbose and self.category in ( 0 , -1 ) else OutputC ( log , True , True ) 

        else    :

            from ostap.logger.utils  import MuteC  , NoContext
            context  = NoContext () if self.verbose and self.category in ( 0 , -1 ) else MuteC     ()
            
        from ostap.logger.utils import NoContext
        from ostap.utils.utils  import ImplicitMT
        
        context2 = ImplicitMT ( True ) if self.multithread else NoContext () 

        with context , context2 :
            
            result = self.__train ()

            ## Training outputs
                    
            rows = [ ( 'Item' , 'Value' ) ]

            vv = set ( self.variables )
            for k in self.signal_vars     : vv.add ( k )
            for k in self.background_vars : vv.add ( k )
            vv = sorted ( vv ) 
            
            row  = 'Variables', ', '.join( vv )
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
            
            if self.workdir and os.path.exists ( self.workdir ) and os.path.isdir ( self.workdir ) :
                row = "Workdir" , self.workdir 
                rows.append ( row )
                
            if self.weights_files : 
                from ostap.utils.basic import commonpath
                wd = commonpath ( self.weights_files ) if 1 < len ( self.weights_files ) else ''
                if wd :
                    row = "Weights directory" , wd
                    rows.append ( row )
                lw = len ( wd )
                for w in self.weights_files :
                    row = "Weights file" , w if not lw else w[lw+1:]
                    rows.append ( row )

            if self.class_files : 
                from ostap.utils.basic import commonpath 
                wd = commonpath ( self.class_files ) if 1 < len ( self.class_files ) else ''
                if wd :
                    row = "Classes directory" , wd
                    rows.append ( row )
                lw = len ( wd )
                for w in self.class_files :
                    row = "Class file" , w if not lw else w[lw+1:]
                    rows.append ( row )
                    
            if self.plots :
                from ostap.utils.basic import commonpath 
                wd = commonpath ( self.plots ) if 1 < len ( self.plots ) else ''
                if wd :
                    row = "Plots directory" , wd
                    rows.append ( row )
                lw = len ( wd )
                for w in self.plots :
                    row = "Plot file" , w if not lw else w[lw+1:]
                    rows.append ( row )
                
            if self.output_file and os.path.exists ( self.output_file  ) : 
                row = "Output ROOT file" , self.output_file
                rows.append ( row )
                
            if self.log_file and os.path.exists ( self.log_file  ) : 
                row = "Log file" , self.log_file 
                rows.append ( row )
                
            if self.tar_file and os.path.exists ( self.tar_file  ) and tarfile.is_tarfile ( self.tar_file ) : 
                row = "Tar file" , self.tar_file 
                rows.append ( row )
                
            import ostap.logger.table as T
            title = "TMVA %s outputs" % self.name 
            table = T.table (  rows , title = title , prefix = "# " , alignment = "lw" )
            self.logger.info ( "%s\n%s" % ( title , table ) ) 
            
            if self.verbose                          and \
                   self.tar_file                     and \
                   os.path.exists ( self.tar_file  ) and \
                   tarfile.is_tarfile ( self.tar_file ) :
                self.logger.info ( 'Content of the tar file %s' % self.tar_file )
                with tarfile.open ( self.tar_file , 'r' ) as tar : tar.list ()
                
            if self.verbose and self.output_file       and \
                   os.path.exists ( self.output_file ) and \
                   os.path.isfile ( self.output_file ) :
                import ostap.io.root_file 
                with ROOT.TFile.open ( self.output_file , 'r' ) as rf :
                    table = rf.as_table()
                    self.logger.info ( "Content of the output ROOT file %s\n%s" % ( self.output_file , table ) )

                   
        if log and os.path.exists ( log ) and os.path.isfile ( log ) :
            try :
                shutil.move ( log , self.dirname )
                log = os.path.join ( self.dirname , os.path.basename ( log ) )
            except :
                pass
            
        if log and os.path.exists ( log ) and os.path.isfile ( log ) :
            self.__log_file = os.path.abspath ( log )
        else : 
            self.__log_file = None 

        return result
                
    # =========================================================================
    ## train TMVA 
    #  @code
    #  trainer.train ()
    #  @endcode  
    #  @return the names of output XML files with the weights 
    def __train ( self ) :
        """Train the TMVA:
        - returns the names of output XML file with the weights 
        >>> trainer.train () 
        """

        rf = []
        import os, glob
        for f in glob.glob ( self.__pattern_xml ) :
            rf.append ( f ) 
            os.remove ( f ) 
        for f in glob.glob ( self.__pattern_C   ) :
            rf.append ( f )
            os.remove ( f )
        if rf : self.logger.debug ( "Remove existing xml/class-files %s" % rf  ) 

        ## ROOT 6/18 crashes here...
        ## if root_version_int < 61800 : 
        ##    ROOT.TMVA.Tools.Instance()
        with ROOT.TFile.Open ( self.output_file, 'RECREATE' )  as outFile :

            self.logger.debug ( 'Output ROOT file: %s ' %  outFile.GetName() )

            # =================================================================
            ## the final adjustment
            # =================================================================
            
            opts = self.__bookingoptions
            from ostap.utils.basic import isatty
            OK1  = self.verbose and self.category in ( 0 , -1 )
            OK2  = OK1 and isatty ()             
            opts = opts_replace ( opts , 'DrawProgressBar:' , OK1 ) 
            opts = opts_replace ( opts , 'Color:'           , OK2 ) 
            self.__bookingoptions = opts 

            # =================================================================
            #
            # =================================================================
            
            all_vars = []
            
            for v in self.variables :
                
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

            ## if self.prefilter         : all_vars.append ( self.prefilter         )                
            ## if self.signal_cuts       : all_vars.append ( self.signal_cuts       )
            ## if self.signal_weight     : all_vars.append ( self.signal_weight     )
            ## if self.background_cuts   : all_vars.append ( self.background_cuts   )
            ## if self.background_weight : all_vars.append ( self.background_weight )

            # =================================================================
            ## prefilter/prescale signal if required 
            if self.prefilter_signal or self.prefilter or 1 != self.prescale_signal or self.signal_vars :
                
                if self.signal_weight     : all_vars.append ( self.signal_weight     )
                
                import ostap.trees.trees
                avars = self.signal.the_variables ( all_vars )
                
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
                                           silent    = silent )
                
                self.__signal      = self.__SigTR
                self.__signal_cuts = ROOT.TCut() 
                
            missvars = [ v for v in self.signal_vars if not v in self.signal ]
            assert not missvars , "Variables %s are not in signal sample!" % missvars                          

            # =================================================================
            ## prefilter/prescale background if required 
            if self.prefilter_background or self.prefilter or 1 != self.prescale_background or self.background_vars :
                
                if self.background_weight : all_vars.append ( self.background_weight )
                
                import ostap.trees.trees
                avars = self.background.the_variables ( all_vars )
                
                bcuts = {}
                if self.prefilter_background : bcuts.update ( { 'PreFilter/Background' : self.prefilter_background } )                  
                if self.prefilter            : bcuts.update ( { 'PreFilter/Common'     : self.prefilter            } )
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
                                           save_vars = avars  ,
                                           new_vars  = self.background_vars     , 
                                           prescale  = self.prescale_background ,  
                                           silent    = silent )
                
                self.__background      = self.__BkgTR
                self.__background_cuts = ROOT.TCut() 

            missvars = [ v for v in self.background_vars if not v in self.background ]
            assert not missvars , "Variables %s are not in background sample!" % missvars                         
                         
            # =====================================================================
            ## check for signal weights
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
                            self.logger.error ( "Method '%s' does not support negative (signal) weights" % m[1] )
                            
            # =================================================================
            ## check for background weights
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
                            self.logger.error ( "Method '%s' does not support negative (background) weights" % m[1] )
                            
                            
            NS = -1
            NB = -1
            SW = None
            BW = None
            
            if 0<= self.signal_train_fraction <1 or 0<= self.background_train_fraction < 1 or self.verbose :
                
                sc = ROOT.TCut ( self.    signal_cuts )
                bc = ROOT.TCut ( self.background_cuts )
                if self.    signal_weight : sc *= self.    signal_weight
                if self.background_weight : sc *= self.background_weight

                if ( 6 , 20 ) <= root_info :
                    from ostap.frames.frames import frame_statVar 
                    ss = frame_statVar ( self.signal     , '1' , sc )
                    sb = frame_statVar ( self.background , '1' , bc )
                else :
                    ss = self.signal    .statVar ( '1' , sc )
                    sb = self.background.statVar ( '1' , bc )
                    
                if isinstance ( ss , WSE ) : ss = ss.values()
                if isinstance ( sb , WSE ) : sb = sb.values()
                
                NS = ss.nEntries ()
                SW = ss.sum      ()            
                NB = sb.nEntries ()
                BW = sb.sum      ()
                
                if 0 < self.signal_train_fraction < 1 :
                    nt = math.ceil ( NS * self.signal_train_fraction )
                    bo = self.configuration.split ( ':' )
                    bo = [ b for b in bo if not b.startswith('nTrain_Signal') ]
                    bo = [ b for b in bo if not b.startswith('nTest_Signal' ) ]
                    nt =  'nTrain_Signal=%s' % nt
                    self.__configuration = ':'.join ( [ nt ] + bo ) 
                    self.logger.info ( "Extend TMVA configuration for '%s'" % nt ) 
                    
                if 0 < self.background_train_fraction < 1 :
                    nt = math.ceil ( NB * self.background_train_fraction )
                    bo = self.configuration.split ( ':' )
                    bo = [ b for b in bo if not b.startswith('nTrain_Background') ]
                    bo = [ b for b in bo if not b.startswith('nTest_Background' ) ]
                    nt =  'nTrain_Background=%s' % nt
                    self.__configuration = ':'.join ( [ nt ] + bo ) 
                    self.logger.info ( "Extend TMVA configuration for '%s'" % nt ) 

            # =================================================================
            # The table
            # =================================================================            

            ## if self.verbose :
            if 1 < 2 : 
                
                rows = [ ( 'Item' , 'Value' ) ]
                
                row  = 'Name'      , self.name
                rows.append ( row )

                vv  = set ( self.variables )
                for k in self.signal_vars     : vv.add ( k )
                for k in self.background_vars : vv.add ( k )
                vv = sorted ( vv )
                
                row = 'Variables'      , ' '.join ( vv )
                rows.append ( row )
                
                if self.spectators :
                    row = 'Spectators' , ' '.join ( self.spectators )
                    rows.append ( row )
                    
                for i , o in enumerate ( [ o for o in self.bookingoptions.split ( ':' ) if not o in trivial_opts ] ) :
                    if 0 == i : row = 'Booking options' , o
                    else      : row = ''                , o 
                    rows.append ( row )
                    
                for i , o in enumerate ( [ o for o in self.configuration.split ( ':' ) if not o in trivial_opts ] ) :
                    if 0 == i : row = 'TMVA configuraton' , o
                    else      : row = ''                  , o 
                    rows.append ( row )
                    
                if self.signal_cuts : 
                    row = 'Signal cuts' , str ( self.signal_cuts ) 
                    rows.append ( row )

                if self.signal_weight : 
                    row = 'Signal weight' , str ( self.signal_weight ) 
                    rows.append ( row )
                    row = 'Signal total'     , '%s' % NS
                    rows.append ( row ) 
                    row = 'Signal weighted'  , '%s' % SW
                    rows.append ( row ) 
                else :
                    row = 'Signal total'     , '%s' % NS
                    rows.append ( row ) 

                if 0 < self.signal_train_fraction < 1 :
                    row = 'Signal train fraction' , '%.1f%%' % ( 100 *  self.signal_train_fraction ) 
                    rows.append ( row )
                    
                if self.background_cuts : 
                    row = 'Background cuts' , str ( self.background_cuts ) 
                    rows.append ( row )
                    
                if self.background_weight : 
                    row = 'Background weight'    , str ( self.background_weight ) 
                    rows.append ( row )
                    row = 'Background total'     , '%s' % NB
                    rows.append ( row ) 
                    row = 'Bacgground weighted'  , '%s' % BW
                    rows.append ( row ) 
                else :
                    row = 'Background total'     , '%s' % NB
                    rows.append ( row ) 
                    
                if 0 < self.background_train_fraction < 1 :
                    row = 'Backgroundtrain fraction' , '%.1f%%' % ( 100 *  self.background_train_fraction ) 
                    rows.append ( row )
                    
                ## for m in self.methods :
                ##    row = 'Method' , '%s #%s' % ( m[1] , m[0] )
                ##    rows.append ( row )
                ##    for i , o in enumerate ( m[2].split(':' ) ) :
                ##        if 0 == i : row = 'Method configruration' , o
                ##        else      : row =   ''                    , o 
                ##        rows.append ( row ) 
                        
                import ostap.logger.table as T
                title = "TMVA Trainer %s start " % self.name 
                table = T.table (  rows , title = title , prefix = "# " , alignment = "lw" )
                self.logger.info ( "%s\n%s" % ( title , table ) ) 
                
                
                
            bo = self.bookingoptions.split (':')
            bo.sort() 
            if self.verbose : self.logger.info  ( 'Book TMVA-factory %s ' % bo ) 
            else            : self.logger.debug ( 'Book TMVA-factory %s ' % bo ) 

            factory = ROOT.TMVA.Factory (                
                self.name             ,
                outFile               ,
                self.bookingoptions   )

            factory.SetVerbose( self.verbose )
        
            ## 
            dataloader = ROOT.TMVA.DataLoader ( self.name )

            
            avars = set ( self.variables )
            for v in self.signal_vars     : avars.add ( v ) 
            for v in self.background_vars : avars.add ( v )
            avars = sorted ( avars )
            
            all_vars = [] 
            ## for v in self.variables :
            for v in avars :
                vv = v
                if isinstance ( vv , str ) : vv = ( vv , 'F' )
                all_vars.append ( vv[0] ) 
                dataloader.AddVariable  ( *vv )    
            self.logger.info ( "Variables           : %s" % str ( self.variables ) ) 

            for v in self.spectators :
                vv = v
                if isinstance ( vv , str ) : vv = ( vv , 'F' )             
                all_vars.append ( vv[0] ) 
                dataloader.AddSpectator ( *vv )
                #

            if self.verbose : self.logger.info ( "Loading 'Signal'     sample" ) 
            dataloader.AddTree ( self.signal     , 'Signal'     , 1.0 , ROOT.TCut ( self.    signal_cuts ) )
            
            if self.verbose : self.logger.info ( "Loading 'Background' sample" )             
            dataloader.AddTree ( self.background , 'Background' , 1.0 , ROOT.TCut ( self.background_cuts ) )
            #
            if self.signal_weight :
                dataloader.SetSignalWeightExpression     ( self.signal_weight     )
                self.logger.info ( "Signal     weight: '%s'" % ( attention ( self.signal_weight     ) ) )
            if self.background_weight :
                dataloader.SetBackgroundWeightExpression ( self.background_weight )
                self.logger.info ( "Background weight: '%s'" % ( attention ( self.background_weight ) ) )
                
            self.logger.info     ( "Configuration    : '%s'" % str ( self.configuration ) )
            dataloader.PrepareTrainingAndTestTree(
                ROOT.TCut ( self.signal_cuts     ) ,
                ROOT.TCut ( self.background_cuts ) ,
                self.configuration            )
            #
            for m in self.methods :
                bo = m[2].split(':')
                bo.sort() 
                if  self.verbose : self.logger.info  ( "Book %11s/%d method %s" % ( m[1] , m[0] , bo ) )
                else             : self.logger.debug ( "Book %11s/%d method %s" % ( m[1] , m[0] , bo ) )
                
                factory.BookMethod ( dataloader , *m )

            # Train MVAs
            ms = self.method_names 
            self.logger.info  ( "Train    all methods %s " % str ( ms ) )
            factory.TrainAllMethods    ()
            ## Test MVAs
            self.logger.info  ( "Test     all methods %s " % str ( ms ) )
            factory.TestAllMethods     ()
            # Evaluate MVAs
            self.logger.info  ( "Evaluate all methods %s " % str ( ms ) )
            factory.EvaluateAllMethods ()

        ## ROC curves 
        if ( self.make_plots or self.verbose ) and ( 6 , 24 ) <= root_info :
            import ostap.plotting.canvas
            from ostap.utils.utils import batch , keepCanvas
            with batch ( ROOT.ROOT.GetROOT().IsBatch() or not self.show_plots ) , keepCanvas() : 
            ## cnv = factory.GetROCCurve ( dataloader )
                cnv = factory.GetROCCurve ( self.name )
                if cnv :
                    cnv.Draw()
                    cnv >> ( "%s/plots/ROC" % self.dirname )
                    del cnv
                    
        ## AUC for ROC curves
        if self.verbose : 
            rows = [ ('Method' , 'AUC' ) ]
            for m in self.methods :
                mname = m[1]
                if m[0] == ROOT.TMVA.Types.kCuts and root_info < ( 6 , 24 ) :
                    self.logger.warning ( 'Skip ROC/AUC for %s' % m[1] ) 
                    continue 
                ## auc = factory.GetROCIntegral ( dataloader , mname )
                auc = factory.GetROCIntegral ( self.name , mname )
                row = mname , '%.6g' % auc
                rows.append ( row ) 
            import ostap.logger.table as T
            title = "ROC/AUC compare"
            table = T.table ( rows , prefix = "# " , title = title , alignment = "ll" )
            self.logger.info ( "%s:\n%s" % ( title , table  ) )
            

        # check the output.
        if os.path.exists ( self.output_file ) :

            if self.output_file == '%s.root' % self.name and os.path.exists ( self.dirname ) and os.path.isdir ( self.dirname ) : 
                import shutil
                try :
                    shutil.move ( self.output_file , self.dirname )
                    noof = os.path.join ( self.dirname , self.output_file )
                    if os.path.exists ( noof ) : self.__output_file = noof
                except :
                    pass  
            try : 
                with ROOT.TFile.Open ( self.output_file, 'READ' ) as outFile : 
                    self.logger.debug ( "Output ROOT file is %s" % outFile.GetName() ) 
                    if self.verbose : outFile.ls()
            except :
                pass


        del dataloader
        del factory 


        if  self.make_plots :
            self.makePlots ()                
        
        import glob, os 
        self.__weights_files = tuple ( [ f for f in glob.glob ( self.__pattern_xml   ) ] )
        self.__class_files   = tuple ( [ f for f in glob.glob ( self.__pattern_C     ) ] ) 
        self.__plots         = tuple ( [ f for f in glob.glob ( self.__pattern_plots ) ] )

        
        tfile = make_tarfile ( output  = '.'.join ( [ self.name , 'tgz'] ) ,
                               files   = self.weights_files + self.class_files + self.plots + ( self.log_file , ) ,
                               verbose = self.verbose  ,
                               tmp     = True          ) 
        
        ## if os.path.exists ( tfile ) :
        ##     self.logger.debug  ( "Remove existing tar-file %s" % tfile ) 
        ##     try :
        ##         os.remove ( tfile )
        ##     except :
        ##         pass 
            
        ## # create temporary tar-file
        
        ## tmptar = CleanUp.tempfile ( prefix = 'ostap-tmp-tarfile-' , suffix = '.tgz' )                
        ## with tarfile.open ( tmptar , 'w:gz' ) as tar :
        ##     for x in self.weights_files : tar.add ( x )
        ##     for x in self.  class_files : tar.add ( x )
        ##     for x in self.  plots       : tar.add ( x )
        ##     if self.log_file and os.path.exists ( self.log_file ) and os.path.isfile ( self.log_file ) :
        ##         tar.add ( self.log_file )
                
        ## assert os.path.exists ( tmptar ) and os.path.isfile ( tmptar ) and tarfile.is_tarfile ( tmptar ) , \
        ##        'Non-existing or invalid temporary tar-file!'

        ## ## copy it
        ## import shutil
        ## shutil.move ( tmptar , tfile )

        ## assert os.path.exists ( tfile ) and os.path.isfile ( tfile ) and tarfile.is_tarfile ( tfile ) , \
        ##        'Non-existing or invalid tar-file!'
          
        self.__weights_files = tuple ( [ os.path.abspath ( f ) for f in self.weights_files ] ) 
        self.__class_files   = tuple ( [ os.path.abspath ( f ) for f in self.class_files   ] ) 
        self.__plots         = tuple ( [ os.path.abspath ( f ) for f in self.__plots       ] ) 
        
        ## finally set tar-file 
        if os.path.exists ( tfile ) and tarfile.is_tarfile( tfile ) :            
            if os.path.exists ( self.dirname ) and os.path.isdir ( self.dirname ) :
                try :
                    shutil.move ( tfile , self.dirname )
                    ntf = os.path.join ( self.dirname , tfile )
                    if os.path.exists ( ntf ) and tarfile.is_tarfile( ntf ) : tfile = os.path.abspath ( ntf )
                except :
                    pass 
            self.__tar_file = os.path.abspath ( tfile ) 
                        
        return self.tar_file 

    # =========================================================================
    ## make selected standard TMVA plots 
    def makePlots ( self , name = None , output = None , ) :
        """Make selected standard TMVA plots"""

        self.logger.warning ( "makePlots: method is disbaled!" )
        return 

        name   = name   if name   else self.name
        output = output if output else self.output_file
        
        if not output :
            self.logger.warning ('No output file is specified!')
            return 
        if not os.path.exists ( output ) or not os.path.isfile ( output ) :
            self.logger.error   ('No output file %s is found !' % output )
            return
        
        try :
            import ostap.io.root_file
            with ROOT.TFile.Open ( output , 'READ' , exception = True ) as o :
                pass   
        except IOError :
            self.logger.error ("Output file %s can't be opened!"   % output )
            return
          
        #
        ## make the plots in TMVA  style
        #
        
        plots = [
            ##
            ( ROOT.TMVA.variables      ,  ( name , output     ) ) ,
            ( ROOT.TMVA.correlations   ,  ( name , output     ) ) ,
            ##
            ( ROOT.TMVA.mvas           ,  ( name , output , 0 ) ) ,
            ( ROOT.TMVA.mvas           ,  ( name , output , 1 ) ) ,
            ( ROOT.TMVA.mvas           ,  ( name , output , 2 ) ) ,
            ( ROOT.TMVA.mvas           ,  ( name , output , 3 ) ) ,
            ##
            ( ROOT.TMVA.efficiencies   ,  ( name , output , 0 ) ) ,
            ( ROOT.TMVA.efficiencies   ,  ( name , output , 1 ) ) ,
            ( ROOT.TMVA.efficiencies   ,  ( name , output , 2 ) ) ,
            ( ROOT.TMVA.efficiencies   ,  ( name , output , 3 ) ) ,
            ##
            ## ( ROOT.TMVA.paracoor       ,  ( name , output     ) ) ,
            ## 
            ]

        ## it crashes if runnig on-line 
        ## if  ( 6 , 24 ) <= root_info :
        ## plots.append ( ( ROOT.TMVA.mvaeffs        ,  ( name , output ) ) ) 

        ## if hasattr ( ROOT.TMVA , 'network'                ) :
        ##     plots.append ( ( ROOT.TMVA.network            , ( name , output ) ) ) 
        ## if hasattr ( ROOT.TMVA , 'nannconvergencetest'    ) :
        ##     plots.append ( ( ROOT.TMVA.annconvergencetest , ( name , output ) ) )

        ## if [ m for m in self.methods if m[0] == ROOT.TMVA.Types.kLikelihood ] :
        ##    plots.append ( ( ROOT.TMVA.likelihoodrefs     , ( name , output ) ) )

        ## if [ m for m in self.methods if m[0] == ROOT.TMVA.Types.kBDT ] :            
        ##     if hasattr ( ROOT.TMVA , 'BDT'                    ) :
        ##         plots.append ( ( ROOT.TMVA.BDT                , ( name , output ) ) )                
        ##     if hasattr ( ROOT.TMVA , 'BDTControlPlots'        ) :
        ##         plots.append ( ( ROOT.TMVA.BDTControlPlots    , ( name , output ) ) )
            
        ## if [ m for m in self.methods if m[0] == ROOT.TMVA.Types.kBoost ] :                
        ##     if hasattr ( ROOT.TMVA , 'BoostControlPlots'      ) :
        ##         plots.append ( ( ROOT.TMVA.BoostControlPlots  , ( name , output ) ) ) 

        ## change to some temporary directory
                
        from ostap.utils.utils import batch, keepCanvas    
        for fun, args  in plots :            
            with batch ( ROOT.ROOT.GetROOT().IsBatch () or not self.show_plots ) , keepCanvas() , rootWarning ()  :            
                logger.info ( 'makePlots: Execute macro ROOT.TMVA.%s%s' % ( fun.__name__ , str ( args ) ) )
                fun ( *args )

# =============================================================================
## make selected standard TMVA plots 
def make_Plots ( name , output , show_plots = True ) :
    """Make selected standard TMVA plots"""

    ## if (6,29) <= root_info :
    ##    logger.warning ("function is disabled")
    ##    return 
    
    if not output :
        logger.warning ('No output file is specified!')
        return 
    if not os.path.exists ( output ) or not os.path.isfile ( output ) :
        logger.error   ('No output file %s is found !' % output )
        return

    try :
        import ostap.io.root_file
        with ROOT.TFile.Open ( output , 'READ' , exception = True ) as o :
            pass   
    except IOError :
        logger.error ("Output file %s can't be opened!"   % output )
        return
        
    output = os.path.abspath ( output )
        
    ## make the plots in TMVA  style
    #
    logger.info ('make_Plots: Making the standard TMVA plots') 
    from ostap.utils.utils import batch , cmd_exists, keepCanvas

    plots = [
        ##
        ( ROOT.TMVA.variables      ,  ( name , output     ) ) ,
        ( ROOT.TMVA.correlations   ,  ( name , output     ) ) ,
        ##
        ( ROOT.TMVA.mvas           ,  ( name , output , 0 ) ) ,
        ( ROOT.TMVA.mvas           ,  ( name , output , 1 ) ) ,
        ( ROOT.TMVA.mvas           ,  ( name , output , 2 ) ) ,
        ( ROOT.TMVA.mvas           ,  ( name , output , 3 ) ) ,
        ##
        ( ROOT.TMVA.efficiencies   ,  ( name , output , 0 ) ) ,
        ( ROOT.TMVA.efficiencies   ,  ( name , output , 1 ) ) ,
        ( ROOT.TMVA.efficiencies   ,  ( name , output , 2 ) ) ,
        ( ROOT.TMVA.efficiencies   ,  ( name , output , 3 ) ) ,
        ##
        ( ROOT.TMVA.paracoor       ,  ( name , output     ) ) ,
        ## 
        ## ( ROOT.TMVA.likelihoodrefs ,  ( name , output     ) ) ,
        ]
    
    ## if (6,24) <= root_info :
    ##    plots.append   ( ( ROOT.TMVA.mvaeffs            , ( name , output ) ) )
    ## else :
    logger.warning ( 'make_Plots: Skip    macro ROOT.TMVA.%s%s' % ( 'mvaeffs'  , str ( ( name , output ) ) ) ) 
            
        
    if hasattr ( ROOT.TMVA , 'network'                ) :
        plots.append ( ( ROOT.TMVA.network            , ( name , output ) ) ) 
    if hasattr ( ROOT.TMVA , 'nannconvergencetest'    ) :
        plots.append ( ( ROOT.TMVA.annconvergencetest , ( name , output ) ) )

    if hasattr ( ROOT.TMVA , 'BDT'                    ) :
        plots.append ( ( ROOT.TMVA.BDT                , ( name , output ) ) )
        
    if hasattr ( ROOT.TMVA , 'BDTControlPlots'        ) :
        plots.append ( ( ROOT.TMVA.BDTControlPlots    , ( name , output ) ) )
        
    if hasattr ( ROOT.TMVA , 'BoostControlPlots'      ) :
        plots.append ( ( ROOT.TMVA.BoostControlPlots  , ( name , output ) ) ) 
        
    workdir = CleanUp.tempdir ( prefix = 'ostap-tmva-%s-plots' % name  )
    from   ostap.utils.utils import keepCWD
    from   ostap.utils.basic import make_dirs
    import glob 
    with keepCWD ( workdir ) :

        make_dirs ( '%s/plots' % name )
        
        logger.info ( "Use temporary working directory:'%s'" % os.getcwd() )
        
        for fun, args  in plots :
            
            with batch ( ROOT.ROOT.GetROOT().IsBatch () or not show_plots ) , keepCanvas() , rootWarning ()  :            
                logger.info ( 'make_Plots: Execute macro ROOT.TMVA.%s%s' % ( fun.__name__ , str ( args ) ) )
                fun ( *args )

        plots = tuple ( [ f for f in glob.glob ( pattern_PLOTS % name ) ] )

        if plots :
            
            ## tarfile with plots 
            tfile = make_tarfile ( output  = '.'.join ( [ '%s_plots' % name , 'tgz' ] ) ,
                                   files   = plots        ,
                                   tmp     = True         ) 
            logger.info ( "Tarfile with plots: '%s'" % tfile )
            return tfile 
                
        return '' 
            
# =============================================================================
## @class Reader
#  Rather generic python interface to TMVA-reader
#  @code 
#  from ostap.tools.tmva import Reader
#  >>> reader = Reader ( 'MyTMVA' ,
#  ... variables = [ ## name      accessor  
#  ...              ( 'pt'   , lambda s : s.pt ) ,
#  ...              ( 'ip'   , lambda s : s.ip ) ,
#  ...                'var1'                     ,   ## s.var1 will be used 
#  ...                'var2'                     ] , ## s.var2 will be used 
#  ... weights_files = 'my_weights.xml' )
#  @endcode
#  <code>weights_fiels</code> can be :
#  - single xml-file with weights for the given method
#  - single tar/tgz/tar.gz-file with weights files
#  - list of xml-files with weights
#  If the xml-filenames follow TMVA Trainer convention, the training method will be
#  extracted from the file name, otherwise it needs to be specified as dictionary   
#  with the key being the method name: 
#  @code
#  weights_files = { 'MPL' : 'my_weights.xml', 'BDTG' : 'weights2.xml' } 
#  @endcode
#  A list of 2-element tuples is also supported :
#  @code
#  weights_files = [ ('MPL','my_weights.xml') , ('BDTG','weights2.xml') ] 
#  @endcode
#
#  Use the reader
#  @code
#  >>> tree =  ....  ## TTree/TChain/RooDataSet with data
#  >>> for entry in tree :
#  ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
#  ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
#  ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   )
#  @endcode
#
#  A bit more efficient form is :
#  @code
#  >>> tree =  ....  ## TTree/TChain/RooDataSet with data
#  >>> mlp_fun  =  reader.MLP
#  >>> bdgt_fun =  reader.BDTG
#  >>> for entry in tree :
#  ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
#  ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
#  ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg) )
#  @endcode
#
#  - It it naturally merges with Ostap's SelectorWithVars utility
#
#  @see TMVA::Reader
#  @date   2013-10-02
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
#  Thanks to Alexander BARANOV
class Reader(object)  :
    """ Rather generic python interface to TMVA-reader

    >>> from ostap.tools.tmva import Reader
    >>> r = Reader ( 'MyTMVA' ,
    ... variables = [ ## name      accessor  
    ...              ( 'pt'   , lambda s : s.pt ) ,
    ...              ( 'ip'   , lambda s : s.ip ) ,
    ...                'var1'                     ,   ## s.var1 will be used 
    ...                'var2'                     ] , ## s.var2 will be used 
    ... weights_files = 'my_weights.xml' )

    - attention: It is *not* CPU-efficient:
    Ugly tricks with arrays are required to bypass some technical limitations

    'weights_files' can be :
    - single xml-file with weights for the given method
    - single tar/tgz/tar.gz-file with weights files  (output of 'Trainer.tar_file')
    - list of xml-files with weights                 (output of 'Trainer.weights_files')

    If the xml-filenames follow TMVA Trainer convention, the training method will be
    extracted from the file name, otherwise it needs to be specified as dictionary
    with the key being the method name, or list of 2-element tuples :
    
    >>> weights_files = {  'MPL':'my_weights.xml'  ,  'BDTG':'weights2.xml'  } 
    >>> weights_files = [ ('MPL','my_weights.xml') , ('BDTG','weights2.xml') ]
    
    - Use the reader
    >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    >>> for entry in tree :
    ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
    ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
    ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg) )

    - A bit more efficient form is :
    >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    >>> mlp_fun  =  reader.MLP
    >>> bdgt_fun =  reader.BDTG
    >>> for entry in tree :
    ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
    ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
    ...     print ( 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)) 
    
    It it naturally merges with Ostap's SelectorWithVars utility 
    """    
    ## constructor
    #  @code 
    #  from ostap.tools.tmva import Reader
    #  >>> r = Reader ( 'MyTMVA' ,
    #  ... variables = [ ## name      accessor  
    #  ...              ( 'pt'   , lambda s : s.pt ) ,
    #  ...              ( 'ip'   , lambda s : s.ip ) ,
    #  ...                'var1'                     ,   ## s.var1 will be used 
    #  ...                'var2'                     ] , ## s.var2 will be used 
    #  ... weights_files = 'my_weights.xml' )
    #  @endcode
    #  <code>weights_fiels</code> can be :
    #  - single xml-file with weights for the given method
    #  - single tar/tgz/tar.gz-file with weights files
    #  - list of xml-files with weights
    #  If the xml-filenames follow TMVA Trainer convention, the training method will be
    #  extracted from the file name, otherwise it needs to be specified as dictionary   
    #  with the key being the method name: 
    #  @code
    #  weights_files = { 'MPL' : 'my_weights.xml', 'BDTG' : 'weights2.xml' } 
    #  @endcode
    #  A list of 2-element tuples is also supported :
    #  @code
    #  weights_files = [ ('MPL','my_weights.xml') , ('BDTG','weights2.xml') ] 
    #  @endcode
    def __init__ ( self                 ,
                   name                 , 
                   variables            ,
                   weights_files        ,
                   options      = ''    ,
                   logger       = None  , 
                   verbose      = True  ) :
        """Construct the TMVA reader         
        >>> from ostap.tools.tmva import Reader
        >>> r = Reader ( 'MyTMVA' ,
        ... variables = [ ## name      accessor  
        ...              ( 'pt'   , lambda s : s.pt ) ,
        ...              ( 'ip'   , lambda s : s.ip ) ,
        ...                'var1'                     ,   ## s.var1 will be used 
        ...                'var2'                     ] , ## s.var2 will be used 
        ... weights_files = 'my_weights.xml' )
        
        - attention: It is *not* CPU-efficient:
        Ugly tricks with arrays are required to bypass some technical limitations
        
        weights_fiels can be :
        - single xml-file with weights for the given method
        - single tar/tgz/tar.gz-file with weights files
        - list of xml-files with weights
        
        If the xml-filenames follow TMVA Trainer convention, the training method will be
        extracted from the file name, otherwise it needs to be specified as dictionary
        with the key being the method name, or list of 2-element tuples :
        
        >>> weights_files = {  'MPL':'my_weights.xml'  ,  'BDTG':'weights2.xml'  } 
        >>> weights_files = [ ('MPL','my_weights.xml') , ('BDTG','weights2.xml') ]
        
        """            

        self.__name   = name
        self.__logger = logger if logger else getLogger ( self.name )

        if root_version_int < 61800 : 
            ROOT.TMVA.Tools.Instance()
        
        verbose = True if verbose else False
        ##
        options = opts_replace ( options , 'V:'      ,     verbose )
        options = opts_replace ( options , 'Silent:' , not verbose )
        ##
        from ostap.utils.basic import isatty
        #
        options = opts_replace ( options , 'Color:'           , verbose and isatty () ) 

        self.__reader = ROOT.TMVA.Reader( options , verbose )

        ## treat the weigths files
        self.__weights = WeightsFiles ( weights_files )

        ##  book the variables:
        #   dirty trick with arrays is needed due to a bit strange reader interface.
        #   [TMVA reader needs the address of 'float' (in C++ sense) variable]
        from array import array

        self.__variables = []
        
        variables = list  ( variables ) ; variables.sort ()
        variables = tuple ( variables )
        
        for v in variables : 

            if   isinstance ( v , str ) :
                
                vname  = v
                vfun   = lambda s : getattr ( s , vname )
                vfield = array ( 'f' , [1] )                  ## NB: note the type 
                
            elif isinstance ( v , tuple ) and 2 == len ( v ) :

                vname  = v[0]
                vfun   = v[1]
                vfield = array ( 'f' , [1] )                  ## NB: note the type here 

            else :
                
                self.logger.error ('Invalid variable description!'     )
                raise AttributeError ( 'Invalid variable description!' )

            ##                     name    accessor   address 
            self.__variables += [ ( vname , vfun     , vfield  ) ] 


        self.__variables = tuple ( self.__variables )
        
        ## declare all variables to TMVA.Reader 
        for v in self.__variables :
            self.__reader.AddVariable ( v[0] , v[2] )            
        
        self.__methods = self.weights.methods

        for method , xml in  items_loop ( self.weights.files ) :
            m = self.__reader.BookMVA ( method , xml )
            assert  m , 'Error in booking %s/%s' % (  method  , xml )
            self.logger.debug ('TMVA Reader is booked for method:%s xml: %s' % ( method , xml ) )
            
        self.logger.info ('TMVA Reader  booked methods are %s' % str ( self.__methods ) )
        self.__methods = tuple ( self.__methods )

    @property
    def name ( self ) :
        """'name' - the name of the reader"""
        return self.__name

    @property
    def logger ( self ) :
        """'logger' : the logger instace for this Trainer"""
        return self.__logger
    
    @property
    def reader ( self ) :
        """'reader' - the  actual TMVA.Reader object"""
        return self.__reader

    @property
    def weights ( self ) :
        """'weigths' : weights-files """
        return self.__weights 
        
    @property
    def methods ( self ) :
        """'methods' - the  list/tuple of booked TMVA methods"""
        return tuple (self.__methods)

    @property
    def variables ( self ) :
        """'variables' - helper structure to access TMVA variables
        >>> variables = [ ## name      accessor  
        ...              ( 'pt'   , lambda s : s.pt ) ,
        ...              ( 'ip'   , lambda s : s.ip ) ,
        ...                'var1'                     ,   ## s.var1 will be used 
        ...                'var2'                     ] , ## s.var2 will be used 
        """        
        return self.__variables 

    # =========================================================================
    ## helper class to get TMVA decision for certain method 
    #  @code
    #  reader = ...
    #  var = reader[ method ]
    #  val = var ( entry )
    #  @endcode 
    class Var (object) :
        """Helper class to get TMVA decision for the certain method
        >>>  reader = ...
        >>>  var = reader[ method ]
        >>>  val = var ( entry )
        """
        def __init__ ( self , reader , method ) :
            self.__reader = reader
            self.__method = str ( method )
            self.__nvars  = len ( reader.variables )
        # =====================================================================
        @property
        def nvars ( self ) :
            """'nvars' : number of TMVA variables"""
            return self.__nvars
        @property
        def reader ( self ) :
            """'reader' : TMVA reader """
            return self.__reader
        @property
        def method ( self ) :
            """'method' : TMVA method name"""
            return self.__method
        # =====================================================================
        ## the main method 
        def __call__ ( self , entry , cut_efficiency = 0.9  ) :
            return self.eval ( entry , cut_efficiency )
        # =====================================================================
        ## Evaluate the method from TTree/TChain/RooAbsData using the accessors, defined  early
        # @code 
        # tree   = ...
        # method = ...
        # print('Response is %s' % method.eval ( tree ) )
        # @endcode 
        def eval ( self , entry , cut_efficiency = 0.9 ) :
            """Evaluate the method fomr TTree/RooAbsData using the
            accessors, defined  early
            >>> tree   = ...
            >>> method = ...
            >>> print('Response is %s'    % method.eval ( tree ) )
            """
            return self.__reader( self.__method , entry , cut_efficiency )        

    # =========================================================================
    ## helper class to get TMVA decision for certain method 
    class Method (Var) :
        """Helper class to get TMVA decision for certain method
        >>>  reader = ...
        >>>  var = reader[ method ]
        >>>  val = var ( entry )
        """
        # =====================================================================
        ## Evaluate the reader
        #  @code
        #  tree   = ...
        #  method = reader.MLP
        #  print('Response %s' % method ( tree ))
        #  @endcode
        #  or using parameters
        #  @code
        #  pt , y, phi = ...
        #  print('Response %s' % method ( pt , y , phi ))
        #  @endcode 
        def __call__ ( self , arg ,  *args ) :
            if   isinstance ( arg , num_types ) : return self.evaluate ( arg , *args ) 
            elif isinstance ( arg , bool      ) : return self.evaluate ( arg , *args ) 
            return self.eval ( arg , *args )
        # =====================================================================
        ## Evaluate the method from parameters 
        # @code 
        # method       = ...
        # pt, eta, phi = 5 ,  3.0 , 0  ## variables 
        # print('Response is %s'    % method.evaluate ( pt ,  eta , phi ) )
        # @endcode 
        def evaluate ( self , *args ) :
            """Evaluate the method from parameters 
            >>> method       = ...
            >>> pt, eta, phi = 5 ,  3.0 , 0  ## variables 
            >>> print('Response is %s'    % method.evaluate ( pt ,  eta , phi ) )
            """
            return self.reader.evaluate ( self.method , *args )        

    # ========================================================================
    ## helper utility to  get the corresponding function from the  reader:
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
        """Helper utility to  get the corresponding function from the  reader:
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
    
    # ========================================================================
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
        ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg) )
        """        
        if not method in self.__methods :
            return AttributeError( 'No method %s is booked!' %  method )
        return Reader.Var  ( self , method ) 

    # =========================================================================
    ## evaluate TMVA
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
    #  (actually TMVA reader needs the address of 'float' (in C++ sense) variable
    def __call__ ( self , method , entry , cut_efficiency = 0.90 ) :
        """Evaluate TMVA
        - Use the reader
        >>> tree =  ....  ## TTree/TChain/RooDataSet with data
        >>> for entry in tree :
        ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
        ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
        ...     print('MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   )
        - It is not CPU efficient :-( 
        - Ugly trick with arrays is needed due to some pure technical problem
        [actually TMVA reader needs the address of 'float' (in C++ sense) variable]
        """
        
        ## loop over all variables 
        for v in self.__variables :
            vfun    = v[1]           ## accessor function 
            v[2][0] = vfun ( entry ) ## fill variable from the tree/chain 
            
        ## evaluate TMVA 
        return self.__reader.EvaluateMVA ( method , cut_efficiency ) 

    # ========================================================================
    ## evaluate TMVA
    #  @code
    #  reader = ...
    #  pt, y  = ...  ##
    #  print ('MLP response is: ', reader.MLP.evaluate ( pt , y ))
    #  @endcode
    def evaluate ( self , method , *args ) :
        """Evaluate TMVA
        >>> reader = ...
        >>> pt, y  = ...  ##
        >>> print('MLP response is: ', reader.MLP.evaluate ( pt , y ))
        """
        l1 = len ( args             )
        l2 = len ( self.__variables )
        
        assert l1 == l2 or l1 == l2 + 1, \
               "TMVA.Reader.evaluate: Invalid length of the argument"

        cut_efficiency = 0.9 
        if l1 == l2 + 1 : cut_efficiency = float ( args[-1] ) 

        ## vector of doubles 
        from ostap.math.base import vDoubles
        
        ## vector of doubles 
        vd = vDoubles ( l2  )  
        for i in  range ( l2 ) :
            vd[i] = float (  args[i] )
            
        ## evaluate TMVA 
        return self.__reader.EvaluateMVA ( vd , method , cut_efficiency ) 

        
_canvas = []
# =============================================================================
## start TMVA gui 
def tmvaGUI ( filename , new_canvas = True ) :
    """Start TMVA-GUI
    """
    ## ROOT.gROOT.LoadMacro('TMVAGui.C')
    if new_canvas :
        from ostap.plotting.canvas import getCanvas
        _c = getCanvas ('glTMVA' , 'TMVA' )
        if not _c in _canvas : _canvas.append ( _c )
    #
    ## start GUI
    return ROOT.TMVA.TMVAGui( filename )

# =============================================================================
## convert input structure to Ostap.TMVA.MAPS
def _inputs2map_ ( inputs ) :
    """Convert input structure to Ostap.TMVA.MAPS
    """
    from ostap.core.core import cpp, std, Ostap 
    MAP     = std.map ( 'std::string', 'std::string' )
    
    _inputs = MAP()
    assert isinstance ( inputs , ( dict , tuple , list ) ) , \
           'Invalid type of "inputs": %s' % inputs
    
    if   isinstance ( inputs , dict  ) :
        for k , v  in items_loop ( inputs ) : _inputs[k] = v
    elif isinstance ( inputs , ( tuple , list ) ) :
        for i in inputs :
            if isinstance ( i , str ) :
                ##
                a , s , b = i.partition ( ':' ) 
                a = a.strip()
                b = b.strip()
                if a and s and b :
                    k , v = a , b
                else : 
                    k , v = i , i
            else                      : k , v = i
            _inputs[k] = v 

    ## 
    assert not _inputs .empty() and _inputs.size() == len ( inputs ), \
           'Invalid MAP size %s for %s' % ( _inputs.size() , inputs ) 

    return _inputs

# =============================================================================
## convert weights structure to Ostap.TMVA.MAP
def _weights2map_ ( weights ) :
    
    from ostap.core.core import cpp, std, Ostap
    MAP = std.map   ( 'std::string', 'std::string' )
    
    if not isinstance ( weights , WeightsFiles ) : 
        weights  = WeightsFiles ( weights  )

    _map = MAP() 
    for method , xml in items_loop ( weights.files ) :
        _map [ method ] = xml

    assert not _map .empty() , \
           'Invalid MAP size %s for' % ( _map.size() , weights )
    
    assert not _map.empty() , "Invalid weights_files: %s"  % weights.files
    return _map , weights  

# =============================================================================
def _add_response_tree  ( tree , verbose , *args ) :
    """Specific action to ROOT.TChain
    """
            
    import ostap.trees.trees
    from   ostap.core.core           import Ostap, ROOTCWD
    from   ostap.io.root_file        import REOPEN
    from   ostap.utils.progress_conf import progress_conf

    vars = set() 
    if verbose :
        vars = set ( tree.branches() ) | set ( tree.leaves () ) 
                
    tdir  = tree.GetDirectory ()
    tname = tree.full_path
    with ROOTCWD () , REOPEN ( tdir ) as tfile : 

        fname = tfile.GetName() 
        logger.debug ( "Procesing %s file" % fname )
        
        if not tfile.IsWritable() :
            logger.error  ( "Can't write TTree back to the file %s" % fname )
            return Ostap.StatusCode ( 1000 ), tree 
        
        tdir.cd()
        
        ## add the progress bar 
        if verbose : sc = Ostap.TMVA.addResponse ( tree , progress_conf () , *args )
        else       : sc = Ostap.TMVA.addResponse ( tree ,                    *args )
        
        if sc.isFailure() : logger.error ( 'Error from Ostap::TMVA::addResponse %s' % sc )
        
        ## tfile.Write() ##  "" ) ## , ROOT.TFile.kOverwrite )
        tfile.Write( "" , ROOT.TFile.kOverwrite )
        
    tree = ROOT.TChain ( tname )
    tree.Add ( fname ) 
    
    if verbose :
        vars  = sorted ( ( set ( tree.branches() ) | set ( tree.leaves () ) ) - vars ) 
        title = "Added TMVA responses [%d entries in '%s'" % ( len ( tree ) , fname  ) 
        table = tree.table ( vars ,  prefix = '# ' , title = title )
        logger.info ( '%s\n%s' % ( title , table ) ) 
        
    return sc , tree 

# =============================================================================
def _add_response_chain ( chain , verbose , *args ) :
    """Specific action to ROOT.TChain
    """
    
    import ostap.trees.trees
    
    files  = chain.files()
    cname  = chain.GetName() 
    
    if not files :
        logger.warning ( 'addTMVAResponse: empty chain (no files)' )
        return Ostap.StatusCode ( 900 ) , chain 

    status = None 

    tree_verbose  = verbose and      len ( files ) < 5
    chain_verbose = verbose and 5 <= len ( files )

    vars = set() 
    if verbose :
        vars = set ( chain.branches() ) | set ( chain.leaves () ) 
        
    from ostap.utils.progress_bar import progress_bar
    nfiles = 0 
    for f in progress_bar ( files , len ( files ) , silent = not chain_verbose  ) :
        
        with ROOT.TFile.Open ( f , 'UPDATE' ,  exception = True ) as rfile :
            ## get the tree
            tt =   rfile.Get ( cname )
            ## treat the tree
            sc , nt = _add_response_tree ( tt  , tree_verbose , *args )
            if status is None or sc.isFailure() : status = sc
            nfiles += 1
            
    newc = ROOT.TChain ( cname )
    for f in  files : newc.Add ( f  )

    if verbose :
        vars  = sorted ( ( set ( newc.branches() ) | set ( newc.leaves () ) ) - vars )
        title = 'Added TMVA responses [%d entries in %d files]' % ( len ( newc ) , nfiles )
        table = newc.table ( vars ,  prefix = '# ' , title = title )
        logger.info ( '%s\n%s' % ( title , table ) )
                      
    return status, newc
        
# =============================================================================
## Helper function to add TMVA response into dataset
#  @code
#  tar_file = trainer.tar_file
#  dataset  = ...
#  inputs = [ 'var1' , 'var2' , 'var2' ]
#  dataset.addTMVAResponse (  inputs , tar_file , prefix = 'tmva_' )
#  @endcode
#  @param dataset  input dataset to be updated 
#  @param inputs   input variables
#  @param weights_files files with TMVA weigths (tar/gz or xml)
#  @param prefix   prefix for TMVA-variable
#  @param suffix   suffix for TMVA-variable
#  @param options  options to be used in TMVA Reader
#  @param verbose  verbose operation?
#  @param aux       obligatory for the cuts method, where it represents the efficiency cutoff
def addTMVAResponse ( dataset                ,   ## input dataset to be updated
                      inputs                 ,   ## input variables 
                      weights_files          ,   ## files with TMVA weigths (tar/gz or xml)
                      prefix   = 'tmva_'     ,   ## prefix for TMVA-variable 
                      suffix   = '_response' ,   ## suffix for TMVA-variable
                      options  = ''          ,   ## TMVA-reader options
                      verbose  = True        ,   ## verbosity flag 
                      aux      = 0.9         ) : ## for Cuts method : efficiency cut-off
    """
    Helper function to add TMVA  response into dataset
    >>> tar_file = trainer.tar_file
    >>> dataset  = ...
    >>> inputs = [ 'var1' , 'var2' , 'var2' ]
    >>> dataset.addTMVAResponse (  inputs , tar_file , prefix = 'tmva_' )
    """
    assert isinstance ( dataset , ( ROOT.TTree , ROOT.RooAbsData ) ),\
           'Invalid dataset type!'

    from ostap.core.core           import cpp, std, Ostap
    from ostap.utils.progress_conf import progress_conf
    from ostap.utils.basic         import isatty

    _inputs       = _inputs2map_  ( inputs        )
    
    weights_files = WeightsFiles  ( weights_files ) 
    _map , _w     = _weights2map_ ( weights_files )
    
    options = opts_replace ( options , 'V:'      ,     verbose    )
    options = opts_replace ( options , 'Silent:' , not verbose    )
    options = opts_replace ( options , 'Color:'  ,     isatty  () )
    
    args = _inputs, _map, options, prefix , suffix , aux
    
    if   isinstance ( dataset , ROOT.TChain     ) :
        
        sc , newdata = _add_response_chain ( dataset , verbose , *args )
        if sc.isFailure() : logger.error ( 'Error from Ostap::TMVA::addResponse %s' % sc )
        return newdata
    
    elif isinstance ( dataset , ROOT.TTree      ) :
        
        sc , newdata = _add_response_tree  ( dataset , verbose , *args )
        if sc.isFailure() : logger.error ( 'Error from Ostap::TMVA::addResponse %s' % sc )
        
    else :

        ## add progress bar 
        if verbose :  sc = Ostap.TMVA.addResponse  ( dataset , progress_conf () , *args )
        else       :  sc = Ostap.TMVA.addResponse  ( dataset                    , *args )
        
        if sc.isFailure() : logger.error ( 'Error from Ostap::TMVA::addResponse %s' % sc )        
        newdata = dataset

    return newdata


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    
# =============================================================================
##                                                                      The END 
# =============================================================================
