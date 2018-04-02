#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file tmva.py
#  Python interface to basic TMVA functionality: Trainer and Reader 
#
#  Actually for the Trainer, it is a bit simplified version of Albert's code 
#   - thanks to Albert PUIG
#  Inspired from
#  @see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py
#
#  This file is a part of 
#  <a href="http://cern.ch/lhcb-comp/Analysis/Bender/index.html">Bender project</a>
#  <b>``Python-based Interactive Environment for Smart and Friendly 
#   Physics Analysis''</b>
#
#  The package has been designed with the kind help from
#  Pere MATO and Andrey TSAREGORODTSEV. 
#  And it is based on the 
#  <a href="http://cern.ch/lhcb-comp/Analysis/LoKi/index.html">LoKi project:</a>
#  ``C++ ToolKit for Smart and Friendly Physics Analysis''
#
#  By usage of this code one clearly states the disagreement 
#  with the smear campaign of Dr.O.Callot et al.: 
#  ``No Vanya's lines are allowed in LHCb/Gaudi software.''
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

This file is a part of BENDER project:
``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the LoKi project:
    ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement 
with the smear campain of Dr.O.Callot et al.: 
    ``No Vanya's lines are allowed in LHCb/Gaudi software.''

"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2013-10-02"
__version__ = '$Revision$'
__all__     = (
    "Trainer"         , ## the basic TMVA trainer 
    "Reader"          , ## the basic TMVA reader
    "addTMVAResponce" , ## add TMVA responce to RooDataSet
    "tmvaGUI"
    )
# =============================================================================
import ROOT, os
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.tmva' )
else                       : logger = getLogger ( __name__           )
# =============================================================================
pattern_XML   = "%s/weights/%s*.weights.xml"
pattern_CLASS = "%s/weights/%s*.class.C" 
# =============================================================================
## @class WeightFiles
#  helper structure  to deal with weights files
import ostap.utils.utils as Utils 
class WeightsFiles(Utils.CleanUp) :
    """Helper structure  to deal with weights files
    """
    def __init__ ( self , weights_files ) :

        ## string ? treat it a a single tar-file 
        if isinstance ( weights_files , str  ) :
            
            import tarfile, os  
            assert os.path.exists  ( weights_files ) , \
                   "Non-existing ``weights_file''  %s"   %  weights_files 
            
            if tarfile.is_tarfile ( weights_files ) :
                def xml_files ( archive ) :
                    for tarinfo in archive:
                        if os.path.splitext(tarinfo.name)[1] == ".xml":
                            yield tarinfo
                            
                with tarfile.open ( weights_files , 'r' ) as tar :
                    ## tar.list() 
                    xmls = [ f for f in xml_files ( tar ) ] 
                    tmpdir = self.tmpdir 
                    tar.extractall ( path = tmpdir , members = xml_files ( tar ) )
                    logger.debug ('Un-tar into temporary directory %s' % tmpdir ) 
                    weights_files  = [ os.path.join ( tmpdir, x.name ) for x  in xmls ]
                    self.tmpfiles += weights_files
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
               "Invalid type of ``weight_files''  %s "    % weights_files

        for method , xml in weights_files.iteritems() :            
            assert os.path.exists ( xml ) and os.path.isfile ( xml ), \
                   "No weights file '%s'" %  xml 
            
        self.__methods       = tuple ( [ i for  i in  weights_files.keys() ] )
        import copy
        self.__weights_files = copy.deepcopy ( weights_files ) 
        
    @property
    def methods ( self ) :
        "``methods'': the known methods from weights-file"
        return self.__methods
    @property
    def files   ( self ) :
        "``files'': the weigths file"
        import copy
        return copy.deepcopy ( self.__weights_files ) 
        
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
#  ...   signal          = treeSignal            ,  ## TTree for ``signal'' sample  
#  ...   background      = treeBackgrund         ,  ## TTree for ``background''  sample 
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
#  @thanks Albert PUIG
class Trainer(object):
    """Helper class to train TMVA:  
    
    >>> from ostap.tools.tmva import Trainer
    >>> t = Trainer(
    ...   ## TMVA methods  
    ...   methods =  [ ( ROOT.TMVA.Types.kMLP ,
    ...                'MLP',
    ...                'H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator' ) ] , 
    ...   variables = [ 'dtfchi2' , 'ctau', 'ptb' , 'vchi2' ] , ## list  of variables 
    ...   signal          = treeSignal            ,  ## TTree for ``signal'' sample  
    ...   background      = treeBackgrund         ,  ## TTree for ``background''  sample 
    ...   signal_cuts     = cuts_for_signal       ,
    ...   background_cuts = cuts_for_background   )

    Use the trainer: 
    >>> t.train()     ## Use the trainer

    Get the results:
    >>> xml      = t.weights_files ## get the weights XML   files 
    >>> classes  = t.class_files   ## get the classes (C++) files
    >>> output   = t.output_file   ## get the output ROOT file  
    >>> tar_file = t.tar_file      ## get the tar-file with XML&C++
    
    Actuially it is a bit simplified version of the original code by Albert PUIG,
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
    #  ...   signal          = treeSignal            ,  ## TTree for ``signal'' sample  
    #  ...   background      = treeBackgrund         ,  ## TTree for ``background''  sample 
    #  ...   signal_cuts     = cuts_for_signal       ,
    #  ...   background_cuts = cuts_for_background   )
    #  @endcode 
    #  For more detailes
    #  @see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
    def __init__(  self                       ,
                   methods                    ,
                   variables                  ,  ## list of variables 
                   signal                     ,  ## signal sample/tree
                   background                 ,  ## background sample/tree 
                   signal_cuts       = ''     ,  ## signal cuts 
                   background_cuts   = ''     ,  ## background cuts 
                   spectators        = []     ,
                   bookingoptions    = "Transformations=I;D;P;G,D" , 
                   configuration     = "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents" ,
                   signal_weight     = None   ,                
                   background_weight = None   ,
                   ##
                   output_file       = ''     ,  ## the name of output file 
                   verbose           = True   ,
                   name              = 'TMVA' ) :
        """Constructor with list of methods
        
        >>> from ostap.tools.tmva import Trainer
        >>> t = Trainer(
        ...   ## TMVA methods  
        ...   methods =  [ ( ROOT.TMVA.Types.kMLP ,
        ...                'MLP',
        ...                'H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator' ) ] , 
        ...   variables = [ 'dtfchi2' , 'ctau', 'ptb' , 'vchi2' ] , ## list  of variables 
        ...   signal          = treeSignal            ,  ## TTree for ``signal'' sample  
        ...   background      = treeBackgrund         ,  ## TTree for ``background''  sample 
        ...   signal_cuts     = cuts_for_signal       ,
        ...   background_cuts = cuts_for_background   )        
        - For more detailes
        see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
        """
        
        self.__methods           = tuple ( methods    )
        self.__variables         = tuple ( variables  )
        
        from ostap.trees.trees import Chain
        if isinstance ( signal     , Chain ) : signal     =     signal.chain ()
        if isinstance ( background , Chain ) : background = background.chain ()
        
        self.__signal            = signal
        self.__signal_cuts       = signal_cuts  
        self.__signal_weight     = signal_weight
        
        self.__background        = background
        self.__background_cuts   = background_cuts 
        self.__background_weight = background_weight
        
        self.__spectators        = [] 

        self.__verbose = True if verbose else False 
        self.__name    = name 

        self.__bookingoptions   = bookingoptions
        self.__configuration    = configuration
        
        ## outputs 
        self.__weights_files = []
        self.__class_files   = []
        self.__output_file   = output_file if output_file else '%s.root' % self.name
        self.__tar_file      = None 
        self.__log_file      = None 
        #
        ## minor adjustment
        #
        opts = self.__bookingoptions
        if      self.verbose and 0 <= opts.find ('!V:')   : opts.replace('!V:', 'V:')
        if not  self.verbose and 0 > opts.find ('Silent') : opts+= ':Silent'
        
        if 0 > opts.find ( "AnalysisType=" ) :
            opts += ":AnalysisType=Classification"
            logger.debug('Trainer(%s): booking options are appended with ":AnalysisType=Classification"' % name )
            
        if  self.verbose :
            
            if   0 <= opts.find ( '!V:'      ) : opts.replace ('!V:','V:')
            elif 0 <= opts.find (  'V:'      ) : pass 
            else                               : opts += ':V:'

            if   0 <= opts.find ( '!Silent:' ) : pass
            elif 0 <= opts.find (  'Silent:' ) : opts.replace ('Silent:' , '!Silent:')
            else                               : opts += ':!Silent:'
            
        else :
            
            if   0 <= opts.find ( '!V:'      ) : pass
            elif 0 <= opts.find (  'V:'      ) : opts.replace ('V:','!V:')
            else                               : opts += ':!V:'
            
            if   0 <= opts.find ( '!Silent:' ) : opts.replace ('!Silent:' , 'Silent:')
            elif 0 <= opts.find (  'Silent:' ) : pass
            else                               : opts += ':Silent:'


        from ostap.utils.basic import isatty 
        if isatty() and self.verbose :
            
            if   0 <= opts.find ( '!DrawProgressBar' ) : opts.replace ('!DrawProgressBar', 'DrawProgressBar' )
            elif 0 <= opts.find (  'DrawProgressBar' ) : pass
            else                                       : opts += ':DrawProgressBar:'
            
            if   0 <= opts.find ( '!Color' )           : opts.replace ('!Color', 'Color' )
            elif 0 <= opts.find (  'Color' )           : pass
            else                                       : opts += ':Color:'
            
        else :
            
            if   0 <= opts.find ( '!DrawProgressBar' ) : pass  
            elif 0 <= opts.find (  'DrawProgressBar' ) : opts.replace ('DrawProgressBar', '!DrawProgressBar' ) 
            else                                       : opts += ':!DrawProgressBar:'
            
            if   0 <= opts.find ( '!Color' )           : pass 
            elif 0 <= opts.find (  'Color' )           : opts.replace ('Color', '!Color' )
            else                                       : opts += ':!Color:'



        self.__bookingoptions = opts

        ## clean-up 
        dirname    = str(self.name)
        for s in ( ' ' , '%' , '!' , '>' , '<' , '\n' , '?' ) :  
            while s in dirname : dirname = dirname.replace ( ' ' , '_' )
            
        pattern_xml = pattern_XML   % ( dirname ,  dirname )
        pattern_C   = pattern_CLASS % ( dirname ,  dirname )

        rf = []
        import glob,os
        for f in glob.glob ( pattern_xml ) :
            rf.append ( f ) 
            os.remove ( f ) 
        for f in glob.glob ( pattern_C   ) :
            rf.append ( f ) 
            os.remove ( f )
        if rf : logger.debug ( "Trainer(%s): remove existing xml/class-files %s" % ( self.name , rf ) ) 
        
        self.__dirname     = dirname
        self.__pattern_xml = pattern_xml 
        self.__pattern_C   = pattern_C 
        
    @property
    def name    ( self ) :
        """``name''    : the name of TMVA trainer"""
        return self.__name
    
    @property
    def methods ( self ) :
        """``methods'' : the list of TMVA methods to be used"""
        return tuple(self.__methods)
    
    @property
    def variables ( self ) :
        """``variables'' : the list of variables  to be used for training"""
        return tuple(self.__variables)

    @property
    def spectators ( self ) :
        """``spectators'' : the list of spectators to be used"""
        return tuple(self.__spectators)

    @property
    def signal ( self ) :
        """``signal'' :  TTree for signal events"""
        return self.__signal
    
    @property
    def signal_cuts ( self ) :
        """``signal_cuts'' :  cuts to be applied for ``signal'' sample"""
        return str(self.__signal_cuts)

    @property
    def signal_weight ( self ) :
        """``signal_weight'' : weight to be applied for ``signal'' sample"""
        return self.__signal_weight
    

    @property
    def background ( self ) :
        """``background'' :  TTree for background events"""
        return self.__background
    
    @property
    def background_cuts ( self ) :
        """``background_cuts'' :  cuts to be applied for ``backgroud'' sample """
        return str(self.__background_cuts)

    @property
    def background_weight ( self ) :
        """``background_weight'' : weight to be applied for ``background'' sample"""
        return self.__background_weight
    
    @property
    def bookingoptions ( self ) :
        """``bookingoptions'' : options used to book TMVA::Factory"""
        return str(self.__bookingoptions)

    @property
    def configuration ( self ) :
        """``configuration'' : options used to book TMVA"""
        return str(self.__configuration)
    
    @property
    def verbose ( self ) :
        """``verbose'' : verbosity  flag"""
        return self.__verbose

    @property
    def dirname  ( self ) :
        """``dirname''  : the output directiory name"""
        return str(self.__dirname) 

    @property
    def weights_files ( self ) :
        """``weights_files'' : the list/tuple of final files with TMVA weights"""
        return tuple(self.__weights_files)
    @property
    def class_files ( self ) :
        """``class_files'' : the list/tuple of final files with TMVA classes"""
        return tuple(self.__class_files)
    @property
    def output_file ( self ) :
        """``output_file''  : the output file (if any)"""
        return str(self.__output_file) if  self.__output_file else None
    
    @property
    def tar_file ( self ) :
        """``tar_file''  : the compressed (gz) tar file with results"""
        return str(self.__tar_file) if self.__tar_file else None 
    @property
    def log_file ( self ) :
        """``log_file''  : the name of log-file """
        return str(self.__log_file) if self.__log_file else None 

    # =========================================================================
    ## train TMVA 
    #  @code
    #  trainer.train ()
    #  @endcode  
    #  @return the name of output XML file with the weights 
    def train ( self ) :
        """Train the TMVA:
        - returns the names of output XML file with the weights 
        >>> trainer.train () 
        """
        
        import glob,os
        rf = [] 
        for f in glob.glob ( self.__pattern_xml ) :
            rf.append ( f ) 
            os.remove ( f ) 
        for f in glob.glob ( self.__pattern_C   ) :
            rf.append ( f )
            os.remove ( f )
        if rf : logger.debug ( "Trainer(%s): remove existing xml/class-files %s" % ( self.name , rf ) ) 
        
        ROOT.TMVA.Tools.Instance()

        outFile = ROOT.TFile.Open   ( self.output_file, 'RECREATE' )        
        logger.debug ( 'Trainer(%s): output ROOT file: %s ' % ( self.name , outFile.GetName() ) )

        factory = ROOT.TMVA.Factory (
            self.name             ,
            outFile               ,
            self.bookingoptions   )
        logger.debug ( 'Trainer(%s): book TMVA-factory %s ' % ( self.name , self.bookingoptions ) )
        factory.SetVerbose( self.verbose )
        
        ## 
        dataloader = ROOT.TMVA.DataLoader ( self.dirname )
        
        #
        for v in self.variables :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )
            dataloader.AddVariable  ( *vv )
            
        for v in self.spectators :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )             
            dataloader.AddSpectator ( *vv )
        #
        if self.signal_cuts :
            logger.info ( "Trainer(%s): Signal       cuts:``%s''" % ( self.name ,     self.signal_cuts ) ) 
        if self.background_cuts :
            logger.info ( "Trainer(%s): Background   cuts:``%s''" % ( self.name , self.background_cuts ) )
        # 
        dataloader.AddTree ( self.signal     , 'Signal'     , 1.0 , ROOT.TCut ( self.    signal_cuts ) )
        dataloader.AddTree ( self.background , 'Background' , 1.0 , ROOT.TCut ( self.background_cuts ) )
        #
        if self.signal_weight :
            dataloader.SetSignalWeightExpression     ( self.signal_weight     )
            logger.info ( "Trainer(%s): Signal     weight:``%s''" % ( self.name ,     self.signal_weight ) )            
        if self.background_weight :
            dataloader.SetBackgroundWeightExpression ( self.background_weight )
            logger.info ( "Trainer(%s): Background weight:``%s''" % ( self.name , self.background_weight ) )
            
        logger.info     ( "Trainer(%s): Configuration    :``%s''" % ( self.name , self.configuration ) )
        dataloader.PrepareTrainingAndTestTree(
            ROOT.TCut ( self.signal_cuts     ) ,
            ROOT.TCut ( self.background_cuts ) ,
            self.configuration            )
        #
        for m in self.methods :
            factory.BookMethod ( dataloader , *m )

        from ostap.logger.utils import tee_cpp , NoContext
        logfile = self.name + '.log'
        if os.path.exists ( logfile ) :
            try    : os.remove ( logfile )
            except : pass
        ##with tee_cpp ( log_file ) if self.verbose else  NoContext() :
        with NoContext() :
            # Train MVAs
            ms = tuple( i[1] for i in  self.methods )
            logger.info  ( "Trainer(%s): Train    all methods %s " % ( self.name , ms ) )
            factory.TrainAllMethods    ()
            # Test MVAs
            logger.info  ( "Trainer(%s): Test     all methods %s " % ( self.name , ms ) )
            factory.TestAllMethods     ()
            # Evaluate MVAs
            logger.info  ( "Trainer(%s): Evaluate all methods %s " % ( self.name , ms ) )
            factory.EvaluateAllMethods ()
            
        # Save the output.
        logger.debug ( "Trainer(%s): Output ROOT file %s is closed" % ( self.name , outFile.GetName() ) )
            
        outFile.Close()

        import glob, os 
        self.__weights_files = [ f for f in glob.glob ( self.__pattern_xml ) ]
        self.__class_files   = [ f for f in glob.glob ( self.__pattern_C   ) ]
        self.__log_file      = logfile if os.path.exists ( logfile ) else None 
        
        del dataloader
        del factory 

        tfile = self.name + '.tgz'
        if os.path.exists ( tfile ) :
            logger.debug  ( "Trainer(%s): Remove existing tar-file %s" % ( self.name , tfile ) )

        import tarfile
        with tarfile.open ( tfile , 'w:gz' ) as tar :
            for x in self.weights_files : tar.add ( x )
            for x in self.  class_files : tar.add ( x )
            logger.info  ( "Trainer(%s): Tar/gz  file  : %s" % ( self.name , tfile ) ) 
            if self.verbose : tar.list ()

        ## finally set tar-file 
        if os.path.exists ( tfile ) and tarfile.is_tarfile( tfile ) :
            self.__tar_file = tfile 

        logger.info  ( "Trainer(%s): Weights files : %s" % ( self.name , self.weights_files ) )
        logger.debug ( "Trainer(%s): Class   files : %s" % ( self.name , self.class_files   ) ) 
        logger.info  ( "Trainer(%s): Output  file  : %s" % ( self.name , self.output_file   ) ) 
            
        return self.weights_files


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
#  extracted frmo the file name, otherwise it needs to be specified as dictionary   
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
#  ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   
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
#  ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   
#  @endcode
#
#  - It it natually merges with Ostap's SelectorWithVars utility
#
#  @see TMVA::Reader
#  @date   2013-10-02
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
#  @thanks Alexander BARANOV
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

    ``weights_files'' can be :
    - single xml-file with weights for the given method
    - single tar/tgz/tar.gz-file with weights files  (output of ``Trainer.tar_file'')
    - list of xml-files with weights                 (output of ``Trainer.weights_files'')

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
    ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   

    - A bit more efficient form is :
    >>> tree =  ....  ## TTree/TChain/RooDataSet with data
    >>> mlp_fun  =  reader.MLP
    >>> bdgt_fun =  reader.BDTG
    >>> for entry in tree :
    ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
    ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
    ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)
    
    It it natually merges with Ostap's SelectorWithVars utility 
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
    def __init__ ( self          ,
                   name          , 
                   variables     ,
                   weights_files ) :
        """Constrct the reader         
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
        
        ROOT.TMVA.Tools.Instance()
        
        self.__reader = ROOT.TMVA.Reader()
        self.__name   = name

        ## treat the weigths files
        weights = WeightsFiles ( weights_files )

        ##  book the variables:
        #   dirty trick with arrays is needed due to a bit strange reader interface.
        #   [TMVA reader needs the address of ``float''(in C++ sense) variable]
        from array import array

        self.__variables = []
        
        for v in variables : 

            if   isinstance ( v , str ) :
                
                vname  = v
                fvun   = lambda s : getattr ( s , vname )
                vfield = array ( 'f' , [1] )                  ## NB: note the type 
                
            elif isinstance ( v , tuple ) and 2 == len ( v ) :

                vname  = v[0]
                vfun   = v[1]
                vfield = array ( 'f' , [1] )                  ## NB: note the type here 

            else :
                
                logger.error ('Reader(%s): Invalid variable description!' % name )
                raise AttributeError, 'Invalid variable description!'

            ##                     name    accessor   address 
            self.__variables += [ ( vname , vfun     , vfield  ) ] 


        self.__variables = tuple ( self.__variables )
        
        ## declare all variables to TMVA.Reader 
        for v in self.__variables :
            self.__reader.AddVariable ( v[0] , v[2] )            

        
        self.__methods = weights.methods

        for method , xml in  weights.files.iteritems() :
            m = self.__reader.BookMVA ( method , xml )
            assert  m , 'Error in booking %s/%s' % (  method  , xml )
            logger.debug ('TMVA Reader(%s) is booked for method:%s xml: %s' % (  self.__name ,
                                                                                 method      ,
                                                                                 xml         ) )
            
        logger.info ('TMVA Reader(%s) booked methods are %s' %  ( self.__name , self.__methods ) )
        self.__methods = tuple ( self.__methods )

    @property
    def name ( self ) :
        """``name'' - the name of the reader"""
        return self.__name

    @property
    def reader ( self ) :
        """``reader'' - the  actual TMVA.Reader object"""
        return self.__reader
    
    @property
    def methods ( self ) :
        """``methods'' - the  list/tuple of booked TMVA methods"""
        return tuple (self.__methods)

    @property
    def variables ( self ) :
        """``variables'' - helper structure to access TMVA variables
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
            """``nvars'' : number of TMVA variables"""
            return self.__nvars
        @property
        def reader ( self ) :
            """``reader'' : TMVA reader """
            return self.__reader
        @property
        def method ( self ) :
            """``method'' : TMVA method name"""
            return self.__method
        # =====================================================================
        ## the main method 
        def __call__ ( self , entry , cut_efficiency = 0.9  ) :
            return self.eval ( entry , cut_efficiency )
        # =====================================================================
        ## Evaluate the method from TTree/TCahin/RooAbsData using the accessors, defined  early
        # @code 
        # tree   = ...
        # method = ...
        # print 'Response is %s' % method.eval ( tree )
        # @endcode 
        def eval ( self , entry , cut_efficiency = 0.9 ) :
            """Evaluate the method fomr TTree/RooAbsData using the
            accessors, defined  early
            >>> tree   = ...
            >>> method = ...
            >>> print 'Response is %s'    % method.eval ( tree ) 
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
        #  print 'Responce %s' % method ( tree )
        #  @endcode
        #  or using parameters
        #  @code
        #  pt , y, phi = ...
        #  print 'Responce %s' % method ( pt , y , phi )
        #  @endcode 
        def __call__ ( self , arg ,  *args ) :
            if isinstance ( arg , ( float , int , long , bool ) ) :
                return self.evaluate ( arg , *args ) 
            return self.eval ( arg , *args )
        # =====================================================================
        ## Evaluate the method from parameters 
        # @code 
        # method       = ...
        # pt, eta, phi = 5 ,  3.0 , 0  ## variables 
        # print 'Response is %s'    % method.evaluate ( pt ,  eta , phi ) 
        # @endcode 
        def evaluate ( self , *args ) :
            """Evaluate the method from parameters 
            >>> method       = ...
            >>> pt, eta, phi = 5 ,  3.0 , 0  ## variables 
            >>> print 'Response is %s'    % method.evaluate ( pt ,  eta , phi ) 
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
    #  ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)
    # @encode        
    def __getitem__ ( self , method ) :
        """Helper utility to  get the correspondig function from the  reader:
        - Use the reader
        >>> tree =  ....  ## TTree/TChain/RooDataSet with data
        >>> mlp_fun  =  reader['MLP']  ## <-- here! 
        >>> bdgt_fun =  reader['BDTG'] ## <-- here!
        >>> for entry in tree :
        ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
        ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
        ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   
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
    #  ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)
    # @encode        
    def __getattr__ ( self , method ) :
        """Helper utility to  get the correspondig function from the  reader:
        - Use the reader
        >>> tree =  ....  ## TTree/TChain/RooDataSet with data
        >>> mlp_fun  =  reader.MLP  ## <-- here! 
        >>> bdgt_fun =  reader.BDTG ## <-- here!
        >>> for entry in tree :
        ...     mlp  = mlp_fun  ( entry )  ## evaluate MLP-TMVA
        ...     bdtg = bdtg_fun ( entry )  ## evalaute BDTG-TMVA
        ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   
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
    #  ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   
    #  @endcode 
    #  @attention it is *not* CPU efficient
    #  Ugly trick with arrays is needed due to some technical problems
    #  (actually TMVA reader needs the address of ``float''(in C++ sense) variable
    def __call__ ( self , method , entry , cut_efficiency = 0.90 ) :
        """Evaluate TMVA
        - Use the reader
        >>> tree =  ....  ## TTree/TChain/RooDataSet with data
        >>> for entry in tree :
        ...     mlp  = reader ( 'MLP'  , entry )  ## evaluate MLP-TMVA
        ...     bdtg = reader ( 'BDTG' , entry )  ## evalaute BDTG-TMVA
        ...     print 'MLP/BDTG for  this event are %s/%s' %  (mlp , bdtg)   
        - It is not CPU efficient :-( 
        - Ugly trick with arrays is needed due to some pure technical problem
        [actually TMVA reader needs the address of ``float''(in C++ sense) variable]
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
    #  print 'MLP response is: ', reader.MLP.evaluate ( pt , y )
    #  @endcode
    def evaluate ( self , method , *args ) :
        """Evaluate TMVA
        >>> reader = ...
        >>> pt, y  = ...  ##
        >>> print 'MLP response is: ', reader.MLP.evaluate ( pt , y )
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
        for k , v  in inputs.itertems() : _inputs[k] = v
    elif isinstance ( inputs , ( tuple , list ) ) :
        for i in inputs :
            if isinstance ( i , str ) : k , v = i , i
            else                      : k , v = i
            _inputs[k] = v 

    ## 
    assert not _inputs .empty() and _inputs.size() == len ( inputs ), \
           'Invalid MAP size %s for %s' % ( _inputs.size() , inputs ) 

    return _inputs 
# =============================================================================
## convert weights structure to Ostap.TMVA.PAIRS 
def _weights2map_ ( weights_files ) :
    
    from ostap.core.core import cpp, std, Ostap
    MAP = std.map   ( 'std::string', 'std::string' )
    
    weights  = WeightsFiles ( weights_files )
    _weights = MAP() 
    for method , xml in weights.files.iteritems() :
        _weights [ method ] = xml

    assert not _weights .empty() , \
           'Invalid MAP size %s for' % ( _weights.size() , weights )
    
    assert not _weights.empty() , "Invalid weights_files: %s"  % weights.files
    return _weights 
    
# =============================================================================
## Helper function to add TMVA response into dataset
#  @code
#  tar_file = trainer.tar_file
#  dataset  = ...
#  inputs = [ 'var1' , 'var2' , 'var2' ]
#  dataset.addTMVAResponse (  inputs , tar_file , prefix = 'tmva_' )
#  @endcode 
def addTMVAResponse ( dataset        ,
                      inputs         ,
                      weights_files  ,
                      prefix   = ''  , 
                      suffix   = ''  ,
                      aux      = 0.9 ) :
    """
    Helper function to add TMVA  responce into dataset
    >>> tar_file = trainer.tar_file
    >>> dataset  = ...
    >>> inputs = [ 'var1' , 'var2' , 'var2' ]
    >>> dataset.addTMVAResponce (  inputs , tar_file , prefix = 'tmva_' )
    """
    from ostap.core.core import cpp, std, Ostap
    PP = std.pair   ( 'std::string', 'std::string' )
    VP = std.vector ( PP )
    
    _inputs  = _inputs2map_  ( inputs        )
    _weights = _weights2map_ ( weights_files )
    
    sc = Ostap.TMVA.addResponse ( dataset  ,
                                  _inputs  ,
                                  _weights ,
                                  prefix   ,
                                  suffix   ,
                                  aux      )
    if sc.isFailure() :
        logger.error ( 'Error from Ostap::TMVA::addResponse %s' % sc )
    return sc 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    
# =============================================================================
# The END 
# =============================================================================
