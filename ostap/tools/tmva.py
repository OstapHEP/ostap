3#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file py.py
#
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
-  Trainer
-  Reader 
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
    "Trainer" ,
    "Reader"  ,
    "tmvaGUI"
    )
# =============================================================================
import ROOT
# =============================================================================
from   ostap.logger.logger  import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.tmva' )
else                       : logger = getLogger( __name__ )
# =============================================================================
    
pattern_XML   = "%s/weights/%s*.weights.xml"
pattern_CLASS = "%s/weights/%s*.class.C" 

# =============================================================================
## @class TMVATrainer
#  Helper class to train TMVA
#
#  @code
#
#  from ostap.tools.tmva import Trainer 
#  t = Trainer( methods =  [
#  ## type, name, configuration 
#  ( ROOT.TMVA.Types.kMLP ,
#    "MLP",
#    "H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator"
#  ) ] )
#
#  varlist = [
#    'dtfchi2' , 
#    'ctau'    , 
#    'ptb'     , 
#    'vchi2'   
#  ]
#
#  t.train ( var_list                              ,
#            signal          = treeSignal          ,
#            background      = treeBackgrund       ,
#            outputfile      = 'output.root'       ,
#            signal_cuts     = cuts_for_signal     ,
#            background_cuts = cuts_for_background ,
#            spectators      = []                  ) 
#
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
    #
    #  from ostap.tools.tmva import Trainer 
    #  t = Trainer( methods =  [
    #  ## type, name, configuration 
    #  ( ROOT.TMVA.Types.kMLP ,
    #    'MLP',
    #    'H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator'
    #  ) ] )
    #
    #  varlist = [
    #    'dtfchi2' , 
    #    'ctau'    , 
    #    'ptb'     , 
    #    'vchi2'   
    #  ]
    #
    #  ## start the actual training 
    #  t.train ( var_list                              ,
    #            signal          = treeSignal          ,
    #            background      = treeBackgrund       ,
    #            outputfile      = 'output.root'       ,
    #            signal_cuts     = cuts_for_signal     ,
    #            background_cuts = cuts_for_background ,
    #            spectators      = []                  ) 
    #
    Actuially it is a bit simplified version of the original code by Albert PUIG,
    inspired from
    http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py
    """
    # =========================================================================
    ## constructor
    #  @code
    # 
    #  from ostap.tools.tmva import Trainer 
    #  t = Trainer( methods = [
    #  ## type                   name   configuration 
    #  ( ROOT.TMVA.Types.kMLP , "MLP", "H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator" ) 
    #  ] )
    #
    #  @endcode
    #  For more detailes
    #  @see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
    def __init__(self, methods , verbose = True ,  name = 'TMVA' ):
        """Constructor with list of methoods
        >>> from ostap.tools.tmva import Trainer
        >>> methods = ....
        >>> trainer = Trainer ( methods )        
        For more detailes
        see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
        """
        self.name    = name 
        self.methods = methods
        self.verbose = verbose 

        ROOT.TMVA.Tools.Instance()
        
    ## define the verbosity 
    def setVerbose(self, verbosity ): self.verbose =  bool( verbosity )

    # =========================================================================
    ## train TMVA 
    #  @code  
    #  varlist = [
    #    'dtfchi2'       , 
    #    'ctau'          ,
    #    'ptb'           , 
    #    'vchi2'         , 
    #  ]
    #
    #  trainer.train ( var_list                       ,
    #          signal          = treeSignal          ,
    #          background      = treeBackgrund       ,
    #          outputfile      = 'output.root'       ,
    #          signal_cuts     = cuts_for_signal     ,
    #          background_cuts = cuts_for_background ,
    #          spectators      = []                  ) 
    #  @endcode  
    #  @return the name of output XML file with weights 
    def train ( self   , var_list               ,
                signal , background             ,
                output_file       = ''          ,
                signal_cuts       = ''          ,
                background_cuts   = ''          ,
                spectators        = []          ,
                bookingoptions    = "Transformations=I;D;P;G,D" , 
                configuration     = "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" ,
                signal_weight     = None        ,                
                background_weight = None        ,
                ) :
        """Train the TMVA:
        >>> trainer.train ( var_list , ... ) 
        """
        #
        ## get the logger
        #
        import os
        name = self.name
        if not output_file : output_file = '%s.root'  % name 
        
        dirname    = str(self.name)
        for s in ( ' ' , '%' , '!' , '>' , '<' , '\n' , '?' ) :  
            while s in dirname : dirname = dirname.replace ( ' ' , '_' )
            
        pattern_xml = pattern_XML   % ( dirname ,  dirname )
        pattern_C   = pattern_CLASS % ( dirname ,  dirname )
        
        import glob,os
        for f in glob.glob ( pattern_xml ) :
            logger.debug ( 'Remove existing weight-file %s' % dirname ) 
            os.remove ( f ) 
        for f in glob.glob ( pattern_C   ) :
            logger.debug ( 'Remove existing class-file  %s' % dirname ) 
            os.remove ( f ) 

        outFile = ROOT.TFile.Open   ( output_file, 'RECREATE' )
        logger.info ( 'Trainer(%s): output ROOT file: %s ' % ( name , output_file  ) )

        if self.verbose and 0 > bookingoptions.find ('Silent') : bookingoptions+= ':Silent'
        #
        if 0 > bookingoptions.find( "AnalysisType=" ) :
            bookingoptions += ":AnalysisType=Classification"
            logger.info('Trainer(%s): booking options are appended with ":AnalysisType=Classification"' % name )
            
        if self.verbose and 0 <= bookingoptions.find ('!V:') :
            bookingoptions.replace('!V:', 'V')

        if   0 <= bookingoptions.find ('!V:') : pass 
        elif 0 <= bookingoptions.find ( 'V:') : pass 
        elif self.verbose : bookingoptions +=  ':V:'
        else              : bookingoptions += ':!V:'


        if   0 <= bookingoptions.find ('!Silent') : pass 
        elif 0 <= bookingoptions.find ( 'Silent') : pass 
        elif self.verbose : bookingoptions += ':!Silent'
        else              : bookingoptions +=  ':Silent'
        
        if   0 <= bookingoptions.find ('!Color') : pass 
        elif 0 <= bookingoptions.find ( 'Color') : pass 
        elif self.verbose : bookingoptions +=  ':Color'
        else              : bookingoptions += ':!Color'

        if   0 <= bookingoptions.find ('!DrawProgressBar') : pass 
        elif 0 <= bookingoptions.find ( 'DrawProgressBar') : pass 
        elif self.verbose : bookingoptions +=  ':DrawProgressBar'
        else              : bookingoptions += ':!DrawProgressBar'


        factory = ROOT.TMVA.Factory (
            self.name             ,
            outFile               ,
            bookingoptions        )
        logger.info ( 'Trainer(%s): book TMVA-factory %s ' % ( name , bookingoptions ) )


        dataloader = ROOT.TMVA.DataLoader ( dirname ) 
        
        factory.SetVerbose(self.verbose)
        #
        for v in var_list :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )
            dataloader.AddVariable  ( *vv )
            
        for v in spectators :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )             
            dataloader.AddSpectator ( *vv )
        #
        signalWeight     = 1.0
        backgroundWeight = 1.0
        #
        if signal_cuts :
            logger.info ( 'Trainer(%s): Signal       cuts: "%s" ' % ( name ,     signal_cuts ) ) 
        if background_cuts :
            logger.info ( 'Trainer(%s): Background   cuts: "%s" ' % ( name , background_cuts ) )
        # 
        dataloader.AddTree ( signal     , 'Signal'     ,     signalWeight ,
                             ROOT.TCut (      signal_cuts ) )
        dataloader.AddTree ( background , 'Background' , backgroundWeight ,
                             ROOT.TCut (  background_cuts ) )
        #
        if signal_weight :
            dataloader.SetSignalWeightExpression     ( signalweight     )
            logger.info ( 'Trainer(%s): Signal     weight: "%s" ' % ( name ,     signal_weight ) )
            
        if background_weight :
            dataloader.SetBackgroundWeightExpression ( backgroundweight )
            logger.info ( 'Trainer(%s): Background weight: "%s" ' % ( name , background_weight ) )
            
        logger.info ( 'Trainer(%s): Configuration  : "%s" ' % ( name , configuration ) )
        dataloader.PrepareTrainingAndTestTree(
            ROOT.TCut ( signal_cuts     ) ,
            ROOT.TCut ( background_cuts ) ,
            configuration                 )
        #
        for m in self.methods :
            factory.BookMethod ( dataloader , *m )

        # Train MVAs
        ms = tuple( i[1] for i in  self.methods )
        logger.info  ( "Trainer(%s): Train    all Methods %s " % ( name , ms ) )
        factory.TrainAllMethods    ()
        # Test MVAs
        logger.info  ( "Trainer(%s): Test     all Methods %s " % ( name , ms ) )
        factory.TestAllMethods     ()
        # Evaluate MVAs
        logger.info  ( "Trainer(%s): Evaluate all Methods %s " % ( name , ms ) )
        factory.EvaluateAllMethods ()
        # Save the output.
        logger.debug ( "Trainer(%s): Output ROOT file %s is closed" % ( name , output_file ) )  
        outFile.Close()
        
        # get the weights files
        import glob, os
        
        self.weight_files = [ f for f in glob.glob ( pattern_xml ) ]
        self.class_files  = [ f for f in glob.glob ( pattern_C   ) ]
        self.output_file  = output_file
        
        logger.info  ( "Trainer(%s): Weights files : %s" % ( name , self.weight_files ) )
        logger.info  ( "Trainer(%s): Output  file  : %s" % ( name , self.output_file  ) ) 
        logger.debug ( "Trainer(%s): Class   files : %s" % ( name , self.class_files  ) ) 
        
        del dataloader
        del factory 

        return self.weight_files[:]
    
# =============================================================================
## @class Reader
#  Rather generic python interface to TMVA-reader
#  @attention It is *not* CPU-efficient
#  Ugly tricks with arrays are required to bypass some technical limitations
#
#  @code
#
#  r = Reader ( 'MyTMVA' ,
#       variables = [
#       ## name      accessor  
#       ( 'pt'   , lambda s : s.pt ) ,
#       ## name      accessor  
#       ( 'ip'   , lambda s : s.ip ) ,
#       ## name     
#         'var1'                     , ## use s.var1 
#       ## name     
#         'var2'                     , ## use s.var2 
#       ] ,
#       weights_file = 'my_weights.xml'
#      )
#  
#  @endcode
#  @see TMVA::Reader
#  @date   2013-10-02
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
#  @thanks Alexander BARANOV
class Reader(object)  :
    """ Rather generic python interface to TMVA-reader
    #
    #  r = Reader ( 'MyTMVA' ,
    #       variables = [
    #       ## name      accessor  
    #       ( 'pt'   , lambda s : s.pt ) ,
    #       ## name      accessor  
    #       ( 'ip'   , lambda s : s.ip ) ,
    #       ## name     
    #         'var1'                     , ## use s.var1 
    #       ## name     
    #         'var2'                     , ## use s.var2 
    #       ] ,
    #       weights_file = 'my_weights.xml'
    #      )
    """
    def __init__ ( self          ,
                   name          , 
                   variables     ,
                   weights_files ) :
        
        ROOT.TMVA.Tools.Instance()
        
        self.reader = ROOT.TMVA.Reader()
        self.name   = name

        ##  book the variables:
        #   dirty trick with arrays is needed due to a bit strange reader interface.
        #   [TMVA reader needs the address of ``float''(in C++ sense) variable]
        from array import array

        self._variables = []
        
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
            self._variables += [ ( vname , vfun     , vfield  ) ] 


        ## declare all variables to TMVA.Reader 
        for v in self._variables :
            self.reader.AddVariable ( v[0] , v[2] )            

        import os
        
        if   isinstance ( weights_files , str  ) : weights_files = [ weights_file ]
        elif isinstance ( weights_files , dict ) :
            weights_files = [ (k,v) for k,v in weights_files.iteritems() ]
        
        self.methods = []
        
        for wf in  weights_files :
            
            if isinstance ( wf , str ) : method , xml = None, wf 
            else                       : method , xml =       wf
            
            if not os.path.exists ( xml ) or not os.path.isfile ( xml ) : 
                raise IOError("No weights file '%s'"  %  xml )

            if not method :
                f      = os.path.split ( xml ) [-1]
                p,s,m  = f.rpartition('_')
                method = m[:m.find('.weights.xml')]

            if not method :
                if not s : raise AttributeError("Can't extract method name from %s" % xml )
                
                
            self.reader.BookMVA ( method , xml )
            self.methods.append ( method ) 
            logger.info ('TMVA Reader(%s) is booked for method:%s xml: %s' % (  self.name ,
                                                                                method    ,
                                                                                xml       ) )
        logger.info ('TMVA Reader(%s) booked methods are %s' %  ( self.name , self.methods ) )
        self.methods = tuple ( self.methods )


    # =========================================================================
    ## helper class to get TMVA decision for certain method 
    class Var (object) :
        """Helper class to get TMVA decision for certain method
        >>>  reader = ...
        >>>  var = reader[ method ]
        >>>  val = var ( entry )
        """
        def __init__ ( self , reader , method ) :
            self.reader = reader
            self.method = method
        def __call__ ( self , entry , cut_efficiency = 0.9 ) :
            return self.reader( self.method , entry , cut_efficiency )
            ##v   =  self.reader( self.method , entry , cut_efficiency )
            ## print v 
            ##return v 

    ## =======================================================================
    def __getitem__ ( self , method ) :
        if not method in self.methods :
            return KeyError( 'No method %s is booked!' %  method )
        return Reader.Var  ( self , method )
    
    ## =======================================================================
    def __getattr__ ( self , method ) :
        if not method in self.methods :
            return AttributeError( 'No method %s is booked!' %  method )
        return Reader.Var  ( self , method ) 

    # =========================================================================
    ## evaluate TMVA
    #  @attention it is *not* CPU efficient
    #  Ugly trick with arrays is needed due to some technical problems
    #  (actually TMVA reader needs the address of ``float''(in C++ sense) variable
    def __call__ ( self , method , entry , cut_efficiency = 0.90 ) :
        """Evaluate TMVA
        - It is not CPU efficient :-( 
        - Ugly trick with arrays is needed due to some pure technical problem
        [actually TMVA reader needs the address of ``float''(in C++ sense) variable]
        """
        
        ## loop over all variables 
        for v in self._variables :
            vfun    = v[1]           ## accessor function 
            v[2][0] = vfun ( entry ) ## fill variable from the tree/chain 
            
        ## evaluate TMVA 
        return self.reader.EvaluateMVA ( method , cut_efficiency ) 


_canvas = []
# =============================================================================
## start TMVA gui 
def tmvaGUI ( filename , new_canvas = True ) :
    """
    Start TMVA-GUI
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
