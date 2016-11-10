#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
# $Id$
# ==========================================================================================
## @file py.py
#
#  Python interface to basic TMVA functionality: Trainer and Reader 
#
#  Actually for the Trainer, it is a bit simplified version of Albert's code 
#  @thanks Albert PUIG
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
#  @thanks Albert PUIG
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
    
# =============================================================================
## @class TMVATrainer
#  Helper class to train TMVA
#
#  @code
#
#  from PyTMVA import Trainer 
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
    #  from PyTMVA import Trainer 
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
    #  from TMVA import Trainer 
    #  t = Trainer( methods = [
    #  ## type                   name   configuration 
    #  ( ROOT.TMVA.Types.kMLP , "MLP", "H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator" ) 
    #  ] )
    #
    #  @endcode
    #  For more detailes
    #  @see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
    def __init__(self, methods , verbose = True ):
        """Constructor with list of methoods
        >>> from PyTMVA import Trainer
        >>> methods = ....
        >>> trainer = Trainer ( methods )        
        For more detailes
        see http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/cvs/tmva/test/TMVAClassification.py.
        """
        self.methods = methods
        self.verbose = verbose 

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
    def train ( self   , var_list             ,
                signal , background           ,
                outputfile      = 'TMVA.root' ,
                signal_cuts     = ''          ,
                background_cuts = ''          ,
                spectators      = []          ,
                bookingoptions  = "Transformations=I;D;P;G,D" , 
                configuration   = "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" ) :
        """Train the TMVA:
        >>> trainer.train ( var_list , ... ) 
        """
        #
        ## get the logger
        #
        import os 
        name = os.path.split    ( outputfile )[-1] 
        name = os.path.splitext ( name       )[ 0]

        outFile = ROOT.TFile.Open   ( outputfile, 'RECREATE' )
        logger.info ( 'Trainer(%s): output ROOT file: %s ' % ( name , outputfile  ) )

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
            "TMVAClassification"  ,
            outFile               ,
            bookingoptions        )
        logger.info ( 'Trainer(%s): book TMVA-factory %s ' % ( name , bookingoptions ) )
        
        factory.SetVerbose(self.verbose)
        #
        for v in var_list :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' ) 
            factory.AddVariable  ( *vv )
            
        for v in spectators :
            vv = v
            if isinstance ( vv , str ) : vv = ( vv , 'F' )             
            factory.AddSpectator ( *vv )
        #
        signalWeight     = 1.0
        backgroundWeight = 1.0
        #
        logger.info ( 'Trainer(%s): Signal     cuts: "%s" ' % ( name ,     signal_cuts ) ) 
        logger.info ( 'Trainer(%s): Background cuts: "%s" ' % ( name , background_cuts ) )
        # 
        factory.AddTree ( signal     , 'Signal'     ,     signalWeight ,
                          ROOT.TCut (      signal_cuts ) )
        factory.AddTree ( background , 'Background' , backgroundWeight ,
                          ROOT.TCut (  background_cuts ) )
        #
        logger.info ( 'Trainer(%s): Configuration  : "%s" ' % ( name , configuration ) )
        factory.PrepareTrainingAndTestTree(
            ROOT.TCut ( signal_cuts     ) ,
            ROOT.TCut ( background_cuts ) ,
            configuration                 )
        #
        for m in self.methods :
            factory.BookMethod( *m )

        # Train MVAs
        logger.info  ( "Trainer(%s): Train    all Methods" % name )  
        factory.TrainAllMethods    ()
        # Test MVAs
        logger.info  ( "Trainer(%s): Test     all Methods" % name )  
        factory.TestAllMethods     ()
        # Evaluate MVAs
        logger.info  ( "Trainer(%s): Evaluate all Methods" % name )  
        factory.EvaluateAllMethods ()
        # Save the output.
        logger.info  ( "Trainer(%s): Output ROOT file %s is closed" % ( name , outputfile ) )  
        outFile.Close()
        # Now move the weights files
        import glob, os, shutil 
        for weightFile in glob.glob("weights/TMVAClassification*.xml"):
            typeWeight = os.path.splitext(weightFile)[0].split("_")[-1].replace('.weights', '')
            outWeight  = os.path.splitext(outputfile)[0] + '_%s_weights.xml' % typeWeight
            shutil.move   ( weightFile, outWeight)
            logger.info   ( "Trainer(%s): Weights file is created: %s" % ( name , outWeight ) ) 
        shutil.rmtree ( 'weights')

        ## return the name of weigth file w
        return outWeight 

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
    def __init__ ( self         ,
                   name         , 
                   variables    ,
                   weights_file ) :
        
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
            
        self.reader.BookMVA( self.name , weights_file )
        logger.info ('TMVA Reader(%s) is booked: %s ' % ( self.name , weights_file ) ) 

    # =========================================================================
    ## evaluate TMVA
    #  @attention it is *not* CPU efficient
    #  Ugly trick with arrays is needed due to some technicla problems
    #  (actually TMVA reader needs the address of ``float''(in C++ sense) variable
    def __call__ ( self , entry ) :
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
        return self.reader.EvaluateMVA( self.name ) 

# =============================================================================
## start TMVA gui 
def tmvaGUI ( filename , new_canvas = True ) :
    """ Start TMVA-GUI
    """
    ## ROOT.gROOT.LoadMacro('TMVAGui.C')
    if new_canvas :
        from ostap.plotting.canvas import getCanvas
        _c = getCanvas ('glTMVA' , 'TMVA' , 1000 , 800 )
        
    ## start GUI
    return ROOT.TMVA.TMVAGui( filename )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
