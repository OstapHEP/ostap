#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_tmva3.py
#
#  Test for TMVA machinery
# 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-10-25 
# =============================================================================
"""Test for TVMA machinery in  Ostap
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-10-26"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   ostap.core.core          import ROOTCWD,SE
from   ostap.stats.counters     import table_counters   
from   ostap.utils.timing       import timing
from   ostap.utils.progress_bar import progress_bar 
from   ostap.utils.cleanup      import CleanUp
from   ostap.tools.tmva         import Reader, addTMVAResponse
from   ostap.utils.root_utils   import batch_env 
import ostap.io.root_file 
import ROOT, array, os, random 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_tmva3' )
else : 
    logger = getLogger ( __name__ )
# ==============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

# ==============================================================================
## prepare traing and testing data for TMVA 
def prepare_data ( nB = 10000 , nS = 10000 ) :
    """ Prepare training and testing data for TMVA"""

    logger = getLogger ( 'test_tmva2:prepare_data' )
    data_file = CleanUp.tempfile ( suffix = '.root' , prefix = 'ostap-test-tools-tmva2-' )
    
    logger.info('Prepare input ROOT file with data  %s' % data_file )
    import ostap.io.root_file 
    with ROOT.TFile ( data_file ,'recreate') as test_file:
        
        treeSignal = ROOT.TTree('S','signal     tree')
        treeBkg    = ROOT.TTree('B','background tree')
        treeSignal.SetDirectory ( test_file ) 
        treeBkg   .SetDirectory ( test_file ) 
        
        var1 = array.array ( 'd', [0] )
        var2 = array.array ( 'd', [0] )
        var3 = array.array ( 'd', [0] )
        var4 = array.array ( 'd', [0] )
        var5 = array.array ( 'd', [0] )
        
        treeSignal.Branch ( 'var1'  , var1 , 'var1/D' )
        treeSignal.Branch ( 'var2'  , var2 , 'var2/D' )
        treeSignal.Branch ( 'var3'  , var3 , 'var3/D' )
        
        treeSignal.Branch ( 'VARS1' , var4 , 'VARS1/D' )
        treeSignal.Branch ( 'VARS2' , var5 , 'VARS2/D' )
        
        treeBkg   .Branch ( 'var1'  , var1 , 'var1/D' )
        treeBkg   .Branch ( 'var2'  , var2 , 'var2/D' )
        treeBkg   .Branch ( 'var3'  , var3 , 'var3/D' )
        
        treeBkg   .Branch ( 'VARB1' , var4 , 'VARB1/D' )
        treeBkg   .Branch ( 'VARB2' , var5 , 'VARB2/D' )
        
        ## fill background tuple: 
        for i in range ( nB ) : 
            
            x = random.uniform ( -2.0  , 2.0 )
            y = random.uniform ( -2.0  , 2.0 )
            z = random.gauss   (   .0  , 0.5 )
            w = random.gauss   (  0.25 , 0.5 )
            u = random.gauss   (  0.25 , 0.5 )
            
            var1[0] =  x + 0.1 * y  
            var2[0] =  x - 0.1 * y  
            var3[0] = -x +       z
            var4[0] =  w
            var5[0] =  u
            
            treeBkg.Fill()
            
        ## fill signal tuple: 
        for i in range ( nS ) : 
            
            x = random.gauss  (  0.0  , 0.1 )
            y = random.gauss  (  0.0  , 0.2 )
            z = random.gauss  (  0.5  , 0.5 )
            w = random.gauss  ( -0.25 , 0.5 )
            u = random.gauss  ( -0.25 , 0.5 )
            
            var1[0] =  x
            var2[0] =  y  
            var3[0] =  z 
            var4[0] =  w
            var5[0] =  u
            
            treeSignal.Fill()
            
        test_file.Write()
        test_file.ls()
        treeSignal = None
        treeBkg    = None
        
    return data_file

# =============================================================================
## Run TMVA test 
def test_tmva3 () :
    """ Run TMVA test """

    logger = getLogger ( 'test_tmva3' )

    nB = 10000
    nS = 10000
    
    data_file = prepare_data ( nB , nS ) 
    
    logger.info('Create and train TMVA')
    
    with ROOT.TFile.Open( data_file ,'READ') as datafile : 
        datafile.ls()
        tSignal  = datafile['S']
        tBkg     = datafile['B']
        
        #
        ## book TMVA trainer
        #
        from ostap.tools.tmva import Trainer 
        trainer = Trainer (
            name    = 'TestTMVA2' ,   
            methods = [ # type               name   configuration
            ( ROOT.TMVA.Types.kMLP        , "MLP1"         , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=500:HiddenLayers=N+1:TestRate=5:!UseRegulator" ) ,                
            ( ROOT.TMVA.Types.kMLP        , "MLP2"         , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=500:HiddenLayers=N+2:TestRate=5:!UseRegulator" ) ,
            ( ROOT.TMVA.Types.kMLP        , "MLP3"         , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=500:HiddenLayers=N+3:TestRate=5:!UseRegulator" ) ,
            ( ROOT.TMVA.Types.kMLP        , "MLP4"         , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=500:HiddenLayers=N+4:TestRate=5:!UseRegulator" ) ,
            ( ROOT.TMVA.Types.kMLP        , "MLP5"         , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=500:HiddenLayers=N+5:TestRate=5:!UseRegulator" ) ,
            ##
            ( ROOT.TMVA.Types.kBDT        , "BDTGD2"       , "H:!V:NTrees=1000:VarTransform=G,D,G,D:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=100:MaxDepth=2" ) ,
            ( ROOT.TMVA.Types.kBDT        , "BDTGD3"       , "H:!V:NTrees=1000:VarTransform=G,D,G,D:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=100:MaxDepth=3" ) ,
            ( ROOT.TMVA.Types.kBDT        , "BDTGD4"       , "H:!V:NTrees=1000:VarTransform=G,D,G,D:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=100:MaxDepth=4" ) ,
            ## 
            ( ROOT.TMVA.Types.kBDT        , "BDTB"        , "H:!V:NTrees=1000:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" )  , 
            ( ROOT.TMVA.Types.kBDT        , "BDTD"        , "H:!V:NTrees=1000:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" ) ,
            ## 
            ( ROOT.TMVA.Types.kCuts       , "Cuts"        , "H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" ) ,
            ##
            ( ROOT.TMVA.Types.kFisher     , "Fisher"      , "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" ),
            ( ROOT.TMVA.Types.kFisher     , "FisherG"     , "H:!V:VarTransform=Gauss"  ),
            ( ROOT.TMVA.Types.kFisher     , "FisherB"     , "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" ),
            ##
            ( ROOT.TMVA.Types.kSVM        , "SVM"         , "H:!V:Gamma=0.25:Tol=0.001:VarTransform=G" ) ,
            ##
            ( ROOT.TMVA.Types.kLikelihood , "Likelihood"  , "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" ) ,
            ( ROOT.TMVA.Types.kLikelihood , "LikelihoodD" , "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" ), 
            ( ROOT.TMVA.Types.kHMatrix    , "HMatrix"     , "H:!V:VarTransform=None" ) ,
            ( ROOT.TMVA.Types.kRuleFit    , "RuleFit"     , "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" ),
            ( ROOT.TMVA.Types.kPDERS      , "PDERS"       , "H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" ) ,
            ( ROOT.TMVA.Types.kKNN        , "KNN"         , "H:!V:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" ) ,
            ] ,
            #
            ## common variablss:
            variables       = [ 'var2' , 'var1' ,  'var3' ] , ## Variables for training
            signal_vars     = { 'VAR1' : 'VARS1' , 'VAR2' : 'VARS2' } , 
            background_vars = { 'VAR1' : 'VARB1' , 'VAR2' : 'VARB2' } , 
            ##                         
            signal         = tSignal                  , ## `Signal' sample
            background     = tBkg                     , ## `Background' sample         
            verbose        = True                     , 
            workdir        = CleanUp.tempdir ( prefix = 'ostap-tmva2-workdir-' ) ) ##  working directory 
        
        with timing ( 'for TMVA training' , logger ) : 
            weights_files = trainer.train ()            
            tar_file      = trainer.tar_file
            trainer_name  = trainer.name
            tmva_output   = trainer.output_file

    # =============================================================================
    ## Use TMVA for classification 
    # =============================================================================
    logger.info('Five ways to use TMVA for classification')

    signal_inputs      = ( 'var1' ,  'var2' , 'var3' , 'VAR1 : VARS1' , 'VAR2 : VARS2' ) 
    background_inputs  = ( 'var1' ,  'var2' , 'var3' , 'VAR1 : VARB1' , 'VAR2 : VARB2' ) 
                              
    # =============================================================================
    ## (1) the most efficient way: add TMVA decisions directly into TTree (fast)
    # =============================================================================
    logger.info ( '(1) Add TMVA decision directly into existing TTree (fast)' ) 
    with timing ( "Add TMVA response to signal TTree" , logger =logger ) : 
        tSignal = ROOT.TChain ( 'S' ) ;  tSignal.Add ( data_file )
        addTMVAResponse ( tSignal ,
                          inputs        = signal_inputs , 
                          weights_files = tar_file      ,
                          prefix        = 'tmva_'       ,
                          suffix        = '_response'   )
        
    with timing ( "Add TMVA response to background TTree" , logger =logger ) :     
        tBkg    = ROOT.TChain ( 'B' ) ;  tBkg.Add    ( data_file )
        addTMVAResponse ( tBkg    ,
                          inputs        = background_inputs , 
                          weights_files = tar_file    ,
                          prefix        = 'tmva_'     ,
                          suffix        = '_response' )
        
# =============================================================================
if '__main__' == __name__ :

    with timing ( 'test_tmva2', logger = logger ) :
        test_tmva3 () 

# =============================================================================
#                                                                       The END
# =============================================================================    
