#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_tmva.py
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
from   ostap.stats.counters     import counters_table  
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
    logger = getLogger ( 'ostap.test_tools_tmva2' )
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
        
        treeSignal.Branch ( 'var1' , var1 , 'var1/D' )
        treeSignal.Branch ( 'var2' , var2 , 'var2/D' )
        treeSignal.Branch ( 'var3' , var3 , 'var3/D' )
        
        treeBkg   .Branch ( 'var1' , var1 , 'var1/D' )
        treeBkg   .Branch ( 'var2' , var2 , 'var2/D' )
        treeBkg   .Branch ( 'var3' , var3 , 'var3/D' )
        
        ## fill background tuple: 
        for i in range ( nB ) : 
            
            x = random.uniform ( -2.0 , 2.0 )
            y = random.uniform ( -2.0 , 2.0 )
            z = random.gauss   (   .0 , 0.5 )
            
            var1[0] =  x + 0.1 * y  
            var2[0] =  x - 0.1 * y  
            var3[0] = -x +       z
            
            treeBkg.Fill()
            
        ## fill signal tuple: 
        for i in range ( nS ) : 
            
            x = random.gauss  (  0.0 , 0.1 )
            y = random.gauss  (  0.0 , 0.2 )
            z = random.gauss  (  0.5 , 0.5 )
            
            var1[0] =  x
            var2[0] =  y  
            var3[0] =  z 
            treeSignal.Fill()
            
        test_file.Write()
        test_file.ls()
        treeSignal = None
        treeBkg    = None
        
    return data_file

# =============================================================================
## Run TMVA test 
def test_tmva2() :
    """ Run TMVA test """

    logger = getLogger ( 'test_tmva2' )

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
            ( ROOT.TMVA.Types.kMLP        , "MLP"         , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=200:HiddenLayers=N+3:TestRate=5:!UseRegulator" ) ,
            ( ROOT.TMVA.Types.kBDT        , "BDTG"        , "H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=100:MaxDepth=2" ) , 
            ( ROOT.TMVA.Types.kBDT        , "BDTB"        , "H:!V:NTrees=1000:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" )  , 
            ( ROOT.TMVA.Types.kBDT        , "BDTD"        , "H:!V:NTrees=1000:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" ) ,        
            ( ROOT.TMVA.Types.kCuts       , "Cuts"        , "H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" ) ,
            ( ROOT.TMVA.Types.kFisher     , "Fisher"      , "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" ),
            ( ROOT.TMVA.Types.kFisher     , "FisherG"     , "H:!V:VarTransform=Gauss"  ),
            ( ROOT.TMVA.Types.kFisher     , "FisherB"     , "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" ),
            ( ROOT.TMVA.Types.kSVM        , "SVM"         , "H:!V:Gamma=0.25:Tol=0.001:VarTransform=Norm" ) ,
            ( ROOT.TMVA.Types.kLikelihood , "Likelihood"  , "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" ) ,
            ( ROOT.TMVA.Types.kLikelihood , "LikelihoodD" , "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" ), 
            ( ROOT.TMVA.Types.kHMatrix    , "HMatrix"     , "H:!V:VarTransform=None" ) ,
            ( ROOT.TMVA.Types.kRuleFit    , "RuleFit"     , "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" ),
            ## ( ROOT.TMVA.Types.kPDERS      , "PDERS"       , "H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" ) ,
            ## ( ROOT.TMVA.Types.kKNN        , "KNN"         , "H:!V:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" ) ,
            ] ,
            variables = [ 'var1' , 'var2' ,  'var3' ] , ## Variables for training 
            signal         = tSignal                  , ## ``Signal'' sample
            background     = tBkg                     , ## ``Background'' sample         
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

    # =============================================================================
    ## (1) the most efficient way: add TMVA decisions directly into TTree (fast)
    # =============================================================================
    logger.info ( '(1) Add TMVA decision directly into existing TTree (fast)' ) 
    with timing ( "Add TMVA response to signal TTree" , logger =logger ) : 
        tSignal = ROOT.TChain ( 'S' ) ;  tSignal.Add ( data_file )
        addTMVAResponse ( tSignal ,
                          inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                          weights_files = tar_file ,
                          prefix        = 'tmva_'     ,
                          suffix        = '_response' )
        
    with timing ( "Add TMVA response to background TTree" , logger =logger ) :     
        tBkg    = ROOT.TChain ( 'B' ) ;  tBkg.Add    ( data_file )
        addTMVAResponse ( tBkg    ,
                          inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                          weights_files = tar_file ,
                          prefix        = 'tmva_'     ,
                          suffix        = '_response' )
        
    # ===============================================================================
    ## (2) Add TMVA decision during TTree -> RooDataSet transformation (can be slow)
    # ===============================================================================
    logger.info ( '(2) Add TMVA decision during TTree -> RooDataSet transformation (can be slow)')
    variables = [ 'var1' , 'var2' , 'var3' ]
    
    ## (2.1) Create TMVA reader
    reader = Reader ( 'MyMLP' ,
                      variables     = [ ('var1' , lambda s : s.var1 )   ,
                                        ('var2' , lambda s : s.var2 )   ,
                                        ('var3' , lambda s : s.var3 ) ] ,
                      weights_files = tar_file   )
    
    from ostap.fitting.pyselectors import Variable     
    for m in reader.methods[:] :
        variables += [ Variable( 'tmva1_%s_response' % m , 'TMVA(%s)' % m , accessor = reader[m] ) ]

    with timing ( "Add TMVA response during TTree->RooDatSet transforomation (signal)"     , logger =logger ) : 
        tSignal   = ROOT.TChain ( 'S' ) ;  tSignal.Add ( data_file )
        ds_S1, _  = tSignal.fill_dataset ( variables )
    with timing ( "Add TMVA response during TTree->RooDatSet transforomation (background)" , logger =logger ) : 
        tBkg      = ROOT.TChain ( 'B' ) ;  tBkg.Add    ( data_file )
        ds_B1, _  = tBkg   .fill_dataset ( variables )
        
    logger.info ( 'Created signal     dataset\n%s' %  ds_S1.table ( prefix = '# ' ) )
    logger.info ( 'Created background dataset\n%s' %  ds_B1.table ( prefix = '# ' ) )

    # ===============================================================================
    ## (3) Add TMVA decision directly into existing RooDataSet (fast)  
    # ===============================================================================
    logger.info ( '(3) Add TMVA decision directly into existing RooDataSet (fast)' ) 
    with timing ( "Add TMVA response to signal RooDataSet" , logger =logger ) : 
        addTMVAResponse ( ds_S1  ,
                          inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                          weights_files = tar_file ,
                          prefix        = 'tmva2_'     ,
                          suffix        = '_response' )
        
    with timing ( "Add TMVA response to background TTree" , logger =logger ) :     
        addTMVAResponse ( ds_B1   ,
                          inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                          weights_files = tar_file ,
                          prefix        = 'tmva2_'    ,
                          suffix        = '_response' )
        
    logger.info ( 'Updated signal     dataset\n%s' %  ds_S1.table ( prefix = '# ' ) )
    logger.info ( 'Updated background dataset\n%s' %  ds_B1.table ( prefix = '# ' ) )

    # ===============================================================================
    ## (4) Calcuate TMVA decision on-fly via the explict loop over TTree entries (slow)
    # ===============================================================================
    logger.info ( '(4) Calcuate TMVA decision on-fly via the explict loop over TTree entries (can be slow)')
    with timing ( "TMVA response via explicit loop over signal TTree" , logger =logger ) :
        tSignal   = ROOT.TChain ( 'S' ) ;  tSignal.Add ( data_file )
        counters  = {}
        methods   = reader.methods[:] 
        for m in methods : counters[m] = SE() 
        for evt in tSignal :
            for method in methods : counters[method] += reader ( method , evt )
        title = 'Signal     response (TTree)'
        table = counters_table ( counters , title = title , prefix = '# ' )
        logger.info ( '%s\n%s' % ( title , table ) )
        
    with timing ( "TMVA response via explicit loop over background TTree" , logger =logger ) :
        tBkg      = ROOT.TChain ( 'B' ) ;  tBkg.Add    ( data_file )
        counters  = {}
        methods   = reader.methods[:] 
        for m in methods : counters[m] = SE() 
        for evt in tBkg :
            for method in methods : counters[method] += reader ( method , evt )
        title = 'Background response (TTree)'
        table = counters_table ( counters , title = title , prefix = '# ' )
        logger.info ( '%s\n%s' % ( title , table ) )

    # ===============================================================================
    ## (5) Calcuate TMVA decision on-fly via the explict loop over RooDataSet entries (slow)
    # ===============================================================================
    logger.info ( '(5) Calcuate TMVA decision on-fly via the explict loop over RooDataSet entries (can be slow)')
    with timing ( "TMVA response via explicit loop over signal     RooDataSet" , logger =logger ) :
        counters  = {}
        methods   = reader.methods[:] 
        for m in methods : counters[m] = SE() 
        for evt, _ in ds_S1 :
            for method in methods : counters[method] += reader ( method , evt )
        title = 'Signal     response (RooDataSet)'
        table = counters_table ( counters , title = title , prefix = '# ' )
        logger.info ( '%s\n%s' % ( title , table ) )
            
    with timing ( "TMVA response via explicit loop over background RooDataSet" , logger = logger ) :        
        counters  = {}
        methods   = reader.methods[:] 
        for m in methods : counters[m] = SE() 
        for evt, _ in ds_B1 :
            for method in methods : counters[method] += reader ( method , evt )
        title = 'Background response (RooDataSet)'
        table = counters_table ( counters , title = title , prefix = '# ' )
        logger.info ( '%s\n%s' % ( title , table ) )
    
# =============================================================================
if '__main__' == __name__ :

    with timing ( 'test_tmva2', logger = logger ) :
        test_tmva2 () 

# =============================================================================
#                                                                       The END
# =============================================================================    
