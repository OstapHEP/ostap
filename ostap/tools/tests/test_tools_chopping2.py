#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_chopping2.py
#  Test for TMVA `chopping' (k-fold cross-validation) machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-10-25 
# =============================================================================
""" Test for TVMA-chopping (k-fold cross-validation) machinery in  Ostap
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-10-26"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   ostap.core.core          import ROOTCWD
from   ostap.utils.progress_bar import progress_bar 
from   ostap.utils.timing       import timing
from   ostap.utils.root_utils   import batch_env 
from   ostap.utils.cleanup      import CleanUp
import ostap.io.root_file
import ROOT, os, array, random  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_chopping2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

data_file = CleanUp.tempfile ( suffix = '.root' , prefix = 'ostap-test-tools-chopping-' ) 

if not os.path.exists( data_file ) :
    
    nB = 20000
    nS = 10000

    s_evt_per_run = 927
    b_evt_per_run = 511
    
    logger.info('Prepare input ROOT file with data %s' % data_file )
    with ROOT.TFile( data_file ,'recreate') as test_file:
        
        treeSignal = ROOT.TTree('S','signal     tree')
        treeBkg    = ROOT.TTree('B','background tree')
        treeSignal.SetDirectory ( test_file ) 
        treeBkg   .SetDirectory ( test_file ) 
        
        var1 = array.array ( 'd',  [ 0 ] )
        var2 = array.array ( 'd',  [ 0 ] )
        var3 = array.array ( 'd',  [ 0 ] )
        var4 = array.array ( 'd',  [ 0 ] )
        vevt = array.array ( 'i',  [ 0 ] )
        vrun = array.array ( 'i',  [ 0 ] )
        
        treeSignal.Branch ( 'var1' , var1 , 'var1/D' )
        treeSignal.Branch ( 'var2' , var2 , 'var2/D' )
        treeSignal.Branch ( 'var3' , var3 , 'var3/D' )
        treeSignal.Branch ( 'VARS' , var4 , 'VARS/D' )
        
        treeSignal.Branch ( 'evt'  , vevt , 'evt/I'  )
        treeSignal.Branch ( 'run'  , vrun , 'run/I'  )
        
        treeBkg   .Branch ( 'var1' , var1 , 'var1/D' )
        treeBkg   .Branch ( 'var2' , var2 , 'var2/D' )
        treeBkg   .Branch ( 'var3' , var3 , 'var3/D' )
        treeBkg   .Branch ( 'VARB' , var4 , 'VARB/D' )
        treeBkg   .Branch ( 'evt'  , vevt , 'evt/I'  )
        treeBkg   .Branch ( 'run'  , vrun , 'run/I'  )
        
        ievt = 0
        irun = 1 
        ## fill background tuple: 
        #for i in progress_bar ( range ( nB ) ) : 
        for i in range ( nB ) : 
            
            x = random.uniform ( -2.0 , 2.0 )
            y = random.uniform ( -2.0 , 2.0 )
            z = random.gauss   (   .0 , 0.5 )
            w = random.gauss   (   .2 , 0.5 )
            
            var1[0] =  x + 0.1 * y  
            var2[0] =  x - 0.1 * y  
            var3[0] = -x +       z
            var4[0] =  w
            
            ievt += 1
            if 0 ==  ( ievt % b_evt_per_run ) :
                irun += 1
                ievt  = 1 

            vevt[0] = ievt
            vrun[0] = irun

            treeBkg.Fill()
            
        ievt = 0
        irun = 1 
        ## fill signal tuple: 
        #for i in progress_bar ( range ( nS ) ) : 
        for i in range ( nS ) : 
            
            x = random.gauss  (  0.0 , 0.1 )
            y = random.gauss  (  0.0 , 0.2 )
            z = random.gauss  (  0.5 , 0.5 )
            w = random.gauss  ( -0.2 , 0.5 )
            
            var1[0] =  x
            var2[0] =  y  
            var3[0] =  z
            var4[0] =  w
            
            ievt += 1
            if 0 == ( ievt % s_evt_per_run ) :
                irun += 1
                ievt  = 1 

            vevt[0] = ievt
            vrun[0] = irun
            
            treeSignal.Fill()
             
            
        test_file.Write()
        test_file.ls()
        treeSignal = None
        treeBkg    = None

# ===========================================================================
##   number of    categories 
N  = 7
logger.info('Create and train TMVA')

# ============================================================================
## Train TMVA
# ============================================================================
cSignal = ROOT.TChain ( 'S' ) ; cSignal.Add ( data_file )
cBkg    = ROOT.TChain ( 'B' ) ; cBkg   .Add ( data_file )
# 
## book TMVA trainer
from ostap.tools.chopping import Trainer 
trainer = Trainer (
    N        = N                 , ## ATTENTION! 
    category = "137*evt+813*run" , ## ATTENTION!
    ## other  arguments as for ``plain'' TMVA
    name    = 'TestChopping' ,   
    methods = [ # type               name   configuration
    ( ROOT.TMVA.Types.kMLP        , "MLP1"       , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=400:HiddenLayers=N+1:TestRate=5:!UseRegulator" ) ,
    ( ROOT.TMVA.Types.kMLP        , "MLP5"       , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=400:HiddenLayers=N+5:TestRate=5:!UseRegulator" ) ,
    ( ROOT.TMVA.Types.kBDT        , "BDTG2"      , "H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" ) , 
    ( ROOT.TMVA.Types.kBDT        , "BDTG4"      , "H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=4" ) , 
    ( ROOT.TMVA.Types.kCuts       , "Cuts"       , "H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" ) ,
    ( ROOT.TMVA.Types.kFisher     , "Fisher"     , "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" ),
    ( ROOT.TMVA.Types.kLikelihood , "Likelihood" , "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )
    ] ,
    ## 
    variables       = [ 'var1' , 'var2' ,  'var3' ] , ## Variables for training
    signal_vars     = { 'VARX' : 'VARS' } , 
    background_vars = { 'VARX' : 'VARB' } ,      
    ##     
    signal         = cSignal                  , ## `Signal' sample
    background     = cBkg                     , ## `Background' sample         
    verbose        = True                     ,
    make_plots     = True                     ,   
    logging        = True                     ,  ## produce  log-files 
    parallel       = True                     ,  ## parallel training
    ## 
    prefilter      = 'var1>-1.8'  ,
    ##
    chop_signal     = True ,
    chop_background = True ,
    ## 
    signal_train_fraction     = 0.85 , 
    background_train_fraction = 0.85 ,
    ##
    workdir        = CleanUp.tempdir ( prefix = 'ostap-chopping-workdir-' ) , ##  working directory 
    ## parallel_conf  = { 'ncpus' : 0 , 'ppservers' : 'auto' }
    )

# train it!  
with timing ( 'for TMVA/Chopping training' , logger ) :
    trainer.train () 
    tar_file      = trainer.tar_file
    
# =============================================================================
# remove unnesessary output files
for f in trainer.output_files :
    if os.path.exists ( f ) and os.path.isfile( f ) :
        try    : os.remove ( f )
        except : pass

# =============================================================================
## Use trained TMVA/Chopping
#  There are four alternatives
#  - add TMVA/Chopping decision to (input) TTrees                        (fast)       
#  - add TMVA/Chopping decision to (existing) RooDataSet      (moderately fast)
#  - add TMVA/Chopping decision during creation of new RooDataSet        (slow) 
#  - use MVA/Chopper Reader directly                              (rather slow)
#    but it is very flexible and powerful with respect to variable transformations
# =============================================================================

# =============================================================================
## A) Add TMVA/Chopping decision to (input) TTrees 
# =============================================================================
with timing ( "Add TMVA/Chopping response to input TTree" , logger = logger ) as time_A :
    
    cSignal = ROOT.TChain ( 'S' ) ; cSignal.Add ( data_file )
    cBkg    = ROOT.TChain ( 'B' ) ; cBkg   .Add ( data_file )
    
    signal_inputs      = ( 'var1' ,  'var2' , 'var3' , 'VARX : VARS' ) 
    background_inputs  = ( 'var1' ,  'var2' , 'var3' , 'VARX : VARB' )   
       
    config = { 'chopping'      : "137*evt+813*run"             ,
               'N'             : N                             , 
               'weights_files' : tar_file                      ,
               'prefix'        : 'tmva_'                       ,
               'suffix'        : '_response'                   }
    
    from ostap.tools.chopping import addChoppingResponse
    
    addChoppingResponse ( cSignal  , inputs = signal_inputs     , **config ) 
    addChoppingResponse ( cBkg     , inputs = background_inputs , **config )

    cSignal = ROOT.TChain ( 'S' ) ; cSignal.Add ( data_file )
    cBkg    = ROOT.TChain ( 'B' ) ; cBkg   .Add ( data_file )

    logger.info ('TTree   SIG (with TMVA decisions):\n%s' % cSignal.table ( prefix = '# ') ) 
    logger.info ('TTree   BKG (with TMVA decisions):\n%s' % cBkg   .table ( prefix = '# ') ) 

                          
# =============================================================================
##                                                                      The END
# =============================================================================    
