#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_chopping.py
#  Test for TMVA ``chopping''(k-fold cross-validation) machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-10-25 
# =============================================================================
"""Test for TVMA-chopping (k-fold cross-validation) machinery in  Ostap
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-10-26"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   ostap.core.meta_info     import root_info
from   ostap.core.core          import ROOTCWD
from   ostap.utils.progress_bar import progress_bar 
from   ostap.utils.timing       import timing
from   ostap.utils.utils        import batch_env 
from ostap.utils.cleanup        import CleanUp
import ostap.io.root_file
import ROOT, os, array, random  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_chopping' )
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

    if root_info <  ( 6 , 15 ) : 
        nB = 2000
        nS = 1000

    s_evt_per_run = 927
    b_evt_per_run = 511
    
    logger.info('Prepare input ROOT file with data %s' % data_file )
    with ROOT.TFile( data_file ,'recreate') as test_file:
        
        treeSignal = ROOT.TTree('S','signal     tree')
        treeBkg    = ROOT.TTree('B','background tree')
        treeSignal.SetDirectory ( test_file ) 
        treeBkg   .SetDirectory ( test_file ) 
        
        var1 = array.array ( 'd', [0] )
        var2 = array.array ( 'd', [0] )
        var3 = array.array ( 'd', [0] )
        vevt = array.array ( 'i', [0] )
        vrun = array.array ( 'i', [0] )
        
        treeSignal.Branch ( 'var1' , var1 , 'var1/D' )
        treeSignal.Branch ( 'var2' , var2 , 'var2/D' )
        treeSignal.Branch ( 'var3' , var3 , 'var3/D' )
        treeSignal.Branch ( 'evt'  , vevt , 'evt/I'  )
        treeSignal.Branch ( 'run'  , vrun , 'run/I'  )
        
        treeBkg   .Branch ( 'var1' , var1 , 'var1/D' )
        treeBkg   .Branch ( 'var2' , var2 , 'var2/D' )
        treeBkg   .Branch ( 'var3' , var3 , 'var3/D' )
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
            
            var1[0] =  x + 0.1 * y  
            var2[0] =  x - 0.1 * y  
            var3[0] = -x +       z
            
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
            
            var1[0] =  x
            var2[0] =  y  
            var3[0] =  z
            
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
    variables = [ 'var1' , 'var2' ,  'var3' ] , ## Variables for training 
    signal         = cSignal                  , ## ``Signal'' sample
    background     = cBkg                     , ## ``Background'' sample         
    verbose        = True     ,
    make_plots     = True     ,   
    logging        = True     ,  ## produce  log-files 
    parallel       = True     ,  ## parallel training
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
    
    config = { 'chopper'       : "137*evt+813*run"             ,
               'N'             : N                             , 
               'inputs'        : ( 'var1' ,  'var2' , 'var3' ) ,
               'weights_files' : tar_file                      ,
               'prefix'        : 'tmva_'                       ,
               'suffix'        : '_response'                   }
    
    from ostap.tools.chopping import addChoppingResponse
    
    addChoppingResponse ( cSignal  , **config ) 
    addChoppingResponse ( cBkg     , **config )

    cSignal = ROOT.TChain ( 'S' ) ; cSignal.Add ( data_file )
    cBkg    = ROOT.TChain ( 'B' ) ; cBkg   .Add ( data_file )

    logger.info ('TTree   SIG (with TMVA decisions):\n%s' % cSignal.table ( prefix = '# ') ) 
    logger.info ('TTree   BKG (with TMVA decisions):\n%s' % cBkg   .table ( prefix = '# ') ) 

# =============================================================================
## B) Add TMVA/Choppnig desision to (existing) RooDataSet 
# =============================================================================  
## B.1) prepare "existing" datasets
if True : 
    ## category function
    category = lambda s :  int ( s.evt*137 + 813*s.run ) % N
    
    ## prepare dataset with TMVA/Chopping result
    
    from ostap.fitting.pyselectors import SelectorWithVars, Variable
    
    variables = [
        Variable ( 'var1' , 'variable#1' ) ,
        Variable ( 'var2' , 'variable#2' ) ,
        Variable ( 'var3' , 'variable#3' ) ,
        ## extra: needed for addChoppingResponse 
        Variable ( 'evt'  , 'event'      ) ,
        Variable ( 'run'  , 'run'        ) ,
        ## extra: needed for cross-checks  
        Variable ( 'cat'  , 'category'   , accessor = category ) ,
        ## Variable ( 'cat'  , 'category'   , accessor = '((137*evt+813*run)%%%d)*1.0' % N ) ,
        ]
    
    cSignal = ROOT.TChain ( 'S' ) ; cSignal.Add ( data_file )
    cBkg    = ROOT.TChain ( 'B' ) ; cBkg   .Add ( data_file )
    
    dsS1, _ = cSignal.fill_dataset ( variables , selection = 'var1<100' , use_frame = -1 )
    dsB1, _ = cBkg   .fill_dataset ( variables , selection = 'var1<100' , use_frame = -1 )
# =============================================================================
## B.2) add response 
with timing ( "Add TMVA/Chopping response to (existing) RooDataSet " , logger = logger ) as time_B :
    
    config = { 'chopper'       : "137*evt+813*run"             ,
               'N'             : N                             , 
               'inputs'        : ( 'var1' ,  'var2' , 'var3' ) ,
               'weights_files' : tar_file                      ,
               'prefix'        : 'tmva_'                       ,
               'suffix'        : '_response'                   }
    
    addChoppingResponse ( dsS1 , **config ) 
    addChoppingResponse ( dsB1 , **config )

    logger.info ('dataset SIG (with TMVA decisions):\n%s' % dsS1.table ( prefix = '# ') ) 
    logger.info ('dataset BKG (with TMVA decisions):\n%s' % dsB1.table ( prefix = '# ') ) 

# =============================================================================
## C) Add TMVA/Choppnig desision to (newly created) RooDataSet 
# =============================================================================  
with timing ( "Add TMVA/Chopping response to (newly created) RooDataSet " , logger = logger ) as time_C :
    
    from ostap.tools.chopping import Reader
    reader = Reader (
        N             = N         , ##  number of   categories
        categoryfunc  = category  , ## category 
        ## other argument  as for plain TMVA     
        name          = 'ChopReader' ,
        variables     = [ ('var1' , lambda s : s.var1 )   ,
                          ('var2' , lambda s : s.var2 )   ,
                          ('var3' , lambda s : s.var3 ) ] ,
        weights_files = tar_file )

    extended_vars = list ( variables )

    for m in reader.methods :
        extended_vars += [ Variable ( 'tmva_%s_response' % m , 'TMVA response (%s)' % m , accessor = reader[m] ) ]
                
    cSignal = ROOT.TChain ( 'S' ) ; cSignal.Add ( data_file )
    cBkg    = ROOT.TChain ( 'B' ) ; cBkg   .Add ( data_file )
    
    dsS2, _ = cSignal.fill_dataset ( extended_vars , selection = 'var1<100' , use_frame = -1 )
    dsB2, _ = cBkg   .fill_dataset ( extended_vars , selection = 'var1<100' , use_frame = -1 )
    
    logger.info ('dataset SIG (with TMVA decisions):\n%s' % dsS2.table ( prefix = '# ') ) 
    logger.info ('dataset BKG (with TMVA decisions):\n%s' % dsB2.table ( prefix = '# ') ) 

    decisions = [ 'tmva_%s_response' % m for m in reader.methods ]

    del extended_vars
    del reader 

# =============================================================================
## Check TMVA results 
# =============================================================================
cSignal = ROOT.TChain ( 'S' ) ; cSignal.Add ( data_file )
cBkg    = ROOT.TChain ( 'B' ) ; cBkg   .Add ( data_file )
    
s_stats   = cSignal.statVars ( decisions )
b_stats   = cBkg   .statVars ( decisions )
s1_stats  = dsS1   .statVars ( decisions )
b1_stats  = dsB1   .statVars ( decisions )
s2_stats  = dsS2   .statVars ( decisions )
b2_stats  = dsB2   .statVars ( decisions )

table = [ ( '' , 'Method ', 'Signal' , 'Background ' ) ]
for t, s,b in ( ( 'A' , s_stats , b_stats  ) ,
                ( 'B' , s1_stats, b1_stats ) ,
                ( 'C' , s2_stats, b2_stats ) ) : 
    for d in decisions :
        
        ms = s [ d ]
        mb = b [ d ]
        
        row = t , d , \
              '%+7.3f +/- %-7.3f' % ( float ( ms.mean() ) , ms.rms() ) ,\
              '%+7.3f +/- %-7.3f' % ( float ( mb.mean() ) , mb.rms() )
        table.append ( row )
        
import ostap.logger.table as T
table  = T.table (  table , title = 'TMVA performance', prefix = '# ' , alignment = 'llcc' )
logger.info ( 'TMVA performance for Signal and Background\n%s' % table )


rows = [ ( 'Method' , 'Time [s]' ) ] 
row  = 'A' , '%.1f' % time_A.delta
rows.append ( row )
row  = 'B' , '%.1f' % time_B.delta
rows.append ( row )
row  = 'C' , '%.1f' % time_C.delta
rows.append ( row )
title = 'Timing'
table = T.table ( rows , title = title , prefix = '# ' , alignment = 'cc' )
logger.info ( '%s\n%s' % ( title , table ) ) 

    
                          
# =============================================================================
##                                                                      The END
# =============================================================================    
