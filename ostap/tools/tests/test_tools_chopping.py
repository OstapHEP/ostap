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
import ROOT, os
import ostap.io.root_file
from   ostap.core.meta_info     import root_info
from   builtins                 import range
from   ostap.core.core          import ROOTCWD
from   ostap.utils.progress_bar import progress_bar 
from   array                    import array
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_chopping' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
from ostap.utils.cleanup import CleanUp
data_file = CleanUp.tempfile ( suffix = '.root' , prefix = 'ostap-test-tools-chopping-' ) 

if not os.path.exists( data_file ) :
    import random
    
    nB = 20000
    nS = 10000

    if root_info < (6,15) : 
        nB = 2000
        nS = 1000

    s_evt_per_run = 927
    b_evt_per_run = 511
    
    logger.info('Prepare input ROOT file with data %s' % data_file )
    with ROOT.TFile.Open( data_file ,'recreate') as test_file:
        ## test_file.cd()
        treeSignal = ROOT.TTree('S','signal     tree')
        treeBkg    = ROOT.TTree('B','background tree')
        treeSignal.SetDirectory ( test_file ) 
        treeBkg   .SetDirectory ( test_file ) 
        
        from array import array 
        var1 = array ( 'd', [0] )
        var2 = array ( 'd', [0] )
        var3 = array ( 'd', [0] )
        vevt = array ( 'i', [0] )
        vrun = array ( 'i', [0] )
        
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

##   number of    categories 
N  = 7
logger.info('Create and train TMVA')
with ROOT.TFile.Open( data_file ,'READ') as datafile : 
    datafile =  ROOT.TFile.Open( data_file ,'READ')
    datafile.ls()
    tSignal  = datafile['S']
    tBkg     = datafile['B']

    #
    ## book TMVA trainer
    #
    from ostap.tools.chopping import Trainer 
    trainer = Trainer (
        N        = N                 , ## ATTENTION! 
        category = "137*evt+813*run" , ## ATTENTION!
        ## other  arguments as for ``plain'' TMVA
        name    = 'TestChopping' ,   
        methods = [ # type               name   configuration
        ( ROOT.TMVA.Types.kMLP        , "MLP"        , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=400:HiddenLayers=N+5:TestRate=5:!UseRegulator" ) ,
        ## ( ROOT.TMVA.Types.kBDT        , "BDTG"       , "H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" ) , 
        ## ( ROOT.TMVA.Types.kCuts       , "Cuts"       , "H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" ) ,
        ## ( ROOT.TMVA.Types.kFisher     , "Fisher"     , "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" ),
        ## ( ROOT.TMVA.Types.kLikelihood , "Likelihood" , "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )
        ] ,
        variables = [ 'var1' , 'var2' ,  'var3' ] , ## Variables for training 
        signal         = tSignal                  , ## ``Signal'' sample
        background     = tBkg                     , ## ``Background'' sample         
        verbose        = True     ,
        ## make_plots     = False    ,   
        logging        = True     ,  ## produce  log-files 
        parallel       = True     ,  ## parallel training
        prefilter      = 'var1>-1.8'  , 
        workdir        = CleanUp.tempdir ( prefix = 'ostap-chopping-workdir-' ) , ##  working directory 
        ## parallel_conf  = { 'ncpus' : 0 , 'ppservers' : 'auto' }
        )

    from ostap.utils.timing import timing
    
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
#  There are two alternatives
#  - usage of TMVA/Chopper Reader : it can be  rather slow,
#    but it is very flexible and powerful with respect to variable transformations
#  - addChoppingResponse function : it is less flexible, but very CPU efficient 
# =============================================================================

## category function
category = lambda s :  int ( s.evt*137 + 813*s.run ) % N

## prepare dataset with TMVA/Chopping result

from ostap.fitting.pyselectors import SelectorWithVars,  Variable     
## 1) Book RooDataset                 
variables = [
    Variable ( 'var1' , 'variable#1' ) ,
    Variable ( 'var2' , 'variable#2' ) ,
    Variable ( 'var3' , 'variable#3' ) ,
    ## extra: needed for addChoppingResponse 
    Variable ( 'evt'  , 'event'      ) ,
    Variable ( 'run'  , 'run'        ) ,
    ## extra: needed for cross-checks  
    Variable ( 'cat'  , 'category'   , accessor = category ) ,
    ]


## 2) create TMVA/Chopping reader
from ostap.tools.chopping import Reader

# =============================================================================
reader = Reader (
    N             = N         , ##  number of   categories
    categoryfunc  = category  , ## category 
    ## other argument  as for plain TMVA     
    name          = 'ChopReader' ,
    variables     = [ ('var1' , lambda s : s.var1 )   ,
                      ('var2' , lambda s : s.var2 )   ,
                      ('var3' , lambda s : s.var3 ) ] ,
    weights_files = tar_file )


methods = reader.methods[:]

## # =============================================================================
## ## A: Use TMVA/Chopping  reader
## #  - It can be slow, but it allows on-flight variables transformation
## #  - much more efficient alternativeis <code>addChoppingResponse</code> function
## # =============================================================================

## # =============================================================================
## ## 2.1) few trivial tests: use the methods/reader as simple function
## for m in methods :
##     method   = reader[m]
##     ## response = [ method ( i  , 1.1 , 0.8 , 0.3 ) for i in  range ( reader.N ) ] 
##     response = method.stat ( 1.1 , 0.8 , 0.3 )
##     logger.info ( 'Simple test: method %10s,response %s' % ( m , response ) )
##     del method 
## # =============================================================================

## ## 2.2) declare/add TMVA/Chopping  variables 
## for m in methods :
##     variables += [ Variable ( 'tmva_%s' % m , 'TMVA(%s)' % m , accessor = reader[m] ) ]

## # =============================================================================
## ## The END of TMVA/Choppiny reader fragment 
## # =============================================================================


## 3)  Run Ostap to   fill   RooDataSet 
from ostap.fitting.pyselectors import SelectorWithVars     
dsS = SelectorWithVars (
    variables = variables    ,
    selection = "var1 < 100" , 
    )
dsB = SelectorWithVars (
    variables = variables    , 
    selection = "var1 < 100" ,
    )

## read input data file 
with ROOT.TFile.Open( data_file ,'READ') as datafile :
    
    datafile.ls()
    tSignal  = datafile['S']
    tBkg     = datafile['B']

    tSignal.process ( dsS )
    logger.info ( 'Reader signal     histo: %s'   %  reader.histo.stat() ) 
    reader.histo.Reset()
    
    tBkg   .process ( dsB )
    logger.info ( 'Reader background histo: %s'   %  reader.histo.stat() )
    
    ds1 = dsS.data
    ds2 = dsB.data

    del variables   ## attention: reader must be deleted explicitely 
    del reader      ## attention: reader must be deleted explicitely 

    # =========================================================================
    ## B: addChoppingResponse
    #  Much better alternative to TMVA/Chopping reader:
    #  - it has much better performance  :-) 
    #  - but it is less flexible with  repsect to varibale  transformation :-(
    # =========================================================================
    from ostap.tools.chopping import addChoppingResponse

    logger.info ('dataset SIG (no TMVA decisions yet):\n%s' % ds1.table ( prefix = '# ') ) 
    logger.info ('dataset BKG (no TMVA decisions yet):\n%s' % ds2.table ( prefix = '# ') ) 

    addChoppingResponse ( ds1  ,
                          chopper       = "137*evt+813*run"  ,
                          N             =  N                 , 
                          inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                          weights_files = tar_file ,
                          prefix        = 'tmva_'     ,
                          suffix        = '_response' )
    
    addChoppingResponse ( ds2    ,
                          chopper       = "137*evt+813*run"  ,
                          N             =  N                 , 
                          inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                          weights_files = tar_file ,
                          prefix        = 'tmva_'     ,
                          suffix        = '_response' )
    
    # =========================================================================
    ## The END of addChoppingResponse  fragment
    # =========================================================================
    
    logger.info ('dataset SIG (with TMVA decisions):\n%s' % ds1.table ( prefix = '# ') ) 
    logger.info ('dataset BKG (with TMVA decisions):\n%s' % ds2.table ( prefix = '# ') ) 


decisions = [ 'tmva_%s_response' % m for m in methods ]
s_stats   = ds1.statVars ( decisions )
b_stats   = ds2.statVars ( decisions )

table = [ ( 'Method ' , 'Signal' , 'Background ' ) ] 
for d in decisions :
    
    ms = s_stats [ d ]
    mb = b_stats [ d ]
    
    row = d , \
          '%+7.3f +/- %-7.3f' % ( float ( ms.mean() ) , ms.rms() ) ,\
          '%+7.3f +/- %-7.3f' % ( float ( mb.mean() ) , mb.rms() )
    table.append ( row )
    
import ostap.logger.table as T
table  = T.table (  table , title = 'TMVA performance', prefix = '# ' , alignment = 'lcc' )
logger.info ( 'TMVA performance for Signal and Background\n%s' % table )

# =============================================================================
##                                                                      The END
# =============================================================================    
