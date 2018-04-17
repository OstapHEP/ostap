#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file test_chopping.py
#
#  Test for TMVA ``chopping'' machinery
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
import ROOT, os
import ostap.io.root_file 
from   ostap.core.core          import ROOTCWD
from   ostap.utils.progress_bar import progress_bar 
from   array                    import array
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'test_chopping' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
from ostap.utils.utils import CleanUp
data_file = CleanUp.tempfile ( suffix = '.root' , prefix = 'test_chopping_' ) 

if not os.path.exists( data_file ) :
    import random 
    nB = 20000
    nS = 10000
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
        #for i in progress_bar ( xrange ( nB ) ) : 
        for i in xrange ( nB ) : 
            
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
        #for i in progress_bar ( xrange ( nS ) ) : 
        for i in xrange ( nS ) : 
            
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

"""
##   number of    categories 
N  = 7
logger.info('Create and train TMVA')
with ROOT.TFile.Open( data_file ,'READ') as datafile : 
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
        name    = 'TestTMVA' ,   
        methods = [ # type               name   configuration
        ( ROOT.TMVA.Types.kMLP        , "MLP"        , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=200:HiddenLayers=N+3:TestRate=5:!UseRegulator" ) ,
        ( ROOT.TMVA.Types.kBDT        , "BDTG"       , "H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" ) , 
        ( ROOT.TMVA.Types.kCuts       , "Cuts"       , "H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" ) ,
        ( ROOT.TMVA.Types.kFisher     , "Fisher"     , "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" ),
        ( ROOT.TMVA.Types.kLikelihood , "Likelihood" , "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )
        ] ,
        variables = [ 'var1' , 'var2' ,  'var3' ] , ## Variables for training 
        signal         = tSignal                  , ## ``Signal'' sample
        background     = tBkg                     , ## ``Background'' sample         
        verbose        = False    )

    from ostap.utils.timing import timing

    # sequential trainig 
    #with timing ( 'for TMVA training' , logger ) : 
    #    weights_files = trainer.train ()
    #    tar_file      = trainer.tar_file
    
    #parallel trainig 
    with timing ( 'for TMVA training' , logger ) : 
        trainer.ptrain() 
        tar_file      = trainer.tar_file
        
# =============================================================================
# remove unnesessary output files
for f in trainer.output_files :
    if os.path.exists ( f ) and os.path.isfile( f ) :
        try    : os.remove ( f )
        except : pass
    
# =============================================================================
## Use trained TMVA
# =============================================================================


## 1) create TMVA reader
from ostap.tools.chopping import Reader
category = lambda s :  int ( s.evt*137 + 813*s.run ) % N
reader = Reader(
    N             = N         , ##  number of   categories
    categoryfunc  = category  , ## category 
    ## other argument  as for plain TMVA     
    name = 'ChopReader' ,
    variables     = [ ('var1' , lambda s : s.var1 )   ,
                      ('var2' , lambda s : s.var2 )   ,
                      ('var3' , lambda s : s.var3 ) ] ,
    weights_files = tar_file )

methods = reader.methods[:]

# =============================================================================
## 1') few trivial tests: use the methods/reader as simple function
for m in methods :
    method   = reader[m]
    ## response = [ method ( i  , 1.1 , 0.8 , 0.3 ) for i in  range ( reader.N ) ] 
    response = method.stat ( 1.1 , 0.8 , 0.3 )
    logger.info ( 'Simple test: method %10s,response %s' % ( m , response ) )
    del method 
# =============================================================================

from ostap.fitting.selectors import SelectorWithVars,  Variable     
## 2) Book RooDataset                 
variables = [
    Variable ( 'var1' , 'variable#1' , accessor = lambda s : s.var1 ) ,
    Variable ( 'var2' , 'variable#2' , accessor = lambda s : s.var2 ) ,
    Variable ( 'var3' , 'variable#3' , accessor = lambda s : s.var3 ) ,
    ## extra: needed for addChoppingResponse 
    Variable ( 'evt'  , 'event'      , accessor = lambda s : s.evt  ) ,
    Variable ( 'run'  , 'run'        , accessor = lambda s : s.run  ) ,
    ## extra: needed for cross-checks  
    Variable ( 'cat'  , 'category'   , accessor = category          ) ,
    ]

## 3) declare/add TMVA  variables 
for m in methods :
    variables += [ Variable ( 'tmva_%s' % m , 'TMVA(%s)' % m , accessor = reader[m] ) ]
    
## 4)  Run Ostap to   fill   RooDataSet 
from ostap.fitting.selectors import SelectorWithVars     
dsS = SelectorWithVars (
    variables = variables + [ Variable ( 'signal' , 'signal' , -1 , 3 , lambda s : 1 ) ] ,
    selection = "var1 < 100" , 
    )
dsB = SelectorWithVars (
    variables = variables + [ Variable ( 'signal' , 'signal' , -1 , 3 , lambda s : 0 ) ] ,
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

    del variables 
    del reader

    from ostap.tools.chopping import addChoppingResponse

    logger.info ('dataset SIG: %s' %  ds1 )
    logger.info ('dataset BKG: %s' %  ds2 )
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
    
    logger.info ('dataset SIG: %s' %  ds1 )
    logger.info ('dataset BKG: %s' %  ds2 )
    
for m in methods :
    
    logger.info('TMVA:%-11s for signal     %s' % ( m, ds1.statVar('tmva_%s' % m ) ) )
    logger.info('TMVA:%-11s for background %s' % ( m, ds2.statVar('tmva_%s' % m ) ) )


"""
# =============================================================================
# The END
# =============================================================================    
