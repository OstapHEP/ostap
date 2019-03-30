#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file test_tmva.py
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
import ROOT, os
import ostap.io.root_file 
from   builtins                 import range  
from   ostap.core.core          import ROOTCWD
from   ostap.utils.progress_bar import progress_bar 
from   array                    import array
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'test_tmva' )
else : 
    logger = getLogger ( __name__ )
# ==============================================================================
from ostap.utils.utils import CleanUp
data_file = CleanUp.tempfile ( suffix = '.root' , prefix = 'test_tmva_' )
if not os.path.exists( data_file ) :
    import random 
    nB = 20000
    nS = 20000
    logger.info('Prepare input ROOT file with data  %s' % data_file )
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
        
        treeSignal.Branch ( 'var1' , var1 , 'var1/D' )
        treeSignal.Branch ( 'var2' , var2 , 'var2/D' )
        treeSignal.Branch ( 'var3' , var3 , 'var3/D' )
        
        treeBkg   .Branch ( 'var1' , var1 , 'var1/D' )
        treeBkg   .Branch ( 'var2' , var2 , 'var2/D' )
        treeBkg   .Branch ( 'var3' , var3 , 'var3/D' )
        
        ## fill background tuple: 
        #for i in progress_bar ( range ( nB ) ) : 
        for i in range ( nB ) : 
            
            x = random.uniform ( -2.0 , 2.0 )
            y = random.uniform ( -2.0 , 2.0 )
            z = random.gauss   (   .0 , 0.5 )
            
            var1[0] =  x + 0.1 * y  
            var2[0] =  x - 0.1 * y  
            var3[0] = -x +       z
            
            treeBkg.Fill()
            
        ## fill signal tuple: 
        #for i in progress_bar ( range ( nS ) ) : 
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
        verbose        = True  )
    
    from ostap.utils.timing import timing
    with timing ( 'for TMVA training' , logger ) : 
        weights_files = trainer.train ( log = True )
        tar_file      = trainer.tar_file

# =============================================================================
if os.path.exists ( trainer.output_file ) :
    os.remove ( trainer.output_file )

# =============================================================================
## Use trained TMVA
# =============================================================================


## 1) create TMVA reader
from ostap.tools.tmva import Reader
reader = Reader( 'MyMLP' ,
                 variables     = [ ('var1' , lambda s : s.var1 )   ,
                                   ('var2' , lambda s : s.var2 )   ,
                                   ('var3' , lambda s : s.var3 ) ] ,
                 weights_files = tar_file   )

methods = reader.methods[:]


# =============================================================================
## 1') few trivial tests: use the methods/reader as simple function
for m in methods :
    method = reader[m]
    logger.info ( 'Method  %12s , response %s' % ( m , method ( 1.1 , 0.8 , 0.3 ) ) )
    del method 
# =============================================================================
                  

from ostap.fitting.selectors import SelectorWithVars, Variable     
## 2) Book RooDataset                 
variables = [
    #Variable( 'var1' , 'variable#1' , accessor = lambda s : s.var1 ) ,
    #Variable( 'var2' , 'variable#2' , accessor = lambda s : s.var2 ) ,
    #Variable( 'var3' , 'variable#3' , accessor = lambda s : s.var3 ) ,
    Variable( 'var1' , 'variable#1' ) ,
    Variable( 'var2' , 'variable#2' ) ,
    Variable( 'var3' , 'variable#3' ) ,
    ]

## ## 3) declare/add TMVA  variables 
## for m in methods :
##     variables += [ Variable( 'tmva_%s' % m , 'TMVA(%s)' % m , accessor = reader[m] ) ]
    
## 4)  Run Ostap to   fill   RooDataSet 
dsS = SelectorWithVars (
    variables = variables , ##   + [ Variable ( 'signal' , 'signal' , -1 , 3 , lambda s : 1 ) ] ,
    selection = "var1 < 100" , 
    )
dsB = SelectorWithVars (
    variables = variables , ## + [ Variable ( 'signal' , 'signal' , -1 , 3 , lambda s : 0 ) ] ,
    selection = "var1 < 100" ,
    )

## read input data file 
with ROOT.TFile.Open( data_file ,'READ') as datafile :
    
    datafile.ls()
    tSignal  = datafile['S']
    tBkg     = datafile['B']

    from ostap.utils.timing import timing
    with timing("Process signal"     ) : tSignal.process ( dsS )
    with timing("Process backrgound" ) : tBkg   .process ( dsB )
    
    ds1 = dsS.data
    ds2 = dsB.data

    del variables 
    del reader

    from ostap.tools.tmva import addTMVAResponse

    logger.info ('dataset SIG: %s' %  ds1 )
    logger.info ('dataset BKG: %s' %  ds2 )
    addTMVAResponse ( ds1 ,
                      inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                      weights_files = tar_file ,
                      prefix        = 'tmva_'     ,
                      suffix        = '_response' )
    addTMVAResponse ( ds2    ,
                      inputs        = ( 'var1' ,  'var2' , 'var3' ) ,
                      weights_files = tar_file ,
                      prefix        = 'tmva_'     ,
                      suffix        = '_response' )
    
    logger.info ('dataset SIG: %s' %  ds1 )
    logger.info ('dataset BKG: %s' %  ds2 )


for m in methods :
    
    logger.info('TMVA:%-11s for signal     %s' % ( m, ds1.statVar('tmva_%s_response' % m ) ) )
    logger.info('TMVA:%-11s for background %s' % ( m, ds2.statVar('tmva_%s_response' % m ) ) )

# =============================================================================
# The END
# =============================================================================    
