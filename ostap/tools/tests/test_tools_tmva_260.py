import ROOT, random, os, tempfile 

def  get_name () : 
    
    the_file= tempfile.NamedTemporaryFile ( delete = False , suffix = '.root' )
    the_file.close()
    os.unlink ( the_file.name )
    return the_file.name 
        
def prepare_data ( nB =  1000 , nS = 1000 ) : 

    data_file = get_name() 
    
    test_file = ROOT.TFile.Open( data_file , 'NEW')
    test_file.cd()
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
    for i in range ( nB ) : 
            
        x = random.uniform ( -2.0 , 2.0 )
        y = random.uniform ( -2.0 , 2.0 )
        z = random.gauss   (   .0 , 0.5 )
            
        var1 [ 0 ] =  x + 0.1 * y  
        var2 [ 0 ] =  x - 0.1 * y  
        var3 [ 0 ] = -x +       z
            
        treeBkg.Fill()
            
    ## fill signal tuple: 
    for i in range ( nS ) : 
            
        x = random.gauss  (  0.0 , 0.1 )
        y = random.gauss  (  0.0 , 0.2 )
        z = random.gauss  (  0.5 , 0.5 )
            
        var1 [ 0 ] =  x
        var2 [ 0 ] =  y  
        var3 [ 0 ] =  z
            
        treeSignal.Fill()
            
    treeSignal.Write ()
    treeBkg   .Write ()
    test_file .Write ()
    
    test_file.ls()
    test_file.Close()
    
    ROOT.gROOT.cd()
    
    return data_file 



##  (1) prepare input data 
data_file = prepare_data ( nB = 10000 , nS = 10000 )    

## input data:
tSignal = ROOT.TChain ( 'S' ) ; tSignal.Add ( data_file )
tBkg    = ROOT.TChain ( 'B' ) ; tBkg   .Add ( data_file )

CONFIGURATION = "nTrain_Background=5000:nTrain_Signal=5000:SplitMode=Random:NormMode=NumEvents"
BOOKING       = "Transformations=I;D;P;G,D,G,D,G,D:AnalysisType=Classification:V:!Silent:DrawProgressBar:Color:"

variables = ( 'var1' , 'var2' , 'var3' )

method = ROOT.TMVA.Types.kMLP , "MLP" , "H:!V:EstimatorType=CE:VarTransform=N:NCycles=200:HiddenLayers=N+5:TestRate=5:!UseRegulator" 

tSignal.Print ( 'vvv' )
tBkg   .Print ( 'vvv' )

import ostap.trees.trees
print ( tSignal.table() )
print ( tSignal.table() )

name = 'TBkg-test'

with ROOT.TFile.Open(  get_name() , 'NEW' ) as out_file :
    
    out_file.cd ()
    
    factory = ROOT.TMVA.Factory ( name , out_file , BOOKING )
    
    factory.SetVerbose( True )
    
    dataloader = ROOT.TMVA.DataLoader ( name )
    for v in variables : dataloader.AddVariable ( v, 'F')

    the_cut = ROOT.TCut() 
    dataloader.AddTree ( tSignal  , 'Signal'     , 1.0 , the_cut ) 
    dataloader.AddTree ( tBkg     , 'Background' , 1.0 , the_cut )
    
    dataloader.PrepareTrainingAndTestTree( the_cut , the_cut , CONFIGURATION )
    
    di = dataloader.DataInput()
    
    print ( "AFTER:", di.GetNTrees("Signal")     , ' #sig-trees' )
    print ( "AFTER:", di.GetNTrees("Background") , ' #bkg-trees' )
    
    
    print ( "AFTER:", di.GetSignalEntries()      ,  ' signal entries/1')
    print ( "AFTER:", di.GetBackgroundEntries()   , ' bkgentries/1')
    
    print ( "AFTER:", di.GetEntries('Signal')     , ' signal entries/2')
    print ( "AFTER:", di.GetEntries('Background') , ' bkgentries/2')
    
    factory.BookMethod ( dataloader , *method )
    
    factory.TrainAllMethods    ()
    factory.TestAllMethods     ()
    factory.EvaluateAllMethods ()

    
## the rest is irrelevannt, crush is above

