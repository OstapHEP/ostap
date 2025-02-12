#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/tools/tests/test_tools_splot.py
# Test sPlotting macheinery (outside of RooFit)
# - fill histogram from the TTree 
# - fit the histogram
# - make sPlot
# - write sPlot results into oroigin Ttree
# ============================================================================= 
""" Test sPlotting macheinery (outside of RooFit)
- fill histoigram from TTree 
- fit the histogram
- make sPlot
- write sPlot results into oroigin Ttree
"""
# ============================================================================= 
from   ostap.core.core            import Ostap, hID
from   ostap.trees.data           import Data
from   ostap.utils.timing         import timing 
from   ostap.utils.progress_bar   import progress_bar
from   ostap.fitting.variables    import FIXVAR 
from   ostap.tools.splot          import sPlot2D
from   ostap.plotting.canvas      import use_canvas 
from   ostap.utils.utils          import wait, batch_env 
import ostap.fitting.models       as     Models
import ostap.io.zipshelve         as     DBASE 
import ostap.fitting.pyselectors 
import ostap.trees.trees
import ostap.histos.histos
import ROOT, math, random, array  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_tools_splot2' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================


xmin , xmax  =    0 ,   10 
ymin , ymax  = xmin , xmax
zmin , zmax  = xmin , xmax
mean , sigma =    5 , 0.50 
xvar = ROOT.RooRealVar ( 'x' , 'x-var' , xmin , xmax )
yvar = ROOT.RooRealVar ( 'y' , 'y-var' , ymin , ymax )
zvar = ROOT.RooRealVar ( 'z' , 'z-var' , zmin , zmax )

# =============================================================================
## create a file with tree 
def create_tree ( fname , nentries = 1000 ) :
    """ Create a file with a tree
    >>> create_tree ( 'file.root' ,  1000 ) 
    """
    
    import ostap.io.root_file
    
    from array import array 
    x , y , z = array ( 'd', [ 0 ] ), array ( 'd', [ 0 ] ) , array ( 'd', [ 0 ] )
    from ostap.core.core import ROOTCWD

    with ROOTCWD() , ROOT.TFile.Open( fname , 'new' ) as root_file:
        
        root_file.cd ()
        
        tree = ROOT.TTree ( 'S','tree' )
        tree.SetDirectory ( root_file  )
        
        tree.Branch ( 'x'  , x , 'x/D' )
        tree.Branch ( 'y'  , y , 'y/D' )
        tree.Branch ( 'z'  , z , 'z/D' )
        
        for i in range ( nentries ) :
            
            xv , yv, zv = -100 , -100 , -100
            
            if 0 == i % 2 :  
                while not xmin < xv < xmax : xv = random.gauss   ( mean , sigma )
                while not ymin < yv < ymax : yv = random.gauss   ( mean , sigma )
                while not zmin < zv < zmax : zv = random.gauss   ( 2.5  , 0.75  ) 
            else :
                while not xmin < xv < xmax : xv = random.uniform ( xmin , xmax  ) 
                while not ymin < yv < ymax : yv = random.uniform ( ymin , ymax  ) 
                while not zmin < zv < zmax : zv = random.gauss   ( 7.5  , 0.75  ) 

            x [ 0 ] = xv
            y [ 0 ] = yv
            z [ 0 ] = zv
            
            tree.Fill()
            
        root_file.Write()
        
# =============================================================================
def prepare_data ( nfiles = 50 ,  nentries = 500  ) :

    from ostap.utils.cleanup import CleanUp    
    files = [ CleanUp.tempfile ( prefix = 'ostap-test-tools-splot2-%d-' % i ,
                                 suffix = '.root' ) for i in range ( nfiles)  ]
    
    for f in progress_bar ( files ) : create_tree ( f , nentries )
    return files

# =============================================================================
## Test sPlotting macheinery (outside of RooFit)
# - fill histogram for Ttree 
# - fit the histogram
# - make sPlot
# - write sPlot results into origin TTree
def test_splotting2 () :
    """ Test sPlotting macheinery (outside of RooFit)
    - fill histogram from TTree 
    - fit the histogram
    - make sPlot
    - write sPlot results into oroigin Ttree
    """
    
    files = prepare_data ( 10 , 100000 )
    
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )
    
    chain  = data.chain

    title = 'Input TTree/TChain'
    logger.info ( '%s:\n%s' % ( title , chain.table ( title = title , prefix = '# ' ) ) ) 

    ## 
    ## dataset , _ = chain.make_dataset ( ( xvar , yvar , zvar ) )
    ## title = 'Dataset form TTree'
    ## logger.info ( '%s:\n%s' % ( title , dataset.table ( title = title , prefix = '# ' ) ) ) 
    ##
    
    histo = ROOT.TH2F ( hID() , 'input histogram' , 200 , xmin , xmax , 200 , ymin , ymax )
    data.chain.project ( histo , 'x, y' ) 

    dataset  = histo 

    signalx  = Models.Gauss_pdf ( 'GX' , xvar = xvar , mean = ( mean ,  xmin , xmax  ) , sigma = ( sigma, 0.5 * sigma , 2.0 * sigma ) )
    signaly  = Models.Gauss_pdf ( 'GY' , xvar = yvar , mean = signalx.mean , sigma = signalx.sigma )
    
    ## 2 models 
    model    = Models.Fit2D     ( signal_x = signalx , signal_y = signaly )
    
    model.SS = 0.5 * len ( data.chain )
    model.BB = 0.5 * len ( data.chain )
    model.SB.setMin ( -100000 )  
    model.BS.setMin ( -100000 )  

    model.SB = 0.0
    model.BS = 0.0
    
    with FIXVAR ( [ model.SB , model.BS ] ) : 
        result , _ = model.fitTo ( dataset , silent = True , refit = 5 )
        
    result , _ = model.fitTo ( dataset , quiet = True , refit = 5 )
    dataset    = model.histo_data.dset
    
    with use_canvas ( 'test_splot2: X' ) : model.draw1 ( dataset , nbins = 50 )
    with use_canvas ( 'test_splot2: Y' ) : model.draw2 ( dataset , nbins = 50 )
    
    with timing ( "Splotting!" , logger = logger ) :
        splot = sPlot2D ( model , dataset , xbins = 100 , ybins = 100  ) ## SPLOT IT!

    with timing ( "Adding sPlot results to TTree/TChain" , logger = logger ) :
        chain = data.chain
        splot.add_to_tree ( chain , 'x' , 'y' , parallel = False ) 
        
    with use_canvas ( 'test_splot2: Z ' , wait = 5 ) :

        chain = data.chain
        
        chain.draw ( 'z' , 'SS_sw' ,          color = 2 )
        chain.draw ( 'z' , 'BB_sw' , 'same' , color = 4 )
        chain.draw ( 'z' , 'SB_sw' , 'same' , color = 6 )
        chain.draw ( 'z' , 'BS_sw' , 'same' , color = 7 )

        title = 'TTree with sPlot-results'
        logger.info ( '%s:\n%s' % ( title , chain.table ( title  = title , prefix = '# ' ) ) ) 

    for var in ( 'SS_sw' , 'BB_sw' , 'SB_sw' , 'BS_sw' ) :        
        logger.info  ('%s variable sum : %s' % ( var , chain.sumVar  ( var ) ) )
        
    with DBASE.tmpdb()  as db :
        db [ 'histo'  ] = histo 
        db [ 'splot2' ] = splot
        db.ls()
            
# =============================================================================
if '__main__' ==  __name__  :

    with timing ( "test_splotting2" , logger = logger ) :
        test_splotting2 ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
