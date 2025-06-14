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
from   ostap.core.core          import Ostap, hID 
from   ostap.trees.data         import Data
from   ostap.utils.basic        import typename 
from   ostap.utils.timing       import timing 
from   ostap.utils.progress_bar import progress_bar
from   ostap.fitting.variables  import FIXVAR 
from   ostap.tools.splot        import COWs, sPLOT, hPlot1D
from   ostap.plotting.canvas    import use_canvas 
from   ostap.utils.root_utils   import batch_env
import ostap.fitting.models     as     Models 
import ostap.io.zipshelve       as     DBASE 
import ostap.trees.trees
import ostap.histos.histos
import ROOT, math, random, array  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_tools_splot' )
else                       : logger = getLogger ( __name__           )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

xmin , xmax  = 3.0 , 3.2
ymin , ymax  = 0   ,  10 
mean , sigma = 3.1 , 0.01
# =============================================================================
## create a file with tree 
def create_tree ( fname , nentries = 1000 ) :
    """ Create a file with a tree
    >>> create_tree ( 'file.root' ,  1000 ) 
    """
    
    import ostap.io.root_file
    
    from array import array 
    var1 = array ( 'd', [ 0 ] )
    var2 = array ( 'd', [ 0 ] )
    
    from ostap.core.core import ROOTCWD

    with ROOTCWD() , ROOT.TFile.Open( fname , 'new' ) as root_file:
        root_file.cd () 
        tree = ROOT.TTree ( 'S','tree' )
        tree.SetDirectory ( root_file  )
        
        tree.Branch ( 'x'  , var1 , 'x/D' )
        tree.Branch ( 'y'  , var2 , 'y/D' )
        
        for i in range ( nentries ) : 

            x, y = -100 , -100
        
            if 0 == i % 2 : 
                while not xmin < x < xmax : x = random.gauss   ( mean , sigma )
                while not ymin < y < ymax : y = random.gauss   ( 2.5  , 0.7   )
            else       :
                while not xmin < x < xmax : x = random.uniform ( xmin , xmax  ) 
                while not ymin < y < ymax : y = random.gauss   ( 7.5  , 0.7   )
                
            var1 [ 0 ] = x 
            var2 [ 0 ] = y  
            
            tree.Fill()
            
        root_file.Write()
        
# =============================================================================
def prepare_data ( nfiles = 50 ,  nentries = 500  ) :

    from ostap.utils.cleanup import CleanUp    
    files = [ CleanUp.tempfile ( prefix = 'ostap-test-tools-splot-%d-' % i ,
                                 suffix = '.root' ) for i in range ( nfiles)  ]
    
    for f in progress_bar ( files ) : create_tree ( f , nentries )
    return files

# =============================================================================
## Test sPlotting macheinery (outside of RooFit)
# - fill histogram for Ttree 
# - fit the histogram
# - make sPlot
# - write sPlot results into origin TTree
def test_splotting  () :
    """ Test sPlotting machinery (outside of RooFit)
    - fill histogram from TTree 
    - fit the histogram
    - make sPlot
    - write sPlot results into oroigin Ttree
    """
    
    files = prepare_data ( 10 , 10000 )
    
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data  ( files , 'S' )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )
    
    chain  = data.chain 
    histo = ROOT.TH1D ( hID() , 'x-distibution' , 200 , xmin , xmax  )

    chain.project ( histo , 'x' )
    
    xvar   = ROOT.RooRealVar ( 'x' , 'c-variable' , xmin , xmax )    
    gauss = Models.Gauss_pdf ( 'G' , xvar = xvar ,
                               mean  = ( mean  , mean - 2 * sigma , mean + 2 * sigma ) ,
                               sigma = ( sigma , 0.5 * sigma , 2.0 * sigma           ) )
    
    gauss.sigma.fix()    
    model = Models.Fit1D ( signal = gauss , background = 0 )
    NS = 0.5 * ( histo.Integral() )
    NB = 0.5 * ( histo.Integral() )
    
    model.S.setMax ( 3 * NS )
    model.B.setMax ( 3 * NS )
    
    model.S = NS 
    model.B = NB 

    with FIXVAR ( ( gauss.mean , gauss.sigma ) ) : 
        model.fitTo ( histo , draw = False , silent = True , refit = 5 )
        
    model.fitTo ( histo , draw = False , silent = True )

    with use_canvas ( 'test_tools_splot: Fit') :
        r, f = model.fitTo ( histo , draw = True  , nbins = 100 , silent = True ) 
        
    title = "Fit results for" 
    logger.info ('%s:\n%s' % ( title  , r.table ( title = title , prefix = '# ') ) ) 

    ## Binned machinery with binnings  
    for fast in ( True , False ) :

        sp  = hPlot1D ( model , dataset = histo  , nbins = 500 , fast = fast ) ## SPLOT IT! 
        sph = sp.hweights['S']
        
        with use_canvas ( 'test_tools_splot: sPlot') :
            sph .draw ()
            value =  float ( sph ( r.mean_G * 1 ) ) 
            assert 0.7 < value < 1.5 , "Something totally wrong here, fast=%s!" % fast 
                
        fnsp = Ostap.Functions.FuncTH1 ( sph , 'x' )
        
        with DBASE.tmpdb()  as db :
            
            db ['histo;fast=%s'   % fast ] = histo 
            db ['splot;fast=%s'   % fast ] = sp    
            db ['splot,h;fast=%s' % fast ] = sph   
            db ['splot,f;fast=%s' % fast ] = fnsp  
            db.ls()

        with timing ( "Adding   BINNED sPlot results to TTree/TChain" , logger = logger ) :
            chain = data.chain
            sp.add_to_tree ( chain , 'x' , prefix = 'f' if fast else 'n' , parallel = False  )


    ds = model.histo_data.dset 
    ## unbineed machinery 
    with timing ( "Adding UNBINNED sPLOT results to TTree/TChain" , logger = logger ) :
        
        chain = data.chain 
        sp    = sPLOT ( model , dataset = ds ) ## SPLOT IT! 
        sp.splot2tree ( chain , 'x' , prefix = 'u' , parallel = True )

    ## COWs
    with timing ( "Adding UNBINNED COWs results to TTree/TChain" , logger = logger ) :
        
        chain = data.chain 
        sp    = COWs ( model , dataset = ds ) ## SPLOT IT!
        names = tuple ( v.name + '_cw' for v in model.alist2 ) 
        sp.cows2tree ( chain , 'x' , names = names  , parallel = True )
          
    chain = data.chain
    
    nS_hw = chain.statVar ( 'nS_hw' )
    fS_hw = chain.statVar ( 'fS_hw' )
    uS_sw = chain.statVar ( 'uS_sw' )
    S_cw  = chain.statVar ( 'S_cw'  )
    B_cw  = chain.statVar ( 'B_cw'  )

    dun   = chain.statVar ( 'uS_sw - nS_hw' )
    duf   = chain.statVar ( 'uS_sw - fS_hw' )
    dnf   = chain.statVar ( 'nS_hw - fS_hw' )
    duc   = chain.statVar ( 'uS_sw -  S_cw' )
    
    logger.info ( 'nS_hw: %s' % nS_hw )
    logger.info ( 'fS_hw: %s' % fS_hw )    
    logger.info ( 'uS_sw: %s' % uS_sw )
    logger.info ( 'S_cw : %s' %  S_cw )
    
    logger.info ( 'u-n  : %s' % dun   )
    logger.info ( 'u-f  : %s' % duf   )
    logger.info ( 'n-f  : %s' % dnf   )
    logger.info ( 'u-c  : %s' % duc   )

    hsn  = ROOT.TH1D( hID() , 'y for signal'     , 100 , ymin , ymax )
    hbn  = ROOT.TH1D( hID() , 'y for background' , 100 , ymin , ymax )
    hsf  = ROOT.TH1D( hID() , 'y for signal'     , 100 , ymin , ymax )
    hbf  = ROOT.TH1D( hID() , 'y for background' , 100 , ymin , ymax )
    hsu  = ROOT.TH1D( hID() , 'y for signal'     , 100 , ymin , ymax )
    hbu  = ROOT.TH1D( hID() , 'y for background' , 100 , ymin , ymax )
    hsc  = ROOT.TH1D( hID() , 'y for signal'     , 100 , ymin , ymax )
    hbc  = ROOT.TH1D( hID() , 'y for background' , 100 , ymin , ymax )
    
    hsn.red  ()
    hbn.blue ()
    hsf.red  ()
    hbf.blue ()
    hsu.red  ()
    hbu.blue ()
    hsc.red  ()
    hbc.blue ()

    with use_canvas ( 'test_tools_splot: y/n (UNBINNED)', wait = 2 ) :
            
        chain.project ( hsu , 'y' , 'uS_sw'   )
        hsu.draw ()

        chain.project ( hbu , 'y' , 'uB_sw'   )        
        hbu.draw ('same')

    with use_canvas ( 'test_tools_splot: y/n (N-BINNED)', wait = 2 ) :
            
        chain.project ( hsn , 'y' , 'nS_hw'   )
        hsn.draw ()

        chain.project ( hbn , 'y' , 'nB_hw'   )        
        hbn.draw ('same')

    with use_canvas ( 'test_tools_splot: y/n (F-BINNED)', wait = 2 ) :
            
        chain.project ( hsf , 'y' , 'fS_hw'   )
        hsf.draw ()

        chain.project ( hbf , 'y' , 'fB_hw'   )        
        hbf.draw ('same')

    with use_canvas ( 'test_tools_splot: y/n (COWs)', wait = 2 ) :
            
        chain.project ( hsc , 'y' , 'S_cw'   )
        hsc.draw ()

        chain.project ( hbc , 'y' , 'B_cw'   )        
        hbc.draw ('same')
        
# =============================================================================
if '__main__' ==  __name__  :

    with timing ( "test_splotting" , logger = logger ) :
        test_splotting ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================

