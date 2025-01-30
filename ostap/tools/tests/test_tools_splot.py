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
from   __future__               import print_function
import ostap.trees.trees
import ostap.histos.histos
from   ostap.core.core          import Ostap, hID 
from   ostap.trees.data         import Data
from   ostap.utils.timing       import timing 
from   ostap.utils.progress_bar import progress_bar
from   ostap.fitting.variables  import FIXVAR 
from   ostap.tools.splot        import sPlot1D
from   ostap.plotting.canvas    import use_canvas 
from   ostap.utils.utils        import wait
import ostap.fitting.models     as     Models 
import ostap.io.zipshelve       as     DBASE 
import ROOT, math, random, array  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_tools_splot' )
else                       : logger = getLogger ( __name__           )
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
    
    files = prepare_data ( 50 , 100000 )
    
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )
    
    chain  = data.chain 
    histo = ROOT.TH1D       ( hID() , 'x-distibution' , 200 , xmin , xmax  )

    chain.project ( histo , 'x' )
    
    xvar   = ROOT.RooRealVar ( 'x' , 'c-variable' , xmin , xmax )
    
    gauss = Models.Gauss_pdf ( 'G' , xvar = xvar ,
                               mean  = ( mean  , mean - 2 * sigma , mean + 2 * sigma ) ,
                               sigma = ( sigma , 0.5 * sigma , 2.0 * sigma           ) )
    
    gauss.sigma.fix()    
    model = Models.Fit1D( signal = gauss , background = 0 )
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
  
    with DBASE.tmpdb()  as db :
        
        for fast in ( True , False ) :
            
            sp  = sPlot1D ( model , histo  , nbins = 500 , fast = fast ) ## SPLOT IT! 
            sph = sp.hweights['S']
            
            with use_canvas ( 'test_tools_splot: sPlot') :
                sph .draw ()
                value =  float ( sph ( r.mean_G * 1 ) ) 
                assert 0.7 < value < 1.5 , "Something totally wrong here, fast=%s!" % fast 
                
            fnsp = Ostap.Functions.FuncTH1 ( sph , 'x' )
            
            db ['histo;fast=%s'   % fast ] = histo 
            db ['splot;fast=%s'   % fast ] = sp    
            db ['splot,h;fast=%s' % fast ] = sph   
            db ['splot,f;fast=%s' % fast ] = fnsp  
            db.ls()
            
            with timing ( "Adding sPlot results to TTree/TChain" , logger = logger ) :
                chain = data.chain
                sp.add_to_tree ( chain , 'x' , prefix = 'f' if fast else 'n' , parallel = False ) 

    chain = data.chain 
    nS_sw = chain.statVar ( 'nS_sw' )
    fS_sw = chain.statVar ( 'fS_sw' )
    diff  = chain.statVar ( 'nS_sw-fS_sw' )
    
    logger.info ( 'nS_sw: %s' % nS_sw )
    logger.info ( 'fS_sw: %s' % fS_sw )    
    logger.info ( 'diff : %s' % diff )


    hs  = ROOT.TH1D( hID() , 'y for signal'     , 100 , ymin , ymax )
    hb  = ROOT.TH1D( hID() , 'y for background' , 100 , ymin , ymax )
    hsf = ROOT.TH1D( hID() , 'y for signal'     , 100 , ymin , ymax )
    hbf = ROOT.TH1D( hID() , 'y for background' , 100 , ymin , ymax )

    hs .red  ()
    hb .blue ()
    hsf.red  ()
    hbf.blue ()
    
    with use_canvas ( 'test_tools_splot: y/n', wait = 5 ) :
            
        chain.project ( hs , 'y' , 'nS_sw'   )
        hs.draw ()

        chain.project ( hb , 'y' , 'nB_sw'   )        
        hb.draw ('same')

    with use_canvas ( 'test_tools_splot: y/f ', wait = 5 ) :

        chain.project ( hsf, 'y' , 'fS_sw'   )
        hsf.draw()

        chain.project ( hbf, 'y' , 'fB_sw'   )        
        hbf.draw('same')
        
    with use_canvas ( 'test_tools_splot: delta', wait = 5 ) :
        chain.draw ( 'nS_sw-fS_sw' )

# =============================================================================
if '__main__' ==  __name__  :

    with timing ( "test_splotting" , logger = logger ) :
        test_splotting ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================

