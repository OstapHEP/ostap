1#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight4.py
#  Test for 3D-reweighting machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2023-01-20
# =============================================================================
"""Test for 3D-reweighting machinery
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2025-09-15"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   ostap.core.pyrouts     import Ostap 
from   ostap.histos.histos    import h1_axis, h2_axes, h3_axes 
from   ostap.utils.timing     import timing
from   ostap.logger.colorized import attention, allright
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.root_utils import batch_env 
from   ostap.utils.cleanup    import CleanUp
import ostap.io.zipshelve     as     DBASE
import ostap.logger.table     as     T 
import ostap.logger.table     as     T 
import ostap.io.root_file 
import ostap.parallel.kisa
import ROOT, random, math, os, time 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_reweight3' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for ND-Reweighting machinery')
# ============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
testdata   = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-test-tools-reweight4-' )

tag_data_r   = 'DATA_R_histogram'
tag_data_r1  = 'DATA_R1_histogram'
tag_data_r2  = 'DATA_R2_histogram'
tag_data_r3  = 'DATA_R3_histogram'

tag_data     = 'DATA_tree'
tag_mc       = 'MC_tree'

dbname       = CleanUp.tempfile ( suffix = '.db' , prefix ='ostap-test-tools-reweight4-'   )

NDATA        =  500000
NMC          = 1000000

rmax = 20 

def prepare_data ( ) : 
    #
        
    seed =  1234567890 
    random.seed ( seed ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )   
    ## prepare "data" histograms:


    ## X,Y,Z histogramss
    
    ir = 100 
    hr1_histo  = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )
    hr2_histo  = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )
    hr3_histo  = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )
    hr_histo   = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )

    with ROOT.TFile.Open ( testdata ,'recreate') as mc_file :
               
        mc_file.cd() 
        
        datatree  = ROOT.TTree ( tag_data , 'data-tree' )
        datatree .SetDirectory ( mc_file ) 
        
        from array import array 
        r1var = array    ( 'f', [0])
        r2var = array    ( 'f', [0])
        r3var = array    ( 'f', [0])
        
        datatree.Branch ( 'r1' , r1var , 'r1/F' )
        datatree.Branch ( 'r2' , r2var , 'r2/F' )
        datatree.Branch ( 'r3' , r3var , 'r3/F' )

        ## Gaussian 3D-component with correlations 
               
        nentries = 0  
        while nentries <  NDATA :
            
            r1 = random.uniform ( 0 , rmax ) 
            r2 = random.uniform ( 0 , rmax ) 
            r3 = random.uniform ( 0 , rmax ) 
            
            if 0.75 < random.uniform ( 0 , 1 ) : 
                rr1 = min ( r1 , r2 , r3 ) 
                rr2 =     ( r1 + r2 + r3 ) / 3 
                rr3 = max ( r1 , r2 , r3 ) 
                r1, r2 , r3 = rr1  , rr2  , rr3 
                
                
            hr1_histo.Fill ( r1 ) 
            hr2_histo.Fill ( r2 ) 
            hr3_histo.Fill ( r3 ) 
                        
            r1var [ 0 ] = rr1 
            r2var [ 0 ] = rr2 
            r3var [ 0 ] = rr3 
            
            datatree.Fill()
            nentries += 1 
            

        ## write the histogram 
    
        mc_file [ tag_data_r1 ] = hr1_histo
        mc_file [ tag_data_r2 ] = hr2_histo
        mc_file [ tag_data_r3 ] = hr3_histo
                
        mc_file [ tag_data_r  ] = hr1_histo + hr2_histo + hr3_histo

        ## 
        mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
        mctree .SetDirectory ( mc_file ) 
        
        from array import array 

        r1var = array  ( 'f', [ 0.0 ] )
        r2var = array  ( 'f', [ 0.0 ] )
        r3var = array  ( 'f', [ 0.0 ] )

        mctree.Branch ( 'r1' , r1var , 'r1/F' )
        mctree.Branch ( 'r2' , r2var , 'r2/F' )
        mctree.Branch ( 'r3' , r3var , 'r3/F' )
        
        nentries = 0 
        while nentries <  NMC :


            r1 = random.gauss ( 0.5  * rmax , 0.40 * rmax )
            if r1 < 0 or  rmax < r1 : continue
            
            r2 = random.gauss ( 0.5  * rmax , 0.40 * rmax )
            if r2 < 0 or  rmax < r2 : continue

            r3 = random.gauss ( 0.5  * rmax , 0.40 * rmax )
            if r3 < 0 or  rmax < r3 : continue
            
            rrr = r1 , r2  , r3 
            
            if 0.75 < random.uniform ( 0 , 1 ) : rrr = sorted ( rrr )
            
            r1var [ 0 ] = rrr[0] 
            r2var [ 0 ] = rrr[1]
            r3var [ 0 ] = rrr[2]
            
            mctree.Fill()
            
            nentries += 1 
            
        mctree.Write()

if not os.path.exists( testdata ) :
    with timing ( "Prepare input data" , logger = logger ) :
        prepare_data ()
        
# =============================================================================
## Read data from DB
# =============================================================================
with ROOT.TFile.open ( testdata , 'r' ) as dbroot :
    
    logger.info ( 'Test data is fetched from DBASE "%s"' % testdata )   
    dbroot.ls_table ()
    
    hr1_data  = dbroot [ tag_data_r1  ] . clone()
    hr2_data  = dbroot [ tag_data_r2  ] . clone()
    hr3_data  = dbroot [ tag_data_r3  ] . clone()
    hr_data   = dbroot [ tag_data_r   ] . clone()

# =============================================================================
## prebook MC histograms
# =============================================================================

ir  = 72  
mc_r1  = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )
mc_r2  = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )
mc_r3  = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )
mc_r   = h1_axis ( [ rmax/ir*i for i in range ( ir + 1 ) ] )

## prepare re-weighting machinery 
maxIter = 15

## check database 
if not os.path.exists( dbname ) :
    logger.info('Create new weights DBASE') 
    with DBASE.open ( dbname , 'c' ) : ##  create new empty db 
        pass 
else :
    logger.info('Existing weights DBASE will be used') 
    with DBASE.open ( dbname , 'r' ) as db : 
        db.ls() 
        
# =============================================================================
## make reweighting iterations
from   ostap.tools.reweight           import Weight, makeWeights,  WeightingPlot, W2Data  
from   ostap.fitting.pyselectors      import Variable 
import ostap.parallel.parallel_fill

# =============================================================================
## configuration of reweighting 
weightings = (
    ## variable          address in DB    
    Weight.Var ( 'r1'  , 'r1-reweight' ) , 
    Weight.Var ( 'r2'  , 'r2-reweight' ) , 
    Weight.Var ( 'r3'  , 'r3-reweight' ) , 
    )

# =============================================================================
## variables to be used in MC-dataset 
variables  = [
    Variable ( 'r1' , 'r1-var'  , 0  , rmax ) , 
    Variable ( 'r2' , 'r2-var'  , 0  , rmax ) ,
    Variable ( 'r3' , 'r3-var'  , 0  , rmax ) ,
    ]

# =============================================================================
datatree   = ROOT.TChain ( tag_data ) ; datatree.Add ( testdata )  
title      = 'Data/target dataset'
logger.info ( '%s:\n%s' % ( title , datatree.table2 ( variables = [ 'r1' , 'r2' , 'r3' ] ,
                                                      title     = title    ,
                                                      prefix    = '# '     ) ) )

vct_data  = datatree.statVct ( 'r1,r2,r3' )
n_data    = len ( datatree ) 
## table of global statistics 
glob_stat = [ ( '#' , 'Mahalanobis' , 'Hotelling' , 'KL/DATA-MC' , 'KL/MC-DATA' , 'KL-sym' ) ] 
# =============================================================================
with timing ( 'Prepare initial MC-dataset:' , logger = logger ) :
    mctree   = ROOT.TChain ( tag_mc   ) ; mctree   .Add ( testdata )
    ## fill dataset from input MC tree
    mcds_ , _ = mctree.make_dataset ( variables = variables ,
                                      selection = '0<=r1 && r1<%s && 0<=r2 && r2<%s && 0<=r3 && r3<%s' % ( rmax , rmax, rmax ) )

# =============================================================================
## Configuration of reweighting plots 
# =============================================================================
plots  = [
    WeightingPlot ( 'r1'       , 'weight' , 'r1-reweight'  , hr1_data , mc_r1  ) , ## ignore = True ) ,  
    WeightingPlot ( 'r2'       , 'weight' , 'r2-reweight'  , hr2_data , mc_r2  ) , ## ignore = True ) ,  
    WeightingPlot ( 'r3'       , 'weight' , 'r3-reweight'  , hr3_data , mc_r3  ) , ## ignore = True ) ,  
    WeightingPlot ( 'r1,r2,r3' , 'weight' , 'r-reweight'   , hr_data  , mc_r   ) ,   
    ]

converged = False 
# =============================================================================
## start reweighting iterations:
for iter in range ( 1 , maxIter + 1 ) :

    tag = 'Reweighting iteration #%d' % iter
    logger.info ( allright ( tag ) ) 
    
    # =========================================================================
    ## 0) The weighter object
    weighter = Weight ( dbname , weightings )
    ## 1a) create new "weighted" mcdataset
    mcds = mcds_.Clone()
    
    with timing ( tag + ': add weight to MC-dataset' , logger = logger ) :
        # =========================================================================
        ## 1b) add  "weight" variable to dataset 
        mcds.add_reweighting ( weighter ,  name = 'weight' , progress = True )
        title = 'Reweighted dataset at iteration #%d' % iter 
        logger.info ( '%s:\n%s' % ( title , mcds.table2 ( variables = [ 'r1' , 'r2' , 'r3' ] ,
                                                          title     = title    ,
                                                          cuts      = 'weight' , 
                                                          prefix    = '# '     ) ) )        
        ## get MC data in a form of vector 
        vct_mc = mcds.statVct ( 'r1,r2,r3', 'weight' )
        n_mc   = int ( mcds.nEff ( 'weight' ) )  
        trow   = ( '%d'   % iter ,
                   '%.4g' % vct_mc  .mahalanobis                 ( vct_data ) ,
                   '%.4g' % Ostap.Math.hotelling ( vct_mc , n_mc , vct_data , n_data ) ,                 
                   '%.4g' % vct_mc  .asymmetric_kullback_leibler ( vct_data ) , 
                   '%.4g' % vct_data.asymmetric_kullback_leibler ( vct_mc   ) , 
                   '%.4g' % vct_data.           kullback_leibler ( vct_mc   ) )
        glob_stat.append ( trow )
        ## title = 'MC-data as vector at #%s' % iter
        ## logger.info ( '%s\n%s' % ( title , vct_mc.table ( title = title , prefix = '# ' , correlations =True ) ) ) 
        
    # =========================================================================
    ## 2) update weights

    ## pmax  = 1.0 if iter <= 4 else 0.66 + 0.34 * math.exp ( - ( iter - 4.0 ) / 4.0 )
    ##  power = lambda n : pmax if n <= 1 else pmax * ( 0.66666 / n + 0.33333 )

    the_plots = plots
    power =  1.0/3 
    
    ## if   iter < 3 : the_plots = plots [ : 1 ]
    ## elif iter < 4 : the_plots = plots [ : 4 ]
    ## else          : the_plots = plots

    ## weigth truncation: avoid very large change of weights for  single iteration 
    wtruncate = ( 0.1 , 10 ) if iter < 4 else ( 0.5 , 2.0 )  
    
    with timing ( tag + ': make actual reweighting:' , logger = logger ) :
        
        # =========================================================================
        ## 2a) the most important line: perform single iteration step  
        active , _ = makeWeights (
            mcds                   , ## what to be reweighted
            the_plots              , ## reweighting plots/setup
            dbname                 , ## DBASE with reweighting constant 
            delta      = 0.02      , ## stopping criteria
            minmax     = 0.08      , ## stopping criteria  
            power      = power     , ## tune: effective power
            wtruncate  = wtruncate , ## truncate weights 
            make_plots = True      , ## make control plots 
            tag        = tag       ) ## tag for printout

    if not active and 5 <= iter : 
        logger.info    ( allright ( 'No more iterations, converged after #%d' % iter ) )
        converged = True 
        break
    
    mcds = Ostap.MoreRooFit.delete_data ( mcds ) 
    del mcds
    
else :
    
    converged = False 
    logger.error ( "No convergency!" )
    
# ===========================================================================
title = 'Weighter object'
logger.info ( '%s:\n%s' % ( title , weighter.table ( prefix = '# ' ) ) )
# ============================================================================
## draw the convergency graphs 
graphs = weighter.graphs ()
for key in graphs : 
    with use_canvas ( "Convergency graph for '%s'" % key ) :
        graph = graphs [ key ]
        graph.draw ( 'a' )
# =============================================================================

# =============================================================================
if converged : # ==============================================================
    # =========================================================================
    with timing ( "Add weight column to initial MC-tree" , logger = logger ) : 
        mctree   = ROOT.TChain ( tag_mc   ) ; mctree   .Add ( testdata )  
        weighter = Weight ( dbname , weightings )
        mctree   = mctree.add_reweighting ( weighter ,  name = 'weight' )
        mctree   = ROOT.TChain ( tag_mc   ) ; mctree   .Add ( testdata )  

    # =======================================================================
    ## compare DATA  and MC before and after reweighting
    # ========================================================================

    datatree   = ROOT.TChain ( tag_data ) ; datatree.Add ( testdata )  
    mctree     = ROOT.TChain ( tag_mc   ) ; mctree  .Add ( testdata )  
    # ========================================================================
    title = 'Data/target dataset'
    logger.info ( '%s:\n%s' % ( title , datatree.table2 ( variables = [ 'r1' , 'r2' , 'r3' ] ,
                                                          title     = title    ,
                                                          prefix    = '# '     ) ) )
    # =============================================================================
    title = 'MC-tree before reweighting' 
    logger.info ( '%s:\n%s' % ( title , mctree.table2   ( variables = [ 'r1' , 'r2' , 'r3' ] ,
                                                        title     = title    ,
                                                          prefix    = '# '     ) ) )
    # =============================================================================
    title = 'MC-tree after reweighting' 
    logger.info ( '%s:\n%s' % ( title , mctree.table2   ( variables = [ 'r1' , 'r2' , 'r3' ] ,
                                                          title     = title    ,
                                                          cuts      = 'weight' , 
                                                          prefix    = '# '     ) ) )
    
    # =============================================================================
    vct_final = mctree.statVct ( 'r1,r2,r3' , 'weight' )
    n_mc = int ( mctree.nEff('weight')  )
    trow = ( '*' ,
             '%.4g'  % vct_final .mahalanobis                 ( vct_data  ) , 
             '%.4g'  % Ostap.Math.hotelling ( vct_final , n_mc ,  vct_data  , n_data ) ,         
             '%+.4g' % vct_final .asymmetric_kullback_leibler ( vct_data  ) , 
             '%+.4g' % vct_data  .asymmetric_kullback_leibler ( vct_final ) , 
             '%+.4g' % vct_data  .           kullback_leibler ( vct_final ) )
    glob_stat.append ( trow )

    title = 'Global DATA/MC similarity'
    table = T.table ( glob_stat , title = title , prefix = '# ' ) 
    logger.info ( '%s\n%s' % ( title , table ) ) 


# ===========================================================================
from   ostap.tools.reweight import backup_to_ROOT, restore_from_ROOT
root_file = backup_to_ROOT    ( dbname     )
new_db    = restore_from_ROOT ( root_file  )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
