#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight2.py
#  Test for 2D-reweighting machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-11
# =============================================================================
"""Test for 2D-reweighting machinery using GReweighter 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-05-10"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   ostap.core.pyrouts     import hID
from   ostap.histos.histos    import h1_axis, h2_axes 
from   ostap.math.base        import axis_range
from   ostap.utils.utils      import vrange, batch_env 
from   ostap.utils.timing     import timing
from   ostap.trees.data_utils import Data 
from   ostap.logger.colorized import attention, allright  
from   ostap.utils.cleanup    import CleanUp
import ostap.io.zipshelve     as     DBASE
import ROOT, array, random, math, os, time 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_gbreweight' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for 2D-Reweighting machinery using GBReweighter')
# ============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

testdata   = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-test-tools-gbreweight-' )
## dbname     = CleanUp.tempfile ( suffix = '.db' , prefix ='ostap-test-tools-gbreweight-'   )
## testdata = 'testdata.root'
## dbname   = 'testdb.db'

tag_data   = 'DATA2_histogram'
tag_datax  = 'DATAX_histogram'
tag_datay  = 'DATAY_histogram'
tag_mc     = 'MC2_tree'
 
## if os.path.exists ( testdata ) : os.remove ( testdata ) 
## if os.path.exists ( dbname   ) : os.remove ( dbname   )

import ostap.parallel.kisa

N1   = 1000000
N2   = 200000

N1   = 100000
N2   = 10000

xmax = 20.0
ymax = 15.0 

def prepare_data ( ) : 
    #
        
    seed =  1234567890 
    random.seed ( seed ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )   
    ## prepare "data" histograms:
    # 1) 2D histograms
    ix , iy  = 30 , 30
    hdata    = h2_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                         [ ymax/iy*i for i in range ( iy + 1 ) ] )
    # 2) non-equal binning 1D histograms for x-component    
    hxdata   = h1_axis ( [    i      for i in range ( 5  ) ] +
                         [  5+i*0.2  for i in range ( 50 ) ] +
                         [ 15+i      for i in range ( 6  ) ] )
    # 2) equal binning 1D histograms for y-component    
    hydata   = h1_axis ( [ i*0.5     for i in range ( 31 ) ] )

    assert hdata .xmax() == xmax , 'XMAX is invalid!'
    assert hdata .ymax() == ymax , 'YMAX is invalid!'
    assert hxdata.xmax() == xmax , 'XMAX is invalid!'
    assert hydata.xmax() == ymax , 'XMAX is invalid!'

    with ROOT.TFile.Open ( testdata ,'recreate') as mc_file:
        mc_file.cd() 
        
        datatree  = ROOT.TTree ( 'DATA_tree', 'data-tree' )
        datatree .SetDirectory ( mc_file ) 
        
        xvar = array.array  ( 'f', [0])
        yvar = array.array  ( 'f', [0])
        datatree.Branch ( 'x' , xvar , 'x/F' )
        datatree.Branch ( 'y' , yvar , 'y/F' )
    
        for i in  range ( N1 ) :

            x , y = -1, -1

            while not 0 <= x < xmax or not 0 <= y < ymax :
                
                v1 = random.gauss ( 0 , 3 )
                v2 = random.gauss ( 0 , 2 )
                
                x  =  5 + v1 + v2   
                y  = 12 + v1 - v2
                
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
 
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()

            
        for i in  range ( N1 ) :

            x , y = -1 , -1

            while not 0 <= x < xmax or not 0 <= y < ymax :
                
                v1 = random.gauss ( 0 , 3 )
                v2 = random.gauss ( 0 , 2 )

                x  = 15 + v1 
                y  =  4 + v2    
                
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
 
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()

            
        for i in  range ( N1 ) :

            x , y = -1 , -1

            while not 0 <= x < xmax or not 0 <= y < ymax :
                
                v1 = random.gauss ( 0 , 3 )
                v2 = random.gauss ( 0 , 2 )

                x  = 15 + v1 - v2  
                y  = 12 + v1 + v2     
                
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
 
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()


        for i in range ( 0 ,  4 * N1 ) :
            
            x = random.uniform ( 0 , xmax ) 
            y = random.uniform ( 0 , ymax )
            
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
            
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()
            
        datatree.Write()
        
        ## write the histogram 
        mc_file [ tag_data  ] = hdata
        mc_file [ tag_datax ] = hxdata
        mc_file [ tag_datay ] = hydata
        
        mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
        mctree .SetDirectory ( mc_file ) 
        
        xvar = array.array  ( 'f', [ 0.0 ] )
        yvar = array.array  ( 'f', [ 0.0 ] )
        mctree.Branch ( 'x' , xvar , 'x/F' )
        mctree.Branch ( 'y' , yvar , 'y/F' )

        for i in  range ( 2 * N2 ) :

            xv = random.uniform ( 0 , xmax ) 
            yv = random.uniform ( 0 , ymax ) 

            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            
            mctree.Fill()
            
        fx = lambda x : ( x/10.0-1 ) ** 2
        fy = lambda y : ( y/ 7.5-1 ) ** 2
        

        for i in  range ( N2 ) :
            
            while True : 
                xv = random.uniform ( 0 , xmax )
                if random.uniform   ( 0 ,  1 ) < fx ( xv ) : break
                
            while True : 
                yv = random.uniform ( 0 , ymax )
                if random.uniform   ( 0 ,  1 ) < fy ( yv ) : break 
                
            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            
            mctree.Fill()

            
        mctree.Write()
        mc_file.ls()


# =============================================================================
def test_gbreweight() :


    logger = getLogger("test_gbreweight")
    
    try :
        
        from ostap.tools.reweighter import Reweighter
        rw = Reweighter ()
        
    except ImportError :
    
        logger.error ('GBReweighter is not available!')
        return 


    if not os.path.exists( testdata ) :
        with timing ( "Prepare input data" , logger = logger ) :
            prepare_data ()
            
    # =========================================================================
    ## Input data/mc samples 
    # =========================================================================
    data = Data ( testdata , 'DATA_tree' )
    mc   = Data ( testdata , tag_mc      ) 
    
    ddata , wdata = data.chain.slice ( 'x , y' , transpose = True )
    dmc   , wmc   = mc  .chain.slice ( 'x , y' , transpose = True )

        
    ## train BDT
    rw.reweight ( original = dmc , target = ddata ) 
    
    ## new weights 
    wnew = rw.weight ( original = dmc )
    ## mc.chain.add_new_branch ( wnew , name = 'w')
    mc.chain.add_new_buffer ( 'w' , wnew )
    
    ## reload data 
    data = Data ( testdata , 'DATA_tree' )
    mc   = Data ( testdata , tag_mc      ) 
    
    wsum = mc.chain.sumVar ( 'w' )
    wvar = '%d*w/%s' % ( len ( data.chain ) , wsum.value() )
    
    nn   = '%s' % ( len(data.chain)*1.0 / len ( mc.chain) ) 
    
    for phi in vrange ( 0 , 2*math.pi , 10 ) :
        
        dvar = '%.5f*x+%.5f*y' % ( math.cos ( phi ) , math.sin ( phi ) )
        
        mn , mx = data.chain.statVar ( dvar ).minmax() 
        h1 = ROOT.TH1D ( hID() , '' , 100 , *axis_range ( mn , mx , delta = 0.05 ) ) 
        h2 = h1.clone()
        h3 = h1.clone()
        
        data.chain.project ( h1 , dvar        ) ## data 
        mc  .chain.project ( h2 , dvar , nn   ) ## original (non-weighted) MC 
        mc  .chain.project ( h3 , dvar , wvar ) ## weighted  MC 
        
        mn , mx = h1.minmax()
        mn , mx = axis_range ( 0 , mx , delta = 0.7 )
        h1.SetMaximum ( mx ) 
        h1.blue()
        h2.green()
        h3.red() 
    
        h1.draw ( '' )
        h2.draw ( 'same hist')
        h3.draw ( 'same')
        
        time.sleep ( 2 ) 
        
# =============================================================================
if '__main__' == __name__ :
    

    test_gbreweight() 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
