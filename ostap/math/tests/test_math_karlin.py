#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_karlin.py
#  Test module for Karlin-Shapley polinomials
#  @see Ostap::Math::KarlinShapley 
#  @see Ostap::Math::KarlinStudden
# ============================================================================= 
""" Test module for Karlin-Shapley polinomials
- see Ostap::Math::KarlinShapley 
- see Ostap::Math::KarlinStudden
"""
# ============================================================================= 
from   ostap.core.core        import Ostap, hID 
from   ostap.utils.timing     import timing
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.ranges     import vrange
from   ostap.utils.root_utils import batch_env 
import ostap.histos.graphs   
import ostap.math.models 
import ostap.math.karlin 
import ROOT, math 
# ============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_karlin' ) 
else                       : logger = getLogger ( __name__                 )
# ============================================================================= 
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

functions = set()
graphs    = set()
    
# =============================================================================
def make_graphs ( obj , index ) :
    
    npars  = obj.npars()

    ## reset 
    for  i in range ( npars ) : obj.setPar ( i , 0 )
    
    troots = obj.troots()
    NP     = 200 

    grs    = []
    while len ( grs ) < len  ( troots ) : grs.append ( ROOT.TGraph() )
    

    for phi in vrange ( 0 , 2 * math.pi , NP ) :
        
        obj.setPar ( index , phi )
        
        for i, g in enumerate ( grs ) :
            
            v = obj.troots()[i]
            g.append ( phi , v ) 

    ## reset 
    for  i in range ( npars ) : obj.setPar ( i , 0 )
    
    return tuple ( grs ) 
        
# =============================================================================
## Evolution of Karlin-Shapley t-roots 
def test_karlin_shapley_troots () :
    """ Evolution of Karlin-Shapley t-roots"""
    
    logger = getLogger ( 'test_karlin_shapley_troots' ) 
    logger.info ( 'Show evolution of Karlin-Shapley t-roots' )
    
    low  = 0
    high = 2 * math.pi 
    
    h1 = ROOT.TH1F ( hID() , '' , 1 , 0 , 2*math.pi ) 
    h1.SetMinimum  ( -0.1 )
    h1.SetMaximum  (  1.1 )
    
    for n in range ( 2 , 10 ) :
        
        ks    = Ostap.Math.KarlinShapley ( n , 0 , 1 )
        
        npars = ks.npars()
        
        logger.info ( 'Karlin-Shapley(n=%s) t-roots %s ' %  ( n , ks.troots() ) ) 

        for  phase in range ( 1 , npars ) :
            
            title = 'KarlinShapley(n=%s) t-roots for par %d' % ( n , phase )
            with use_canvas ( title , wait = 1 ) :
                
                h1.draw()
                grs = make_graphs ( ks , phase )
                
                for i,g in enumerate ( grs , start = 1 ) :
                    g.SetLineColor   ( i )
                    g.SetMarkerColor ( i )
                    g.SetLineWidth   ( 2 )
                    g.draw ('c')
                    
                graphs.add ( grs ) 
                
        functions.add ( ks )

# =============================================================================
## Evolution of Karlin-Studden t-roots 
def test_karlin_studden_troots () :
    """Evolution of Karlin-Studden t-roots"""
    
    logger = getLogger ( 'test_karlin_studden_troots' ) 
    logger.info ( 'Show evolution of Karlin-Studden t-roots' )
    
    low  = 0
    high = 2 * math.pi 
    
    h1 = ROOT.TH1F ( hID() , '' , 1 , 0 , 2*math.pi ) 
    h1.SetMinimum  ( -0.1 )
    h1.SetMaximum  (  1.1 )
    
    for n in range ( 2 , 10 ) :
        
        ks    = Ostap.Math.KarlinStudden ( n , 0 , 1 )
        
        npars = ks.npars()
        
        logger.info ( 'Karlin-Studden(n=%s) t-roots %s ' %  ( n , ks.troots() ) ) 
        logger.info ( 'Karlin-Studden(n=%s) z-roots %s ' %  ( n , ks.zroots() ) ) 

        for  phase in range ( 1 , npars ) :
            
            title = 'KarlinStudden(n=%s) t-roots for par %d' % ( n , phase )
            with use_canvas ( title , wait = 1 ) :
                
                h1.draw()
                grs = make_graphs ( ks , phase )
                
                for i,g in enumerate ( grs , start = 1 ) :
                    g.SetLineColor   ( i )
                    g.SetMarkerColor ( i )
                    g.SetLineWidth   ( 2 )
                    g.draw ('c')
                    
                graphs.add ( grs ) 
                
        functions.add ( ks )
        
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :
    
    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for i , f in enumerate ( functions , start = 1 ) :
            db [ '%03d:%s' % ( i , f.__class__.__name__ ) ] = f 
        db ['functions'   ] = functions 
        db ['graphs'      ] = graphs 
        db.ls() 
            
# =============================================================================
if '__main__' == __name__ :

    with timing ( 'Evolution of Karlin-Shapley t-roots' , logger ) : 
        test_karlin_shapley_troots ()
        
    with timing ( 'Evolution of Karlin-Studden t-roots' , logger ) : 
        test_karlin_studden_troots ()
    
    ## check finally that everything is serializeable:
    with timing ('test_db' , logger ) :
        test_db ()
        

# =============================================================================
##                                                                      The END 
# ============================================================================= 
