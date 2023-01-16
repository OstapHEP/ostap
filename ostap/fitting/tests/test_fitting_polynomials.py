#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_models.py
# Test module for varipsu polynomial PDFs
# ============================================================================= 
""" Test module for varipsu polynomial PDFs
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, hID , rooSilent, Ostap 
from   ostap.utils.timing       import timing
from   builtins                 import range
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.fitting.background import make_bkg
from   ostap.logger.colorized   import attention
from   ostap.utils.utils        import vrange
import ostap.histos.graphs 
import ostap.logger.table       as     T
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_polynomials' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

## observable 
x       = ROOT.RooRealVar ( 'x' , 'Some test variable' , 0, 10)
varset  = ROOT.RooArgSet  ( x )
dataset = ROOT.RooDataSet ( dsID() , 'Test data' , varset )

NS      = 10000
NB      = 10000 

for i in range ( NS ) :
    v = random.gammavariate ( 3 , 1/0.8 )
    while not v in x :
        v = random.gammavariate ( 3 , 1/0.8 )
    x.value = v
    dataset.add ( varset )
    
for i in range ( NB ) :
    v = random.uniform ( 0 , 10 ) 
    x.value = v
    dataset.add ( varset )    
    
logger.info ('DATASET\n%s' % dataset )


models  = {} 
results = {}
plots   = {}
graphs  = set()
    

Nmin , Nmax = 2 , 8

# =============================================================================
## Test Bernstein positive polynomials PolyPos_pdf
def test_poly_positive() :
    """Test Bernstein positive polynomials PolyPos_pdf"""

    logger = getLogger( "test_poly_positive" )

    logger.info ( "Test Bernstein positive polynomials PolyPos_pdf" )

    for n in range ( Nmin , Nmax ) :
        
        model = Models.PolyPos_pdf ( 'B%s' % n  , xvar  = x  , power = n  )

        title = 'PolyPos_pdf(n=%d)' % n 
        model.fitTo ( dataset , silent = True )
        with use_canvas ( title , wait = 1 ) :
            result , plot = model.fitTo ( dataset , silent = True , draw = True , nbins = 100 )
            logger.info ( '%s\n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
        models  [ title ] = model  
        results [ title ] = result 
        plots   [ title ] = plot  

# =============================================================================
## Test Karlin-Shapley positive polynomials KarlinShapley_pdf
def test_poly_karlin_shapley () :
    """Test Karlin-Shapley positive polynomials KarlinShapley_pdf"""

    logger = getLogger( "test_poly_karlin_shapley" )

    logger.info ( "Test Karlin-Shapley positive polynomials KarlinShapley_pdf" ) 
    
    for n in range ( Nmin , Nmax ) :
        
        model = Models.KarlinShapley_pdf ( 'KSh%s' % n  , xvar  = x  , power = n  , xmin = 0 , xmax = 10 )
        
        title = 'KarlinShapley_pdf(n=%d)' % n 
        model.fitTo ( dataset , silent = True , refit = 5 )
        with use_canvas ( title , wait = 1 ) :
            result , plot = model.fitTo ( dataset , silent = True , draw = True , nbins = 100 )
            logger.info ( '%s\n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
            
        models  [ title ] = model  
        results [ title ] = result 
        plots   [ title ] = plot  


# =============================================================================
## Test Karlin-Studden positive polynomials KarlinStudden_pdf
def test_poly_karlin_studden () :
    """Test Karlin-Shapley positive polynomials KarlinStudden_pdf"""

    logger = getLogger( "test_poly_karlin_studden" )

    logger.info ( "Test Karlin-Studden positive polynomials KarlinStudden_pdf" ) 
    
    for n in range ( Nmin , Nmax ) :
        
        model = Models.KarlinStudden_pdf ( 'KSt%s' % n  , xvar  = x  , power = n  , xmin = 0 , scale = 10 )
        
        title = 'KarlinStudden_pdf(n=%d)' % n 
        model.fitTo ( dataset , silent = True , refit = 5 )
        with use_canvas ( title , wait = 1 ) :
            result , plot = model.fitTo ( dataset , silent = True , draw = True , nbins = 100 )
            logger.info ( '%s\n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            model.pdf.setPars()
        models  [ title ] = model  
        results [ title ] = result 
        plots   [ title ] = plot  


# =============================================================================
## Test native RooFit polynomials RooPoly_pdf
def test_poly_roopoly () :
    """Test native RooFit polynomials RooPoly_pdf"""

    logger = getLogger( "test_poly_roopoly" )

    logger.info ( "Test native RooFit polynomials RooPoly_pdf" ) 
    
    for n in range ( 2 , 4 ) :
        
        model = Models.RooPoly_pdf ( 'R%s' % n  , xvar  = x  , power = n  )
        
        title = 'RooPoly_pdf(n=%d)' % n 
        model.fitTo ( dataset , silent = True , refit = 5 )
        with use_canvas ( title , wait = 1 ) :
            result , plot = model.fitTo ( dataset , silent = True , draw = True , nbins = 100 )
            logger.info ( '%s\n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
        models  [ title ] = model  
        results [ title ] = result 
        plots   [ title ] = plot  


# =============================================================================
## Test RooFit Chebyshev polynomials RooCheb_pdf
def test_poly_roocheb () :
    """Test RooFit Chebyshev polynomials RooCheb_pdf"""

    logger = getLogger( "test_poly_roocgeb" )

    logger.info ( "Test RooFit Chebyshev polynomials RooCheb_pdf" ) 
    
    for n in range ( Nmin , Nmax ) :
        
        model = Models.RooCheb_pdf ( 'C%s' % n  , xvar  = x  , power = n  )
        
        title = 'RooCheb_pdf(n=%d)' % n 
        model.fitTo ( dataset , silent = True , refit = 5 )
        with use_canvas ( title , wait = 1 ) :
            result , plot = model.fitTo ( dataset , silent = True , draw = True , nbins = 100 )
            logger.info ( '%s\n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
        models  [ title ] = model  
        results [ title ] = result 
        plots   [ title ] = plot  


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

    print ( 'LEN', [ len(g) for g in grs ] )
    
    ## reset 
    for  i in range ( npars ) : obj.setPar ( i , 0 )
    
    return tuple ( grs ) 
        
# =============================================================================
## Evolution of Karlin-Shapley t-roots 
def test_karlin_shapley_troots () :
    """Evolution of Karlin-Shapley t-roots"""
    
    logger = getLogger ( 'test_karlin_shapley_troots' ) 

    low  = 0
    high = 2 * math.pi 
    
    h1 = ROOT.TH1F ( hID() , '' , 1 , 0 , 2*math.pi ) 
    h1.SetMinimum  ( 0 )
    h1.SetMaximum  ( 1 )
    
    for n in range ( 2 , 6 ) :
        
        ks    = Ostap.Math.KarlinShapley ( n , 0 , 1 )
        
        npars = ks.npars()
        
        for  phase in range ( 1 , npars ) :
            
            title = 'KarlinShapley(n=%s) t-roots for par %d' % ( n , phase )
            with use_canvas ( title , wait = 2 ) :
                
                h1.draw()
                grs = make_graphs ( ks , phase )
                
                for i,g in enumerate ( grs , start = 1 ) :
                    g.SetLineColor   ( i )
                    g.SetMarkerColor ( i )
                    g.draw('c')
                    
                graphs.add ( grs ) 

                
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        db['x'        ] = x
        db['varset'   ] = varset
        db['dataset'  ] = dataset        
        for m in models :
            db['model:%s'      % m ] = models [ m ] 
            db['models_pdf:%s' % m ] = models [ m ].pdf            
        db['models'   ] = models
        for r in results :
            db ['result:%s' % r ] = results [ r ] 
        db['results'   ] = results
        for p in plots :
            db [ 'plot:%s'  % p ] = plots   [ p ] 
        db['plots'     ] = plots 
        db['graphs'    ] = graphs 
        db.ls() 


# =============================================================================
if '__main__' == __name__ :
    
    with timing ( 'Bernstein positive polynomials: PolyPos_pdf' , logger ) : 
        test_poly_positive() 
        
    with timing ( 'KarlinShapley positive polynomials: KarlinShapley_pdf' , logger ) : 
        test_poly_karlin_shapley () 
    
    with timing ( 'KarlinStudden positive polynomials: KarlinStudden_pdf' , logger ) : 
        test_poly_karlin_studden () 

    with timing ( 'Native RooFit polynomials: RooPoly_pdf' , logger ) : 
        test_poly_roopoly() 

    with timing ( 'RooFit Chebyshev polynomials: RooCheb_pdf' , logger ) : 
        test_poly_roocheb() 

    test_karlin_shapley_troots() 
       
    test_db()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
