#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_bwi.py
# test for some Breit-Wigner models 
# ============================================================================= 
""" Test module for Breit-Wigner wirth interferenc e 
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models    as     Models 
from   ostap.core.core         import Ostap
import ostap.io.zipshelve      as     DBASE
import ostap.logger.table      as     T 
from   ostap.utils.timing      import timing
from   ostap.utils.utils       import wait
from   ostap.plotting.canvas   import use_canvas
from   ostap.fitting.variables import SETVAR 
from   ostap.utils.utils       import vrange 
from   builtins                import range
import ROOT, time, math  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_bwi' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================


models  = set() 
results = [] 

# =============================================================================
def test_fitting_bwi1 () :
    
    
    logger = getLogger ( 'test_fitting_bwi1' )
    logger.info ( 'Test simple "Breit-Wigner with interference" model' )
    
    m1    = 1.0 
    m2    = 1.0
    m0    = 5.0
    g0    = 1.0 
    
    ## create BreitWigner function 
    breitwigner = Ostap.Math.BreitWigner ( m0 , g0 , m1 , m2  )
    phasespace  = Ostap.Math.PhaseSpace2 ( m1 , m2 )
    
    ## mass observable 
    mass        = ROOT.RooRealVar ( 'mass' , 'mass-observable' , 0 , 15 )  
    
    ## create BW function
    bw          = Models.BreitWigner_pdf ( 'BW'  ,
                                           m0     = ( m0 , m0 - 1 , m0 + 1 ) , 
                                           gamma  = ( g0 , g0/2   , g0 * 2 ) , 
                                           breitwigner = breitwigner , xvar = mass )
    
    ## create BWI function
    bwi         = Models.BWI_pdf        ( 'BWI' ,
                                          m0          = bw.m0    ,
                                          gamma       = bw.gamma ,                                     
                                          breitwigner = bw.breitwigner , xvar = mass )
    bwi.magnitude.phis =  0.05 , 
    bwi.phase    .phis = -0.7  , 
    bwi.scale1         =  0.3
    bwi.scale2         = 11.0 
    
    ## model without interference
    bkg1   = Models.PSLeftExpoPol_pdf ( 'B1' , xvar = mass       , 
                                        phasespace  = phasespace ,
                                        scale       = ROOT.RooFit.RooConst ( 1 ) ,
                                        tau         = ROOT.RooFit.RooConst ( 0 ) ,
                                        power       = 5 )
    model1   = Models.Fit1D ( signal = bw , background = bkg1 , suffix = '1' ) 
    model1.S = 100
    model1.B =  50
    
    ## model with interference
    bkg2   = Models.PSLeftExpoPol_pdf ( 'B2' , xvar = mass       ,
                                        phasespace  = phasespace ,
                                        scale       = ROOT.RooFit.RooConst ( 1 ) ,
                                        tau         = ROOT.RooFit.RooConst ( 0 ) ,
                                        power       = 1 )
    model2   = Models.Fit1D ( signal = bwi , background = bkg2 , suffix = '2' ) 
    model2.S = 630
    model2.B = 330
    
    ds1 = model1.generate ( 1000 )
    ds2 = model2.generate ( 1000 )
    
    models.add  ( model1 ) 
    models.add  ( model2 ) 
    models.add  ( model1.signal ) 
    models.add  ( model2.signal ) 
    models.add  ( model1.background ) 
    models.add  ( model2.background ) 
    

    # ==============================================================================
    with use_canvas ( 'Fit ds1 with model1' , wait = 1 ) : 
        r11 , f = model1.fitTo ( ds1 , silent = True )
        r11 , f = model1.fitTo ( ds1 , silent = True , draw = True )
        title = 'Fit ds1 with model1'
        logger.info ( '%s\n%s' % ( title , r11.table ( title , prefix = '# ' ) ) )
        results.append ( r11 ) 
        
    with use_canvas ( 'Fit ds2 with model2' , wait = 1 ) : 
        r22 , f = model2.fitTo ( ds2 , silent = True )
        r22 , f = model2.fitTo ( ds2 , silent = True , draw = True )
        title = 'Fit ds2 with model2'
        logger.info ( '%s\n%s' % ( title , r22.table ( title , prefix = '# ' ) ) ) 
        results.append ( r22 ) 
        
    with use_canvas ( 'Fit ds1 with model2' , wait = 1 ) : 
        r12 , f = model2.fitTo ( ds1 , silent = True )
        r12 , f = model2.fitTo ( ds1 , silent = True , draw = True )
        title = 'Fit ds1 with model2'
        logger.info ( '%s\n%s' % ( title , r12.table ( title , prefix = '# ' ) ) ) 
        results.append ( r12 ) 
        
    with use_canvas ( 'Fit ds2 with model1' , wait = 1 ) : 
        r21 , f = model1.fitTo ( ds2 , silent = True )
        r21 , f = model1.fitTo ( ds2 , silent = True , draw = True )
        title = 'Fit ds2 with model1'
        logger.info ( '%s\n%s' % ( title , r21.table ( title , prefix = '# ' ) ) ) 
        results.append ( r21 ) 


data  = {} 
rs    = {}
plots = {}

# =============================================================================
def test_fitting_bwi2 () :
## if 1 < 2 :
    
    logger = getLogger ( 'test_fitting_bwi2' )
    logger.info ( 'Test simple "Breit-Wigner with interference" model' )

    GeV = 1.0
    MeV = 0.001 * GeV
    
    xmin, xmax = 3.85 * GeV , 3.90 * GeV  
    mass = ROOT.RooRealVar( 'mass' , '' , xmin , xmax )

    m0   = 3.872 * GeV 
    g0   = 0.9   * MeV 

    m1   = 3.1 * GeV
    m2   = 0.5 * GeV 

    breit_wigner = Ostap.Math.BreitWigner ( m0 , g0 , m1 , m2 )

    ## create BW function
    bw   = Models.BreitWigner_pdf ( 'BW2'  ,
                                    xvar        = mass , 
                                    m0          = ( m0 , m0 - 3 * MeV , m0 + 3 * MeV  ) , 
                                    gamma       = ( g0 , g0/2   , g0 * 2 ) , 
                                    breitwigner = breit_wigner )
    
    bw.m0   .fix()
    bw.gamma.fix()
    
    b     = ROOT.RooRealVar( 'b'     , 'magnitude of the coherent background' , 1 , 0 , 200 ) 
    phi_b = ROOT.RooRealVar( 'phi_b' , 'phase of the coherent background'     , 0 , -4 * math.pi , 4 * math.pi )
    
    ## create BWI function
    bwi   = Models.BWI_pdf        ( 'BWI2'                        ,
                                    xvar        = mass            , 
                                    m0          = bw.m0           ,
                                    gamma       = bw.gamma        ,                                     
                                    breitwigner = bw.breitwigner  ,
                                    magnitude   = b               ,
                                    phase       = phi_b           )
    
    ## apply mass-resolution
    resolution = Models.ResoGauss ( 'Gauss'           ,
                                    xvar  = mass      ,
                                    sigma = 2.3 * MeV )
    
    cnv_conf = { 'nbins' : 10000 , 'buffer' : 0.3 } 
    signal   = Models.Convolution_pdf ( name       = 'X'        ,
                                        pdf        = bwi        ,
                                        resolution = resolution , **cnv_conf)
    
    model    = Models.Fit1D ( signal = signal , background = -1 )
    model.S  = 1
    model.B  = 1 

    models.add ( bw     )
    models.add ( bwi    )
    models.add ( signal )
    models.add ( model  )

    NN       = 1000 
    for vb in vrange ( 1 , 100 , 4 ) :
        with SETVAR ( b ) :
            b.setVal ( vb )
            for vphi in vrange ( 0 , math.pi , 4 ) :
                with SETVAR ( b ) :                    
                    phi_b.setVal ( phi_b ) 
                    ds = model.generate ( 1000 )
                    data [ ( vb, vphi ) ] = ds


    rows = [ ( 'b[gen]' ,  'phi/pi[gen]' , 'b/fit' , 'phi/pi[fit]' , 'status' , 'cov2' , 'fS' , 'fB' , 'fI' , 'nSfit' , 'nBfit' , 'nS*' , 'nB*' ) ]
    
    for key in data :
        vb , vphi = key
        ds        = data [ key ]

        with SETVAR ( phi_b ) , SETVAR ( b ) :
            
            b    .setVal ( vb   )
            phi_b.setVal ( vphi )
            with use_canvas ( 'tets_fititng_bwi2: b=%.2f phi/pi=%+.2f' % ( vb , vphi / math.pi ) ) : 
                r , f    = model.fitTo ( ds , silent = True , draw =True , nbins = 50 )
            
            iS,iB,iI = bwi.ffs() 
            
            nS       = r.S * 1
            nB       = r.B * 1
            
            b_fit    = r.b * 1.0 
            phi_fit  = r.phi_b / math.pi 
            
            nSt      = nS * ( iS + iI )
            nBt      = nB + nS * iB

            
            row      = '%.1f'  % vb                           , \
                       '%+.3f' % ( vphi / math.pi )           , \
                       b_fit  .toString ( '%+.2f +/- %-.2f' ) , \
                       phi_fit.toString ( '%+.2f +/- %-.2f' ) , \
                       '%s' % r.status () , \
                       '%s' % r.covQual() , \
                       '%+.2f' % iS , \
                       '%+.2f' % iB , \
                       '%+.2f' % iI , \
                       nS .toString ( '%.0f +/- %-.0f' ) , \
                       nB .toString ( '%.0f +/- %-.0f' ) , \
                       nSt.toString ( '%.0f +/- %-.0f' ) , \
                       nBt.toString ( '%.0f +/- %-.0f' ) 
                               
              
        rows.append ( row ) 
        
        plots [ key ] = f
        rs    [ key ] = r        
        
        print ( key , r ) 

    title = 'Fits witn interference' 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rrcccc') 
    logger.info ( '%s\n%s' % ( title , table ) )

    
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for m in models :
            db['model:' + m.name ] = m
            db['roo:%s' % m.name ] = m.pdf
        db['models'   ] = models
        for i, r in enumerate ( results ) :
            db ['result:%s' % r.name ] = r
        db['results'   ] = results

        for key in data :
            db ['d_%s' % str  (key ) ] = data  [ key ]
            db ['r_%s' % str  (key ) ] = rs    [ key ]
            db ['p_%s' % str  (key ) ] = plots [ key ]
            
        db ['d'] = data 
        db ['r'] = rs
        db ['p'] = plots 
        
        db.ls() 


# =============================================================================
if '__main__' == __name__ :

    with timing ('Test BWI-1' , logger ) :
        test_fitting_bwi1 ()

    with timing ('Test BWI-2' , logger ) :
        test_fitting_bwi2 ()

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# =============================================================================
