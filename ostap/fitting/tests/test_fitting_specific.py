#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_speficis.py
# Test module for soem specific fit models 
# ============================================================================= 
""" Test module for soem specific fit models 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
from   ostap.fitting.specific   import * 
import ostap.io.zipshelve       as     DBASE
import ostap.fitting.models     as     Models
import ostap.fitting.roofit 
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_specific' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

mB  = ROOT.RooRealVar ( 'mB'  , 'some beauty mass' , 5.0 , 5.5 )
mBc = ROOT.RooRealVar ( 'mBc' , 'some beauty mass' , 6.1 , 6.5 )

mC  = ROOT.RooRealVar ( 'mC' , 'some charm   mass' , 1.6 ,  2.1 )
mY  = ROOT.RooRealVar ( 'mY' , 'some upsilon mass' , 8.5 , 12.5 )


results = {}
plots   = {}
models  = {}

# =============================================================================
## test Bd_pdf 
def test_Bd () :
    """Test `Bd_pdf`
    - see `Bd_pdf`
    """

    logger = getLogger ( 'test_Bd' )
    logger.info ( 'Test Bd_pdf' )
    
    pdf  = Models.Fit1D ( signal      = Bd_pdf  ( xvar = mB ) ,
                          background  = 'convex decreasing 2'  )
    pdf.S = 1000
    pdf.B = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'Bd_pdf' , wait = 2 ) :
        
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
        
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )

        results [ 'test_Bd'] = r
        plots   [ 'test_Bd'] = f
        models  [ 'test_Bd'] = pdf 
    
# =============================================================================
## test Bds_pdf 
def test_Bs () :
    """Test `Bs_pdf`
    - see `Bs_pdf`
    """

    logger = getLogger ( 'test_Bs' )
    logger.info ( 'Test Bs_pdf' )

    pdf  = Models.Fit1D ( signal      = Bs_pdf  ( xvar = mB ) ,
                          background  = 'convex decreasing 2'  )
    pdf.S = 1000
    pdf.B = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'Bs_pdf' , wait = 2 ) :
        
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
        
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )

        results [ 'test_Bs'] = r
        plots   [ 'test_Bs'] = f
        models  [ 'test_Bs'] = pdf 

# =============================================================================
## test Bu_pdf 
def test_Bu () :
    """Test `Bu_pdf`
    - see `Bu_pdf`
    """

    logger = getLogger ( 'test_Bu' )
    logger.info ( 'Test Bu_pdf' )

    pdf  = Models.Fit1D ( signal      = Bu_pdf  ( xvar = mB ) ,
                          background  = 'convex decreasing 2'  )
    pdf.S = 1000
    pdf.B = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'Bu_pdf' , wait = 2 ) :
        
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
        
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )

        results [ 'test_Bu'] = r
        plots   [ 'test_Bu'] = f
        models  [ 'test_Bu'] = pdf 

# =============================================================================
## test Bc_pdf 
def test_Bc () :
    """Test `Bc_pdf`
    - see `Bc_pdf`
    """

    logger = getLogger ( 'test_Bc' )
    logger.info ( 'Test Bc_pdf' )

    pdf  = Models.Fit1D ( signal      = Bc_pdf  ( xvar = mBc ) ,
                          background  = 'convex decreasing 2'  )
    pdf.S = 1000
    pdf.B = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'Bc_pdf' , wait = 2 ) :
        
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )

        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )

        results [ 'test_Bc'] = r
        plots   [ 'test_Bc'] = f
        models  [ 'test_Bc'] = pdf 
        
# =============================================================================
## test D0_pdf 
def test_D0 () :
    """Test `D0_pdf`
    - see `D0_pdf`
    """

    logger = getLogger ( 'test_D0' )
    logger.info ( 'Test D0_pdf' )

    pdf  = Models.Fit1D ( signal      = D0_pdf  ( xvar = mC  ) ,
                          background  = 'convex decreasing 2'  )
    pdf.S = 1000
    pdf.B = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'D0_pdf' , wait = 2 ) :
        
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
        
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )

        results [ 'test_D0'] = r
        plots   [ 'test_D0'] = f
        models  [ 'test_D0'] = pdf 
        

# =============================================================================
## test Dp_pdf 
def test_Dp () :
    """Test `Dp_pdf`
    - see `Dp_pdf`
    """

    logger = getLogger ( 'test_Dp' )
    logger.info ( 'Test Dp_pdf' )

    pdf  = Models.Fit1D ( signal      = Dp_pdf  ( xvar = mC  ) ,
                          background  = 'convex decreasing 2'  )
    pdf.S = 1000
    pdf.B = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'Dp_pdf' , wait = 2 ) :
        
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
        
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )

        results [ 'test_Dp'] = r
        plots   [ 'test_Dp'] = f
        models  [ 'test_Dp'] = pdf 
        
# =============================================================================
## test Ds_pdf 
def test_Ds () :
    """Test `Ds_pdf`
    - see `Ds_pdf`
    """

    logger = getLogger ( 'test_Ds' )
    logger.info ( 'Test Ds_pdf' )

    pdf  = Models.Fit1D ( signal      = Ds_pdf  ( xvar = mC  ) ,
                          background  = 'convex decreasing 2'  )
    pdf.S = 1000
    pdf.B = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'Ds_pdf' , wait = 2 ) :
        
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )

        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )

        results [ 'test_Ds'] = r
        plots   [ 'test_Ds'] = f
        models  [ 'test_Ds'] = pdf 
        
# =============================================================================
## test DpDs_pdf 
def test_DpDs () :
    """Test `DpDs_pdf`
    - see `DpDs_pdf`
    """

    logger = getLogger ( 'test_DpDs' )
    logger.info ( 'Test DpDs_pdf' )

    pdf  = DpDs_pdf  ( xvar       = mC ,  name = 'DpDs'     , 
                       background  = 'convex decreasing 2'  )
    pdf.NDp = 700
    pdf.NDs = 300
    pdf.B   = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'DpDs_pdf' , wait = 2 ) :
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
                
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )
        
        results [ 'test_DpDs'] = r
        plots   [ 'test_DpDs'] = f
        models  [ 'test_DpDs'] = pdf 


# =============================================================================
## test BdBs_pdf 
def test_BdBs () :
    """Test `BdBs_pdf`
    - see `BdBs_pdf`
    """

    logger = getLogger ( 'test_BdBs' )
    logger.info ( 'Test BdBs_pdf' )

    pdf  = BdBs_pdf  ( xvar       = mB ,  name = 'BdBs'     , 
                       background  = 'convex decreasing 2'  )
    pdf.NBd = 700
    pdf.NBs = 300
    pdf.B   = 200

    ds    = pdf.generate ( 1200 )

    with use_canvas ( 'BdBs_pdf' , wait = 2 ) :
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
                
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )
        
        results [ 'test_BdBs'] = r
        plots   [ 'test_BsBs'] = f
        models  [ 'test_BdBs'] = pdf 

        

# =============================================================================
## test Manca_pdf 
def test_Manca () :
    """Test `Manca_pdf`
    - see `Manca_pdf`
    """

    logger = getLogger ( 'test_Manca' )
    logger.info ( 'Test Manca_pdf' )

    pdf  = Manca_pdf  ( xvar       = mY       ,
                        name       = 'Manca'  , 
                        background = 'expo+'  ,
                        m1s        = 9.460    , 
                        s1s        = 0.040    )
    
    pdf.N1S = 700
    pdf.N2S = 300 
    pdf.N3S = 100 
    pdf.B   = 500

    ds    = pdf.generate ( 1600 )

    with use_canvas ( 'Manca_pdf' , wait = 2 ) :
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
                
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )
        
        results [ 'test_Manca'] = r
        plots   [ 'test_Manca'] = f
        models  [ 'test_Manca'] = pdf 

# =============================================================================
## test Manca2_pdf 
def test_Manca2 () :
    """Test `Manca2_pdf`
    - see `Manca2_pdf`
    """

    logger = getLogger ( 'test_Manca2' )
    logger.info ( 'Test Manca2_pdf' )

    pdf  = Manca2_pdf  ( xvar       = mY        ,
                         name       = 'Manca2'  , 
                         background = 'expo+'   ,
                         m1s        = 9.460     , 
                         s1s        = 0.040     )
    
    pdf.N1S = 700
    pdf.N2S = 300 
    pdf.N3S = 100 
    pdf.B   = 500

    ds    = pdf.generate ( 1600 )

    with use_canvas ( 'Manca2_pdf' , wait = 2 ) :
        r , _ = pdf.fitTo ( ds , silent = True )
        r , f = pdf.fitTo ( ds , silent = True , draw = True , nbins = 100 )
                
        if 0 != r.status() or 3 != r.covQual() :
            logger.warning ('Fit result\n%s' % r.table ( prefix = '# ' ) ) 
        else :
            logger.info    ('Fit result\n%s' % r.table ( prefix = '# ' ) )
        
        results [ 'test_Manca2'] = r
        plots   [ 'test_Manca2'] = f
        models  [ 'test_Manca2'] = pdf 
        
        
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing
    
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        
        db [ 'mB'      ] = mB
        db [ 'mBc'     ] = mBc
        db [ 'mC'      ] = mC 
        db [ 'mY'      ] = mY
        
        for key in models  : db [ '%s:model'  % key ] = models  [ key ]
        for key in results : db [ '%s:result' % key ] = results [ key ]
        for key in plots   : db [ '%s:plot'   % key ] = plots   [ key ]
        
        db.ls()
        
        
# =============================================================================
if '__main__' == __name__ :

    """
    test_Bd     ()
    test_Bs     ()
    test_Bu     ()
    test_Bc     ()
    
    test_D0     ()
    test_Dp     ()
    test_Ds     ()
    """

    test_BdBs   ()
    test_DpDs   ()

    """
    test_Manca  ()
    test_Manca2 ()
    """
    
    test_db     () 


# =============================================================================
##                                                                      The END 
# =============================================================================
