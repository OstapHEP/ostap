# -*- coding: utf-8 -*-
# ==============================================================================
# @file ostap/stats/tests/test_stats_gof1d.py 
# Test Goddness-of-fits 1D 
# Copyright (c) Ostap developpers.
# ==============================================================================
""" # Test averages for inconsistend data 
"""
# ==============================================================================
import ostap.fitting.models  as     M 
import ostap.stats.gofnd     as     GnD 
from   ostap.plotting.canvas import use_canvas
from   ostap.utils.timing    import timing
from   ostap.logger.pretty   import pretty_float 
import ROOT, random   
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_gofnd' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
xvar   = ROOT.RooRealVar ( 'x', '', 0, 10)
yvar   = ROOT.RooRealVar ( 'y', '', 0, 10)
varset = ROOT.RooArgSet  ( xvar , yvar   )

xgauss = M.Gauss_pdf     ( 'GX' , xvar = xvar , mean = ( 5 , 4 , 6 ) , sigma = ( 0.5 , 0.2 , 2.5 ) )
ygauss = M.Gauss_pdf     ( 'GY' , xvar = yvar , mean = ( 5 , 4 , 6 ) , sigma = ( 0.5 , 0.2 , 2.5 ) )
gauss2 = xgauss*ygauss

NG        = 100
NG2       =   2
data_good = gauss2.generate ( NG + NG2 , sample = False )
data_bad  = gauss2.generate ( NG       , sample = False )
for i in range ( NG2 ) :

    x, y  = -10, -10 
    while not 0 < x < 10 : x = random.gauss ( 9 , 1.0 )
    while not 0 < y < 10 : y = random.gauss ( 9 , 1.0 )
    xvar.setVal ( x )
    yvar.setVal ( y )
    data_bad.add ( varset )

# ===============================================================================
def test_ppd () :
    
    logger = getLogger ("test_ppd")
    
    ppd    = GnD.PPD ( sigma = 0.1 )
    pdf    = gauss2 

    ## presumably good fit
    with timing ( "Good fit" , logger = logger ) : 
        rgood, _     = pdf.fitTo  ( data_good , quiet = True , refit = 5 )
        tgood        = ppd        ( pdf , data_good )
        tgood, pgood = ppd.pvalue ( pdf , data_good )

    ## presumably bad fit 
    with timing ( "Bad  fit" , logger = logger ) : 
        rbad, _     = pdf.fitTo  ( data_bad , quiet = True , refit = 5 )
        tbad        = ppd        ( pdf , data_bad )
        tbad, pbad  = ppd.pvalue ( pdf , data_bad )

        
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )

    gp = pgood * 100 
    bp = pbad  * 100
    
    if ge : logger.info ( 'Good fit T=%s [10^%+d] p-value %+.1f +/- %.1f [%%]' % ( gt , ge , gp.value() , gp.error () ) )
    else  : logger.info ( 'Good fit T=%s          p-value %+.1f +/- %.1f [%%]' % ( gt ,      gp.value() , gp.error () ) )
    if be : logger.info ( 'Bad  fit T=%s [10^%+d] p-value %+.1f +/- %.1f [%%]' % ( bt , be , bp.value() , bp.error () ) )
    else  : logger.info ( 'Bad  fit T=%s          p-value %+.1f +/- %.1f [%%]' % ( bt ,      bp.value() , bp.error () ) )

    
# ===============================================================================
if '__main__' == __name__ :

    test_ppd ()

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
    
