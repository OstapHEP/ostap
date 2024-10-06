#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==============================================================================
# @file ostap/stats/tests/test_stats_gof1d.py 
# Test Goddness-of-fits 1D 
# Copyright (c) Ostap developpers.
# ==============================================================================
""" # Test averages for inconsistend data 
"""
# ==============================================================================
from   ostap.plotting.canvas import use_canvas
from   ostap.utils.timing    import timing
from   ostap.logger.pretty   import pretty_float
from   ostap.plotting.canvas import use_canvas
import ostap.fitting.models  as     M 
import ostap.stats.gofnd     as     GnD 
import ostap.logger.table    as     T 
import ROOT, random   
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_gofnd' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
xvar   = ROOT.RooRealVar ( 'x', '', 0, 10)
yvar   = ROOT.RooRealVar ( 'y', '', 0, 10)
varset = ROOT.RooArgSet  ( xvar , yvar   )

xgauss = M.Gauss_pdf     ( 'GX' , xvar = xvar , mean = ( 5 , 4 , 6 ) , sigma = ( 1.0 , 0.5 , 2.5 ) )
ygauss = M.Gauss_pdf     ( 'GY' , xvar = yvar , mean = ( 5 , 4 , 6 ) , sigma = ( 1.0 , 0.5 , 2.5 ) )
gauss2 = xgauss*ygauss

NG        = 100
NG2       =  10
data_good = gauss2.generate ( NG + NG2 , sample = False )
data_bad  = gauss2.generate ( NG       , sample = False )
for i in range ( NG2 ) :
    x = random.uniform ( 0 , 10 )
    y = random.uniform ( 0 , 10 )    
    xvar.setVal ( x )
    yvar.setVal ( y )
    data_bad.add ( varset )

# ===============================================================================
def test_ppd () :
    
    logger = getLogger ("test_ppd")
    from ostap.stats.gof_np import np,sp,s2u,cdist
    if not np or not sp or not s2u or not cdist:
        logger.warning ('No numpy/scipy/s4u/cdist: skip the test!')
        return         

    pdf    = gauss2 

    rgood , _ = pdf.fitTo  ( data_good , quiet = True , refit = 5 )
    with use_canvas ( 'Good fit x-projection' ) : pdf.draw1 ( data_good , nbins = 50 ) 
    with use_canvas ( 'Good fit y-projection' ) : pdf.draw2 ( data_good , nbins = 50 ) 
    
    rbad  , _ = pdf.fitTo  ( data_bad , quiet = True , refit = 5 )
    with use_canvas ( 'Bad  fit x-projection' ) : pdf.draw1 ( data_bad  , nbins = 50 ) 
    with use_canvas ( 'Bad  fit y-projection' ) : pdf.draw2 ( data_bad  , nbins = 50 ) 

    ## 't/good', 'x[..]' ,
    ## 't/bad' , 'x[..]' ,
    rows  = [ ( 'Distance' , 'sigma' , 'p-value/good[%]' , 'p-value/bad[%]' ) ]
    
    sigma = '' 
    for conf in ( { 'psi' : 'linear'      } ,
                  { 'psi' : 'squared'     } ,
                  { 'psi' : 'logarithm'   } ,
                  { 'psi' : 'coulomb'     } ,   
                  { 'psi' : 'cityblock'   } ,   
                  { 'psi' : 'seuclidean'  } ,   
                  { 'psi' : 'cosine'      } ,    
                  { 'psi' : 'chebyshev'   } ,  
                  { 'psi' : 'canberra'    } ,  
                  { 'psi' : 'braycurtis'  } , 
                  { 'psi' : 'mahalanobis' } ,
                  { 'psi' : 'gaussian' , 'sigma' : 1.00 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.50 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.20 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.10 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.05 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.02 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.01 } ) : 
        
        ppd = GnD.PPD ( Nperm = 1000 , **conf )
        
        ## presumably good fit
        with timing ( "Good fit %s" % conf [ 'psi' ] , logger = logger ) :
            pdf.load_params ( rgood , silent = True ) 
            tgood        = ppd        ( pdf , data_good )
            tgood, pgood = ppd.pvalue ( pdf , data_good )
        
        ## presumably bad fit 
        with timing ( "Bad  fit %s" % conf [ 'psi' ] , logger = logger ) : 
            pdf.load_params ( rbad , silent = True ) 
            tbad        = ppd        ( pdf , data_bad )
            tbad, pbad  = ppd.pvalue ( pdf , data_bad )
            
        gp = pgood * 100 
        bp = pbad  * 100

        gt , ge = pretty_float ( tgood )
        bt , be = pretty_float ( tbad  )

        sigma = conf.get ( 'sigma' , '' )
        sigma = '' if not sigma else '%.2f' % sigma 
        row = conf ['psi'] , sigma ,\
            '%4.1f +/- %.1f' % ( gp.value() , gp.error () ) , \
            '%4.1f +/- %.1f' % ( bp.value() , bp.error () ) 
        rows.append ( row )
        
        
    title= 'Goodness-of-fit test'
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# ===============================================================================
if '__main__' == __name__ :

    test_ppd ()

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
    
