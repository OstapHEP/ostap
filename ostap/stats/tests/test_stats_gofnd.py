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
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.timing     import timing
from   ostap.logger.pretty    import pretty_float
from   ostap.plotting.canvas  import use_canvas
from   ostap.math.math_ve     import significance
from   ostap.utils.root_utils import batch_env
from ostap.utils.basic        import numcpu
from   ostap.logger.symbols   import plus_minus , greek_lower_sigma 
import ostap.fitting.models   as     M 
import ostap.stats.gofnd      as     GnD
from   ostap.stats.gof_utils  import clip_pvalue 
import ostap.logger.table     as     T 
import ROOT, random   
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_gofnd' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
batch_env ( logger ) 
# =============================================================================

xvar   = ROOT.RooRealVar ( 'x', '', 0, 10)
yvar   = ROOT.RooRealVar ( 'y', '', 0, 10)
varset = ROOT.RooArgSet  ( xvar , yvar   )

xgauss = M.Gauss_pdf     ( 'GX' , xvar = xvar , mean = ( 5 , 4 , 6 ) , sigma = ( 0.5 , 0.1 , 2.5 ) )
ygauss = M.Gauss_pdf     ( 'GY' , xvar = yvar , mean = ( 5 , 4 , 6 ) , sigma = ( 0.5 , 0.1 , 2.5 ) )
gauss2 = xgauss*ygauss

NG        = 100 
NG2       = 100
data_good = gauss2.generate ( NG + NG2 , sample = False )
data_bad  = gauss2.generate ( NG       , sample = False )
for i in range ( NG2 ) :
    x = random.uniform ( 0 , 10 )
    y = random.uniform ( 0 , 10 )    
    xvar.setVal ( x )
    yvar.setVal ( y )
    data_bad.add ( varset )

xgauss.sigma.fix()
ygauss.sigma.fix()
    
pdf       = gauss2 
rgood , _ = pdf.fitTo  ( data_good , quiet = True , refit = 5 )
## with use_canvas ( 'Good fit x-projection' ) : pdf.draw1 ( data_good , nbins = 50 ) 
## with use_canvas ( 'Good fit y-projection' ) : pdf.draw2 ( data_good , nbins = 50 ) 
rbad  , _ = pdf.fitTo  ( data_bad , quiet = True , refit = 5 )
## with use_canvas ( 'Bad  fit x-projection' ) : pdf.draw1 ( data_bad  , nbins = 50 ) 
## with use_canvas ( 'Bad  fit y-projection' ) : pdf.draw2 ( data_bad  , nbins = 50 ) 

# ===============================================================================
def test_PPD () :
    
    logger = getLogger ("test_PPD")

    rows  = [ ( 'Distance' , 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma ) ]
    
    sigma = ''

    if 10 <= numcpu () : tconf = { 'nToys' : 1000 , 'parallel' : True }
    else               : tconf = { 'nToys' :   50 }
        
    for conf in ( { 'psi' : 'linear'     } ,
                  { 'psi' : 'logarithm'  } ,
                  { 'psi' : 'chebyshev'  } ,   
                  { 'psi' : 'coulomb'    } ,   
                  ## { 'psi' : 'inverse2'   } ,   
                  ## { 'psi' : 'squared'    } ,   
                  ## { 'psi' : 'cosine'     } ,   
                  ## { 'psi' : 'canberra'   } ,   
                  { 'psi' : 'braycurtis' } ,   
                  { 'psi' : 'cityblock'  } ,   
                  { 'psi' : 'gaussian' , 'sigma' : 5.00 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 2.00 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 1.00 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.50 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.10 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.05 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.04 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.03 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.02 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.01 } ) : 

        conf.update ( tconf ) 
        ppd = GnD.PPD ( **conf )
        
        ## presumably good fit
        with timing ( "Good fit PPD distance %s %s" % ( conf [ 'psi' ] , conf.get('sigma','') ) , logger = logger ) :
            pdf.load_params ( rgood , silent = True ) 
            tgood        = ppd        ( pdf , data_good )
            tgood, pgood = ppd.pvalue ( pdf , data_good )

        with use_canvas ( title = 'PPD  good %s' % conf  ) : ppd.draw( tvalue = tgood ) 
            
        ## presumably bad fit 
        with timing ( "Bad  fit PPD distance %s %s" % ( conf [ 'psi' ] , conf.get('sigma','') ) , logger = logger ) : 
            pdf.load_params ( rbad  , silent = True ) 
            tbad        = ppd        ( pdf , data_bad )
            tbad, pbad  = ppd.pvalue ( pdf , data_bad )

        with use_canvas ( title = 'PPD  bad %s ' % conf ) as cnv :
            ppd.draw ( tvalue = tbad  )
            
        gp = pgood * 100 
        bp = pbad  * 100

        gt , ge = pretty_float ( tgood )
        bt , be = pretty_float ( tbad  )

        sigma  = conf.get ( 'sigma' , '' )
        sigma  = '' if not sigma else '/%.2f' % sigma
        
        pvg    = clip_pvalue  ( pgood )
        pvb    = clip_pvalue  ( pbad  )        
        nsg    = significance ( pvg   )
        nsb    = significance ( pvb   )
        
        nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
        nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
            
        row = '%s%s' % ( conf ['psi'] , sigma ) ,\
            '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
            '%5.2f %s %.2f' % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb        
        rows.append ( row )
                
    title = 'Goodness-of-Fit PPD tests'
    rows  = T.remove_empty_columns ( rows ) 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 10*'c' )
    logger.info ( '%s:\n%s' % ( title , table ) )

# ===============================================================================
def test_DNN () :
    
    logger = getLogger ("test_DNN")

    rows  = [ ( 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma ) ]
    
    dnn = GnD.DNN ( nToys = 200 , histo = 50  )
    
    ## presumably good fit
    with timing ( "Good fit DNN" , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood        = dnn        ( pdf , data_good )
        tgood, pgood = dnn.pvalue ( pdf , data_good )
        
    ## presumably bad fit 
    with timing ( "Bad  fit DNN" , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad        = dnn        ( pdf , data_bad )
        tbad, pbad  = dnn.pvalue ( pdf , data_bad )
        
    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )
    
    pvg    = clip_pvalue  ( pgood )
    pvb    = clip_pvalue  ( pbad  )        
    nsg    = significance ( pvg   )
    nsb    = significance ( pvb   )
    
    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    row = '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb 
    rows.append ( row )
            
    title= 'Goodness-of-Fit DNN test'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

# ===============================================================================
def test_USTAT () :
    
    logger = getLogger ("test_USTAT")

    
    from ostap.stats.ustat import USTAT 
    rows  = [ ( 'p-value/good[%]' , 'p-value/bad[%]' , '#sigma/good' , '#sigma/bad') ]
    
    ust = USTAT ( nToys = 200  , histo = 50 )
    
    ## presumably good fit
    with timing ( "Good fit USTAT" , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood        = ust        ( pdf , data_good )
        tgood, pgood = ust.pvalue ( pdf , data_good )
        
    ## presumably bad fit 
    with timing ( "Bad  fit USTAT" , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad        = ust        ( pdf , data_bad )
        tbad, pbad  = ust.pvalue ( pdf , data_bad )
        
    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )

    pvg    = clip_pvalue  ( pgood )
    pvb    = clip_pvalue  ( pbad  )        
    nsg    = significance ( pvg   )
    nsb    = significance ( pvb   )

    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    row = '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb 
    rows.append ( row )
            
    title = 'Goodness-of-Fit USTAT test'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# ===============================================================================
if '__main__' == __name__ :

    test_PPD   ()
    test_DNN   ()
    test_USTAT ()

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
    
