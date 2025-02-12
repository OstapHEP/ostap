#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_parallel_toys.py
# Test module for ostap/fitting/toys.py
# - make some fitting toys 
# ============================================================================= 
""" Test module for ostap/fitting/toys.py
- make some fitting toys 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.meta_info         import root_info
from   ostap.utils.timing           import timing
from   ostap.utils.memory           import memory 
from   ostap.core.core              import cpp, VE, dsID, hID, rooSilent
from   ostap.plotting.canvas        import use_canvas
from   ostap.utils.utils            import wait, batch_env 
from   ostap.utils.basic            import numcpu 
from   ostap.fitting.toys           import pull_var
from   ostap.math.base              import num_range
from   ostap.utils.progress_bar     import progress_bar
import ostap.parallel.parallel_toys as     Toys
import ostap.logger.table           as     T
import ostap.fitting.models         as     Models
import ostap.fitting.roofit 
import ROOT, sys, os, random, time 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_toys' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
batch_env ( logger ) 
# =============================================================================
if sys.version_info < ( 3, 0 ) and ( 6 , 18 ) <= root_info < ( 6 , 20 ) :
    try : 
        import dill
        if dill.__version__ < '0.3' :
            os.environ['OSTAP_PARALLEL'] = 'GAUDIMP'
            logger.warning ( "Redefined os.environ['OSTAP_PARALLEL']='%s'" % os.environ.get('OSTAP_PARALLEL','' ) ) 
    except ImportError :
        pass


nominal_mean    = 0.4
nominal_sigma   = 0.1
nominal_S       = 1000
nominal_B       = 100
nominal_F       = float ( nominal_S ) / ( nominal_S + nominal_B )

mass            = ROOT.RooRealVar ( 'mass' , '', 0 , 1 )  
gen_gauss       = Models.Gauss_pdf ( 'GG' , xvar = mass )
fit_gauss       = Models.Gauss_pdf ( 'FG' , xvar = mass )
gen_gauss.mean  = nominal_mean
gen_gauss.sigma = nominal_sigma 

model           = Models.Fit1D ( signal = gen_gauss , background = 'flat'                    )
model_NE        = Models.Fit1D ( signal = gen_gauss , background = 'flat' , extended = False )

model.S         = nominal_S 
model.B         = nominal_B
model_NE.F      = nominal_F

toy_results  = {}

fit_config      = { 'silent'   : True  , 'refit'  : 5 }

NCPU = numcpu()
NT   = max ( 4 , min ( NCPU , 12 ) )
logger.info ( 'Parallel toys will use %d/%s split ' % ( NT , NCPU ) ) 
# ==============================================================================
## Perform toy-study for possible fit bias and correct uncertainty evaluation
#  - generate <code>nToys</code> pseudoexperiments with some PDF <code>pdf</code>
#  - fit teach experiment with the same PDF
#  - store  fit results
#  - calculate staistics of pulls
#  - fill distributions for fit results
#  - fill distribution of pulls 
def test_parallel_toys ( ) :
    """ Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `pdf`
    - fit teach experiment with the same PDF
    - store  fit results
    - calculate statistics of pulls
    - fill distributions of fit results
    - fill distributions of pulls 
    """

    logger = getLogger ( 'test_parallel#%d_toys' % NT ) 
    
    more_vars   = { 'pull:mean_GG'  : pull_var ( 'mean_GG'  , nominal_mean  ) ,  
                    'pull:sigma_GG' : pull_var ( 'sigma_GG' , nominal_sigma ) }
    
    params      = { 'S'        : nominal_S      , 
                    'B'        : nominal_B      ,
                    'mean_GG'  : nominal_mean   ,
                    'sigma_GG' : nominal_sigma  }

    ## reset parameters 
    model.load_params ( params )    
    with timing ( 'Toys      analysis' , logger = logger )  :        
        results , stats = Toys.parallel_toys  (
            pdf         = gen_gauss   ,
            nToys       = 500 * NT    ,
            nSplit      = NT          ,
            data        = [ mass ]    , 
            gen_config  = { 'nEvents' : ( nominal_S + nominal_B )  , 'sample' : True } ,
            fit_config  = fit_config  ,
            init_pars   = params      ,
            more_vars   = more_vars   , 
            silent      = True        , 
            progress    = True        ,
            logger      = logger      )
        
    logger.info( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) ) 

    ## make histos:    
    h_mean       = ROOT.TH1F ( hID ()  , 'mean  of Gauss ' , 200 ,  0.75 * nominal_mean  , nominal_mean  * 1.25 ) 
    h_sigma      = ROOT.TH1F ( hID ()  , 'sigma of Gauss'  , 200 ,  0.50 * nominal_sigma , nominal_sigma * 1.50  )
    
    for r in results [ 'mean_GG'  ] : h_mean  .Fill ( r ) 
    for r in results [ 'sigma_GG' ] : h_sigma .Fill ( r )

    for h in ( h_mean , h_sigma ) :
        with use_canvas ( 'test_paralllel#%d_toys  %s' % ( NT , h.title ) , wait = 1 ) : h.draw()

    toy_results['test_parallel_toys'] = results , stats 

# =============================================================================
## Perform toy-study for possible fit bias and correct uncertainty evaluation
#  - generate <code>nToys</code> pseudoexperiments with some PDF <code>gen_pdf</code>
#  - fit teach experiment with the PDF <code>fit_pdf</code>
#  - store  fit results
#  - fill distributions for fit results
def test_parallel_toys2 ( ) :    
    """ Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `gen_pdf`
    - fit teach experiment with the PDF `fit_pdf`
    - store  fit results
    - fill distributions of fit results
    """    

    logger = getLogger ( 'test_parallel#%d_toys2' % NT ) 

    more_vars   = { 'pull:mean_FG'  : pull_var ( 'mean_FG' , nominal_mean  ) ,  
                    'pull:sigma_FG' : pull_var ( 'sigma_FG', nominal_sigma ) }
    
    gen_pars    = { 'mean_GG' : nominal_mean , 'sigma_GG' : nominal_sigma  } 
    fit_pars    = { 'mean_FG' : nominal_mean , 'sigma_FG' : nominal_sigma  }
    
    ## reset parameters 
    gen_gauss.load_params ( gen_pars )    
    fit_gauss.load_params ( fit_pars )    

    with timing ( 'Toys2     analysis' , logger = logger )  :        
        results , stats = Toys.parallel_toys2 (
            gen_pdf     = gen_gauss   ,
            fit_pdf     = fit_gauss   ,
            nToys       = NT * 500    ,
            nSplit      = NT          ,
            data        = [ mass ]    , 
            gen_config  = { 'nEvents' : 100 , 'sample' : True } ,
            fit_config  = fit_config  ,
            gen_pars    = gen_pars    ,
            fit_pars    = fit_pars    ,
            more_vars   = more_vars   , 
            silent      = True        , 
            progress    = True        ,
            logger      = logger      )
        
    logger.info( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) ) 
                    
    ## make histos
        
    h_mean       = ROOT.TH1F ( hID () , 'mean  of Gauss ' , 200 ,  0.75 * nominal_mean  , nominal_mean  * 1.25 ) 
    h_sigma      = ROOT.TH1F ( hID () , 'sigma of Gauss'  , 200 ,  0.50 * nominal_sigma , nominal_sigma * 1.50  )

    for r in results ['mean_FG'  ] : h_mean .Fill ( r ) 
    for r in results ['sigma_FG' ] : h_sigma.Fill ( r )

    for h in ( h_mean , h_sigma ) :
        with use_canvas ( 'test_parallel#%d_toys2 %s' % ( NT , h.title ) , wait = 1 ) : h.draw()

    toy_results['test_parallel_toys2'] = results , stats 

# ============================================================================
def toys3_action ( fit_result , fit_pdf , dataset ) :
    return float ( fit_result.mean_FG ) 
    
# =============================================================================
## Perform toy-study for possible fit bias and correct uncertainty evaluation
#  - generate <code>nToys</code> pseudoexperiments with some PDF <code>gen_pdf</code>
#  - fit teach experiment with the PDF <code>fit_pdf</code>
#  - store  fit results
#  - fill distributions for fit results
def test_parallel_toys3 ( ) :    
    """ Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `gen_pdf`
    - fit teach experiment with the PDF `fit_pdf`
    - store  fit results
    - fill distributions of fit results
    """    

    logger = getLogger ( 'test_parallel#%d_toys3' % NT ) 

    gen_pars    = { 'mean_GG' : nominal_mean , 'sigma_GG' : nominal_sigma  } 
    fit_pars    = { 'mean_FG' : nominal_mean , 'sigma_FG' : nominal_sigma  }
    
    ## reset parameters 
    gen_gauss.load_params ( gen_pars )    
    fit_gauss.load_params ( fit_pars )    

    with timing ( 'Toys3     analysis' , logger = logger )  :        
        results , stats = Toys.parallel_toys3 (
            gen_pdf     = gen_gauss    ,
            fit_pdf     = fit_gauss    ,
            nToys       = NT * 500     ,
            nSplit      = NT           ,
            data        = [ mass ]     ,
            action      = toys3_action , 
            gen_config  = { 'nEvents' : 100 , 'sample' : True } ,
            fit_config  = fit_config   ,
            gen_pars    = gen_pars     ,
            fit_pars    = fit_pars     ,
            silent      = True         , 
            progress    = True         ,
            logger      = logger       )
        
    logger.info( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) ) 
                    
    ## make histos
        
    h_mean       = ROOT.TH1F ( hID () , 'mean  of Gauss ' , 200 ,  0.75 * nominal_mean  , nominal_mean  * 1.25 )

    for r in results [''  ] : h_mean .Fill ( r ) 

    for h in ( h_mean , ) :
        with use_canvas ( 'test_parallel#%d_toys3 %s' % ( NT , h.title ) , wait = 1 ) : h.draw()

    toy_results['test_parallel_toys3'] = results , stats 



# ==============================================================================
## Perform toy-study for Jackknife (non-extended model) 
def test_parallel_jackknife_NE1 ( ) :
    """ Perform toys-study for Jackknife (non-extended model) 
    """

    logger = getLogger ( 'test_parallel#%d_jackknife_NE1' % NT ) 
    logger.info ( 'Jackknife analysis for non-extended model' ) 
                  
    pdf = model_NE
    
    fit_pars  = { 'F'        : nominal_F     ,
                  'sigma_GG' : nominal_sigma ,
                  'mean_GG'  : nominal_mean  }

    ## reset all parameters 
    pdf.load_params ( fit_pars , silent = True )     
    with timing ( 'Jackknife analysis (non-extended)' , logger = logger )  :
        
        data = pdf.generate ( nominal_S + nominal_B , sample = False )
        with use_canvas ( "Jackknife parallel#%d analysis (non-extended)" % NT , wait = 3 ) :
            r , f = pdf.fitTo ( data , draw = True , nbins = 100 , **fit_config )
            
            ## reset all parameters 
            pdf.load_params ( fit_pars , silent = True )
            
            results, stats = Toys.parallel_jackknife ( pdf        = pdf        ,
                                                       nSplit     = NT         , 
                                                       data       = data       ,
                                                       fit_config = fit_config ,
                                                       fit_pars   = fit_pars   , 
                                                       progress   = True       ,
                                                       silent     = True       ,
                                                       logger      = logger    )

            
    logger.info ( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info ( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) )
    
    toy_results [ 'test_parallel_jackknife_NE1' ] = results , stats
    
# ==============================================================================
## Perform toy-study for Jackknife (non-extended model) 
def test_parallel_jackknife_NE2 ( ) :
    """ Perform toys-study for Jackknife (non-extended model) 
    """

    logger = getLogger ( 'test_parallel#%d_jackknife_NE2' % NT ) 
    logger.info ( 'Jackknife analysis for non-extended model' ) 
                  
    pdf = fit_gauss
    
    fit_pars  = { 'sigma_FG' : nominal_sigma ,
                  'mean_FG'  : nominal_mean  }
    
    ## reset all parameters 
    pdf.load_params ( fit_pars , silent = True ) 
    with timing ( 'Jackknife analysis (non-extended)' , logger = logger )  :
        
        data = pdf.generate ( nominal_S , sample = False )
        with use_canvas ( "Jackknife parallel#%d analysis (non-extended)" % NT , wait = 3 ) :
            r , f = pdf.fitTo ( data , draw = True , nbins = 100 , **fit_config )
            
            ## reset all parameters 
            pdf.load_params ( fit_pars , silent = True )             
            results, stats = Toys.parallel_jackknife ( pdf        = pdf        ,
                                                       nSplit     = NT         , 
                                                       data       = data       ,
                                                       fit_config = fit_config ,
                                                       progress   = True       ,
                                                       silent     = True       ,
                                                       logger     = logger     )
            
    logger.info ( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info ( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) )
    
    toy_results [ 'test_parallel_jackknife_NE2' ] = results , stats
            
# ==============================================================================
## Perform toy-study for Jackknife (extended model) 
def test_parallel_jackknife_EXT ( ) :
    """ Perform toys-study for Jackknife (extended model) 
    """

    logger = getLogger ( 'test_parallel#%d_jackknife_EXT' % NT ) 
    logger.info ( 'Jackknife analysis for extended model' ) 
                  
    pdf = model
    
    fit_pars = { 'S'        : nominal_S      , 
                 'B'        : nominal_B      ,
                 'mean_GG'  : nominal_mean   ,
                 'sigma_GG' : nominal_sigma  }
    
    ## reset all parameters 
    pdf.load_params ( fit_pars , silent = True ) 
    with timing ( 'Jackknife analysis (extended)' , logger = logger )  :        

        data = pdf.generate ( nominal_S + nominal_B , sample = False )        
        with use_canvas ( "Jackknife parallel#%d analysis (extended)" % NT , wait = 3 ) :
            r , f = pdf.fitTo ( data , draw = True , nbins = 100 , **fit_config )

            ## reset all parameters 
            pdf.load_params ( fit_pars , silent = True )             
            results, stats = Toys.parallel_jackknife ( pdf        = pdf        ,
                                                       nSplit     = NT         , 
                                                       data       = data       ,
                                                       fit_config = fit_config ,
                                                       fit_pars   = fit_pars   , 
                                                       progress   = True       ,
                                                       silent     = True       ,
                                                       logger     = logger     )
            
    logger.info ( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info ( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) )

    toy_results [ 'test_parallel_jackknife_EXT' ] = results , stats
        
# ==============================================================================
## Perform toy-study for Bootstrap (non-extended)
def test_parallel_bootstrap_NE1 ( ) :
    """ Perform toys-study for Bootstrap (non-extended)
    """

    logger = getLogger ( 'test_parallel#%d_bootstrap_NE1' % NT  ) 
    logger.info ( 'Bootstrap analysis for non-extended model' ) 

    pdf = model_NE
    
    fit_pars  = { 'F'        : nominal_F     ,
                  'sigma_GG' : nominal_sigma ,
                  'mean_GG'  : nominal_mean  }
    
    ## reset all parameters 
    pdf.load_params ( fit_pars , silent = True )       
    with timing ( 'Bootstrap analysis (non-extended)' , logger = logger )  :        

        data    = pdf.generate ( nominal_S + nominal_B , sample = False ) 
        with use_canvas ( "Booststrap parallel#%d analysis" % NT , wait = 3 ) :
            r , f = pdf.fitTo ( data , draw = True , nbins = 100 , **fit_config )
            
            ## reset all parameters 
            pdf.load_params ( fit_pars , silent = True )       
            results, stats = Toys.parallel_bootstrap ( pdf        = pdf        ,
                                                       nSplit     = NT         , 
                                                       size       = 400        , 
                                                       data       = data       ,
                                                       fit_config = fit_config ,
                                                       fit_pars   = fit_pars   , 
                                                       extended   = False      , 
                                                       progress   = True       ,
                                                       silent     = True       ,
                                                       logger     = logger     )
            
    logger.info ( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info ( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) )
    
    toy_results [ 'test_parallel_bootstrap_NE1' ] = results , stats

# ==============================================================================
## Perform toy-study for Bootstrap (non-extended)
def test_parallel_bootstrap_NE2 ( ) :
    """ Perform toys-study for Bootstrap (non-extended)
    """

    logger = getLogger ( 'test_parallel#%d_bootstrap_NE2' % NT ) 
    logger.info ( 'Bootstrap analysis for non-extended model' ) 

    pdf = fit_gauss 

    fit_pars  = { 'sigma_FG' : nominal_sigma ,
                  'mean_FG'  : nominal_mean  }
    
    ## reset all parameters 
    pdf.load_params ( fit_pars , silent = True )     
    with timing ( 'Bootstrap analysis (non-extended)' , logger = logger )  :        

        data    = pdf.generate ( nominal_S + nominal_B ) 
        with use_canvas ( "Booststrap parallel#%d analysis" % NT , wait = 3 ) :
            r , f = pdf.fitTo ( data , draw = True , nbins = 100 , **fit_config )
            
            ## reset all parameters 
            pdf.load_params ( fit_pars , silent = True )     
            results, stats = Toys.parallel_bootstrap ( pdf        = pdf        ,
                                                       nSplit     = NT         , 
                                                       size       = 400        , 
                                                       data       = data       ,
                                                       fit_config = fit_config ,
                                                       fit_pars   = fit_pars   , 
                                                       extended   = False      , 
                                                       progress   = True       ,
                                                       silent     = True       ,
                                                       logger     = logger     )
            
    logger.info ( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info ( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) )
    
    toy_results [ 'test_parallel_bootstrap_NE2' ] = results , stats

# ==============================================================================
## Perform toy-study for Bootstrap (extended)
def test_parallel_bootstrap_EXT ( ) :
    """ Perform toys-study for Bootstrap (extended)
    """
    
    logger = getLogger ( 'test_parallel#%d_bootstrap_EXT' % NT ) 
    logger.info ( 'Bootstrap analysis for extended model' ) 

    pdf = model

    fit_pars = { 'S'        : nominal_S      , 
                 'B'        : nominal_B      ,
                 'mean_GG'  : nominal_mean   ,
                 'sigma_GG' : nominal_sigma  }
    
    ## reset all parameters 
    pdf.load_params ( fit_pars , silent = True ) 
    with timing ( 'Bootstrap analysis (extended)' , logger = logger )  :        

        data    = pdf.generate ( nominal_S + nominal_B ) 
        with use_canvas ( "(Extended) Booststrap parallel#%d analysis" % NT , wait = 3 ) :
            r , f = pdf.fitTo ( data , draw = True , nbins = 100 , **fit_config )

            ## reset all parameters 
            pdf.load_params ( fit_pars , silent = True )             
            results, stats = Toys.parallel_bootstrap ( pdf        = pdf        ,
                                                       nSplit     = NT         , 
                                                       size       = 400        , 
                                                       data       = data       ,
                                                       fit_config = fit_config ,
                                                       fit_pars   = fit_pars   , 
                                                       extended   = True       , 
                                                       progress   = True       ,
                                                       silent     = True       ,
                                                       logger     = logger     )
            
    logger.info ( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info ( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) )
    
    toy_results [ 'test_parallel_bootstrap_EXT' ] = results , stats
        
# =============================================================================
## Perform toy-study for significance of the signal 
#  - generate <code>nToys</code> pseudoexperiments using background-only hypothesis 
#  - fit teach experiment with "signal+background" hypothesis
#  - store  fit results
#  - fill distributions for fit results
def test_parallel_significance ( ) :
    """ Perform toy-study for significance of the signal 
    - generate `nToys` pseudoexperiments using background-only hypothesis 
    - fit each experiment with signal+background hypothesis
    - store  fit results
    - fill distributions for fit results
    """
    
    logger = getLogger ( 'test_parallel#%d_significance' % NT ) 
    logger.info ( 'Make (parallel) Significance analysis' ) 

    with timing ( 'Significance analysis' , logger = logger )  :        

        ## only background hypothesis
        bkg_only = Models.Bkg_pdf    ( "BKG" , xvar =  mass , power = 0 , tau = 0      )
               
        signal   = Models.Gauss_pdf  ( 'S'   , xvar = mass , mean = 0.5 , sigma = 0.1 )
        
        signal.mean .fix ( 0.4 ) 
        signal.sigma.fix ( 0.1 ) 
        
        ## signal + background hypothesis
        model    = Models.Fit1D      ( signal = signal , background = 1 )
        model.background.tau.fix ( 0 )
        
        NB = 100
        NS =  10
        gen_pars  = { 'tau_BKG'  : 0.   }  ## initial values for generation 
        fit_pars  = { 'B'          : NB  ,
                      'S'          : NS  ,
                      'phi0_Bkg_S' : 0.0 }  ## initial fit values for parameters
        
        results , stats = Toys.parallel_toys2 (
            gen_pdf     = bkg_only    ,
            fit_pdf     = model       ,
            ## nToys       = 1000000     ,
            ## nSplit      = 36          , 
            nToys       = NT * 1000   ,
            nSplit      = NT          , 
            data        = [ mass ]    , 
            gen_config  = { 'nEvents' : NB , 'sample' : True } ,
            fit_config  = fit_config  ,
            gen_pars    = gen_pars    , ## initial values for generation 
            fit_pars    = fit_pars    , ## initial fit values for parameters
            silent      = True        , 
            progress    = True        ,
            logger      = logger      )

    logger.info( 'Results : [%s]' % ( ', '.join ( key for key in results ) ) ) 
    logger.info( 'Stats   : [%s]' % ( ', '.join ( key for key in stats   ) ) )

    # =========================================================================
    ## yields 
    h_Y      = ROOT.TH1F ( hID() , '#S' , 200 , 0 , 50 )    
    for r in results [ 'S'  ] : h_Y .Fill ( r )

    ## get p-value and significance histograms 
    h_P , h_S = h_Y.significance ()

    h_S.red  ()
    h_P.blue ()

    minv = h_P.min_positive () 
    if 0 < minv : h_P.SetMinimum ( minv )
    h_S.SetMinimum ( 0    )
    
    with use_canvas ( 'test_parallel#%d_significance: yields'       % NT , wait = 1 ) : h_Y.draw ( )
    with use_canvas ( 'test_parallel#%d_significance: p-value'      % NT , wait = 1 ) : h_P.draw ( logy = True )
    with use_canvas ( 'test_parallel#%d_significance: significance' % NT , wait = 3 ) : h_S.draw ()
    
    toy_results [ 'test_parallel_significance' ] = results , stats

# =============================================================================
## Save resuls of toys to DBASE 
def test_db () :
    """ Save resuls of toys to DBASE 
    """
    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all toys results into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for key in progress_bar ( toy_results ) :
            r , s = toy_results [ key ] 
            key1 = '%s:results' % key
            key2 = '%s:stats'   % key
            db [ key1 ] = r
            db [ key2 ] = s 
        db.ls()
    while toy_results : toy_results.popitem() 
# =============================================================================

    

# =============================================================================
if '__main__' == __name__ :

    with memory ( 'test_parallel_toys'          , logger = logger ) as mtt : test_parallel_toys          ()
    with memory ( 'test_parallel_toys2'         , logger = logger ) as mt2 : test_parallel_toys2         ()
    with memory ( 'test_parallel_toys3'         , logger = logger ) as mt3 : test_parallel_toys3         ()

    with memory ( 'test_parallel_jackknife_NE1' , logger = logger ) as mj1 : test_parallel_jackknife_NE1 ()
    with memory ( 'test_parallel_jackknife_NE2' , logger = logger ) as mj2 : test_parallel_jackknife_NE2 ()
    with memory ( 'test_parallel_jackknife_EXT' , logger = logger ) as mje : test_parallel_jackknife_EXT ()
    with memory ( 'test_parallel_bootstrap_NE1' , logger = logger ) as mb1 : test_parallel_bootstrap_NE1 ()
    with memory ( 'test_parallel_bootstrap_NE2' , logger = logger ) as mb2 : test_parallel_bootstrap_NE2 ()
    with memory ( 'test_parallel_bootstrap_EXT' , logger = logger ) as mbe : test_parallel_bootstrap_EXT ()
    
    with memory ( 'test_parallel_significance'  , logger = logger ) as mts : test_parallel_significance  ()
    with memory ( 'test_db'                     , logger = logger ) as mdb : test_db                     ()

    
    rows = [ ( 'Test' , 'Memory [MB]' ) ]

    row  = 'Toys'          , '%.1f' % mtt.delta 
    rows.append ( row )

    row  = 'Toys2'         , '%.1f' % mt2.delta 
    rows.append ( row )
    
    row  = 'Toys3'         , '%.1f' % mt3.delta 
    rows.append ( row )

    row  = 'Jackknife_NE1' , '%.1f' % mj1.delta 
    rows.append ( row )

    row  = 'Jackknife_NE2' , '%.1f' % mj2.delta 
    rows.append ( row )

    row  = 'Jackknife_EXT' , '%.1f' % mje.delta 
    rows.append ( row )
    
    row  = 'Bootstrap_NE1' , '%.1f' % mb1.delta 
    rows.append ( row )

    row  = 'Bootstrap_NE2' , '%.1f' % mb2.delta 
    rows.append ( row )

    row  = 'Bootstrap_EXT' , '%.1f' % mbe.delta 
    rows.append ( row )
    
    row  = 'Significance'  , '%.1f' % mts.delta 
    rows.append ( row )
    
    row  = 'test_db'       , '%.1f' % mdb.delta 
    rows.append ( row )
    
    title = 'Memory usage for various toys'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lc' ) 
    logger.info ( '%s:\n%s' % ( title , table ) )
    

# =============================================================================
##                                                                      The END 
# =============================================================================
