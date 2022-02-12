#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_fill.py
# Test module for filling dataset
# - make some fitting toys 
# ============================================================================= 
""" Test module for filling datasets 
- fill some datasets 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random
from   builtins                     import range
import ostap.core.pyrouts 
import ostap.fitting.roofit
import ostap.io.root_file 
from   ostap.utils.progress_bar     import progress_bar
from   ostap.utils.timing           import timing  
from   ostap.trees.data             import Data
from   ostap.fitting.pyselectors    import Variable, SelectorWithVars
from   ostap.logger.colorized       import attention
import ostap.logger.table           as     T
import ostap.parallel.parallel_fill
from   ostap.parallel.parallel      import DILL_PY3_issue 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_fill' )
else : 
    logger = getLogger ( __name__            )
# =============================================================================

# =============================================================================
## create a file with tree 
def create_tree ( fname , nentries = 1000 ) :
    """Create a file with a tree
    >>> create_tree ( 'file.root' ,  1000 ) 
    """
    
    import ROOT, random 
    import ostap.io.root_file
    
    from array import array 
    var0   = array ( 'i', [ 0 ] )
    var1   = array ( 'd', [ 0 ] )
    var2   = array ( 'd', [ 0 ] )
    var3   = array ( 'd', [ 0 ] )

    NV     = 200
    evars  = [ array ( 'd' , [0]   ) for i in range ( NV )  ]
    
    NM     = 100
    mvars  = [ array ( 'd' , 4*[0] ) for i in range ( NM )  ]
    
    from ostap.core.core import ROOTCWD

    with ROOTCWD() , ROOT.TFile.Open ( fname , 'new' ) as root_file:
        
        root_file.cd () 
        tree = ROOT.TTree ( 'S','tree' )
        tree.SetDirectory ( root_file  ) 
        tree.Branch ( 'evt'   , var0 , 'evt/I'   )
        tree.Branch ( 'mass'  , var1 , 'mass/D'  )
        tree.Branch ( 'pt'    , var2 , 'pt/D'    )
        tree.Branch ( 'eta'   , var3 , 'eta/D'   )

        for i in range ( 1 , NV ) :
            tree.Branch ( "sv%d" % i , evars[i] , 'sv%d/D'    % i )
            
        for i in range ( 1 , NM ) :
            tree.Branch ( "vv%d" % i , mvars[i] , 'vv%d[4]/D' % i )
            
                          
        for i in range ( nentries ) : 
            
            m   = random.gauss        ( 3.1 ,  0.015 )
            pt  = random.uniform      ( 0   ,  8     )
            eta = random.uniform      ( 2   ,  4     )
            
            var0 [ 0 ] = i 
            var1 [ 0 ] = m
            var2 [ 0 ] = pt
            var3 [ 0 ] = eta

            for j in range ( NV )  : evars[j] = j 
            for j in range ( NM )  :
                for k in range(4) :
                    mvars[j][k] = 4 * j + k 

            tree.Fill()
            
        root_file.Write()
            
# =============================================================================
## prepare data for tests 
def prepare_data ( nfiles = 50 ,  nentries = 500  ) :
    """ prepare data for tests
    """
    from ostap.utils.cleanup import CleanUp    
    files = [ CleanUp.tempfile ( prefix = 'ostap-test-fitting-fill-%d-' % i ,
                                 suffix = '.root' ) for i in range ( nfiles)  ]

    for f in progress_bar ( files ) : create_tree ( f , nentries )
    return files

GeV = 1.0
MeV = GeV/1000


def ptcut ( s ) : return 3 < s.pt
def xvar  ( s ) : return (s.mass+s.pt+s.eta)/s.eta 

# =============================================================================
def test_fitting_fill_1 () :
    
    logger = getLogger ('test_fitting_fill_1' ) 

    ## prepare data
    with timing ( "Prepare test data" , logger = logger ) : 
        files = prepare_data ( 4 , 5000 )
        data  = Data ( 'S' ,  files )

    mJPsi = ROOT.RooRealVar ( 'mJPsi' , 'mass(J/Psi) [GeV]' , 3.0 * GeV , 3.2 * GeV )

    # =========================================================================
    logger.info ( attention( 'All trivial variables' ) ) 
    # =========================================================================
 
    variables = [
        Variable ( mJPsi     , accessor = 'mass' ) ,
        Variable ( 'massMeV' , 'mass in MeV' , 3000 , 3200 , 'mass*1000.0' ) , 
        Variable ( 'vv102'   , 'vv10[2]'     , -1   ,  100 , '1.0*vv10[2]' ) ,
        ( 'evt' , ) , 
        ( 'pt'  , ) ,
        ( 'eta' , ) ,
        ( 'x'   , 'some variable'  , 0 , 5000 , '(mass+pt+eta)/eta' ) 
        ]
    
    config = { 'variables' : variables  , 'selection' : "pt>7 && eta<3"  }
    
    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 :
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = False )
        ds1_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 :
        logger.info ( attention ( t2.name ) )        
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = False )        
        ds1_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )        
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = True  )
        ds1_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )        
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = True  )
        ds1_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title1 = "All trivial variables"
    table1 = T.table ( table , title = title1 , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title1 , table1 ) ) 
    
    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 :
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = False , max_files = 1 )
        ds1p_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 :
        logger.info ( attention ( t2.name ) )        
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = False , max_files = 1 )        
        ds1p_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )        
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = True  , max_files = 1 )
        ds1p_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )        
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = True  , max_files = 1 )
        ds1p_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title1p = "All trivial variables (parallel)"
    table1p = T.table ( table , title = title1p , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title1p , table1p ) ) 


    # =========================================================================
    logger.info ( attention( 'Trivial variables + CUT' ) ) 
    # =========================================================================
 
    variables = [
        Variable ( mJPsi     , accessor = 'mass' ) ,
        Variable ( 'massMeV' , 'mass in MeV' , 3000 , 3200 , 'mass*1000'   ) , 
        Variable ( 'vv102'   , 'vv10[2]'     , -1   ,  100 , '1.0*vv10[2]' ) , 
        ( 'evt' , ) , 
        ( 'pt'  , ) ,
        ( 'eta' , ) ,
        ( 'x'   , 'some variable'  , 0 , 5000 , '(mass+pt+eta)/eta' ) 
        ] 

    if not DILL_PY3_issue :
        
        config = { 'variables' : variables          ,
                   'selection' : "pt>7 && eta<3"    ,
                   'cuts'   :lambda s : s.pt > 3    } ## ATTENTION: no trivial cuts!
        
    else :
        
        logger.warning ( 'There is an issue with dill+python3: avoid lambda!' ) 
        config = { 'variables' : variables          ,
                   'selection' : "pt>7 && eta<3"    ,
                   'cuts'      : ptcut              } ## ATTENTION: no trivial cuts!
    
    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 : 
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = False )
        ds2_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 : 
        logger.info ( attention ( t2.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = False )
        ds2_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = True  )
        ds2_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = True  )
        ds2_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title2 = "Trivial variables + CUT"
    table2 = T.table ( table , title = title2 , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title2 , table2 ) ) 

    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 : 
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = False , maX_files = 1 )
        ds2p_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 : 
        logger.info ( attention ( t2.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = False , maX_files = 1 )
        ds2p_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = True  , max_files = 1 )
        ds2p_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = True  , max_files = 1 )
        ds2p_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title2p = "Trivial variables + CUT (parallel)"
    table2p = T.table ( table , title = title2p , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title2p , table2p ) ) 


    # =========================================================================
    logger.info ( attention ( 'Non-trivial variables' ) ) 
    # =========================================================================
     
    if not DILL_PY3_issue : 
       
        variables = [
            Variable ( mJPsi     , accessor = 'mass' ) ,
            Variable ( 'massMeV' , 'mass in MeV' , 3000 , 3200  , 'mass*1000'   ) , 
            Variable ( 'vv102'   , 'vv10[2]'     , -1   ,  100  , '1.0*vv10[2]' ) , 
            ( 'evt' , ) , 
            ( 'pt'  , ) ,
            ( 'eta' , ) ,
            ( 'x'   , 'some variable'  , 0 , 5000 , lambda s : (s.mass+s.pt+s.eta)/s.eta ) 
            ]
        
    else :
        
        logger.warning ( 'There is an issue with dill+python3: avoid lambda!' ) 
        variables = [
            Variable ( mJPsi     , accessor = 'mass' ) ,
            Variable ( 'massMeV' , 'mass in MeV' , 3000 , 3200  , 'mass*1000'   ) , 
            Variable ( 'vv102'   , 'vv10[2]'     , -1   ,  100  , '1.0*vv10[2]' ) , 
            ( 'evt' , ) , 
            ( 'pt'  , ) ,
            ( 'eta' , ) ,
            ( 'x'   , 'some variable'  , 0 , 5000 , xvar ) 
            ]


    config = { 'variables' : variables          ,
               'selection' : "pt>7 && eta<3"    }

    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 : 
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = False )
        ds3_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 : 
        logger.info ( attention ( t2.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = False )
        ds3_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = True  )
        ds3_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = True  )
        ds3_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title3 = "Non-trivial variables"
    table3 = T.table ( table , title = title3 , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title3 , table3 ) ) 


    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 : 
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = False , max_files = 1 )
        ds3p_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 : 
        logger.info ( attention ( t2.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = False , max_files = 1 )
        ds3p_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = True  , max_files = 1 )
        ds3p_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = True  , max_files = 1)
        ds3p_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title3p = "Non-trivial variables (parallel)"
    table3p = T.table ( table , title = title3p , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title3p , table3p ) ) 


    # =========================================================================
    logger.info ( attention ( 'Non-trivial variables + CUT' ) ) 
    # =========================================================================
 
    if not DILL_PY3_issue : 
        
        variables = [
            Variable ( mJPsi     , accessor = 'mass' ) ,
            Variable ( 'massMeV' , 'mass in MeV' , 3000 , 3200  , 'mass*1000'   ) , 
            Variable ( 'vv102'   , 'vv10[2]'     , -1   ,  100  , '1.0*vv10[2]' ) , 
            ( 'evt' , ) , 
            ( 'pt'  , ) ,
            ( 'eta' , ) ,
            ( 'x'   , 'some variable'  , 0 , 5000 , lambda s : (s.mass+s.pt+s.eta)/s.eta ) 
            ]
        
    else :

        logger.warning ( 'There is an issue with dill+python3: avoid lambda!' ) 
        variables = [
            Variable ( mJPsi     , accessor = 'mass' ) ,
            Variable ( 'massMeV' , 'mass in MeV' , 3000 , 3200  , 'mass*1000'   ) , 
            Variable ( 'vv102'   , 'vv10[2]'     , -1   ,  100  , '1.0*vv10[2]' ) , 
            ( 'evt' , ) , 
            ( 'pt'  , ) ,
            ( 'eta' , ) ,
            ( 'x'   , 'some variable'  , 0 , 5000 , xvar ) 
            ]
            
    
    if not DILL_PY3_issue :
        
        config = { 'variables' : variables          ,
                   'selection' : "pt>7 && eta<3"    ,
                   'cuts'      :lambda s : s.pt > 3 } ## ATTENTION: no trivial cuts! 
    else :
        
        logger.warning ( 'There is an issue with dill+python3: avoid lambda!' ) 
        config = { 'variables' : variables          ,
                   'selection' : "pt>7 && eta<3"    ,
                   'cuts'      : ptcut              } ## ATTENTION: no trivial cuts! 
        
    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 : 
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = False )
        ds4_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 : 
        logger.info ( attention ( t2.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = False )
        ds4_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = False , use_frame = True  )
        ds4_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.process ( selector , shortcut = True  , use_frame = True  )
        ds4_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title4 = "Non-trivial variables + CUT"
    table4 = T.table ( table , title = title4 , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title4 , table4 ) ) 


    with timing ( "No SHORTCUT, no FRAME" , logger = None ) as t1 : 
        logger.info ( attention ( t1.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = False , max_files = 1 )
        ds4p_1 = selector.data 
    with timing ( "   SHORTCUT, no FRAME" , logger = None ) as t2 : 
        logger.info ( attention ( t2.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = False , max_files = 1 )
        ds4p_2 = selector.data 
    with timing ( "No SHORTCUT,    FRAME" , logger = None ) as t3 : 
        logger.info ( attention ( t3.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = False , use_frame = True  , max_files = 1 )
        ds4p_3 = selector.data 
    with timing ( "   SHORTCUT,    FRAME" , logger = None ) as t4 : 
        logger.info ( attention ( t4.name ) )
        selector = SelectorWithVars ( **config ) 
        data.chain.pprocess ( selector , shortcut = True  , use_frame = True  , max_files = 1 )
        ds4p_4 = selector.data 

    table = [ ('Configuration' , 'CPU' ) ] 

    table.append ( ( t1.name , '%.3fs' % t1.delta ) )
    table.append ( ( t2.name , '%.3fs' % t2.delta ) )
    table.append ( ( t3.name , '%.3fs' % t3.delta ) )
    table.append ( ( t4.name , '%.3fs' % t4.delta ) )

    title4p = "Non-trivial variables + CUT (parallel)"
    table4p = T.table ( table , title = title4p , prefix = '# ' , alignment = 'rr' )
    logger.info ( '%s\n%s' % ( title4p , table4p ) ) 



    logger.info ( '%s\n%s' % ( title1  , table1  ) )
    logger.info ( '%s\n%s' % ( title1p , table1p ) )
    
    logger.info ( '%s\n%s' % ( title2  , table2  ) ) 
    logger.info ( '%s\n%s' % ( title2p , table2p ) ) 

    logger.info ( '%s\n%s' % ( title3  , table3  ) ) 
    logger.info ( '%s\n%s' % ( title3p , table3p ) )
    
    logger.info ( '%s\n%s' % ( title4  , table4  ) ) 
    logger.info ( '%s\n%s' % ( title4p , table4p ) ) 


    
# =============================================================================
if '__main__' == __name__ :

    
    test_fitting_fill_1 ()
    pass

# =============================================================================
##                                                                      The END 
# =============================================================================
