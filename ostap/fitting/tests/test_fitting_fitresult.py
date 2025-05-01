#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_fitresult.py
# Test module for <code>RooFitResuls</code>
# ============================================================================= 
""" Test module for `ROOT.RooFitResult`
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing       import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
from   ostap.logger.colorized   import attention
import ostap.fitting.models     as     Models 
import ostap.logger.table       as     T
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_fitresult' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

## make simple test mass 
mass        = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 30 )

## book very simple data set
varset      = ROOT.RooArgSet  ( mass )
dataset     = ROOT.RooDataSet ( dsID() , 'Test Data set' , varset )  

mmin , mmax = mass.minmax()

NS1 = 10000
NS2 = 10000
NS3 = 10000
NB  = 10000

m1  = VE(10,1.5**2)
m2  = VE(15,2.5**2)
m3  = VE(20,1.5**2)

for i in range(0,NS1) :
    mass.value = m1.gauss () 
    dataset.add ( varset )

for i in range(0,NS2) :
    mass.value = m2.gauss () 
    dataset.add ( varset )

for i in range(0,NS3) :
    mass.value = m3.gauss () 
    dataset.add ( varset )

for i in range(0,NB) :
    mass.value = random.uniform ( mmin , mmax ) 
    dataset.add ( varset   )

logger.info ('DATASET\n%s' % dataset )



def test_fitresult () :


    logger = getLogger('test_fitresult')
    
    signal1 = Models.Gauss_pdf ('G1' , xvar = mass , mean = ( 10 ,  8 , 12 ) , sigma = ( 1.5 , 1.0 , 2.0 ) )
    signal2 = Models.Gauss_pdf ('G2' , xvar = mass , mean = ( 15 , 13 , 17 ) , sigma = ( 2.5 , 2.0 , 3.0 ) )
    signal3 = Models.Gauss_pdf ('G3' , xvar = mass , mean = ( 20 , 18 , 22 ) , sigma = ( 1.5 , 1.0 , 2.0 ) )
    
    model   = Models.Fit1D ( signals    = ( signal1 , signal2 , signal3 ) ,
                             background = -2  )
    model.S = NS1 , NS2 , NS3
    model.B = NB

    result, frame = model. fitTo ( dataset , silent = True )
    result, frame = model. fitTo ( dataset , silent = True )
    result, frame = model. fitTo ( dataset , silent = True , draw = True , nbins = 100  )

    logger.info ('Fit result\n%s' % result.table ( prefix = '# ' ) ) 

    s0 = result.S_0 * 1
    s1 = result.S_1 * 1
    s2 = result.S_2 * 1

    rows = [ ( 'expression' , 'naive' , 'native' , 'generic' ) ]

    ## S_0 + S_1
    fun1 = lambda x,y : x+y 
    row = 'S_0+S_1'                                  , \
          '%s' % ( s0 + s1 )                         , \
          '%s' % result.sum      ( 'S_0' , 'S_1' )   , \
          '%s' % result.evaluate ( fun1  , ( 'S_0' , 'S_1' ) )
    rows.append(row)

    ## S_0 + S_2
    fun2 = lambda x,y : x+y 
    row = 'S_0+S_2'                                 , \
          '%s' % ( s0 + s2 )                        , \
          '%s' % result.sum      ( 'S_0' , 'S_2' )  , \
          '%s' % result.evaluate ( fun2  , ( 'S_0' , 'S_2' ) )
    rows.append(row)

    ## S_1 + S_2
    fun3 = lambda x,y : x+y 
    row = 'S_1+S_2'                                  , \
          '%s' % ( s1 + s2 )                         , \
          '%s' % result.sum      ( 'S_1' , 'S_2' )   , \
          '%s' % result.evaluate ( fun3  , ( 'S_1' , 'S_2' ) )
    rows.append(row)
        
    ## S_0 + S_1 + S_2
    fun4 = lambda x,y,z : x+y+z 
    row = 'S_0+S_1+S_2'                                      , \
          '%s' % ( s0 + s1 + s2 )                            , \
          '%s' % result.sum      ( 'S_0' , 'S_1' , 'S_2' )   , \
          '%s' % result.evaluate ( fun4  , ( 'S_0' , 'S_1' , 'S_2' ) )
    rows.append(row)

    ## S_0 * S_1
    fun5 = lambda x,y : x*y 
    row = 'S_0*S_1'                                  , \
          '%s' % ( s0 * s1 )                         , \
          '%s' % result.multiply ( 'S_0' , 'S_1' )   , \
          '%s' % result.evaluate ( fun5  , ( 'S_0' , 'S_1' ) )
    rows.append(row)

    ## S_0 * S_2
    fun6 = lambda x,y : x*y 
    row = 'S_0*S_2'                                  , \
          '%s' % ( s0 * s2 )                         , \
          '%s' % result.multiply ( 'S_0' , 'S_2' )   , \
          '%s' % result.evaluate ( fun6  , ( 'S_0' , 'S_2' ) )
    rows.append(row)

    ## S_1 * S_2
    fun7 = lambda x,y : x*y 
    row = 'S_1+S_2'                                  , \
          '%s' % ( s1 * s2 )                         , \
          '%s' % result.multiply ( 'S_1' , 'S_2' )   , \
          '%s' % result.evaluate ( fun7  , ( 'S_1' , 'S_2' ) )
    rows.append(row)
        
    ## S_0 * S_1 * S_2
    fun8 = lambda x,y,z : x*y*z 
    row = 'S_0*S_1*S_2'                                       , \
          '%s' % ( s0 * s1 * s2 )                             , \
          '%s' % result.multiply ( 'S_0' , 'S_1' , 'S_2' )    , \
          '%s' % result.evaluate ( fun8  , ( 'S_0' , 'S_1' , 'S_2' ) )
    rows.append(row)


    ## S_1 / ( S_0 + S_2 ) 
    fun9 = lambda x,y,z : y/(x+z)
    row = 'S_1/(S_0+S_2)'             , \
          '%s' % ( s1 / ( s0 + s2 ) ) , \
          ' '                         , \
          '%s' % result.evaluate ( fun9  , ( 'S_0' , 'S_1' , 'S_2' ) )
    rows.append(row)
    
    ## (S_0 + S_1) / ( S_1 + S_2 ) 
    fun10 = lambda x,y,z : (x+y)/(y+z)
    row = '(S_0+S_1)/(S_1+S_2)' , \
          '%s' % ( ( s0 +  s1 ) / ( s1 + s2 ) ) , \
          ' '                                   , \
          '%s' % result.evaluate ( fun10  , ( 'S_0' , 'S_1' , 'S_2' ) )
    rows.append(row)


    for v in ( 0 , 1 , 100 , 500 , 1000 , 5000 , 10000 ) :
        
        ## S_0/(S_0+v) +  S_1/(S_1+v)  + S_2/(S_2+v)
        fun11 = lambda x,y,z : x/(x+v)+y/(y+v)+z/(z+v) 
        row = 'S_0/(S_0+%s)+S_1/(S_1+%s)+S_2/(S_2+%s)' % ( v , v , v ) , \
              '%s' % (  s0/(s0+v) + s1/(s1+v) + s2/(s2+v) ) , \
              ' ' ,  \
              '%s' % result.evaluate ( fun11  , ( 'S_0' , 'S_1' , 'S_2' ) )
        rows.append(row)
        
        
    title = 'Expressions & their uncertainties'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccc' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 
    
    




# =============================================================================
if '__main__' == __name__ :

    with timing ('test_fitresult' , logger ) :
        test_fitresult() 


# =============================================================================
##                                                                      The END 
# ============================================================================= 
