#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_efficiency2.py
# Test module for ostap/fitting/efficiency.py
# ============================================================================= 
""" Test module for ostap/fitting/efficiency.py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import cpp, VE, dsID, Ostap, rooSilent  
from   ostap.fitting.efficiency import Efficiency1D
from   ostap.utils.utils        import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env  
from   ostap.fitting.roofit     import FIXVAR 
from   ostap.fitting.pyvar      import PyVAR, PyVARLite 
import ostap.fitting.models     as     Models 
import ostap.logger.table       as     T
import ROOT, random, math, time  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_efficiency2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

## make
x           = ROOT.RooRealVar ( 'x',  'test' , 0 , 10 )
xmin , xmax = x.minmax()

acc = ROOT.RooCategory( 'cut','cut')
acc.defineType('accept',1)
acc.defineType('reject',0)
varset  = ROOT.RooArgSet  ( x , acc )
ds      = ROOT.RooDataSet ( dsID() , 'test data' ,  varset )

A  = 0.1
B  = 0.8
C  = 0.5
X0 = 6.0

# =============================================================================
## Simple python function to   parameterize efficiency
def eff  ( x , a , b ,  c , x0 ) :
    """ Simple python function to  parameterize efficiency
    """    
    return a + 0.5 * b * ( 1.0 + math.tanh ( c * 1.0 * ( x - x0 ) ) )

# =============================================================================
## simple python function to parameterize efficiency
def eff0 ( x ) :
    """ Simple python function to parameterize efficiency
    """
    return eff ( x , A , B , C , X0 ) 

emax       = 1.0
N          = 1000

# ==============================================================================
for i in range ( N  ) :    
    xv = random.uniform ( xmin , xmax )    
    x.setVal ( xv )    
    ev = random.uniform ( 0 , emax )    
    if eff0( xv )  > ev : acc.setIndex(1)
    else                : acc.setIndex(0)     
    ds.add ( varset ) 

title = 'Input dataset'
logger.info ( '%s:\n%s' % ( title , ds.table ( title = title , prefix = '# ' ) ) ) 
              
np     = 20
dx     = (xmax-xmin)/np 
points = [ dx * i for i in range ( np + 1 ) ]

a       = ROOT.RooRealVar    ( 'a'  , 'a'  , A  , 0.0     , 2 * A  )
b       = ROOT.RooRealVar    ( 'b'  , 'b'  , B  , 0.5 * B , 2 * B  )
c       = ROOT.RooRealVar    ( 'c'  , 'c'  , C  , 0.1 * C , 5 * C  )
x0      = ROOT.RooRealVar    ( 'x0' , 'x0' , X0 , X0 - 2  , X0 + 2 )

vars    = ROOT.RooArgList ( x , a , b , c , x0 ) 
vars2   = ROOT.RooArgList ( x , a , b , c , x0 ) 

# =================================================================================
## make comparison table 
def make_table ( func , title , prefix = "# ") :
    rows = [ ( 'x' , 'fitted eff [%]' , 'true eff [%]' , 'delta [%]' ) ]
    for p in points :

        e1  = 100 * func ( p , error = True ) 
        e2  = 100 * eff0 ( p )
        d   = e1 - e2    
        row = "%4.2f" % p , \
              "%s"    %  e1.toString ( '(%5.2f+-%4.2f)'  ) ,\
              "%.2f"  %  e2  ,\
              "%s"    %  d .toString ( '(%5.2f+-%4.2f)'  )        
        rows.append ( row )
        
    return T.table ( rows , title = title , prefix = prefix ) 
    
# =============================================================================
## use RooFormulaVar to parameterise the efficiency:
def test_formula () :
    """ Use RooFormulaVar to parameterise efficiency:
    """
    logger = getLogger ( "test_formula" ) 
    
    a.value  = A
    b.value  = B         
    c .value = C
    x0.value = X0
    
    ## create RooFormularVar 
    myEff = ROOT.RooFormulaVar ( "effFunc" , "a+0.5*b*(1+tanh((x-x0)*c))" , vars )
    
    with timing ("Using-RooFormularVar" , logger ) :
        
        eff_obj = Efficiency1D( 'E1' , myEff , cut  = acc , xvar = x )
        
        with FIXVAR ( [ a , b  ] ) : eff_obj.fitTo ( ds , silent = True )            
        with FIXVAR ( [ c , x0 ] ) : eff_obj.fitTo ( ds , silent = True )

        r1   = eff_obj.fitTo ( ds , silent = True )        
        logger.info ( "Fit result using-RooFormularVar \n%s" % r1.table ( prefix = "# ") )        

    with use_canvas ("test_formula" , wait = 2 ) :

        eff_obj.draw ( ds , nbins = 25 )
        title = "Compare with true efficiency using RooFormula"             
        logger.info ( "%s\n%s" % ( title , make_table ( eff_obj , title = 'using RooFormula' ) ) ) 
        
# =============================================================================
## use PyVAR stuff
#  @attention For *OLD* PyROOT only!
def test_pyVAR () :
    """Yse PyVAR stuff
    - For *old* PyROOT only!
    """
    logger = getLogger("test_pyVAR") 
    
    # =========================================================================
    ## @class MyEff2
    #  Local "pythonic" variable
    #  @see PyVAR
    class MyEff(PyVAR) :
        """ Local ``pythonic'' variable
        """
        
        def __init__ ( self , name , title = '' , variables = () , clone = None ) :
            
            assert not clone or isinstance ( clone , MyEff ), \
                "MyEff: invalid `clone` type:%s " % typename ( clone )

            ## initializw the base 
            super ( MyEff , self ).__init__ ( name = name , title = title , variables = variables , clone = clone )

        def clone    ( self , newname = '' ) :
            cloned  = MyEff ( newname if newname else self.name , title = self.title , clone = self )
            ROOT.SetOwnership ( cloned  , False )
            return cloned

        ## evaluate the efficiency 
        def evaluate ( self ) :
            """ Evaluate the efficiency
            - call external function `eff`
            """
            return eff ( *self.values )

    a.value  = A
    b.value  = B         
    c .value = C
    x0.value = X0
    
    ## create the efficiency function (RooAbsReal) 
    myEff   = MyEff ( 'myEff2' , variables = vars , title = 'title' )
    
    with timing ("Using PyVAR" , logger ) :

        ## efficiency fitter 
        eff_obj = Efficiency1D( 'E2' , myEff , cut  = acc , xvar = x )
        
        with FIXVAR ( [ a , b  ] ) : eff_obj.fitTo ( ds , silent = True )            
        with FIXVAR ( [ c , x0 ] ) : eff_obj.fitTo ( ds , silent = True )

        r1   = eff_obj.fitTo ( ds , silent = True )        
        logger.info ( "Fit result using PyVAR\n%s" % r1.table ( prefix = "# ") )        

    with use_canvas ("test_PyVAR" , wait = 2 ) :

        eff_obj.draw ( ds , nbins = 25 )
        title = "Compare with true efficiency using PyVAR"             
        logger.info ( "%s\n%s" % ( title , make_table ( eff_obj , title = 'using PyVAR' ) ) ) 
        
# =============================================================================
## use PyVARLite 
def test_pyVARLite  () :
    """ Use PyVARLite
    """
    
    logger = getLogger("test_pyVARLite") 

    myEff  = PyVARLite ( name = 'myEff' , function  = eff , variables = vars ) 
    
    with timing ("Using-PyVARLite" , logger ) :

        a .value = A 
        b .value = B 
        c .value = C
        x0.value = X0
        
        eff_obj = Efficiency1D( 'E3' , myEff  , cut  = acc , xvar = x )


        with FIXVAR ( [ a , b  ] ) : eff_obj.fitTo ( ds , silent = True )            
        with FIXVAR ( [ c , x0 ] ) : eff_obj.fitTo ( ds , silent = True )

        r1   = eff_obj.fitTo ( ds , silent = True )        
        logger.info ( "Fit result using PyVARLite\n%s" % r1.table ( prefix = "# ") )        

    with use_canvas ("test_PyVAR" , wait = 2 ) :

        eff_obj.draw ( ds , nbins = 25 )
        title = "Compare with true efficiency using PyVARLite"             
        logger.info ( "%s\n%s" % ( title , make_table ( eff_obj , title = 'using PyVARLite' ) ) ) 
        

# =============================================================================
if '__main__' == __name__ :
    
    test_formula   ()
    test_pyVAR     ()
    test_pyVARLite ()
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 
