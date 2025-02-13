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
from   ostap.core.meta_info     import old_PyROOT
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env  
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit
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
    """Simple python function to   parameterize efficiency"""
    
    return a + 0.5 * b * ( 1.0 + math.tanh ( c * 1.0 * ( x - x0 ) ) ) 

# =============================================================================
## simple python function to parameterize efficiency
def eff0 ( x ) :
    """Simple python function to parameterize efficiency"""
    return eff ( x , A , B , C , X0 ) 

emax       = 1.0

N          = 1000
for i in range ( N  ) :
    
    xv = random.uniform ( xmin , xmax )
    
    x.setVal ( xv )
    
    ev = random.uniform ( 0 , emax )
    
    if eff0( xv )  > ev : acc.setIndex(1)
    else                : acc.setIndex(0) 
    
    ds.add ( varset ) 

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

        e1  = 100 * func (  p , error = True ) 
        e2  = 100 * eff0 ( p )
        d   = e1 - e2    
        row = "%4.2f" % p , \
              "%s"    %  e1.toString ( '(%5.2f+-%4.2f)'  ) ,\
              "%.2f"  %  e2  ,\
              "%s"    %  d .toString ( '(%5.2f+-%4.2f)'  )
        
        rows.append ( row )
    from ostap.logger.table import table
    return table ( rows , title = title , prefix = prefix ) 
    
# =============================================================================
## use RooFormulaVar to parameterise efficiency:
def test_formula () :
    """Use RooFormulaVar to parameterise efficiency:
    """
    logger = getLogger("test_formula") 
    

    ## create RooFormularVar 
    effFunc = ROOT.RooFormulaVar ( "effFunc" , "a+0.5*b*(1+tanh((x-x0)*c))" , vars )
    
    with timing ("Using-RooFormularVar" , logger ) :
        
        eff1 = Efficiency1D( 'E1' , effFunc , cut  = acc , xvar = x )

        a.fix ( A )
        b.fix ( B )
        c .value = C
        x0.value = X0
        
        r1   = eff1.fitTo ( ds , silent = True )
        
        a.release ()
        b.release ()
        
        c .fix ()
        x0.fix ()
        
        r1   = eff1.fitTo ( ds , silent = True )
        
        c .release()
        x0.release()
        
        r1   = eff1.fitTo ( ds , silent = True )
        
        
        logger.info ( "Fit result using-RooFormularVar \n%s" % r1.table ( prefix = "# ") )
        logger.info ( "Compare with true efficiency (using  RooFormulaVar)\n%s" % make_table (
            eff1 , title = 'using RooFormulaVer') )

    with wait ( 2 ) , use_canvas ("test_formula") : 
        eff1.draw ( ds , nbins = 25 )

    
# =============================================================================
## use PyVAR stuff
#  @attention For *OLD* PyROOT only!
def test_pyVAR () :
    """Yse PyVAR stuff
    - For *old* PyROOT only!
    """
    logger = getLogger("test_pyVAR") 
    
    if not old_PyROOT :
        logger.warning ("test is enabled only for *(very)OLD* PyROOT!")
        return

    from ostap.fitting.pyvar import PyVAR
    
    # =========================================================================
    ## @class MyEff2
    #  Local "pythonic" variable
    #  @see PyVAR
    class MyEff2(PyVAR) :
        """Local ``pythonic'' variable
        """
        def  evaluate ( self ) :

            vlist = self.varlist
            
            _x  = float ( vlist [ 0 ] )
            _a  = float ( vlist [ 1 ] ) 
            _b  = float ( vlist [ 2 ] ) 
            _c  = float ( vlist [ 3 ] ) 
            _x0 = float ( vlist [ 4 ] ) 
            
            return eff ( _x , _a , _b , _c , _x0 )
        
    myEff2  = MyEff2 ( 'myEff2' , vars = vars , title = 'title' )
    the_fun = myEff2.var
    
    with timing ("Using-PyVAR" , logger ) :
        
        
        eff2 = Efficiency1D( 'E2' , the_fun , cut  = acc , xvar = x )
        
        a .fix ( A )
        b .fix ( B )
        c .value = C
        x0.value = X0

        r2   = eff2.fitTo ( ds , silent = True )
        
        a.release ()
        b.release ()
        
        c .fix ()
        x0.fix ()
        
        r2   = eff2.fitTo ( ds , silent = True )
        
        c .release()
        x0.release()
        
        r2   = eff2.fitTo ( ds , silent = True )
        
        logger.info ( "Fit result using-PyVAR \n%s" % r2.table ( prefix = "# ") )
        logger.info ( "Compare with true efficiency (using PyVAR)\n%s" % make_table (
            eff2 , title = 'using PyVAR') )
    
    with wait ( 2 ) , use_canvas ("test_pyVAR") : 
        eff2.draw ( ds , nbins = 25 )
    
# =============================================================================
## use PyVAR2
def test_pyVAR2 () :
    """Use PyVAR2
    """
    
    logger = getLogger("test_pyVAR2") 
    
    if (6,31) <= root_info :
        logger.warning ( 'Test is TEMPORARILY disabled for %s' % ROOT.gROOT.GetVersion() )
        return
    
    from ostap.fitting.pyvar import PyVAR2
    
    myEff3 = PyVAR2 ( name = 'myEff3' , vars = vars , function  = eff )
    
    with timing ("Using-PyVAR2" , logger ) :
        
        eff3 = Efficiency1D( 'E3' , myEff3.var , cut  = acc , xvar = x )
        
        a.fix ( A )
        b.fix ( B )
        c .value = C
        x0.value = X0
        
        r3   = eff3.fitTo ( ds , silent = True )
        
        a.release ()
        b.release ()
        
        c .fix ()
        x0.fix ()
        
        r3   = eff3.fitTo ( ds , silent = True )
        
        c .release()
        x0.release()
        
        r3   = eff3.fitTo ( ds , silent = True )

    
        logger.info ( "Fit result using-PyVAR2 \n%s"     % r3.table ( prefix = "# ") )
        logger.info ( "Compare with true efficiency (using PyVAR2)\n%s" % make_table (
            eff3 , title = 'using PyVAR2') )

    with wait ( 2 ) , use_canvas ("test_pyVAR2") : 
        eff3.draw ( ds , nbins = 25 )

# =============================================================================
## use PyVar stuff
#  @attention For *NEW* PyROOT only!
def test_pyVar () :
    """use PyVar stuff
    - For *NEW* PyROOT only!
    """
    
    logger = getLogger("test_pyVar") 
    
    if old_PyROOT :
        logger.warning ("test is enabled only for *NEW* PyROOT!")
        return

    if (6,31) <= root_info :
        logger.warning ( 'Test is TEMPORARILY disabled for %s' % ROOT.gROOT.GetVersion() )
        return
    
    # =========================================================================
    ## @class MyEff4
    #  Local ``pythonic'' variable
    class MyEff4 (Ostap.Functions.PyVar) :
        
        def __init__ ( self , name , title ,  variables ) :

            vlist = ROOT.RooArgList()
            for v in   variables : vlist.add ( v ) 
            super(MyEff4,self).__init__ ( name , title , vlist )
                        
        def clone ( self , newname ) :
            
            name = newname if newname else self.name             
            nv = MyEff4( name ,  self.title  , self.variables()  )
            ROOT.SetOwnership ( nv ,  False )             
            return nv 
        
        def evaluate ( self ) :
            
            vlist = self.variables() 
            
            _x  = float ( vlist [ 0 ] )
            _a  = float ( vlist [ 1 ] ) 
            _b  = float ( vlist [ 2 ] ) 
            _c  = float ( vlist [ 3 ] ) 
            _x0 = float ( vlist [ 4 ] ) 
            
            return eff ( _x , _a , _b , _c , _x0 )
        
    myEff4  = MyEff4 ( 'myEff4' , variables = vars , title = 'title' )
    the_fun = myEff4
        
    with timing ("Using-PyVar" , logger ) :
        
        eff4 = Efficiency1D( 'E4' , the_fun , cut  = acc , xvar = x )

        a.fix ( A )
        b.fix ( B )
        c .value = C
        x0.value = X0
        
        r4   = eff4.fitTo ( ds , silent = True )
        
        a.release ()
        b.release ()
        
        c .fix ()
        x0.fix ()
        
        r4   = eff4.fitTo ( ds , silent = True )
        
        c .release()
        x0.release()
        
        r4   = eff4.fitTo ( ds , silent = True )


        logger.info ( "Fit result using-PyVar \n%s"      % r4.table ( prefix = "# ") )
        logger.info ( "Compare with true efficiency (using  PyVar)\n%s" % make_table (
            eff4 , title = 'using PyVar') )

    with wait ( 2 ) , use_canvas ("test_pyVar") : 
        eff4.draw ( ds , nbins = 25 )

# =============================================================================
if '__main__' == __name__ :
    
    test_formula ()
    test_pyVAR   ()
    test_pyVAR2  ()
    test_pyVar   ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
