#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
from   ostap.utils.timing      import timing
from   builtins                import range
from   ostap.fitting.ds2numpy  import ds2numpy
import ostap.fitting.models    as     Models
from   ostap.plotting.canvas   import use_canvas
import ostap.math.models 
import ostap.fitting.roofit 
import ROOT, random 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_ds2numpy' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================


# =============================================================================
def test_small_ds():
    
    logger = getLogger ( 'test_ds2numpy_small_ds' ) 
    
    NN = 1000 
    # Создаем переменные RooRealVar
    x      = ROOT.RooRealVar  ("x", "x"        , 0, 10)
    y      = ROOT.RooRealVar  ("y", "y"        , 0, 10)
    i      = ROOT.RooCategory ("i", 'category' , 'odd' , 'even' ) 
    varset = ROOT.RooArgSet   ( x , y , i ) 
    
    # Создаем RooDataSet
    data = ROOT.RooDataSet("data", "data", varset )

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed
    for e in range ( NN ):
        x_val = random_generator.Uniform(0, 10)
        y_val = random_generator.Uniform(0, 10)
        x.setVal   ( x_val )
        y.setVal   ( y_val )
        f = random.randint ( 0 , 100 )
        i.setVal   ( f % 2  )  
        data.add   ( varset )


    g1 = Models.Gauss_pdf ( 'G1' , xvar = x  , mean = 5  , sigma = 1 )
    g2 = Models.Gauss_pdf ( 'G2' , xvar = y  , mean = 5  , sigma = 1 )
    
    ws = ds2numpy ( data, ['x', 'y' , 'i' ] , more_vars = { 'gaus1' : g1 ,
                                                            'gaus2' : g2 } )
    
    cdfs = data.cdfs ( 'x,y' )
    for key in cdfs :
        with use_canvas ( 'test_small_ds: CDF(%s)' % key , wait = 2 ) : 
            fun = cdfs [ key ]
            fun.draw() 
    
    print ( ws  )

# =============================================================================
def test_small_ds_with_weights():

    logger = getLogger ( 'test_ds2numpy_small_dsw' ) 

    NN = 1000 
    # Создаем переменные RooRealVar
    x      = ROOT.RooRealVar("x", "x", 0, 10)
    y      = ROOT.RooRealVar("y", "y", 0, 10)
    i      = ROOT.RooCategory ("i", 'category' , 'a' , 'b' , 'c' , 'd' , 'e' ) 
    varset = ROOT.RooArgSet(x, y , i ) 
    
    # Создаем RooDataSet
    data = ROOT.RooDataSet("data", "data", varset )

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed
    for e in range(NN):
        x_val = random_generator.Uniform(0, 10)
        y_val = random_generator.Uniform(0, 10)
        x.setVal   ( x_val  )
        y.setVal   ( y_val  )
        f = random.randint ( 0 , 100 ) 
        i.setVal   ( f % 5  ) 
        data.add   ( varset )

    ds = data.makeWeighted('x+y')

    ws = ds2numpy ( ds, ['x', 'y' , 'i' ] )
    print ( ws )
        
# =============================================================================
def test_ds_with_weights():
    
    logger = getLogger ( 'test_ds2numpy_ds_with_weights' )
    
    # Создаем переменные RooRealVar
    
    x = ROOT.RooRealVar("x", "x", -101, 101)
    y = ROOT.RooRealVar("y", "y", -101, 101)
    z = ROOT.RooRealVar("z", "z", -101, 101)
    q = ROOT.RooRealVar("q", "q", -101, 101)
    varset = ROOT.RooArgSet(x, y, z,q) 
    
    # Создаем RooDataSet
    data = ROOT.RooDataSet("data", "data", varset ) 

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed

    NN  = 10000
    for _ in range ( NN ):
        x_val = random.gauss   (    0 ,   1 )
        y_val = random.gauss   (   10 ,   1 )
        z_val = random.uniform ( -100 , 100 ) 
        q_val = random.uniform (  -10 ,  10 ) 

        x.setVal ( x_val )
        y.setVal ( y_val )
        z.setVal ( z_val )
        q.setVal ( q_val )
        data.add( varset )

    ds = data.makeWeighted('x+y')

    ws = ds2numpy ( ds , ['x', 'y' , 'q' ] )

# =============================================================================
def test_large_ds_with_weights():
    
    logger = getLogger ( 'test_ds2numpy_large_ds_with_weights' )

    N  = 100 
    NN = 10000
    # Создаем RooDataSet
    variables = []
    for i in range ( N ):
        var_name = "x{}".format(i)
        var = ROOT.RooRealVar(var_name, var_name, 0, 10)
        variables.append(var)
    varset = ROOT.RooArgSet()
    for v in variables : varset.add ( v )
    data   = ROOT.RooDataSet("data", "data", varset )

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed
    for _ in range( NN ):
        values = [random_generator.Uniform(0, 10) for _ in range(100)]
        for i, var in enumerate(variables):
            var.setVal(values[i])
        data.add( varset )

    ds = data.makeWeighted('x1+x2')

    var_lst = list ( set( "x{}".format( random.randint ( 0 , N - 1  ) ) for i in range ( 50 ) ) ) 

    ws = ds2numpy(ds, var_lst )


# ============================================================================
def test_large_ds_without_weights():
    
    logger = getLogger ( 'test_ds2numpy_large_ds_no_weights' )

    N  = 100
    NN = 10000
    # Создаем RooDataSet
    variables = []
    for i in range ( N ):
        var_name = "x{}".format(i)
        var = ROOT.RooRealVar(var_name, var_name, 0, 10)
        variables.append(var)
    varset = ROOT.RooArgSet()
    for v in variables : varset.add ( v )
    data   = ROOT.RooDataSet("data", "data", varset )

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed
    for _ in range ( NN ):
        values = [ random_generator.Uniform(0, 10) for _ in range ( N )]
        for i, var in enumerate(variables):
            var.setVal(values[i])
        data.add ( varset )


    var_lst = list ( set( "x{}".format( random.randint ( 0 , N - 1 ) ) for i in range ( 50 ) ) ) 
    
    ws = ds2numpy(data, var_lst )
    
    cdfs = data.cdfs ( 'x0,x1,x2,x3' )
    for key in cdfs :
        with use_canvas ( 'test_large_ds: CDF(%s)' % key , wait = 2 ) : 
            fun = cdfs [ key ]
            fun.draw() 
    

# ============================================================================

if '__main__' == __name__ :

    ## pass

    with timing ('Test small ds' , logger ) :
        test_small_ds()

    with timing ('Test small dataset with    weights', logger ) :
        test_small_ds_with_weights()

    with timing ('Test large dataset with    weights', logger ) :
       test_ds_with_weights()
        
    with timing ('Test large dataset with    weights', logger ) :
        test_large_ds_with_weights()
    
    with timing ('Test large dataset without weights', logger ) :
        test_large_ds_without_weights()


# =============================================================================
##                                                                      The END 
# =============================================================================
