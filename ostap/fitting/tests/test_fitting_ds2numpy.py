#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
from   ostap.utils.timing      import timing
from   builtins                import range
from   ostap.fitting.ds2numpy  import ds2numpy
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
    x = ROOT.RooRealVar("x", "x", 0, 10)
    y = ROOT.RooRealVar("y", "y", 0, 10)

    # Создаем RooDataSet
    data = ROOT.RooDataSet("data", "data", ROOT.RooArgSet(x, y))

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed
    for _ in range( NN ):
        x_val = random_generator.Uniform(0, 10)
        y_val = random_generator.Uniform(0, 10)
        x.setVal(x_val)
        y.setVal(y_val)
        data.add(ROOT.RooArgSet(x, y))

    with timing ('Test small ds' , logger ) :
        ws = ds2numpy ( data, ['x', 'y'] )


# =============================================================================
def test_small_ds_with_weights():

    logger = getLogger ( 'test_ds2numpy_small_dsw' ) 

    NN = 1000 
    # Создаем переменные RooRealVar
    x = ROOT.RooRealVar("x", "x", 0, 10)
    y = ROOT.RooRealVar("y", "y", 0, 10)
    varset = ROOT.RooArgSet(x, y ) 
    
    # Создаем RooDataSet
    data = ROOT.RooDataSet("data", "data", varset )

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed
    for _ in range(NN):
        x_val = random_generator.Uniform(0, 10)
        y_val = random_generator.Uniform(0, 10)
        x.setVal(x_val)
        y.setVal(y_val)
        data.add( varset )

    ds = data.makeWeighted('x+y')

    with timing ('Test small ds with weights' , logger ) :        
        ws = ds2numpy ( ds, ['x', 'y' ] )
        
# =============================================================================
## def test_ds_with_weights():
if 1 < 2 :
    
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

    NN  = 10
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

    with timing ('Test ds with weights ' , logger ) :        
        ws = ds2numpy ( ds , ['x', 'y' , 'q' ] )

# =============================================================================
def test_large_ds_with_weights():

    
    logger = getLogger ( 'test_ds2numpy_large_ds_with_weights' )

    N  = 100 
    NN = 50000
    # Создаем RooDataSet
    variables = []
    for i in range ( N ):
        var_name = "x{}".format(i)
        var = ROOT.RooRealVar(var_name, var_name, 0, 10)
        variables.append(var)
    varset = ROOT.RooArgSet(*variables)        
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

    with timing ('Test large ds with weights ' , logger ) :        
        ws = ds2numpy(ds, var_lst )


# ============================================================================
def test_large_ds_without_weights():

    logger = getLogger ( 'test_ds2numpy_large_ds_without_weights' )

    N  = 100
    NN = 50000
    # Создаем RooDataSet
    variables = []
    for i in range ( N ):
        var_name = "x{}".format(i)
        var = ROOT.RooRealVar(var_name, var_name, 0, 10)
        variables.append(var)
    varset = ROOT.RooArgSet(*variables)
    data   = ROOT.RooDataSet("data", "data", varset )

    # Заполняем датасет случайными данными
    random_generator = ROOT.TRandom3(42)  # устанавливаем seed
    for _ in range ( NN ):
        values = [ random_generator.Uniform(0, 10) for _ in range ( N )]
        for i, var in enumerate(variables):
            var.setVal(values[i])
        data.add ( varset )


    var_lst = list ( set( "x{}".format( random.randint ( 0 , N - 1 ) ) for i in range ( 50 ) ) ) 
    
    with timing ('Test large ds without weights ' , logger ) :        
        ws = ds2numpy(data, var_lst )

# ============================================================================

if '__main__' == __name__ :

    pass

 ##    with timing ('Test small ds' , logger ) :
##         test_small_ds()

##     with timing ('Test small dataset with    weights', logger ) :
##         test_small_ds_with_weights()

##     with timing ('Test large dataset with    weights', logger ) :
##        test_ds_with_weights()
        
##     with timing ('Test large dataset with    weights', logger ) :
##         test_large_ds_with_weights()
    
##     with timing ('Test large dataset without weights', logger ) :
##         test_large_ds_without_weights()


# =============================================================================
##                                                                      The END 
# =============================================================================
