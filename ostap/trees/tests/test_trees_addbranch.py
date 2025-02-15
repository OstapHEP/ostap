#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/trees/tests/test_trees_addbranch.py
# Four ways to add branch into TTree/Tchain
# - using string formula (TTreeFormula-based)
# - using pure python function
# - using histogram/function
# - using histogram sampling
#
# There are more specific similar actions
# - adding to TTree/TChain TMVA decision
# - adding to TTree/TChain TMVA/Choppnig decision
# - adding to TTree/TChain reweigting information
# Copyright (c) Ostap developers.
# ============================================================================= 
""" Test module for ostap/trees/cuts.py.
Four ways to add branch into TTree/Tchain
- using string formula (TTreeFormula-based)
- using pure python function
- using histogram/function
- using histogram sampling
"""
# ============================================================================= 
from   ostap.core.meta_info     import root_info
from   ostap.core.pyrouts       import hID , Ostap 
from   ostap.trees.data         import Data
from   ostap.math.make_fun      import make_fun1, make_fun2, make_fun3 
from   ostap.utils.timing       import timing 
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import batch_env
import ostap.logger.table       as     T 
import ostap.trees.trees
import ostap.histos.histos
import ROOT, cppyy, math, random, array  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_trees_addbranch' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
batch_env ( logger ) 
# =============================================================================
## define three C++ fnuctions 
cppyy.cppdef( """
namespace MyTest 
{
   double cpp_fun1 ( double x                      ) { return x         ; } 
   double cpp_fun2 ( double x, double y            ) { return x + y     ; } 
   double cpp_fun3 ( double x, double y , double z ) { return x + y + z ; } 
}""" )


## create a file with tree 
def create_tree ( fname , nentries = 1000 ) :
    """Create a file with a tree
    >>> create_tree ( 'file.root' ,  1000 ) 
    """
    
    import ROOT, random 
    import ostap.io.root_file
    
    from array import array 
    var1 = array ( 'd', [ 0 ] )
    var2 = array ( 'd', [ 0 ] )
    var3 = array ( 'd', [ 0 ] )
    
    from ostap.core.core import ROOTCWD

    with ROOTCWD() , ROOT.TFile.Open( fname , 'new' ) as root_file:
        root_file.cd () 
        tree = ROOT.TTree ( 'S','tree' )
        tree.SetDirectory ( root_file  ) 
        tree.Branch ( 'mass'  , var1 , 'mass/D'  )
        tree.Branch ( 'c2dtf' , var2 , 'c2dtf/D' )
        tree.Branch ( 'pt'    , var3 , 'pt/D'    )
        
        for i in range ( nentries ) : 
            
            m  = random.gauss        ( 3.1 ,  0.015 )
            c2 = random.gammavariate ( 2.5 , 0.5    ) / 5 
            pt = random.uniform      ( 0   , 10     )
            
            var1[0] = m
            var2[0] = c2 
            var3[0] = pt
            
            tree.Fill()
            
        root_file.Write()
        
# =============================================================================
def prepare_data ( nfiles = 50 ,  nentries = 500  ) :

    from ostap.utils.cleanup import CleanUp    
    files = [ CleanUp.tempfile ( prefix = 'ostap-test-trees-addbranch-%d-' % i ,
                                 suffix = '.root' ) for i in range ( nfiles)  ]
    
    for f in progress_bar ( files ) : create_tree ( f , nentries )
    return files

# =============================================================================
# Many ways to add branch into TTree/Tchain
# - using string formula (TTreeFormula-based)
# - using pure python function
# - using histogram/function
# - using histogram sampling
def test_addbranch() :
    """ Many ways to add branch into TTree/Tchain
    - using string formula (TTreeFormula-based)
    - using pure python function
    - using histogram/function
    - using histogram sampling
    """
    
    files = prepare_data ( 3 , 100 )
    
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )

    rows = [ ( 'Method' , 'CPU [s]' ) ]

    # =========================================================================
    ## 1) add new branch as TTree-formula:
    # =========================================================================
    with timing ('expression-1' , logger = logger ) as timer :          
        chain = data.chain 
        chain = chain.add_new_branch ( 'sqrt(pt*pt+mass*mass)' , name  = 'et1' )
    rows.append ( ( timer.name , '%.3f' % timer.delta ) ) 
    ## reload the chain and check: 
    assert 'et1' in chain , "Branch `et1' is  not here!"

    """ 

    # =========================================================================
    ## 1) add new branch as TTree-formula:
    # =========================================================================
    with timing ('expression-2' , logger = logger ) as timer :          
        chain = data.chain 
        chain = chain.add_new_branch ( 'et2' , formula = 'sqrt(pt*pt+mass*mass)'  )        
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) ) 
    ## reload the chain and check: 
    assert 'et2' in chain , "Branch `et2' is  not here!"

    # =========================================================================
    ## 2) add several new branches as TTree-formula:
    # =========================================================================
    with timing ('simultaneous' , logger = logger ) as timer :          
        chain = data.chain  
        chain = chain.add_new_branch ( { 'et3' : 'sqrt(pt*pt+mass*mass)' ,
                                         'et4' : 'sqrt(pt*pt+mass*mass)' ,
                                         'et5' : 'sqrt(pt*pt+mass*mass)' ,
                                         'et6' : 'sqrt(pt*pt+mass*mass)' ,
                                         'et7' : 'sqrt(pt*pt+mass*mass)' ,
                                         'et8' : 'sqrt(pt*pt+mass*mass)' ,
                                         'et9' : 'sqrt(pt*pt+mass*mass)' } ) 
    rows.append ( ( timer.name , '%.3f' % timer.delta ) ) 
    ## reload the chain and check: 
    assert 'et3' in chain , "Branch `et3' is  not here!"
    assert 'et4' in chain , "Branch `et4' is  not here!"
    assert 'et5' in chain , "Branch `et5' is  not here!"
    assert 'et6' in chain , "Branch `et6' is  not here!"
    assert 'et7' in chain , "Branch `et7' is  not here!"
    assert 'et8' in chain , "Branch `et8' is  not here!"
    assert 'et9' in chain , "Branch `et9' is  not here!"

    # =========================================================================
    ## 3) add new branch as pure python function 
    # =========================================================================
    with timing ('pyfunc-1' , logger = logger ) as timer :          
        pt2   = lambda tree : tree.pt**2 + tree.mass**2        
        chain = data.chain
        chain = chain.add_new_branch ( 'pt2', function = pt2 )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) ) 
    ## reload the chain and check: 
    assert 'pt2' in data.chain , "Branch `pt2' is  not here!"

    # =========================================================================
    ## 4) add new branch as pure python function 
    # =========================================================================
    with timing ('pyfunc-2' , logger = logger ) as timer :          
        pt3   = lambda tree : tree.pt**2 + tree.mass**2        
        chain = data.chain
        chain = chain.add_new_branch ( pt2 , name = 'pt3' )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) ) 
    ## reload the chain and check: 
    assert 'pt3' in data.chain , "Branch `pt3' is  not here!"

    # =========================================================================
    ## 5) add new branch as histogram-function 
    # =========================================================================
    with timing ('fun-hist-1' , logger = logger ) as timer :          
        h1  = ROOT.TH1D ( hID () , 'some pt-correction' , 100 , 0 , 10 )
        h1 += lambda x :  1.0 + math.tanh( 0.2* ( x - 5 ) )         
        from   ostap.trees.funcs  import FuncTH1
        ptw = FuncTH1 ( h1 , 'pt' )
        chain = data.chain 
        chain = chain.add_new_branch ( ptw , name = 'ptw1' )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'ptw1' in chain , "Branch `ptw1' is  not here!"

    # =========================================================================
    ## 6) add new branch as histogram-function 
    # =========================================================================
    with timing ( 'fun-hist-2' , logger = logger ) as timer :          
        h1  = ROOT.TH1D ( hID () , 'some pt-correction' , 100 , 0 , 10 )
        h1 += lambda x :  1.0 + math.tanh( 0.2* ( x - 5 ) )         
        from   ostap.trees.funcs  import FuncTH1
        ptw   = FuncTH1 ( h1 , 'pt' )
        chain = data.chain 
        chain = chain.add_new_branch ( 'ptw2' , function = ptw )     
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'ptw2' in chain , "Branch `ptw2' is  not here!"

    # =========================================================================
    ## 7) add several functions simultaneously 
    # =========================================================================
    with timing ('sim-funcs' , logger = logger ) as timer :          
        from   ostap.trees.funcs  import FuncTH1
        hh    = ROOT.TH1D ( hID() , 'some pt-correction' , 100 , 0 , 10 )
        h1    = hh + ( lambda x :  1.0 + math.tanh ( 0.1 * ( x - 5 ) ) )
        ptw3  = FuncTH1 ( h1 , 'pt' )
        h2    = h1 + ( lambda x :  1.0 + math.tanh ( 0.2 * ( x - 5 ) ) )
        ptw4  = FuncTH1 ( h2 , 'pt' )
        h3    = h1 + ( lambda x :  1.0 + math.tanh ( 0.3 * ( x - 5 ) ) ) 
        ptw5  = FuncTH1 ( h3 , 'pt' )
        chain = data.chain
        brs   = { 'ptw3' : ptw3 , 'ptw4' : ptw4 , 'ptw5' : ptw4 } 
        chain = chain.add_new_branch ( brs )     
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'ptw3' in chain , "Branch `ptw3' is  not here!"
    assert 'ptw4' in chain , "Branch `ptw4' is  not here!"
    assert 'ptw5' in chain , "Branch `ptw5' is  not here!"
    
    h1 = ROOT.TH1D ( hID() , 'Gauss1' , 120 , -6 , 6 )
    h2 = ROOT.TH2D ( hID() , 'Gauss2' ,  50 , -6 , 6 , 50 , -6 , 6 )
    h3 = ROOT.TH3D ( hID() , 'Gauss3' ,  20 , -6 , 6 , 20 , -6 , 6 , 20 , -6 , 6 )
    for i in range ( 100000 ) :
        g1 = random.gauss ( 0, 1 )
        g2 = random.gauss ( 0, 1 )
        g3 = random.gauss ( 0, 1 )
        h1.Fill ( g1 ) ; h1.Fill ( g2 ) ; h1.Fill ( g3 ) ;
        h2.Fill ( g1 , g2 )
        h3.Fill ( g1 , g2 , g3 )
        
    # =========================================================================
    ## 8) add the variable sampled from the 1D histogram
    # =========================================================================
    with timing ('sample-h1' , logger = logger ) as timer :          
        chain = data.chain 
        chain = chain.add_new_branch ( h1 , xname = 'xh1' ) 
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'xh1' in chain , "Branch `xh1' is  not here!"

    # =========================================================================
    ## 9) add the variable sampled from the 1D histogram
    # =========================================================================
    with timing ('sample-h2' , logger = logger ) as timer :          
        chain = data.chain 
        chain = chain.add_new_branch ( h2 , xname = 'xh2' , yname = 'yh2' ) 
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'xh2' in chain , "Branch `xh2' is  not here!"
    assert 'yh2' in chain , "Branch `yh2' is  not here!"
    
    # =========================================================================
    ## 10) add the variable sampled from the 1D histogram
    # =========================================================================
    with timing ('sample-h3' , logger = logger ) as timer :          
        chain = data.chain 
        chain = chain.add_new_branch ( h3 , xname = 'xh3' , yname = 'yh3' , zname = 'zh3' ) 
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'xh3' in chain , "Branch `xh2' is  not here!"
    assert 'yh3' in chain , "Branch `yh2' is  not here!"
    assert 'zh3' in chain , "Branch `zh2' is  not here!"
    
    # =========================================================================
    ## 11) python function again 
    # =========================================================================
    with timing ('gauss-1' , logger = logger ) as timer :          
        def gauss ( *_ ) : return random.gauss(0,1)    
        chain = data.chain 
        chain = chain.add_new_branch ( gauss , name = 'gauss1'  )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gauss1' in chain , "Branch `gauss1' is  not here!"
    
    # =========================================================================
    ## 12) python function again 
    # =========================================================================
    with timing ('gauss-2' , logger = logger ) as timer :          
        def gauss ( *_ ) : return random.gauss(0,1)    
        chain = data.chain 
        chain = chain.add_new_branch ( 'gauss2' , function = gauss )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gauss2' in chain , "Branch `gauss2' is  not here!"

    """
    
    # =========================================================================
    ## 13) generic fucction with 1 argument 
    # =========================================================================
    with timing ('generic-1D-py-1' , logger = logger ) as timer :          
        def fun1  ( x ) : return x  
        chain = data.chain 
        chain = chain.add_new_branch ( fun1 , name = 'gf1py_1' , arguments = 'pt*10' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf1py_1' in chain , "Branch `gf1py_1' is  not here!"

    """
    
    # =========================================================================
    ## 14) generic function with 1 argument 
    # =========================================================================
    with timing ('generic-1D-py-2' , logger = logger ) as timer :          
        def fun1  ( x ) : return x 
        chain = data.chain 
        chain = chain.add_new_branch ( 'gf1py_2' , function = fun1 , arguments = 'pt*10' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf1py_2' in chain , "Branch `gf1py_2' is  not here!"

    # =========================================================================
    ## 15) generic function with 1 argument 
    # =========================================================================
    with timing ('generic-1D-cxx-1' , logger = logger ) as timer :          
        fun1  = ROOT.MyTest.cpp_fun1 
        chain = data.chain 
        chain = chain.add_new_branch ( fun1 , name = 'gf1cxx_1' , arguments = 'pt*10' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf1cxx_1' in chain , "Branch `gf1cxx_1' is  not here!"

    # =========================================================================
    ## 16) generic function with 1 argument 
    # =========================================================================
    with timing ('generic-1D-cxx-2' , logger = logger ) as timer :          
        fun1  = ROOT.MyTest.cpp_fun1 
        chain = data.chain 
        chain = chain.add_new_branch ( 'gf1cxx_2' , function = fun1 , arguments = 'pt*10' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf1cxx_2' in chain , "Branch `gf1cxx_2' is  not here!"

    # =========================================================================
    ## 17) generic function with 2 arguments 
    # =========================================================================
    with timing ('generic-2D-py-1' , logger = logger ) as timer :          
        def fun2  ( x , y ) : return x + y  
        chain = data.chain 
        chain = chain.add_new_branch ( fun2 , name = 'gf2py_1' , arguments = 'et1,et2' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf2py_1' in chain , "Branch `gf2py_1' is  not here!"

    # =========================================================================
    ## 18) generic function with 2 arguments 
    # =========================================================================
    with timing ('generic-2D-py-2' , logger = logger ) as timer :          
        def fun2  ( x , y ) : return x + y 
        chain = data.chain 
        chain = chain.add_new_branch ( 'gf2py_2' , function = fun2 , arguments = ( 'et1',  'et2' ) ) 
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf2py_2' in chain , "Branch `gf2py_2' is  not here!"

    # =========================================================================
    ## 19) generic function with 2 arguments 
    # =========================================================================
    with timing ('generic-2D-cxx-1' , logger = logger ) as timer :          
        fun2  = ROOT.MyTest.cpp_fun2 
        chain = data.chain 
        chain = chain.add_new_branch ( fun2 , name = 'gf2cxx_1' , arguments = 'et1,et2' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf2cxx_1' in chain , "Branch `gf2cxx_1' is  not here!"

    # =========================================================================
    ## 20) generic function with 2 arguments 
    # =========================================================================
    with timing ('generic-2D-cxx-2' , logger = logger ) as timer :          
        fun2  = ROOT.MyTest.cpp_fun2 
        chain = data.chain 
        chain = chain.add_new_branch ( 'gf2cxx_2' , function = fun2 , arguments = ( 'et1',  'et2' ) ) 
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf2cxx_2' in chain , "Branch `gf2cxx_2' is  not here!"
    
    # =========================================================================
    ## 21) generic function with 3 arguments 
    # =========================================================================
    with timing ('generic-3D-py-1' , logger = logger ) as timer :          
        def fun3  ( x , y , z ) : return x + y + z  
        chain = data.chain 
        chain = chain.add_new_branch ( fun3 , name = 'gf3py_1' , arguments = 'et1,et2,et3' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf3py_1' in chain , "Branch `gf3py_1' is  not here!"

    # =========================================================================
    ## 22) generic fnuction with 3 arguments 
    # =========================================================================
    with timing ('generic-3D-py-2' , logger = logger ) as timer :          
        def fun3  ( x , y , z ) : return x + y + z 
        chain = data.chain 
        chain = chain.add_new_branch ( 'gf3py_2' , function = fun3 , arguments = ( 'et1',  'et2' , 'et3' ) ) 
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf3py_2' in chain , "Branch `gf3py_2' is  not here!"

    # =========================================================================
    ## 23) generic fnuction with 3 arguments 
    # =========================================================================
    with timing ('generic-3D-cxx-1' , logger = logger ) as timer :          
        fun3  = ROOT.MyTest.cpp_fun3 
        chain = data.chain 
        chain = chain.add_new_branch ( fun3 , name = 'gf3cxx_1' , arguments = 'et1,et2,et3' )         
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf3cxx_1' in chain , "Branch `gf3cxx_!' is  not here!"

    # =========================================================================
    ## 24) generic fnuction with 3 arguments 
    # =========================================================================
    with timing ('generic-3D-cxx-2' , logger = logger ) as timer :          
        fun3  = ROOT.MyTest.cpp_fun3 
        chain = data.chain 
        chain = chain.add_new_branch ( 'gf3cxx_2' , function = fun3 , arguments = ( 'et1',  'et2' , 'et3' ) ) 
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )         
    ## reload the chain and check: 
    assert 'gf3cxx_2' in chain , "Branch `gf3cxx_2' is  not here!"

    """
    
    title = 'With ALL variables'
    logger.info ( '%s:\n%s' %  (title , chain.table ( title = title , prefix = '# ' ) ) )
            
    
    title = 'CPU performance'
    logger.info ( '%s:\n%s' % ( title ,  T.table ( rows , title = title , prefix = '# ' , alignment = 'lr' ) ) ) 
        

# =============================================================================
# Many ways to add branch into TTree/Tchain
# - using string formula (TTreeFormula-based)
# - using pure python function
# - using histogram/function
# - using histogram sampling
def test_addbuffer() :
    """ Many ways to add branch into TTree/Tchain
    - using string formula (TTreeFormula-based)
    - using pure python function
    - using histogram/function
    - using histogram sampling
    """
    
    files = prepare_data ( 10 , 1000 )
    
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )

    rows = [ ( 'Method' , 'CPU [s]' ) ]

    tlen = len ( data.chain )
    dlen = tlen // 2

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import numpy
        # =====================================================================
        checked = set () 
        from ostap.trees.trees import numpy_buffer_types    
        for btype in  numpy_buffer_types :            
            if btype in checked  : continue            
            bname  = 'numpy_%s'% btype.__name__
            buffer = numpy.full ( 200 , 10 , dtype = btype )            
            with timing ('numpy %s' % btype.__name__  , logger = logger ) as timer :
                chain  = data.chain
                chain  = chain.add_new_buffer ( bname , buffer )            
            rows.append ( ( timer.name  , '%.3f' % timer.delta ) )            
            assert bname in chain , "Branch `%s' is  not here!" % bname 
            checked.add ( btype )

        title = 'With numpy buffers'
        logger.info ( '%s:\n%s' %  (title , chain.table ( title = title , prefix = '# ' ) ) )
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.warning ( "No numpy arrays!" ) 

        
    from ostap.trees.trees import array_buffer_types    
    for atype in array_buffer_types :
        
        aname  = 'array_%s'% atype
        buffer = array.array ( atype , ( 10 for i in range ( 200 ) ) )         
        with timing ('array %s' % atype , logger = logger ) as timer :
            chain  = data.chain
            chain  = chain.add_new_buffer ( aname , buffer )            
        rows.append ( ( timer.name  , '%.3f' % timer.delta ) )            
        assert aname in chain , "Branch `%s' is  not here!" % bname 

    title = 'With array buffers'
    logger.info ( '%s:\n%s' %  (title , chain.table ( title = title , prefix = '# ' ) ) )
        
    buffer = bytearray ( 200 * [10] )
    with timing ('bytes_array' , logger = logger ) as timer :
        bname  = 'bytes_array'
        chain  = data.chain
        chain  = chain.add_new_buffer ( bname , buffer )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )            
    assert bname in chain , "Branch `%s' is  not here!" % bname 

    buffer = bytes ( buffer  )
    with timing ('bytes' , logger = logger ) as timer :
        bname  = 'bytes'
        chain  = data.chain
        chain  = chain.add_new_buffer ( bname , buffer )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )            
    assert bname in chain , "Branch `%s' is  not here!" % bname 
    
    buffer = memoryview ( buffer )
    with timing ( 'memory_view' , logger = logger ) as timer :
        bname  = 'memory_view'
        chain  = data.chain
        chain  = chain.add_new_buffer ( bname , buffer )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )            
    assert bname in chain , "Branch `%s' is  not here!" % bname 

    buffer = list ( buffer  )
    with timing ( 'list' , logger = logger ) as timer :
        bname  = 'list'
        chain  = data.chain
        chain  = chain.add_new_buffer ( bname , buffer )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )            
    assert bname in chain , "Branch `%s' is  not here!" % bname 

    buffer = tuple ( buffer )
    with timing ( 'tuple' , logger = logger ) as timer :
        bname  = 'tuple'
        chain  = data.chain
        chain  = chain.add_new_buffer ( bname , buffer )
    rows.append ( ( timer.name  , '%.3f' % timer.delta ) )            
    assert bname in chain , "Branch `%s' is  not here!" % bname 

    title = 'With ALL buffers'
    logger.info ( '%s:\n%s' %  (title , chain.table ( title = title , prefix = '# ' ) ) )
            
    title = 'CPU performance'
    logger.info ( '%s:\n%s' % ( title ,  T.table ( rows , title = title , prefix = '# ' , alignment = 'lr' ) ) ) 

# =============================================================================
if '__main__' ==  __name__  :

    test_addbranch ()
    test_addbuffer ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
