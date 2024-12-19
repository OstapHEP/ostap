#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_carlson.py
# Test module for Carlon symmetric forms 
# ============================================================================= 
""" Test module for Carlson symmetric forms 
"""
# ============================================================================= 
from   ostap.core.pyrouts     import Ostap, SE 
from   ostap.utils.gsl        import gslCount
from   ostap.logger.colorized import attention 
import ostap.logger.table     as     T
from   ostap.utils.utils      import wait
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.utils      import vrange 
import ROOT, math, random  
# ============================================================================
from   ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_carlson' ) 
else                       : logger = getLogger ( __name__                  )
# ============================================================================

RF     = Ostap.Math.carlson_RF
RC     = Ostap.Math.carlson_RC
RJ     = Ostap.Math.carlson_RJ
RD     = Ostap.Math.carlson_RD
RG     = Ostap.Math.carlson_RG

RJ_gsl = Ostap.Math.carlson_RJ_gsl 
RD_gsl = Ostap.Math.carlson_RD_gsl 
RF_gsl = Ostap.Math.carlson_RF_gsl 
RC_gsl = Ostap.Math.carlson_RC_gsl 
RG_gsl = Ostap.Math.carlson_RG_gsl 

RG_int = Ostap.Math.carlson_RG_int
RF_int = Ostap.Math.carlson_RF_int
RD_int = Ostap.Math.carlson_RD_int
RC_int = Ostap.Math.carlson_RC_int
RJ_int = Ostap.Math.carlson_RJ_int

prec_TIGHT = 1.e-13
prec_LOOSE = 1.e-9

NTEST      = 100000

# =============================================================================
## Check reference tabulated values,
def test_carlson_values ( ) :
    """ Test predefined values, (section 3), arXiv:math/9409227
    """

    logger = getLogger( 'test_carlson_values')
    logger.info ( 'Test predefined values, (section 3), arXiv:math/9409227' ) 
    
    test_RF = [
        ( 'RF'      , RF     , (1,2,0) , 1.3110287771461    ) ,
        ( 'RF'      , RF     , (2,3,4) , 0.58408284167715   ) ,
        ( 'RF'      , RF     , (1,2,4) , 0.6850858166334359 ) , ## extra  
        ##
        ( 'RF_int'  , RF_int , (1,2,0) , 1.3110287771461  ) ,
        ( 'RF_int'  , RF_int , (2,3,4) , 0.58408284167715 ) ,
        ( 'RF_int'  , RF_int , (1,2,4) , 0.6850858166334359 ) , ## extra  
        
        ##
        ( 'RF_gsl'  , RF_gsl , (1,2,0) , 1.3110287771461  ) ,
        ( 'RF_gsl'  , RF_gsl , (2,3,4) , 0.58408284167715 ) ,
        ( 'RF_gsl'  , RF_gsl , (1,2,4) , 0.6850858166334359 ) , ## extra  
        ## 
        ( 'RF2'     , RF     , (1,2  ) , 1.3110287771461  ) , ## 2-argument form 
        ]
    
    test_RC = [
        ( 'RC'      , RC     , ( 0        , 0.25 ) , math.pi              ) ,
        ( 'RC'      , RC     , ( 0.25 * 9 , 2    ) , math.log(2.0)        ) ,
        ( 'RC'      , RC     , ( 0.25     , -2   ) , math.log(2.0) / 3.0  ) ,
        ##
        ( 'RC_gsl'  , RC_gsl , ( 0        , 0.25 ) , math.pi              ) ,
        ( 'RC_gsl'  , RC_gsl , ( 0.25 * 9 , 2    ) , math.log(2.0)        ) ,
        ( 'RC_gsl'  , RC_gsl , ( 0.25     , -2   ) , math.log(2.0) / 3.0  ) ,
        ##
        ( 'RC_int'  , RC_int , ( 0        , 0.25 ) , math.pi              ) ,
        ( 'RC_int'  , RC_int , ( 0.25 * 9 , 2    ) , math.log(2.0)        ) ,
        ( 'RC_int'  , RC_int , ( 0.25     , -2   ) , math.log(2.0) / 3.0  ) ,
        ]
    
    test_RJ = [
        ( 'RJ'      , RJ     , ( 0 , 1 , 2, 3 ) , 0.77688623778582 ) ,
        ( 'RJ'      , RJ     , ( 2 , 3 , 4, 5 ) , 0.14297579667157 ) ,
        ##
        ( 'RJ_gsl'  , RJ_gsl , ( 0 , 1 , 2, 3 ) , 0.77688623778582 ) ,
        ( 'RJ_gsl'  , RJ_gsl , ( 2 , 3 , 4, 5 ) , 0.14297579667157 ) ,
        ##
        ( 'RJ_int'  , RJ_int , ( 0 , 1 , 2, 3 ) , 0.77688623778582 ) ,
        ( 'RJ_int'  , RJ_int , ( 2 , 3 , 4, 5 ) , 0.14297579667157 ) ,
        ]
    
    test_RD = [
        ( 'RD'     , RD     , ( 0 , 2 , 1 ) ,  1.7972103521034  ) ,
        ( 'RD'     , RD     , ( 2 , 3 , 4 ) ,  0.16510527294261 ) ,
        ##
        ( 'RD_gsl' , RD_gsl , ( 0 , 2 , 1 ) ,  1.7972103521034  ) ,
        ( 'RD_gsl' , RD_gsl , ( 2 , 3 , 4 ) ,  0.16510527294261 ) ,
        ##
        ( 'RD_int' , RD_int , ( 0 , 2 , 1 ) ,  1.7972103521034  ) ,
        ( 'RD_int' , RD_int , ( 2 , 3 , 4 ) ,  0.16510527294261 ) ,
        ]

    test_RG = [
        ( 'RG'     , RG     , ( 0 , 16      , 16 ) ,  math.pi         ) ,
        ( 'RG'     , RG     , ( 2 , 3       , 4  ) ,  1.7255030280692 ) ,
        ( 'RG'     , RG     , ( 0 , 0.0796  , 4  ) ,  1.0284758090288 ) ,
        ##
        ( 'RG_gsl' , RG_gsl , ( 0 , 16      , 16 ) ,  math.pi         ) ,
        ( 'RG_gsl' , RG_gsl , ( 2 , 3       , 4  ) ,  1.7255030280692 ) ,
        ( 'RG_gsl' , RG_gsl , ( 0 , 0.0796  , 4  ) ,  1.0284758090288 ) ,
        ##
        ( 'RG_int' , RG_int , ( 0 , 16      , 16 ) ,  math.pi         ) ,
        ( 'RG_int' , RG_int , ( 2 , 3       , 4  ) ,  1.7255030280692 ) ,
        ( 'RG_int' , RG_int , ( 0 , 0.0796  , 4  ) ,  1.0284758090288 ) ,
        ##
        ( 'RG2'    , RG , ( 16      , 16 ) ,  math.pi         ) , ## 2-argument form 
        ( 'RG2'    , RG , ( 0.0796  , 4  ) ,  1.0284758090288 ) , ## 2-argument form         
        ]

    rows = [ ( 'Function' , 'Arguments' , 'Result' , 'Expected' , 'abs-delta' , 'rel-delta' ) ]
    
    ad_max = -1
    rd_max = -1

    with gslCount () : 
        for test in test_RF + test_RC + test_RJ + test_RD + test_RG :
            
            name , fun , args, r = test 
            result = fun ( *args )
            
            ad = abs ( result - r     )
            rd = abs ( result / r - 1 ) 

            at , rt = '%.4g' % ad , '%.4g' % rd
            
            if  prec_TIGHT < ad : at = attention ( at ) 
            if  prec_TIGHT < rd : rt = attention ( rt ) 
            
            row = name , str(args) , '%+.12f' % result , '%+.12f' % r , at , rt 
            rows.append ( row )
            ad_max = max ( ad_max , ad ) 
            rd_max = max ( rd_max , rd ) 

    table = T.table ( rows ,
                      title = 'Test of Carlson forms'  ,
                      prefix = '# ' , alignment = 'llllll' )
    
    logger.info ( 'Test Carlson forms:\n%s' % table ) 

    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )


# =============================================================================
## Compare local/GSL and plain integration methods 
def test_carlson_cmp ( ) :
    """ Compare local/GSL and plain integration methods
    """

    logger = getLogger( 'test_carlson_1')
    logger.info ( 'Compare local/GSL and plain integration methods  ' ) 
    
    cF1 = SE()
    cF2 = SE()
    
    cD1 = SE()
    cD2 = SE()
    
    cC1 = SE()
    cC2 = SE()

    cJ1 = SE()
    cJ2 = SE()

    cG1 = SE()
    cG2 = SE()

    with gslCount () :
        
        for i in range ( 10000 ) :
            
            x = random.uniform ( 0 , 1 )  
            y = random.uniform ( 0 , 1 ) 
            z = random.uniform ( 0 , 1 ) 
            p = random.uniform ( 0 , 1 )
            
            ## RF
            
            v1   = RF     ( x , y , z )
            v2   = RF_gsl ( x , y , z )
            v3   = RF_int ( x , y , z )
            cF1 += ( v1 / v2 ) - 1 
            cF2 += ( v1 / v3 ) - 1 
            
            v1   = RF     ( x , y , p )
            v2   = RF_gsl ( x , y , p )
            v3   = RF_int ( x , y , p )
            cF1 += ( v1 / v2 ) - 1 
            cF2 += ( v1 / v3 ) - 1 
            
            v1   = RF     ( x , p , z )
            v2   = RF_gsl ( x , p , z )
            v3   = RF_int ( x , p , z )
            cF1 += ( v1 / v2 ) - 1 
            cF2 += ( v1 / v3 ) - 1 
            
            v1   = RF     ( p , x , z )
            v2   = RF_gsl ( p , x , z )
            v3   = RF_int ( p , x , z )
            cF1 += ( v1 / v2 ) - 1 
            cF2 += ( v1 / v3 ) - 1 
            
            ## RJ 
            
            v1   = RJ     ( x , y , z , p )
            v2   = RJ_gsl ( x , y , z , p )
            v3   = RJ_int ( x , y , z , p )
            cJ1 += ( v1 / v2 ) - 1
            cJ2 += ( v1 / v3 ) - 1
            
            v1   = RJ     ( x , y , p , z )
            v2   = RJ_gsl ( x , y , p , z )
            v3   = RJ_int ( x , y , p , z )
            cJ1 += ( v1 / v2 ) - 1
            cJ2 += ( v1 / v3 ) - 1
            
            v1   = RJ     ( x , p , z , y )
            v2   = RJ_gsl ( x , p , z , y )
            v3   = RJ_int ( x , p , z , y )
            cJ1 += ( v1 / v2 ) - 1
            cJ2 += ( v1 / v3 ) - 1
            
            v1   = RJ     ( p , y , z , x )
            v2   = RJ_gsl ( p , y , z , x )
            v3   = RJ_int ( p , y , z , x )
            cJ1 += ( v1 / v2 ) - 1
            cJ2 += ( v1 / v3 ) - 1
            
            ## RD 
            
            v1   = RD     ( x , y , z )
            v2   = RD_gsl ( x , y , z )
            v3   = RD_int ( x , y , z )
            cD1 += ( v1 / v2 ) - 1
            cD2 += ( v1 / v3 ) - 1
            
            v1   = RD     ( x , y , p )
            v2   = RD_gsl ( x , y , p )
            v3   = RD_int ( x , y , p )
            cD1 += ( v1 / v2 ) - 1
            cD2 += ( v1 / v3 ) - 1
            
            v1   = RD     ( x , p , z )
            v2   = RD_gsl ( x , p , z )
            v3   = RD_int ( x , p , z )
            cD1 += ( v1 / v2 ) - 1
            cD2 += ( v1 / v3 ) - 1
            
            v1   = RD     ( p , y , z )
            v2   = RD_gsl ( p , y , z )
            v3   = RD_int ( p , y , z )
            cD1 += ( v1 / v2 ) - 1
            cD2 += ( v1 / v3 ) - 1
            
            ## RC 
            
            v1   = RC     ( x , y )
            v2   = RC_gsl ( x , y )
            v3   = RC_int ( x , y )
            cC1 += ( v1 / v2 ) - 1
            cC2 += ( v1 / v2 ) - 1
            
            v1   = RC     ( x , z )
            v2   = RC_gsl ( x , z )
            v3   = RC_int ( x , z )
            cC1 += ( v1 / v2 ) - 1
            cC2 += ( v1 / v2 ) - 1

            v1   = RC     ( x , p )
            v2   = RC_gsl ( x , p )
            v3   = RC_int ( x , p )
            cC1 += ( v1 / v2 ) - 1
            cC2 += ( v1 / v2 ) - 1
            
            v1   = RC     ( y , z )
            v2   = RC_gsl ( y , z )
            v3   = RC_int ( y , z )
            cC1 += ( v1 / v2 ) - 1
            cC2 += ( v1 / v2 ) - 1
            
            v1   = RC     ( y , p )
            v2   = RC_gsl ( y , p )
            v3   = RC_int ( y , p )
            cC1 += ( v1 / v2 ) - 1
            cC2 += ( v1 / v2 ) - 1
            
            v1   = RC     ( z , p )
            v2   = RC_gsl ( z , p )
            v3   = RC_int ( z , p )
            cC1 += ( v1 / v2 ) - 1
            cC2 += ( v1 / v2 ) - 1
            
            ## RG
            
            v1   = RG     ( x , y , z )
            v2   = RG_gsl ( x , y , z )
            v3   = RG_int ( x , y , z )
            cG1 += ( v1 / v2 ) - 1
            cG2 += ( v1 / v2 ) - 1
            
            v1   = RG     ( x , y , p )
            v2   = RG_gsl ( x , y , p )
            v3   = RG_int ( x , y , p )
            cG1 += ( v1 / v2 ) - 1
            cG2 += ( v1 / v2 ) - 1
            
            v1   = RG     ( x , p , z )
            v2   = RG_gsl ( x , p , z )
            v3   = RG_int ( x , p , z )
            cG1 += ( v1 / v2 ) - 1
            cG2 += ( v1 / v2 ) - 1
            
            v1   = RG     ( p , y , z )
            v2   = RG_gsl ( p , y , z )
            v3   = RG_int ( p , y , z )
            cG1 += ( v1 / v2 ) - 1
            cG2 += ( v1 / v2 ) - 1
            
    rows = [ ( 'Name' , '#' , 'Mean' , 'rms' , 'min' , 'max' ) ]

    for n , c in [ ( 'RF' , cF1 ) ,
                   ( 'RJ' , cJ1 ) ,
                   ( 'RD' , cD1 ) ,
                   ( 'RC' , cC1 ) ,
                   ( 'RG' , cG1 ) ] :
        
        mean  = c.mean ()
        tmean = '%+.5g' % mean
        if prec_TIGHT < abs ( mean ) : tmean = attention ( tmean )

        rms   = c.rms  ()
        trms  = '%+.5g' % rms 
        if prec_TIGHT < abs ( rms  ) : trms  = attention ( trms )

        vmin  = c.min  ()
        tmin  = '%+.5g' % vmin
        if prec_TIGHT < abs ( vmin ) : tmin  = attention ( tmin )

        vmax  = c.max  ()
        tmax  = '%+.5g' % vmax
        if prec_TIGHT < abs ( vmax ) : tmax  = attention ( tmax )

        row =  n , '%d' % c.nEntries() , tmean , trms , tmin , tmax 
        rows.append ( row )

    table = T.table ( rows ,
                      title = 'Test of Carlson forms, Local vs GSL '  ,
                      prefix = '# ' , alignment = 'lrllll' )
    
    logger.info ( 'Test Carlson forms, Local vs GSL:\n%s' % table ) 

    rows = [ ( 'Name' , '#' , 'Mean' , 'rms' , 'min' , 'max' ) ]
    for n , c in [ ( 'RF' , cF2 ) ,
                   ( 'RJ' , cJ2 ) ,
                   ( 'RD' , cD2 ) ,
                   ( 'RC' , cC2 ) ,
                   ( 'RG' , cG2 ) ] :
        
        mean  = c.mean ()
        tmean = '%+.5g' % mean
        if prec_TIGHT < abs ( mean ) : tmean = attention ( tmean )

        rms   = c.rms  ()
        trms  = '%+.5g' % rms 
        if prec_TIGHT < abs ( rms  ) : trms  = attention ( trms )

        vmin  = c.min  ()
        tmin  = '%+.5g' % vmin
        if prec_TIGHT < abs ( vmin ) : tmin  = attention ( tmin )

        vmax  = c.max  ()
        tmax  = '%+.5g' % vmax
        if prec_TIGHT < abs ( vmax ) : tmax  = attention ( tmax )

        row =  n , '%d' % c.nEntries() , tmean , trms , tmin , tmax 
        rows.append ( row )

    table = T.table ( rows ,
                      title = 'Test of Carlson forms, Local vs plain integration'  ,
                      prefix = '# ' , alignment = 'lrllll' )
    
    logger.info ( 'Test Carlson forms, Local vs plain integrtation:\n%s' % table ) 

# =============================================================================
## Test expression for complete elliptic integral \f$ K(k) \f$ via Carlson's symmetric forms
#  @see Ostap::Math::elliptic_K 
#  @see Ostap::Math::carlson_RF
#  @see Eq. (55) in arXiv:math/9409227
def test_carlson_K ( ) :
    """ Test expression for complete elliptic integral E(k)` via Carlson's symmetric forms
    - see Ostap.Math.elliptic_K
    - see Ostap.Math.carlson_RF 
    - see Eq. (55) in arXiv:math/9409227
    """
    
    logger = getLogger( 'test_carlson_K')
    logger.info ( 'Test expression for K(k) via Carlson Forms' ) 
    
    from ostap.math. models import f1_draw
    
    def k1 ( k ) : return Ostap.Math.elliptic_K     ( k ) 
    def k2 ( k ) : return Ostap.Math.carlson_RF     ( 0 , 1-k*k , 1 ) 
    def k3 ( k ) : return Ostap.Math.carlson_RF_gsl ( 0 , 1-k*k , 1 ) 
    def k4 ( k ) : return Ostap.Math.carlson_RF_int ( 0 , 1-k*k , 1 ) 
    def k5 ( k ) : return Ostap.Math.elliptic_K_gsl ( k ) 

    with wait ( 3 ), use_canvas( 'test_carlson_K' ) :
        f1_draw ( k1 ,          xmin = 0 , xmax = 1-1.e-7 , min = 0 , linecolor = 2 , linewidth = 2 )
        f1_draw ( k2 , 'same' , xmin = 0 , xmax = 1-1.e-7 , min = 0 , linecolor = 4 , linewidth = 2 , linestyle = 9 )
        f1_draw ( k3 , 'same' , xmin = 0 , xmax = 1-1.e-7 , min = 0 , linecolor = 8 , linewidth = 2 , linestyle = 9 )
        f1_draw ( k4 , 'same' , xmin = 0 , xmax = 1-1.e-7 , min = 0 , linecolor = 5 , linewidth = 2 , linestyle = 9 )
        f1_draw ( k5 , 'same' , xmin = 0 , xmax = 1-1.e-7 , min = 0 , linecolor = 6 , linewidth = 2 , linestyle = 9 )

        logger.info ( 'Red     line : K (k) complete elliptic integral' ) 
        logger.info ( "Blue    line : K (k) expressed via symmetric Carlson's RF function"     ) 
        logger.info ( "Green   line : K (k) expressed via symmetric Carlson's RF function/GSL" ) 
        logger.info ( "Yellow  line : K (k) expressed via symmetric Carlson's RF function/int" ) 
        logger.info ( "Magenta line : K (k) expressed using GSL" ) 

# =============================================================================
## Test expression for complete elliptic integral \f$ E(k) \f$ via Carlson's symmetric forms
#  @see Ostap::Math::elliptic_E 
#  @see Ostap::Math::carlson_RG
def test_carlson_E ( ) :
    """ Test expression for complete elliptic integral E(k)` via Carlson's symmetric forms
    - see Ostap.Math.elliptic_E 
    - see Ostap.Math.carlson_RG
    """
    
    logger = getLogger( 'test_carlson_E')
    logger.info ( 'Test expression for E(k) via Carlson Forms' ) 

    from ostap.math. models import f1_draw
    
    def e1 ( k ) : return   Ostap.Math.elliptic_E ( k ) 
    def e2 ( k ) : return 2*Ostap.Math.carlson_RG     ( 1-k*k , 1 , 0 ) 
    def e3 ( k ) : return 2*Ostap.Math.carlson_RG_gsl ( 1-k*k , 1 , 0 ) 
    def e4 ( k ) : return 2*Ostap.Math.carlson_RG_int ( 1-k*k , 1 , 0 ) 
    def e5 ( k ) : return   Ostap.Math.elliptic_E ( k ) 

    with wait ( 3 ), use_canvas( 'test_carlson_E' ) :
        f1_draw ( e1 ,          xmin = 0 , xmax = 1 , min = 0 , linecolor = 2 , linewidth = 2 )
        f1_draw ( e2 , 'same' , xmin = 0 , xmax = 1 , min = 0 , linecolor = 4 , linewidth = 2 , linestyle = 9 )
        f1_draw ( e3 , 'same' , xmin = 0 , xmax = 1 , min = 0 , linecolor = 8 , linewidth = 2 , linestyle = 9 )
        f1_draw ( e4 , 'same' , xmin = 0 , xmax = 1 , min = 0 , linecolor = 5 , linewidth = 2 , linestyle = 9 )
        f1_draw ( e5 , 'same' , xmin = 0 , xmax = 1 , min = 0 , linecolor = 6 , linewidth = 2 , linestyle = 9 )

        logger.info ( 'Red     line : E (k) complete elliptic integral' ) 
        logger.info ( "Blue    line : E (k) expressed via symmetric Carlson's RG function" ) 
        logger.info ( "Green   line : E (k) expressed via symmetric Carlson's RG function/GSL" ) 
        logger.info ( "Yellow  line : E (k) expressed via symmetric Carlson's RG function/int" ) 
        logger.info ( "Magenta line : E (k) expressed using GSL" ) 

# =============================================================================
## Test expressions for complete elliptic integral \f$ K(k)-E(k) \f$ via Carlson's symmetric forms
#  @see Ostap::Math::elliptic_K
#  @see Ostap::Math::elliptic_E
#  @see Ostap::Math::elliptic_KmE
#  @see Ostap::Math::elliptic_RD
#  @see Ostap::Math::elliptic_RF
#  @see Ostap::Math::elliptic_RG
def test_carlson_KmE ( ) :
    """ Test expressions for complete elliptic integral `K(k)-E(k)` via Carlson's symmetric forms 
    - see Ostap.Math.elliptic_K
    - see Ostap.Math.elliptic_E
    - see Ostap.Math.elliptic_KmE
    - see Ostap.Math.elliptic_RD
    - see Ostap.Math.elliptic_RF
    - see Ostap.Math.elliptic_RG
    """
    
    logger = getLogger( 'test_carlson_KmE')
    logger.info ( 'Test expression for K(k)-E(k) via Carlson Forms' ) 
    
    from ostap.math. models import f1_draw
    
    def e1 ( k ) : return Ostap.Math.elliptic_K ( k ) - Ostap.Math.elliptic_E ( k ) 
    def e2 ( k ) : return k * k * Ostap.Math.carlson_RD ( 0 , 1 - k * k,1) / 3.0  
    def e3 ( k ) : return Ostap.Math.carlson_RF ( 1-k*k,1) - 2*Ostap.Math.carlson_RG(1-k*k,1) 
    def e4 ( k ) : return Ostap.Math.elliptic_KmE ( k )  

    with wait ( 3 ), use_canvas( 'test_carlson_E' ) :
        f1_draw ( e1 ,          xmin = 0 , xmax = 1-1.e-10 , min = 0 , linecolor = 2 , linewidth = 2 )
        logger.info ( 'Red    line : K(k) - E(k) as they are ' ) 
        f1_draw ( e2 , 'same' , xmin = 0 , xmax = 1-1.e-10 , min = 0 , linecolor = 4 , linewidth = 2 , linestyle = 9 )
        logger.info ( "Blue   line : K(k) - E(k) expression via Carlson's RD" ) 
        f1_draw ( e3 , 'same' , xmin = 0 , xmax = 1-1.e-10 , min = 0 , linecolor = 8 , linewidth = 2 , linestyle = 3 )
        logger.info ( "Green  line : K(k) - E(k) expression via Carlson's RF and RG" ) 
        f1_draw ( e4 , 'same' , xmin = 0 , xmax = 1-1.e-10 , min = 0 , linecolor = 5 , linewidth = 2 , linestyle = 4 )
        logger.info ( "Yellow line : K(k) - E(k) expression via Carlson's RD " ) 
        
        
# =============================================================================
## Test 3-body phase space calculation via elliptic integrals
#  @see Ostap::Math::PhaseSpace3
#  @see Ostap::Math::PhaseSpace3s
#  @see Ostap::Kinematics::phasespace3
#  @see https://indico.cern.ch/event/368497/contributions/1786992/attachments/1134067/1621999/davydychev.PDF
#  @see http://cds.cern.ch/record/583358/files/0209233.pdf
#  @see https://www.researchgate.net/publication/2054534_Three-body_phase_space_symmetrical_treatments
#
#  @see A.Davydychev and R.Delbourgo,
#       "Explicitly symmetrical treatment of three body phase space",
#       J.Phys. A37 (2004) 4871, arXiv:hep-th/0311075",
#       doi = 10.1088/0305-4470/37/17/016
#  @see https://arxiv.org/abs/hep-th/0311075
#  @see https://iopscience.iop.org/article/10.1088/0305-4470/37/17/016
def test_carlson_PS3 ( ) :
    """ Test 3-body phase space calculation via elliptic integrals
    
    - see Ostap.Math.PhaseSpace3
    - see Ostap.Math.PhaseSpace3s
    
    - see Ostap.Kinematics.phasespace3
    - see https://indico.cern.ch/event/368497/contributions/1786992/attachments/1134067/1621999/davydychev.PDF
    - see http://cds.cern.ch/record/583358/files/0209233.pdf
    - see https://www.researchgate.net/publication/2054534_Three-body_phase_space_symmetrical_treatments
    
    - see A.Davydychev and R.Delbourgo, ``Explicitly symmetrical treatment of three body phase space'',
    J.Phys. A37 (2004) 4871, arXiv:hep-th/0311075,
    doi = 10.1088/0305-4470/37/17/016
    - see https://arxiv.org/abs/hep-th/0311075
    - see https://iopscience.iop.org/article/10.1088/0305-4470/37/17/016    
    """
    logger = getLogger( 'test_carlson_PS3')
    logger.info ( 'Test 3-body phase space calculation via elliptic integrals' ) 
    
    ps1 = Ostap.Math.PhaseSpace3  ( 3 , 1 , 0.1 ) 
    ps2 = Ostap.Math.PhaseSpace3s ( 3 , 1 , 0.1 ) ## <--- HERE

    with wait ( 3 ), use_canvas( 'test_carlson_PS3' ) :
        ps1.draw (          xmin = ps1.threshold() , xmax = 50 , linecolor=2 , linewidth = 2 )
        logger.info ( 'Red  line - 3-body phase space via numerical integration' ) 
        ps2.draw ( 'same' , xmin = ps2.threshold() , xmax = 50 , linecolor=4 , linewidth = 2 )
        logger.info ( 'Blue line - analytic expression of 3-body phase space via elliptic integrals' ) 

    xmin = min ( ps1.threshold() , ps2.threshold() )
    xmax = 50
    
    def fun1 ( x ) : return abs ( ps1 ( x ) - ps2 ( x ) )
    def fun2 ( x ) : return       ps1 ( x ) 
    from ostap.math.integral import integral
    i1 = integral ( fun1 , xmin = xmin , xmax = xmax , err = True ) 
    i2 = integral ( fun2 , xmin = xmin , xmax = xmax , err = True ) 
    i1 = i1 ** 0.5
    i2 = i2 ** 0.5
    ii = i1 / i2     
    logger.error ( 'Difference is %s [Analytic expression is wrong!]' % ii ) 
    
# =============================================================================
## Test identity Eq.(49) arXiv:math/9409227' )
#  @see Ostap::Math::carlson_RF 
def test_carlson_Eq49 ( ) :
    """ Test identity Eq.(49) arXiv:math/9409227' )
    - see Ostap::Math.carlson_RF 
    """

    logger = getLogger( 'test_carlson_Eq49') 
    logger.info ( 'Test identity Eq.(49) from arXiv:math/9409227' ) 

    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
        
        x = random.uniform ( 0 , 100 )
        y = random.uniform ( 0 , 100 )
        l = random.uniform ( 0 , 100 )
        m = x*y/l
        
        left  = RF ( x + l , y + l , l ) + RF ( x + m , y + m , m )
        right = RF ( x , y , 0  )
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            
            row = str ( (x,y,l,m) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                  '%.5g' % ad       , \
                  '%.5g' % rd 
            rows.append ( row )
            
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 
                
    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(49) from arXiv:math/9409227; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )

# =============================================================================
## Test identity Eq (51) arXiv:math/9409227' )
#  @see Ostap::Math::carlson_RJ 
def test_carlson_Eq51 ( ) :
    """ Test identity Eq (49) arXiv:math/9409227' )
    - see Ostap::Math.carlson_RF 
    """

    logger = getLogger( 'test_carlson_Eq51') 
    logger.info ( 'Test identity Eq.(451) from arXiv:math/9409227' ) 


    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
            
        x = random.uniform ( 0 , 100 )
        y = random.uniform ( 0 , 100 )
        p = random.uniform ( 0 , 100 )
        l = random.uniform ( 0 , 100 )
        m = x*y/l
        
        a = p*p*(l+m+x+y)
        b = p*(p+l)*(p+m)
        
        left  = RJ ( x + l , y + l , l , p + l ) + RJ ( x + m , y + m , m , p + m )
        right = RJ ( x , y , 0 , p ) - 3*RC ( a , b ) 
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            row = str ( (x,y,p,l,m) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                  '%.5g' % ad       , \
                  '%.5g' % rd 
            rows.append ( row )
        
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 

            
    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(51) from arXiv:math/9409227; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )


# =============================================================================
## Test identity Eq (53) arXiv:math/9409227' )
#  @see Ostap::Math::carlson_RD 
def test_carlson_Eq53 ( ) :
    """ Test identity Eq.(53) arXiv:math/9409227' )
    - see Ostap::Math.carlson_RD 
    """

    logger = getLogger( 'test_carlson_Eq53') 
    logger.info ( 'Test identity Eq.(53) from arXiv:math/9409227' ) 



    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
        
        x = random.uniform ( 0 , 100 )
        y = random.uniform ( 0 , 100 )
        p = random.uniform ( 0 , 100 )
        l = random.uniform ( 0 , 100 )
        m = x*y/l
        
        left  = RD ( l , x + l , y + l ) + RD ( m , x+m , y + m )
        right = RD ( 0 , x , y ) - 3.0 / ( y * math.sqrt( x + y + l + m ))
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            row = str ( (x,y,p,l,m) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                  '%.5g' % ad       , \
                  '%.5g' % rd 
            rows.append ( row )
            
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 

            
    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(53) from arXiv:math/9409227; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )


# =============================================================================
## Test identity Eq.(54) arXiv:math/9409227' )
#  @see Ostap::Math::carlson_RD
def test_carlson_Eq54 ( ) :
    """ Test identity Eq.(54) arXiv:math/9409227' )
    - see Ostap::Math.carlson_RD 
    """

    logger = getLogger( 'test_carlson_Eq54') 
    logger.info ( 'Test identity Eq.(54) from arXiv:math/9409227' ) 


    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
        
        x = random.uniform ( 0 , 100 )
        y = random.uniform ( 0 , 100 )
        z = random.uniform ( 0 , 100 )
        
        left  = RD ( x , y, z ) + RD ( y , z , x ) + RD ( z , x , y ) 
        right = 3.0 / math.sqrt ( x * y * z ) 
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            row = str ( (x,y,z) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                      '%.5g' % ad       , \
                      '%.5g' % rd 
            rows.append ( row )
            
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 


    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(54) from arXiv:math/9409227; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )



# =============================================================================
## Test identity Eq.(19.21.1) from https://dlmf.nist.gov/19 
#  @see Ostap::Math::carlson_RD
#  @see Ostap::Math::carlson_RF
def test_carlson_eq19211 ( ) :
    """ Test identity Eq.(19.21.1) https://dlmf.nist.gov/19' 
    - see Ostap::Math.carlson_RD 
    - see Ostap::Math.carlson_RF 
    """

    logger = getLogger( 'test_carlson_eq19211') 
    logger.info ( 'Test identity Eq.(19.21.1) https://dlmf.nist.gov/19' ) 


    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
        
        z = random.uniform ( 0 , 100 )
        
        left  = RF ( 0 , 1 + z , z ) * RD ( 0 , 1 + z , 1 ) + RD ( 0 , 1 + z , z ) * RF ( 0 , 1 + z , 1 ) 
        right = 3*math.pi /(2*z) 
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            row = str ( (z,) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                      '%.5g' % ad       , \
                      '%.5g' % rd 
            rows.append ( row )
            
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 


    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(19.21.1) from https://dlmf.nist.gov/19; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )


# =============================================================================
## Test identity Eq.(19.21.2) from https://dlmf.nist.gov/19 
#  @see Ostap::Math::carlson_RF
#  @see Ostap::Math::carlson_RD
def test_carlson_eq19212 ( ) :
    """ Test identity Eq.(19.21.2) https://dlmf.nist.gov/19' 
    - see Ostap::Math.carlson_RD 
    - see Ostap::Math.carlson_RF 
    """

    logger = getLogger( 'test_carlson_eq19212') 
    logger.info ( 'Test identity Eq.(19.21.2) https://dlmf.nist.gov/19' ) 


    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
        
        y = random.uniform ( 0 , 100 )
        z = random.uniform ( 0 , 100 )
        
        left  = 3*RF (  0 , y , z ) 
        right = z * RD ( 0 , y , z )  + y * RD ( 0 , z , y ) 
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            row = str ( (z,) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                      '%.5g' % ad       , \
                      '%.5g' % rd 
            rows.append ( row )
            
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 


    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(19.21.2) from https://dlmf.nist.gov/19; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )

# =============================================================================
## Test identity Eq.(19.21.3) from https://dlmf.nist.gov/19 
#  @see Ostap::Math::carlson_RF
#  @see Ostap::Math::carlson_RD
def test_carlson_eq19213 ( ) :
    """ Test identity Eq.(19.21.3) https://dlmf.nist.gov/19' 
    - see Ostap::Math.carlson_RD 
    - see Ostap::Math.carlson_RF 
    - see Ostap::Math.carlson_RG 
    """

    logger = getLogger( 'test_carlson_eq19213') 
    logger.info ( 'Test identity Eq.(19.21.3) https://dlmf.nist.gov/19' ) 


    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
        
        y = random.uniform ( 0 , 100 )
        z = random.uniform ( 0 , 100 )
        
        left  = 6 * RG (  0 , y , z ) 
        right = 3 * z * RF( 0 , y , z ) + z * ( y - z ) * RD ( 0 , y , z ) 
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            row = str ( (z,) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                      '%.5g' % ad       , \
                      '%.5g' % rd 
            rows.append ( row )
            
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 


    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(19.21.3) from https://dlmf.nist.gov/19; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )



# =============================================================================
## Test identity Eq.(19.21.7) from https://dlmf.nist.gov/19 
#  @see Ostap::Math::carlson_RF
#  @see Ostap::Math::carlson_RD
def test_carlson_eq19217 ( ) :
    """ Test identity Eq.(19.21.7) https://dlmf.nist.gov/19' 
    - see Ostap::Math.carlson_RD 
    - see Ostap::Math.carlson_RF 
    """

    logger = getLogger( 'test_carlson_eq19217') 
    logger.info ( 'Test identity Eq.(19.21.7) https://dlmf.nist.gov/19' ) 


    ad_max = -1
    rd_max = -1
    
    rows = [ ( 'Arguments' , 'Left', 'Right' , 'abs-delta'  , 'rel-delta' ) ]
    
    for i in range ( NTEST ) : 
        
        x = random.uniform ( 0 , 100 )
        y = random.uniform ( 0 , 100 )
        z = random.uniform ( 0 , 100 )
        
        left  = ( x - y ) * RD ( y , z , x )  +  ( z - y )  * RD ( x , y , z )
        right = 3 * RF ( x , y , z ) - 3 * math.sqrt ( y / ( x* z ) )
        
        ad = abs ( left  - right      )
        rd = abs ( left  / right  - 1 )
        
        if prec_TIGHT < ad or prec_TIGHT < rd : 
            row = str ( (z,) ) , \
                  '%+.12f' % left   , \
                  '%+.12f' % right  , \
                      '%.5g' % ad       , \
                      '%.5g' % rd 
            rows.append ( row )
            
        ad_max = max ( ad_max , ad ) 
        rd_max = max ( rd_max , rd ) 


    row = '' , '' , '' , '%.5g' % ad_max , '%.5g' % rd_max 
    rows.append ( row )

    title = 'Test Eq.(19.21.7) from https://dlmf.nist.gov/19; bad values from %d tests ' % NTEST
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    if   max ( ad_max , rd_max ) < prec_TIGHT :
        logger.info     ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    elif max ( ad_max , rd_max ) < prec_LOOSE  :
        logger.warning  ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )
    else : 
        logger.error    ('Maximal differences are %.5g/%.5g (abs/rel)' % ( ad_max , rd_max ) )

        
# =============================================================================
## Test incomplete Elliptic integral \f$ F(\phi,k) \f$ 
#  @see Ostap::Math::elliptic_F
#  @see Ostap::Math::elliptic_F_gsl
def test_elliptic_F () :
    """ Test incomplete Elliptic integral F(phi,k) 
    - see Ostap.Math.elliptic_F
    - see Ostap.Math.elliptic_F_gsl
    """
    
    from ostap.math.integral import integral
    rows = [ ('k' , 'diff' ) ]
    for k in vrange ( 0.01 , 0.99 , 20 ) :
        def myfun ( x ) : return abs ( Ostap.Math.elliptic_F     ( x , k ) -
                                       Ostap.Math.elliptic_F_gsl ( x , k ) )
        
        d = integral ( myfun , xmin = 0 , xmax = 20 )
        item = '%.5g' % d 
        if abs ( d ) > prec_TIGHT : item = attention ( item )        
        row = '%.5g' % k , item        
        rows.append ( row )

    title = 'Tesf for Elliptic_F, local vs GSL'
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'll' )    
    logger.info ( '%s\n%s' % ( title , table ) )
    

# =============================================================================
## Test incomplete Elliptic integral \f$ E(\phi,k) \f$ 
#  @see Ostap::Math::elliptic_E
#  @see Ostap::Math::elliptic_E_gsl
def test_elliptic_E () :
    """ Test incomplete Elliptic integral  E(phi,k)  
    - see Ostap.Math.elliptic_E
    - see Ostap.Math.elliptic_E_gsl
    """
    
    from ostap.math.integral import integral
    rows = [ ('k' , 'diff' ) ]
    for k in vrange ( 0.01 , 0.99 , 20 ) :
        def myfun ( x ) : return abs ( Ostap.Math.elliptic_E     ( x , k ) -
                                       Ostap.Math.elliptic_E_gsl ( x , k ) )


        d = integral ( myfun , xmin = 0 , xmax = 20 )
        item = '%.5g' % d 
        if abs ( d ) > prec_TIGHT : item = attention ( item )        
        row = '%.5g' % k , item        
        rows.append ( row )

    title = 'Tesf for Elliptic_E, local vs GSL'
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'll' )    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    

# =============================================================================
## Test elliptic integral \f$ Fm ( am ( u , m ) , m ) \equiv u \f$ 
#  @see Ostap::Math::elliptic_Fm
#  @see Ostap::Math::am
def test_elliptic_Fm () :
    """ Test elliptic integral Fm ( am ( u , m ) , m ) ==  u 
    - see Ostap.Math.elliptic_Fm
    - see Ostap.Math.am
    """
    
    from ostap.math.integral import integral
    rows = [ ( 'm' , 'diff' ) ]
    
    for m in vrange ( 0.01 , 0.99 , 20 ) :
        def myfun ( x ) : return abs ( Ostap.Math.elliptic_Fm ( Ostap.Math.am ( x , m ) , m ) - x )

        d = integral ( myfun , xmin = 0 , xmax = 20 )
        item = '%.5g' % d 
        if abs ( d ) > prec_TIGHT : item = attention ( item )        
        row = '%.5g' % m , item        
        rows.append ( row )


    title = 'Tesf for Fm(am(x,M),m)==x'
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'll' )    
    logger.info ( '%s\n%s' % ( title , table ) )

# =============================================================================
## Test elliptic functions \f$ sn == \sin am \f$
#  @see Ostap::Math::sn
#  @see Ostap::Math::am
def test_elliptic_sn_am () :
    """ Test elliptic functions sn == sin am 
    - see Ostap.Math.sn
    - see Ostap.Math.am
    """
    
    from ostap.math.integral import integral
    rows = [ ( 'm' , 'diff' ) ]
    
    for m in vrange ( 0.01 , 0.99 , 20 ) :
        def myfun ( x ) : return abs ( Ostap.Math.sn ( x , m ) - math.sin ( Ostap.Math.am ( x , m ) ) )
        
        d = integral ( myfun , xmin = 0 , xmax = 20 )
        item = '%.5g' % d 
        if abs ( d ) > prec_TIGHT : item = attention ( item )        
        row = '%.5g' % m , item        
        rows.append ( row )

    title = 'Tesf for sn = sin am '
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'll' )    
    logger.info ( '%s\n%s' % ( title , table ) )

# =============================================================================
## Test elliptic functions \f$ sn^2+cn^2==1f$
#  @see Ostap::Math::sn
#  @see Ostap::Math::cn
def test_elliptic_sn_cn () :
    """ Test elliptic functions sn^2+cn^2==1
    - see Ostap.Math.sn
    - see Ostap.Math.cn
    """
    
    from ostap.math.integral import integral
    rows = [ ( 'm' , 'diff' ) ]
    
    for m in vrange ( 0.01 , 0.99 , 20 ) :
        def myfun ( x ) : return abs ( Ostap.Math.sn ( x , m )**2 + Ostap.Math.cn ( x , m )**2 - 1 )
                                               
        d = integral ( myfun , xmin = 0 , xmax = 20 )
        item = '%.5g' % d 
        if abs ( d ) > prec_TIGHT : item = attention ( item )        
        row = '%.5g' % m , item        
        rows.append ( row )

    title = 'Tesf for sn^2+cn^2==1 '
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'll' )    
    logger.info ( '%s\n%s' % ( title , table ) )

# =============================================================================
## Test elliptic functions \f$ m sn^2 + dn^2==1f$
#  @see Ostap::Math::sn
#  @see Ostap::Math::dn
def test_elliptic_sn_dn () :
    """ Test elliptic functions sn^2+ m * dn^2==1 
    - see Ostap.Math.sn
    - see Ostap.Math.dn
    """
    from ostap.math.integral import integral
    rows = [ ( 'm' , 'diff' ) ]
    
    for m in vrange ( 0.01 , 0.99 , 20 ) :
        def myfun ( x ) : return abs ( m * Ostap.Math.sn ( x , m )**2 + Ostap.Math.dn ( x , m )**2 - 1 )

        d = integral ( myfun , xmin = 0 , xmax = 20 )
        item = '%.5g' % d 
        if abs ( d ) > prec_TIGHT : item = attention ( item )        
        row = '%.5g' % m , item        
        rows.append ( row )

    title = 'Tesf for m*sn^2+ dn^2==1 '
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'll' )    
    logger.info ( '%s\n%s' % ( title , table ) )

# =============================================================================
## Test elliptic functions \f$ dn = d/du am f$
#  @see Ostap::Math::dn
#  @see Ostap::Math::am
def test_elliptic_dn_am () :
    """ Test elliptic functions dn = d/du am 
    - see Ostap.Math.dm
    - see Ostap.Math.am
    """
    from ostap.math.integral   import integral
    from ostap.math.derivative import Derivative 
    rows = [ ( 'm' , 'diff' ) ]

    for m in vrange ( 0.01 , 0.99 , 20 ) :

        def am ( x ) : return Ostap.Math.am ( x , m )
        der_am = Derivative ( am ) 
        
        def myfun ( x ) : return abs ( der_am ( x ) - Ostap.Math.dn ( x , m ) )
                                               
        d = integral ( myfun , xmin = 0 , xmax = 20 )
        item = '%.5g' % d 
        if abs ( d ) > prec_TIGHT * 30 : item = attention ( item )        
        row = '%.5g' % m , item        
        rows.append ( row )

    title = 'Tesf for dn = d/du am '
    table = T.table ( rows , title  = title , prefix = '# ' , alignment = 'll' )    
    logger.info ( '%s\n%s' % ( title , table ) )
    
    
# =============================================================================
if '__main__' == __name__ :

    test_carlson_values  () 
    test_carlson_cmp     () 
    test_carlson_K       ()
    test_carlson_E       ()
    test_carlson_KmE     ()    
    test_carlson_PS3     ()

    test_carlson_Eq49    ()
    test_carlson_Eq51    ()
    test_carlson_Eq53    ()
    test_carlson_Eq54    ()

    test_carlson_eq19211 ()
    test_carlson_eq19212 ()
    test_carlson_eq19213 ()
    test_carlson_eq19217 ()
    
    test_elliptic_F      ()
    test_elliptic_E      ()
    test_elliptic_Fm     ()
    test_elliptic_sn_am  ()
    test_elliptic_sn_cn  ()
    test_elliptic_sn_dn  ()
    test_elliptic_dn_am  ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================


