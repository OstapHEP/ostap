#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_phasespace3.py
# Test module for 3-body phase space 
# ============================================================================= 
""" Test module for 3-boidy phase space 
"""
# ============================================================================= 
from   ostap.core.pyrouts     import Ostap, SE 
from   ostap.utils.gsl        import gslCount
from   ostap.logger.colorized import attention 
from   ostap.utils.root_utils import batch_env 
from   ostap.plotting.canvas  import use_canvas
from   ostap.math.models      import f1_draw
from   ostap.math.minimize    import minimize_scalar
from   ostap.utils.timing     import timing 
import ostap.logger.table     as     T
import ROOT, math, random, itertools   
# ============================================================================
from   ostap.logger.logger    import getLogger
logger = getLogger ( 'ostap.test_math_phasespace3' ) 
# ============================================================================
## set batch from environment 
batch_env ( logger )
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
#
def test_phasespace3 ( ) :
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
    logger = getLogger( 'test_phasespace3')
    logger.info ( 'Test 3-body phase space calculation via elliptic integrals' ) 
    
    masses = ( 3 , 1 , 0.1 )

    ps2 = Ostap.Math.PhaseSpace3s ( *masses ) ## analytical     symemtric 
    ps3 = Ostap.Math.PhaseSpace3a ( *masses ) ## analytical non-symmetric 
    ps4 = Ostap.Math.PhaseSpace3i ( *masses ) ## numerical    integration
    
    with use_canvas( 'test_phasespace3' , wait = 1 ) :
        
        ps2.draw    ( 'same' , xmin = ps2.threshold () , xmax = 50 , linecolor = 2 , linewidth = 2 , linestyle = 9 )
        logger.info ( 'Red   line -     symmetric expression of 3-body phase space via elliptic integrals' )
        
        ps3.draw    ( 'same' , xmin = ps2.threshold () , xmax = 50 , linecolor = 4 , linewidth = 2 , linestyle = 9 )
        logger.info ( 'Blue  line - non-symmetric expression of 3-body phase space via elliptic integrals' )
        
        ps4.draw    ( 'same' , xmin = ps2.threshold () , xmax = 50 , linecolor = 3 , linewidth = 2 , linestyle = 9 )
        logger.info ( 'Green line - explicit numerical integration' )

keep = []
# =============================================================================
def test_phasespace3i_permutations () :
    
    logger = getLogger( 'test_phasespace3i_permutations')

    masses = ( 3 , 1 , 0.1 )
    funcs  = []
    
    for p in itertools.permutations ( masses ) :
        f = lambda x : Ostap.Kinematics.phasespace3i ( x , *p  )
        funcs.append ( f )
        keep.append  ( f )
        
    from ostap.math.models import f1_draw
    
    xmin = sum ( masses ) 
    with use_canvas( 'test_phasespace3i_permutations' , wait = 1 ) :
        for i, f  in enumerate ( funcs ) :
            color = i + 1 
            if i  == 0 : f1_draw ( f ,          line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 
            else       : f1_draw ( f , 'same' , line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 


# =============================================================================
def test_phasespace3a_permutations () :
    
    logger = getLogger( 'test_phasespace3a_permutations')
    
    masses = ( 3 , 1 , 0.1 )
    funcs  = []
    
    for p in itertools.permutations ( masses ) :
        f = lambda x : Ostap.Kinematics.phasespace3a ( x , *p  )
        funcs.append ( f )
        keep.append  ( f )
        
    from ostap.math.models import f1_draw
    
    xmin = sum ( masses ) 
    with use_canvas( 'test_phasespace3i_permutations' , wait = 1 ) :
        for i, f  in enumerate ( funcs ) :
            color = i + 1 
            if i  == 0 : f1_draw ( f ,          line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 
            else       : f1_draw ( f , 'same' , line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 


# =============================================================================
def test_phasespace3s_permutations () :
    
    logger = getLogger( 'test_phasespace3_permutations')

    masses = ( 3 , 1 , 0.1 )
    funcs  = []
    
    for p in itertools.permutations ( masses ) :
        f = lambda x : Ostap.Kinematics.phasespace3s ( x , *p  )
        funcs.append ( f )
        keep.append  ( f )
                
    xmin = sum ( masses ) 
    with use_canvas( 'test_phasespace3s_permutations' , wait = 1 ) :
        for i, f  in enumerate ( funcs ) :
            color = i + 1 
            if i  == 0 : f1_draw ( f ,          line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 
            else       : f1_draw ( f , 'same' , line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 


# =============================================================================
def test_phasespace3_compare () :

    logger = getLogger( 'test_phasespace3a_compare')

    masses = ( 3 , 1 , 0.1 )

    fun1 = lambda x : Ostap.Kinematics.phasespace3   ( x , *masses ) ## default 
    fun2 = lambda x : Ostap.Kinematics.phasespace3s  ( x , *masses ) ##     symmetric form 
    fun3 = lambda x : Ostap.Kinematics.phasespace3a  ( x , *masses ) ## non-symmetric form
    fun4 = lambda x : Ostap.Kinematics.phasespace3i  ( x , *masses ) ## numerical integration 
    
    xmin = sum ( masses ) 

    with use_canvas( 'test_phasespace3_compare' , wait = 1 , keep = True ) as cnv :

        keep.append ( fun1 )
        keep.append ( fun2 )
        keep.append ( fun3 )
        keep.append ( fun4 )

        logger.info ( 'Black line -       default expression of 3-body phase' )
        logger.info ( 'Red   line -     symmetric expression of 3-body phase space via elliptic integrals' )
        logger.info ( 'Blue  line - non-symmetric expression of 3-body phase space via elliptic integrals' )
        logger.info ( 'Green line - explicit numerical integration' )

        f1_draw ( fun1                            , linewidth = 2 , xmin = xmin , xmax = 40 )
        f1_draw ( fun2 , 'same' , line_color = 2  , linewidth = 2 , xmin = xmin , xmax = 40 )
        f1_draw ( fun3 , 'same' , line_color = 4  , linewidth = 2 , xmin = xmin , xmax = 40 )
        f1_draw ( fun4 , 'same' , line_color = 3  , linewidth = 2 , xmin = xmin , xmax = 40 )

# =============================================================================
def test_phasespace3_numbers ( seed = 42 ) :
    
    logger = getLogger( 'test_phasespace3a_numbers' )
    
    masses = (3.0, 1.0, 0.1)
    xmin   = 0.999 * sum ( masses )
    xmax   = 5.000 * sum ( masses )

    rng    = random.Random(seed)
        
    for _ in range ( 100 ) :
        
        x     = rng.uniform ( xmin , xmax )
        
        v_d   = Ostap.Kinematics.phasespace3  ( x , *masses )
        v_s   = Ostap.Kinematics.phasespace3s ( x , *masses )
        v_a   = Ostap.Kinematics.phasespace3a ( x , *masses )
        v_i   = Ostap.Kinematics.phasespace3i ( x , *masses )

        values = v_d , v_s , v_a , v_i        
        vmin   = min ( values )
        vmax   = max ( values )        
        diff   = vmax - vmin 
        
        perms = { Ostap.Kinematics.phasespace3s ( x, *p ) for p in itertools.permutations ( masses) }
        
        logger.info ( f"x={x:5.1f} | i={v_i:.6e} | s={v_s:.6e} | a={v_a:.6e} | perm_spread={max(perms)-min(perms):.2e} | max-min={diff:e}" )

# =============================================================================
def test_phasespace3_timing ( n_cases = 100 , n_x = 500 , seed = 42 ):

    logger = getLogger('test_phasespace3a_timing')
    rng    = random.Random(seed)

    m1_min , m1_max = 0.0 , 5.0
    m2_min , m2_max = 0.0 , 5.0
    m3_min , m3_max = 0.0 , 5.0
    
    # Generate random physical configurations above threshold (m1 + m2 + m3).
    cases = []
    for _ in range(n_cases):
        m1        = rng.uniform ( m1_min , m1_max )
        m2        = rng.uniform ( m2_min , m2_max )
        m3        = rng.uniform ( m3_min , m3_max )
        threshold = m1 + m2 + m3
        for q in range ( n_x ) :
            x     = rng.uniform ( 0.999 * threshold , 10 * threshold ) 
            cases.append ( ( x, ( m1 , m2 , m3 ) ) )

    # Warmup before timing so one-time initialization does not skew measurements.
    for x, masses in cases:

        y     = rng.uniform ( 1.001 * threshold , 9.999 * threshold ) 
        
        Ostap.Kinematics.phasespace3  ( y , *masses )
        Ostap.Kinematics.phasespace3s ( y , *masses )
        Ostap.Kinematics.phasespace3a ( y , *masses )
        Ostap.Kinematics.phasespace3i ( y , *masses )

    with timing('phasespace3'  , logger = logger ) as timer_d:
        for x, masses in cases: Ostap.Kinematics.phasespace3s(x, *masses)
                
    with timing('phasespace3s' , logger = logger ) as timer_s:
        for x, masses in cases: Ostap.Kinematics.phasespace3s(x, *masses)

    with timing('phasespace3a' , logger = logger ) as timer_a:
        for x, masses in cases: Ostap.Kinematics.phasespace3a(x, *masses)

    with timing('phasespace3i' , logger = logger ) as timer_i:
        for x, masses in cases: Ostap.Kinematics.phasespace3i(x, *masses)

    us_d = timer_d.delta 
    us_s = timer_s.delta 
    us_a = timer_a.delta 
    us_i = timer_i.delta 

    logger.info(
        f"Timing ({len(cases)} random cases): "
        f"d={us_d:.3f}s | "
        f"s={us_s:.3f}s | "
        f"a={us_a:.3f}s | "
        f"i={us_i:.3f}s )"
    )
    
# =============================================================================
if '__main__' == __name__ :


    test_phasespace3               ()
    test_phasespace3s_permutations ()
    test_phasespace3i_permutations ()
    test_phasespace3a_permutations ()
    test_phasespace3_compare       ()
    test_phasespace3_numbers       ()
    test_phasespace3_timing        ()

# =============================================================================
##                                                                      The END 
# =============================================================================

