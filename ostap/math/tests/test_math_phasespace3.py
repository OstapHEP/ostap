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
import ostap.logger.table     as     T
from   ostap.utils.utils      import wait
from   ostap.plotting.canvas  import use_canvas
from   ostap.math.models      import f1_draw
from   ostap.math.minimize    import minimize_scalar 
import ROOT, math, random, itertools   
# ============================================================================
from   ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_phasespace3' ) 
else                       : logger = getLogger ( __name__                      )
# ============================================================================

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
    """Test 3-body phase space calculation via elliptic integrals
    
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
    
    masses = ( 3 , 1 , 0.1 )

    ps1 = Ostap.Math.PhaseSpace3  ( *masses ) 
    ps2 = Ostap.Math.PhaseSpace3s ( *masses ) ## <--- HERE
    
    ps3 = lambda x : Ostap.Kinematics.phasespace3a  ( x , *masses ) ## non-symmetric form 

    with wait ( 3 ), use_canvas( 'test_phasespace3' ) :
        ps1.draw (          xmin = ps1.threshold()     , xmax = 50 , linecolor=2 , linewidth = 2 )
        logger.info ( 'Red   line - 3-body phase space via numerical integration' ) 
        ps2.draw ( 'same' , xmin = ps1.threshold()     , xmax = 50 , linecolor=4 , linewidth = 2 , linestyle = 9 )
        logger.info ( 'Blue  line - symmetric     expression of 3-body phase space via elliptic integrals' )        
        f1_draw  ( ps3 , 'same' , xmin = ps1.threshold ()  , xmax = 50 , linecolor=8 , linewidth = 2 , linestyle = 11 )
        logger.info ( 'Green line - non-symmetric expression of 3-body phase space via elliptic integrals' ) 
        

# =============================================================================
def test_phasespace3i_permutations () :
    
    masses = ( 3 , 1 , 0.1 )
    funcs  = []
    
    for p in itertools.permutations ( masses ) :
        f = lambda x : Ostap.Kinematics.phasespace3i ( x , *p  )
        funcs.append ( f )
        
    from ostap.math.models import f1_draw
    
    xmin = sum ( masses ) 
    with wait ( 3 ), use_canvas( 'test_phasespace3i_permutations' ) :
        for i, f  in enumerate ( funcs ) :
            color = i + 1 
            if i  == 0 :
                f1_draw ( f ,          line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 
            else: 
                f1_draw ( f , 'same' , line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 


# =============================================================================
def test_phasespace3a_permutations () :
    
    masses = ( 3 , 1 , 0.1 )
    funcs  = []
    
    for p in itertools.permutations ( masses ) :
        f = lambda x : Ostap.Kinematics.phasespace3a ( x , *p  )
        funcs.append ( f )
        
    from ostap.math.models import f1_draw
    
    xmin = sum ( masses ) 
    with wait ( 3 ), use_canvas( 'test_phasespace3i_permutations' ) :
        for i, f  in enumerate ( funcs ) :
            color = i + 1 
            if i  == 0 :
                f1_draw ( f ,          line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 
            else: 
                f1_draw ( f , 'same' , line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 


# =============================================================================
def test_phasespace3s_permutations () :
    
    masses = ( 3 , 1 , 0.1 )
    funcs  = []
    
    for p in itertools.permutations ( masses ) :
        f = lambda x : Ostap.Kinematics.phasespace3s ( x , *p  )
        funcs.append ( f )
        
    xmin = sum ( masses ) 
    with wait ( 3 ), use_canvas( 'test_phasespace3s_permutations' ) :
        for i, f  in enumerate ( funcs ) :
            color = i + 1 
            if i  == 0 :
                f1_draw ( f ,          line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 
            else: 
                f1_draw ( f , 'same' , line_color = color , linewidth = 2 , xmin = xmin , xmax = 50 ) 

# =============================================================================
def test_phasespace3_compare () :
## if 1 < 2 : 

    masses = ( 3 , 1 , 0.1 )

    fun1 = lambda x : Ostap.Kinematics.phasespace3i  ( x , *masses ) ## numerical integration 
    fun2 = lambda x : Ostap.Kinematics.phasespace3s  ( x , *masses ) ## symmetric form 
    fun3 = lambda x : Ostap.Kinematics.phasespace3a  ( x , *masses ) ## non-symmetric form 
    fun4 = lambda x : Ostap.Kinematics.phasespace3nr ( x , *masses ) ## non-relativistic limit 
    
    xmin = sum ( masses ) 
        
    with wait ( 3 ), use_canvas( 'test_phasespace3_compare' ) :

        for i, f in enumerate ( ( fun1 , fun2 , fun3 , fun4 ) ) :

            color = i + 2
            
            if i  == 0 :
                f1_draw ( f ,          line_color = color , linewidth = 2 , xmin = xmin , xmax = 40 ) 
            else: 
                f1_draw ( f , 'same' , line_color = color , linewidth = 2 , xmin = xmin , xmax = 40 )

    
# =============================================================================
if '__main__' == __name__ :


    test_phasespace3               ()
    test_phasespace3s_permutations ()
    test_phasespace3i_permutations ()
    test_phasespace3a_permutations ()
    test_phasespace3_compare       ()



# =============================================================================
##                                                                      The END 
# =============================================================================

