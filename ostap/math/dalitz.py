#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/dalitz.py
#  Set of useful "kinematic" utilities for dealing with Dalitz plot
#  @see Ostap::Kinematics::Dalitz
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-07-15
# =============================================================================
""" Set of useful ``kinematic'' utilities for Dalitz plot 
- see Ostap.Kinematics.Dalitz
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2009-09-12"
__version__ = "Version $Revision:$"
# =============================================================================
__all__     = (
    )
# =============================================================================
import ROOT, math
from   builtins             import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.dalitz' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
import ostap.math.kinematic
from   ostap.core.core      import Ostap
import ostap.math.base 

# =============================================================================
## make a graph of Dalitz plot
#  @code
#  d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
#  graph = d.graph ()
#  graph.draw ('al')
#  @endcode
def _dp_points_ ( dp , npoints = 250 ) :
    """Make a graph of Dalitz plot
    >>> d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph ()
    >>> graph.draw ('al')
    """
    
    m1  = dp.m1  ()
    m2  = dp.m2  ()
    m3  = dp.m3  ()
    M   = dp.M   () 
    s   = dp.s   ()
    sqs = dp.sqs () 

    s1_min  = ( m1 + m2 ) ** 2
    s1_max  = ( M  - m3 ) ** 2
    
    s2_min  = ( m2 + m3 ) ** 2
    s2_max  = ( M  - m1 ) ** 2

    s3_min  = ( m3 + m1 ) ** 2
    s3_max  = ( M  - m2 ) ** 2
    
    ## four reference minmax points ..
    P = [ ( s1_min , dp.s2_minmax_for_s1 ( s1_min ).first  ) ,
          ( dp.s1_minmax_for_s2 ( s2_max ).second , s2_max ) , 
          ( s1_max , dp.s2_minmax_for_s1 ( s1_max ).second ) ,  
          ( dp.s1_minmax_for_s2 ( s2_min ).first  , s2_min ) ]

    P = [ p for p in P if s1_min <= p[0] <= s1_max and s2_min <= p[1] <= s2_max ]
    
    pnts = []

    from   ostap.utils.utils    import vrange 

    ## fill branches 1 and 3 
    for v in vrange ( s1_min , s1_max , npoints ) :

        s1 = v
        s2 = dp.s2_minmax_for_s1  ( s1 )
        
        s2min  = s2.first
        s2max  = s2.second
        if not s2min < s2max : continue 
        
        x       = s1
        y1 , y2 = s2min , s2max
        
        if s1_min <= x <= s1_max :
            if s2_min <= y1 <= s2_max : pnts.append ( ( x , y1 ) )
            if s2_min <= y2 <= s2_max : pnts.append ( ( x , y2 ) )

    ## fill branches 2 and 4 :
    for v in vrange ( s2_min , s2_max , npoints ) :
                
        s2 = v 
        s1 = dp.s1_minmax_for_s2 ( s2 )
        
        s1min  = s1.first
        s1max  = s1.second
        if s1min < s1max :
            
            y       = s2
            x1, x2  = s1min , s1max 

            if s2_min <= y <= s2_max :
                if s1_min <= x1 <= s1_max : pnts.append ( ( x1 , y ) )
                if s1_min <= x2 <= s1_max : pnts.append ( ( x2 , y ) )

    pnts = P + pnts 


    ## find some point "inside" Dalitz plot
    ## first guess 
    x0 = 0.5 * ( s1_min + s1_max )
    y0 = 0.5 * ( s2_min + s2_max )
    ## find another point if needed 
    while not dp.inside ( x0 , y0 ) :
        import random 
        x0 = random.uniform ( s1_min , s1_max ) 
        y0 = random.uniform ( s2_min , s2_max ) 


    ## collect, eliminate the duplicates and ogranize all points according to phi
        
    from   ostap.math.base      import isequal
    
    points = set ()
    for p in pnts :
        x , y = p 
        dx    = x - x0
        dy    = y - y0        
        phi   = math.atan2 ( dy , dx )

        in_list = False        
        for entry in points  :
            if isequal  ( phi , entry[0] ) :
                in_list = True 
                break
            
        if in_list : continue
        
        point = phi , x , y 
        points.add ( point )

    ## convert set to the list 
    points = list ( points )        
    ## sort the list according to phi 
    points.sort ()
    ## make a closed loop  from the points
    if points and points [0][1:] != points [-1][1:] : 
        points.append ( points[0] ) 

    return  tuple( [ p[1:] for p in points ] )


# =============================================================================
## make a graph of Dalitz plot:  \f$ s_2~vs~s_1\f$ or \f$ m_{23}~vs~m_{12}\f$ 
#  @code
#  d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
#  graph = d.graph ()
#  graph.draw ('al')
#  @endcode
def _dp_graph21_ ( dp , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot: s2 vs s1 
    >>> d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph   ()
    >>> graph = d.graph21 () ## ditto
    >>> graph.draw ('al')
    """
    
    points = _dp_points_ (  dp , npoints )
    
    import ostap.histos.graphs
    
    graph = ROOT.TGraph ( len ( points ) )
    
    for i , point in enumerate ( points ) :
        x , y = point
        if masses  :
            x , y = math.sqrt ( x ) , math.sqrt ( y ) 
        graph[i] =  x, y 

    return graph 
        
# =============================================================================
## make a graph of Dalitz plot:  \f$ s_3~vs~s_1\f$ or \f$ m_{31}~vs~m_{12}\f$ 
#  @code
#  d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
#  graph = d.graph31 ()
#  graph.draw ('al')
#  @endcode
def _dp_graph31_ ( dp , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot:  s3 vs s1 
    >>> d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph31 ()
    >>> graph.draw ('al')
    """
    
    points = _dp_points_ (  dp , npoints )

    pnts   = [ ( p[0], dp.s3( p[0]  , p[1] ) ) for p in points ] 
    
    
    import ostap.histos.graphs
    
    graph = ROOT.TGraph ( len ( pnts ) )
    
    for i , point in enumerate ( pnts ) :
        x , y = point
        if masses  :
            x , y = math.sqrt ( x ) , math.sqrt ( y ) 
        graph[i] =  x, y 

    return graph 

# =============================================================================
## make a graph of Dalitz plot:  \f$ s_3~vs~s_2\f$ or \f$ m_{31}~vs~m_{23}\f$ 
#  @code
#  d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
#  graph = d.graph31 ()
#  graph.draw ('al')
#  @endcode
def _dp_graph32_ ( dp , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot:  s3 vs s2 
    >>> d = Dalitz (5 , 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph31 ()
    >>> graph.draw ('al')
    """
    
    points = _dp_points_ (  dp , npoints )

    pnts   = [ ( p[1], dp.s3( p[0]  , p[1] ) ) for p in points ] 
    
    import ostap.histos.graphs
    
    graph = ROOT.TGraph ( len ( pnts ) )
    
    for i , point in enumerate ( pnts ) :
        x , y = point
        if masses  :
            x , y = math.sqrt ( x ) , math.sqrt ( y ) 
        graph[i] =  x, y 

    return graph 
        

# ============================================================================
Ostap.Kinematics.Dalitz.graph   = _dp_graph21_ 
Ostap.Kinematics.Dalitz.graph21 = _dp_graph21_ 
Ostap.Kinematics.Dalitz.graph31 = _dp_graph31_ 
Ostap.Kinematics.Dalitz.graph32 = _dp_graph32_

# =============================================================================   
## decorated classes 
_decorated_classes_  = (
    Ostap.Kinematics.Dalitz , 
    )
## new methdods 
_new_methods_       = (
    # 
    Ostap.Kinematics.Dalitz.graph   , 
    Ostap.Kinematics.Dalitz.graph21 , 
    Ostap.Kinematics.Dalitz.graph31 , 
    Ostap.Kinematics.Dalitz.graph32 ,
    Ostap.Kinematics.Dalitz , 
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================
##                                                                      The END 
# =============================================================================
