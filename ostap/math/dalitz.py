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
    'DSwap'   , ## swap variables for Dalitz density 
    'DPlot'   , ## vizialize Dalitz density as function \f$(s_1, s_2)\f$ 
    'DPlotM'  , ## vizialize Dalitz density as functuon of \f$(m_{12}, m_{23})\f$   
    'DPlotR'  , ## vizialize Dalitz density as "rectangular plot" 
    'DPlotRM' , ## vizialize Dalitz density as "rectangular mass plot" 
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
from   ostap.core.core      import Ostap, fID 
import ostap.math.base
import ostap.histos.graphs  

# =============================================================================
## make a graph of Dalitz plot
#  @code
#  d     = Dalitz0 ( 0.1   , 0.2 , 0.3 )
#  points = d.points ( M = 5 )
#  @endcode
def _dp0_points_ ( dp , M , npoints = 250 ) :
    """Make a graph of Dalitz plot
    >>> d = Dalitz ( 0.1   , 0.2 , 0.3 )
    >>> points = d.points ( M = 5 )
    """

    m1  = dp.m1 ()
    m2  = dp.m2 ()
    m3  = dp.m3 ()

    assert m1 + m2 + m2 <= M, \
           'Dalitz0.points: Invalid mass %s>%s+%s+%s' %  ( M , m1 , m2 , m3 )
    
    s   = M * M
    
    s1_min  = dp.s1_min (   ) 
    s1_max  = dp.s1_max ( M )
    
    s2_min  = dp.s2_min (   ) 
    s2_max  = dp.s2_max ( M ) 

    s3_min  = dp.s3_min (   ) 
    s3_max  = dp.s3_max ( M ) 
    
    ## four reference minmax points ..
    P = [ ( s1_min                                      , dp.s2_minmax_for_s_s1 ( s , s1_min ).first  ) ,
          ( dp.s1_minmax_for_s_s2 ( s , s2_max ).second , s2_max ) , 
          ( s1_max                                      , dp.s2_minmax_for_s_s1 ( s , s1_max ).second ) ,  
          ( dp.s1_minmax_for_s_s2 ( s , s2_min ).first  , s2_min ) ]

    P = [ p for p in P if s1_min <= p[0] <= s1_max and s2_min <= p[1] <= s2_max ]
    
    pnts = []

    from   ostap.utils.utils    import vrange 

    ## fill branches 1 and 3 
    for v in vrange ( s1_min , s1_max , npoints ) :

        s1 = v
        s2 = dp.s2_minmax_for_s_s1  ( s , s1 )
        
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
        s1 = dp.s1_minmax_for_s_s2 ( s , s2 )
        
        s1min  = s1.first
        s1max  = s1.second
        if s1min < s1max :
            
            y       = s2
            x1, x2  = s1min , s1max 

            if s2_min <= y <= s2_max :
                if s1_min <= x1 <= s1_max : pnts.append ( ( x1 , y ) )
                if s1_min <= x2 <= s1_max : pnts.append ( ( x2 , y ) )

    pnts = P + pnts 
    pnts  = [ p for p in pnts if dp.inside ( s , p[0] , p[1] ) or abs ( dp.distance ( s , p[0] , p[1] ) ) < 1.e-8 * s ]

    ## find some point "inside" Dalitz plot
    ## first guess 
    x0 = 0.5 * ( s1_min + s1_max )
    y0 = 0.5 * ( s2_min + s2_max )
    ## find another point if needed 
    while not dp.inside ( s , x0 , y0 ) :
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

    return  tuple ( [ p[1:] for p in points ] )


# =============================================================================
## make a graph of Dalitz plot:  \f$ s_2~vs~s_1\f$ or \f$ m_{23}~vs~m_{12}\f$ 
#  @code
#  d = Dalitz0  ( 0.1   , 0.2 , 0.3 )
#  graph = d.graph   ( M = 5 )
#  graph = d.graph21 ( M = 5 ) ## ditto 
#  graph.draw ('al')
#  @endcode
def _dp0_graph21_ ( dp , M , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot: s2 vs s1 
    >>> d = Dalitz0 ( 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph   ( M = 5 )
    >>> graph = d.graph21 () ## ditto
    >>> graph.draw ('al')
    """
    assert 0 <= M and dp.s_min() <= M * M ,\
           'Dalitz.graph21: Invalid mass %s' % M  
    
    points = _dp0_points_ (  dp , M , npoints )
    
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
#  d = Dalitz0 ( 0.1   , 0.2 , 0.3 )
#  graph = d.graph31 ( M = 5 )
#  graph.draw ('al')
#  @endcode
def _dp0_graph31_ ( dp , M , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot:  s3 vs s1 
    >>> d = Dalitz0 ( 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph31 ( M = 5 )
    >>> graph.draw ('al')
    """
    
    points = _dp0_points_ (  dp , M , npoints )

    s      = M * M 
    pnts   = [ ( p[0], dp.s3( s , p[0]  , p[1] ) ) for p in points ] 
    
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
#  d = Dalitz0 ( 0.1   , 0.2 , 0.3 )
#  graph = d.graph32 ( M = 5 )
#  graph.draw ('al')
#  @endcode
def _dp0_graph32_ ( dp , M , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot:  s3 vs s2 
    >>> d = Dalitz0 (5 , 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph31 ( M = 5  )
    >>> graph.draw ('al')
    """
    
    points = _dp0_points_ (  dp , M , npoints )

    s      = M * M 
    pnts   = [ ( p[1], dp.s3( s , p[0]  , p[1] ) ) for p in points ] 
    
    import ostap.histos.graphs
    
    graph = ROOT.TGraph ( len ( pnts ) )
    
    for i , point in enumerate ( pnts ) :
        x , y = point
        if masses  :
            x , y = math.sqrt ( x ) , math.sqrt ( y ) 
        graph[i] =  x, y 

    return graph 

        

# =============================================================================
## make a graph of Dalitz plot:  \f$ s_2~vs~s_1\f$ or \f$ m_{23}~vs~m_{12}\f$ 
#  @code
#  d = Dalitz  ( 5 , 0.1   , 0.2 , 0.3 )
#  graph = d.graph   ()
#  graph = d.graph21 () ## ditto 
#  graph.draw ('al')
#  @endcode
def _dp_graph21_ ( dp , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot: s2 vs s1 
    >>> d = Dalitz0 ( 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph   ( M = 5 )
    >>> graph = d.graph21 () ## ditto
    >>> graph.draw ('al')
    """
    dp0 = Ostap.Kinematics.Dalitz0 ( dp ) 
    return _dp0_graph21_ ( dp0 , M = dp.M () , npoints = npoints , masses = masses  )


# =============================================================================
## make a graph of Dalitz plot:  \f$ s_3~vs~s_1\f$ or \f$ m_{31}~vs~m_{12}\f$ 
#  @code
#  d = Dalitz ( 5.0 , 0.1   , 0.2 , 0.3 )
#  graph = d.graph31 ()
#  graph.draw ('al')
#  @endcode
def _dp_graph31_ ( dp , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot:  s3 vs s1 
    >>> d = Dalitz ( 5.0 , 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph31 ()
    >>> graph.draw ('al')
    """
    dp0 = Ostap.Kinematics.Dalitz0 ( dp ) 
    return _dp0_graph31_ ( dp0 , M = dp.M () , npoints = npoints , masses = masses )


# =============================================================================
## make a graph of Dalitz plot:  \f$ s_3~vs~s_2\f$ or \f$ m_{31}~vs~m_{23}\f$ 
#  @code
#  d = Dalitz ( 5.0 , 0.1   , 0.2 , 0.3 )
#  graph = d.graph31 ()
#  graph.draw ('al')
#  @endcode
def _dp_graph32_ ( dp , npoints = 250 , masses = False ) :
    """Make a graph of Dalitz plot:  s3 vs s2 
    >>> d = Dalitz ( 5 , 0.1   , 0.2 , 0.3 )
    >>> graph = d.graph31 ()
    >>> graph.draw ('al')
    """
    dp0 = Ostap.Kinematics.Dalitz0 ( dp ) 
    return _dp0_graph32_ ( dp0 , M = dp.M () , npoints = npoints , masses = masses )


# ============================================================================
Ostap.Kinematics.Dalitz0.graph   = _dp0_graph21_ 
Ostap.Kinematics.Dalitz0.graph21 = _dp0_graph21_ 
Ostap.Kinematics.Dalitz0.graph31 = _dp0_graph31_ 
Ostap.Kinematics.Dalitz0.graph32 = _dp0_graph32_

Ostap.Kinematics.Dalitz .graph   = _dp_graph21_ 
Ostap.Kinematics.Dalitz .graph21 = _dp_graph21_ 
Ostap.Kinematics.Dalitz .graph31 = _dp_graph31_ 
Ostap.Kinematics.Dalitz .graph32 = _dp_graph32_

# =============================================================================
## @class DSwap
#  Swap variables for the Dalitz plot density
class DSwap( object ) :
    """Swap variables for the Dalitz plot density
    """
    def __init__ ( self   ,
                   i1     ,
                   i2     , 
                   f3     , ## f2(s,s1,s2)
                   dalitz ) :
        
        assert isinstance ( i1 , int ) and 1 <= i1 <= 3 , 'Invalid i1 %s' % i1 
        assert isinstance ( i2 , int ) and 1 <= i2 <= 3 , 'Invalid i2 %s' % i2
        assert i1 != i2 , 'Invalid i1/i2 %s/%s' %  ( i1 , i2 )     
        
        if   1 == i1 :
            if   2 == i2 : swap = lambda d , x1 , x2 , s : ( x1 , x2  )
            elif 3 == i2 : swap = lambda d , x1 , x3 , s : ( x1 , d.s2  ( s , x1 , x3 ) ) 
        elif 2 == i1 :
            if   1 == i2 : swap = lambda d , x2 , x1 , s : ( x1 , x2  ) 
            elif 3 == i2 : swap = lambda d , x2 , x3 , s : ( d.s1  ( s , x2 , x3 ) , x2 ) 
        elif 3 == i1 :
            if   1 == i2 : swap = lambda d , x3 , x1 , s : ( x1 , d.s2  ( s , x1 , x3 ) ) 
            elif 2 == i2 : swap = lambda d , x3 , x2 , s : ( d.s1 ( s , x2 , x3 ) , x2  )

        self.__swap   = swap
        self.__dalitz = dalitz
        self.__f3     = f3 

    def __call__ ( self , s , x1 , x2 ) :

        d = self.__dalitz
        
        s1 , s2 = self.__swap ( d , x1 , x2 , s )
        
        if not d.inside ( s , s1 , s2 ) : return 0 ##  needed?
        
        return self.__f3 ( s , s1 , s2 )
    
# =============================================================================
## @class DPlot
#  Helper class to vizualize the Dalitz density via conversion to ROOT.TF2
#  @code
#  dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
#  density =  ...
#  dp = DPlot ( density , dalitz , s = 5 )
#  f2 = dp.tf2 () 
#  @endcode 
class DPlot(object) :
    """Helper class to vizualize the Dalitz density via conversion to ROOT.TF2
    >>> dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
    >>> density =  ...
    >>> dp = DPlot ( density , dalitz , s = 5 )
    >>> f2 = dp.tf2 () 
    """    
    
    def __init__ ( self        ,
                   f3          ,
                   dalitz      ,
                   s           ) : 
        self.__dalitz = dalitz
        self.__f3     = f3

        assert dalitz.s_min () <= s , 'DPlot: Invalid value of s %s' % s
        
        self.__s      = s
        self.__M      = s**0.5 

    ## helper method for using with ROOT.TF2 
    def eval ( self , X ) : return self ( X [0] , X [1] )

    ## the main method 
    def __call__ ( self , s1 , s2 ) :
        """The main method"""
        s  = self.__s 
        if not self.__dalitz.inside ( s , s1 , s2 ) : return 0
        return self.__f3 ( s , s1 , s2 ) 

    ## convert to TF2 
    def tf2  ( self        ,
               xmin = None ,
               xmax = None ,
               ymin = None ,
               ymax = None ) :

        d = self.__dalitz
        s = self.__s  
        M = self.__M  
        
        xmin = xmin if not xmin is None else d.s1_min (   )
        xmax = xmax if not xmax is None else d.s1_max ( M )
        ymin = ymin if not ymin is None else d.s2_min (   )
        ymax = ymax if not ymax is None else d.s2_max ( M )

        assert xmin < xmax, 'Invalid xmin/xmax %s/%s' %  ( xmin , xmax )
        assert ymin < ymax, 'Invalid ymin/ymax %s/%s' %  ( ymin , ymax )
        
        self.__tf = ROOT.TF2 ( fID () , self.eval , xmin , xmax , ymin , ymax )
        return self.__tf
        
    ## get the dalitz plot boundary
    def graph ( self ) :
        """Get Dalitz plot boundary
        """
        gr = _dp0_graph21_ ( self.dalitz , M = self.M )
        self.__gr =  gr
        gr.SetLineWidth ( 2 ) 
        return gr
    
    @property
    def s (  self ) :
        """``s'' : s-parameter"""
        return self.__s
    @s.setter
    def s ( self , value ) :
        v = float ( value ) 
        assert self.dalitz.s_min () <= v , 'Invalid value for s %s' % v 
        self.__s = v
        self.__M = v ** 0.5
        
    @property
    def M ( self ) :
        """``sM'' : <-parameter"""
        return self.__M        
    @property
    def dalitz (  self ) :
        """``dalitz'' : Dalitz configuration"""
        return self.__dalitz
    @property
    def f3 (  self ) :
        """``f3'' : the function f3(s,s1,s2)"""
        return self.__f3 

    # =========================================================================
    ## "transpose it!"
    def transpose  ( self , i1 , i2 ) :
        """Transpose it!"""

        f = DSwap                 ( i1 , i2 , self.f3 , self.dalitz )
        d = self.dalitz.transpose ( i1 , i2 ) 
        ##
        T = type ( self )
        return T ( f , d , self.s )

        


# =============================================================================
## @class DPlotM
#  Visualize the Dalitz density as function of\f$(m_{12},m_{23})\f$ masses 
#  @code
#  dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
#  density =  ...
#  dp = DPlotM ( density , dalitz , s = 5 )
#  f2 = dp.tf2 () 
#  @endcode 
class DPlotM(DPlot) :
    """Helper class to vizualize the Dalitz density as function of (m12,m23) via conversion to ROOT.TF2
    >>> dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
    >>> density =  ...
    >>> dp = DPlotM ( density , dalitz , s = 5 )
    >>> f2 = dp.tf2 () 
    """    

    ## the main method 
    def __call__ ( self , m12 , m23 ) :
        """The main method"""
        s  = self.s
        d  = self.dalitz 
        s1 = m12 * m12
        s2 = m23 * m23
        if not d.inside ( s , s1 , s2 ) : return 0
        jacob = 0.25 / ( m12 * m23 ) 
        return jacob * self.f3 ( s , s1 , s2 ) 

    ## convert to TF2 
    def tf2  ( self        ,
               xmin = None ,
               xmax = None ,
               ymin = None ,
               ymax = None ) :

        d = self.dalitz
        s = self.s  
        M = self.M
        
        xmin = xmin if not xmin is None else d.s1_min (   ) ** 0.5 
        xmax = xmax if not xmax is None else d.s1_max ( M ) ** 0.5 
        ymin = ymin if not ymin is None else d.s2_min (   ) ** 0.5 
        ymax = ymax if not ymax is None else d.s2_max ( M ) ** 0.5 

        assert xmin < xmax, 'Invalid xmin/xmax %s/%s' %  ( xmin , xmax )
        assert ymin < ymax, 'Invalid ymin/ymax %s/%s' %  ( ymin , ymax )
        
        self.__tf = ROOT.TF2 ( fID () , self.eval , xmin , xmax , ymin , ymax )        
        return self.__tf

    ## get the dalitz plot boundary
    def graph ( self ) :
        """Get Dalitz plot boundary
        """
        gr = _dp0_graph21_ ( self.dalitz , M = self.M , masses = True )
        self.__gr =  gr
        gr.SetLineWidth ( 2 ) 
        return gr

        
# =============================================================================
## @class DPlotR
#  Visualize Dalitz density as function of ``rectangular variables''
#  \f[ \left( \begin{array}{l} x_1 \\ x_2 \end{array} \right) \equiv 
#  \left( \begin{array}{l} \cos \theta^{12}_{R23} \\ s_2 \end{array} \right) \f]
#  @code
#  dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
#  density =  ...
#  dp = DPlotR ( density , dalitz , s = 5 )
#  f2 = dp.tf2 () 
#  @endcode
#  @see Ostap::Kinematics::Dalitz0::x1
#  @see Ostap::Kinematics::Dalitz0::x2
#  @see Ostap::Kinematics::Dalitz0::x2s
#  @see Ostap::Kinematics::Dalitz0::J
class DPlotR(DPlot) :
    """Helper class to vizualize the Dalitz density using 'rectangular variables'
    -  x1 : cos theta^{12}_{R23}
    -  x2 : s_2 
    >>> dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
    >>> density =  ...
    >>> dp = DPlotR ( density , dalitz , s = 5 )
    >>> f2 = dp.tf2 ()
    - see Ostap.Kinematics.Dalitz0.x1
    - see Ostap.Kinematics.Dalitz0.x2
    - see Ostap.Kinematics.Dalitz0.x2s
    - see Ostap.Kinematics.Dalitz0.J    
    """    
    def __init__ ( self        ,
                   f3          ,
                   dalitz      ,
                   s           ,
                   T = True    ) :

        DPlot.__init__ ( self , f3 , dalitz , s )
        self.__T = True if T else False 

    ## the main method 
    def __call__ ( self , x1 , x2 ) :
        """The main method"""

        ## swap variables 
        if self.__T :
            x2 , x2 = x1 , x2
            
        s       = self.s
        d       = self.dalitz 
        s1 , s2 = d.x2s ( s , x1 , x2 ) 
        if not d.inside ( s , s1 , s2 ) : return 0
        jacob   = d.J   ( s , s1 , s2 )  
        return jacob * self.f3 ( s , s1 , s2 ) 

    ## convert to TF2 
    def tf2  ( self        ,
               xmin = None ,
               xmax = None  ,
               ymin = None ,
               ymax = None ) :

        d = self.dalitz
        s = self.s  
        M = self.M

        if not self.__T : 
            xmin = xmin if not xmin is None else -1.0 
            xmax = xmax if not xmax is None else +1.0
            ymin = ymin if not ymin is None else d.s2_min (   ) 
            ymax = ymax if not ymax is None else d.s2_max ( M )
        else :
            xmin = xmin if not xmin is None else d.s2_min (   ) 
            xmax = xmax if not xmax is None else d.s2_max ( M )
            ymin = ymin if not ymin is None else -1.0 
            ymax = ymax if not ymax is None else +1.0
            
        assert xmin < xmax, 'Invalid xmin/xmax %s/%s' %  ( xmin , xmax )
        assert ymin < ymax, 'Invalid ymin/ymax %s/%s' %  ( ymin , ymax )
        
        self.__tf = ROOT.TF2 ( fID () , self.eval , xmin , xmax , ymin , ymax )
        return self.__tf

    @property
    def T  ( self ) :
        """``T'' : swap variables/transpose ? """
        return  self.__T 

    ## get the dalitz plot boundary
    def graph ( self ) :
        """Get Dalitz plot boundary
        """
        
        gr = ROOT.TGraph ( 5 )
        d  = self.dalitz
        M  = self.M
        
        gr [ 0 ] = -1.0 , d.s2_min (   )
        gr [ 1 ] = -1.0 , d.s2_max ( M )
        gr [ 2 ] =  1.0 , d.s2_max ( M )
        gr [ 3 ] =  1.0 , d.s2_min (   )
        gr [ 4 ] = -1.0 , d.s2_min (   )

        if self.T : gr = gr.T ()
        
        gr.SetLineWidth ( 2 ) 
        self.__gr =  gr
        return gr

# =============================================================================
## @class DPlotRM
#  Visualize Dalitz density as function of ``rectangular variables''
#  \f[ \left( \begin{array}{l} x_1 \\ x_2 \end{array} \right) \equiv 
#  \left( \begin{array}{l} \cos \theta^{12}_{R23} \\ m_{23} \end{array} \right) \f]
#  @code
#  dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
#  density =  ...
#  dp = DPlotR ( density , dalitz , s = 5 )
#  f2 = dp.tf2 () 
#  @endcode
class DPlotRM(DPlotR) :
    """Helper class to vizualize the Dalitz density using 'rectangular variables'
    -  x1 : cos theta^{12}_{R23}
    -  x2 : m_23 
    >>> dalitz  = Ostap.Kinematics.Dalitz0 ( 0.1 , 0.2 , 0.3 ) 
    >>> density =  ...
    >>> dp = DPlotRM ( density , dalitz , s = 5 )
    >>> f2 = dp.tf2 ()
    """    

    ## the main method 
    def __call__ ( self , x1 , x2 ) :
        """The main method"""
        
        ## swap variables 
        if self.T :
            m23 , x1 = x1 , x2 

        s       = self.s
        d       = self.dalitz
        s2      = m23 * m23 
        s1 , _  = d.x2s ( s , x1 , s2 ) 
        if not d.inside ( s , s1 , s2 ) : return 0
        jacob   = d.J   ( s , s1 , s2 ) / ( 2 * m23  )
        return jacob * self.f3 ( s , s1 , s2 )

    ## convert to TF2 
    def tf2  ( self        ,
               xmin = None ,
               xmax = None ,
               ymin = None ,
               ymax = None ) :

        d = self.dalitz
        s = self.s  
        M = self.M
        
        if not self.T : 
            xmin = xmin if not xmin is None else -1.0 
            xmax = xmax if not xmax is None else +1.0
            ymin = ymin if not ymin is None else d.s2_min (   ) ** 0.5 
            ymax = ymax if not ymax is None else d.s2_max ( M ) ** 0.5 
        else :
            xmin = xmin if not xmin is None else d.s2_min (   ) ** 0.5 
            xmax = xmax if not xmax is None else d.s2_max ( M ) ** 0.5 
            ymin = ymin if not ymin is None else -1.0 
            ymax = ymax if not ymax is None else +1.0

        assert xmin < xmax, 'Invalid xmin/xmax %s/%s' %  ( xmin , xmax )
        assert ymin < ymax, 'Invalid ymin/ymax %s/%s' %  ( ymin , ymax )
        
        self.__tf = ROOT.TF2 ( fID () , self.eval , xmin , xmax , ymin , ymax )
        return self.__tf

    ## get the dalitz plot boundary
    def graph ( self ) :
        """Get Dalitz plot boundary
        """
        
        gr = ROOT.TGraph ( 5 )
        d  = self.dalitz
        M  = self.M
        
        gr [ 0 ] = -1.0 , d.s2_min (   ) ** 0.5
        gr [ 1 ] = -1.0 , d.s2_max ( M ) ** 0.5 
        gr [ 2 ] =  1.0 , d.s2_max ( M ) ** 0.5 
        gr [ 3 ] =  1.0 , d.s2_min (   ) ** 0.5 
        gr [ 4 ] = -1.0 , d.s2_min (   ) ** 0.5 

        if self.T : gr = gr.T ()
        
        gr.SetLineWidth ( 2 ) 
        self.__gr =  gr
        return gr

# =============================================================================   
## decorated classes 
_decorated_classes_  = (
    Ostap.Kinematics.Dalitz  , 
    Ostap.Kinematics.Dalitz0 , 
    )
## new methdods 
_new_methods_       = (
    # 
    Ostap.Kinematics.Dalitz0.graph   , 
    Ostap.Kinematics.Dalitz0.graph21 , 
    Ostap.Kinematics.Dalitz0.graph31 , 
    Ostap.Kinematics.Dalitz0.graph32 ,
    Ostap.Kinematics.Dalitz .graph   , 
    Ostap.Kinematics.Dalitz .graph21 , 
    Ostap.Kinematics.Dalitz .graph31 , 
    Ostap.Kinematics.Dalitz .graph32 ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
