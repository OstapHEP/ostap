#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/plotting/graph_summary.py
# Prepare "summary" plot 
# =============================================================================
"""Prepare ``summary'' plot 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-04-11"
# =============================================================================
__all__     = (
    'DrawConfig'    , ## helper base class  to keep draw attributes
    'Limit'         , ## graphical representation of the upper/lower limit (arrow) 
    'Record'        , ## graphical representation of data poins with several unceratities
    'Average'       , ## graphical representation of average  (box)
    'Average'       , ## graphical representation of average  (box)
    'graph_summary' , ## prepare summary graph 
    'draw_summary'  , ## draw summary graph 
    )
# =============================================================================
import ROOT
from   ostap.core.ostap_types import num_types, string_types  
from   ostap.core.core        import VE, hID
from   ostap.utils.cidict     import cidict 
from   ostap.utils.utils      import vrange  
import ostap.histos.graphs
# =============================================================================
##  @class DrawConfig
#   Helper base class to keep the configuration of drawing attributes 
class DrawConfig(object) :
    """Helper base class to keep the configuration of drawingattributes 
    """
    
    def __init__ ( self , **config ) :
        self.__config = cidict ( transform = lambda k : k.lower().replace('_','' ) )
        self.__config.update   ( config )
         
    @property
    def config ( self ) :
        """Configuration: line/fill/marker color/style/size/..."""
        return self.__config

# =============================================================================
## @class  Limit
#  A graphical represenation of upper/lower limit as
#  arrow from <code>limit</endcode> to <code>zero</endcode>
#  @code
#  l = Limit ( 0.25 , line_color = 3  , arrow_size = 0.15 , ,arrow_style = '->' )  
#  @endcode 
class Limit(DrawConfig) :
    """A graphical represenation of upper/lower limit as
    arrow from `limit`to `zero`
    >>> l = Limit ( 0.25 , line_color = 3  , arrow_size = 0.15 , ,arrow_style = '->' )  
    """
    def __init__  ( self , limit , zero = 0. , **config ) :
        super(Limit,self).__init__  ( **config )
        
        self.__limit = 1.0 * float ( limit ) 
        self.__zero  = 1.0 * float ( zero  )
        
        if not 'marker_style' in self.config : self.config [ 'marker_style' ] = 1
        
        if 'color' in self.config :
            if not 'marker_color' in self.config : self.config [ 'marker_color' ] = self.config [ 'color' ]
            if not 'line_color'   in self.config : self.config [ 'line_color'   ] = self.config [ 'color' ]
            if not 'fill_color'   in self.config : self.config [ 'fill_color'   ] = self.config [ 'color' ]
            
    @property
    def limit ( self ) :
        """The upper/lower limit"""
        return self.__limit
    @property
    def zero  ( self ) :
        """Zero value"""
        return self.__zero

    ## create an arrow for the upper limit at certain level 
    def arrow ( self , level ) :
        """ Create the arrow at certain level """
        
        size  = self.config.get ( 'arrow_size'  , 0.05   )
        style = self.config.get ( 'arrow_style' , '|-|>' )
        angle = self.config.get ( 'arrow_angle' , 30     )
        
        arrow = ROOT.TArrow ( self.limit , level , self.zero , level , size , style )
        arrow.SetAngle ( angle )
        
        arrow.set_line_attributes ( **self.config )
        arrow.set_fill_attributes ( **self.config )

        point = ROOT.TGraph(1)
        point[0] = self.limit , 1. * level
        name = self.config.get('name','')
        if name : point.SetName ( name ) 
        
        point.set_line_attributes   ( **self.config )
        point.set_fill_attributes   ( **self.config )
        point.set_marker_attributes ( **self.config )

        return arrow , point 

# ==============================================================================
##  @class Record
#   Graphical representation of the data point/measurement with several
#   unceratinties, that can be asymmetric
#   @code
#   p1 = Record ( 15 , 0.3 , 0.4 , (-0.1,0.2) , 0.46 , marker_style = 20 , marker_size = 4 )
#   p2 = Record ( VE(15,1**2) , 0.4 , (-0.1,0.2) , marker_style = 20 , marker_size = 4 )
#   @endcode
class Record(DrawConfig) :
    """Graphical representation of the data point/measurement with several
    unceratinties, that can be asymmetric
    >>> p1 = Record ( 15 , 0.3 , 0.4 , (-0.1,0.2) , 0.46 , marker_style = 20 , marker_size = 4 )
    >>> p2 = Record ( VE(15,1**2) , 0.4 , (-0.1,0.2) , marker_style = 20 , marker_size = 4 )
    """
    def  __init__ ( self  , value , *errors , **config ) :

        ## initialize  the base  
        super(Record,self).__init__ ( **config )

        if not 'marker_style' in self.config : self.config [ 'marker_style' ] = 20
        if not 'marker_size'  in self.config : self.config [ 'marker_size'  ] =  2

        if 'color' in self.config :
            if not 'marker_color' in self.config : self.config [ 'marker_color' ] = self.config [ 'color' ]
            if not 'line_color'   in self.config : self.config [ 'line_color'   ] = self.config [ 'color' ] 
            
        covp  = 0.0
        covn  = 0.0
        self.__errsp = []
        self.__errsn = []
        
        if   isinstance ( value , num_types  ) :
            
            self.__value = 1.0 * value
            
        elif isinstance ( value , VE ) and 0 <= value.cov2() :
            
            self.__value = value.value()
            covp += value.cov2()
            covn += value.cov2()
            self.__errsp.append ( covp **0.5 ) 
            self.__errsn.append ( covn **0.5 )
            
        else :
            raise TypeError( 'Invalid value %s/%s ' % ( value , type ( value ) ) ) 


        for i, e in enumerate ( errors ) :
            
            if   isinstance ( e , num_types ) and 0 <= e    : ep, en =  e   , -e  
            elif 2 == len ( e ) and e[0] <= 0 and e[1] >= 0 : ep , en = e[1] , e[0] 
            elif 2 == len ( e ) and e[1] <= 0 and e[0] >= 0 : ep , en = e[0] , e[1] 
            else :
                raise TypeError( 'Invalid errors[%d]=%s ' % ( i , str ( e ) ) ) 

            covp += ep * ep
            covn += en * en
            
            self.__errsp.append ( covp ** 0.5 ) 
            self.__errsn.append ( covn ** 0.5 ) 

        self.__errsp = tuple ( self.__errsp )
        self.__errsn = tuple ( self.__errsn )
        
    @property
    def value ( self ) :
        """``value'' : the value/measurement/data point"""
        return self.__value
    
    @property
    def positive_errors ( self ) :
        """``positive errors'' : positive uncуrtainties (cumulative sum in quadrature)"""
        return self.__errsp

    @property
    def negative_errors ( self ) :
        """``positive errors'' : positive uncertainties (cumulative sum in quadrature)"""
        return self.__errsn

    ## construct a graph objejct for this data point/measurement  at givel level 
    def point ( self , level ):
        """Contruct a graph objejct for this data point/measurement  at givel level"""

        epos , eneg = self.positive_errors, self.negative_errors

        if not epos :
            
            graph = ROOT.TGraph ( 1  )
            graph [ 0 ] = self.value , 1.0 *  level 
            
        elif epos == eneg :

            graph = ROOT.TGraphErrors ( len ( epos ) )
            for  i , e in enumerate  ( epos ) :
                graph [ i ] = VE ( self.value , e*e ) , VE ( 1.0 * level , 0. )

        else : 

            graph = ROOT.TGraphAsymmErrors ( len ( epos ) )
            for  i , e in enumerate  ( zip ( epos , eneg ) ) :
                ep , en = e 
                graph [ i ] = self.value , en , ep , 1. * level , 0. , 0.  
            
        graph.set_line_attributes   ( **self.config )
        graph.set_fill_attributes   ( **self.config )
        graph.set_marker_attributes ( **self.config )
        
        name = self.config.get('name','')
        if name : graph.SetName ( name ) 

        return graph

# ==============================================================================
##  @class Average 
#   Graphical representation of the data point/measurement with several
#   unceratinties, that can be asymmetric
#   @code
#   p1 = Record ( 15 , 0.3 , 0.4 , (-0.1,0.2) , 0.46 , marker_style = 20 , marker_size = 4 )
#   p2 = Record ( VE(15,1**2) , 0.4 , (-0.1,0.2) , marker_style = 20 , marker_size = 4 )
#   @endcode
class Average(Record) :
    """Graphical representation of the data point/measurement with several
    unceratinties, that can be asymmetric
    >>> p1 = Record ( 15 , 0.3 , 0.4 , (-0.1,0.2) , 0.46 , marker_style = 20 , marker_size = 4 )
    >>> p2 = Record ( VE(15,1**2) , 0.4 , (-0.1,0.2) , marker_style = 20 , marker_size = 4 )
    """
    def  __init__ ( self  , value , *errors  , **config ) :

        assert isinstance ( value , VE ) and 0 < value.cov2() , 'Invalid average %s/%s' % ( value , type(value) )
        
        ## initialize  the base  
        super(Average,self).__init__ ( value , *errors  , **config )
        
        if not 'fill_style' in self.config : self.config [ 'fill_style' ] = 1001 
        if not 'fill_color' in self.config : self.config [ 'fill_color' ] = ROOT.kOrange   
        if not 'line_color' in self.config : self.config [ 'line_color' ] = ROOT.kOrange + len ( self.positive_errors )    
        if not 'average_width' in self.config : self.config [ 'average_width' ] = 3 
        
    ## construct a box object for the average 
    def boxes ( self , level ) :
        """Construct a graph objejct for this ``average'' at givel level
        """
        
        epos , eneg  = reversed ( self.positive_errors ) , \
                       reversed ( self.negative_errors )

        boxes  = []
        for i , e in enumerate ( zip ( epos, eneg ) ) :
            ep , en = e 
            b = ROOT.TBox  ( self.value - en , 0., self.value + ep , level )
            
            b.set_fill_attributes ( **self.config )
            ##  adjust  the color
            b.SetFillColor ( self.config ['fill_color'] + i )
            boxes.append ( b )

        ##  mean value 
        line  = ROOT.TLine( self.value , 0 , self.value , level )
        line.set_line_attributes ( **self.config )
        line.SetLineWidth ( self.config[ 'average_width' ] )
        
        boxes.append (  line  )   
        
        return tuple ( boxes ) 

# =============================================================================================
##  make summary (multi) graph
def graph_summary ( data  , average = None  , transpose = False  ) :  
    
    np      = len ( data)
    grpahs  = []
    points  = ROOT.TMultiGraph()

    limits  = []

    ldata = reversed ( data ) if transpose else data 

    for i , record in enumerate ( ldata ) :

        iv = np - i - 0.5

        if   isinstance ( record  , Record ) :
            
            point  = record.point ( iv )
            
            option = record.config.get ( 'option' , '' )
            if option : points.Add ( point , option )
            else      : points.Add ( point          )

        elif isinstance ( record , Limit ) :

            arrow , point = record. arrow ( iv ) 
            
            option = record.config.get ( 'option' , '' )
            if option : points.Add ( point , option )
            else      : points.Add ( point          )
        
            limits.append ( arrow  )
            
        else :
            raise TypeError ( 'Invalid record #%d %s/%s' % ( i , str ( record ) , type ( record ) ) )

        
    assert not average or isinstance ( average , ( VE , Average ) ) ,\
           'Invalid average %s/%s'% ( average , type(average) )   


    boxes = () 
    if isinstance( average , VE ) : average  = Average ( average  )
    if average  : 
        boxes = average.boxes ( np ) 
    
    if transpose  :
        
        points = points.T ()
        limits = [ l.T()  for l in limits  ]
        boxes  = [ b.T()  for b in boxes   ] 
        
        points.GetXaxis().SetNdivisions(0)

    else :

        points.GetYaxis().SetNdivisions(0)
        points.SetMinimum ( 0  )
        points.SetMaximum ( np )
    
    return points, limits, boxes 


# =============================================================================================
##  make summary (multi) graph
def draw_summary ( data      = []     ,
                   transpose = False  ,
                   average   = None   , 
                   vmin      = None   ,
                   vmax      = None   ) : 

    points , limits , boxes  = graph_summary ( data , average  , transpose )

    if transpose :
        
        if vmin is None :
            vmin = points.ymin() 
            for l in limits  : vmin = min ( vmin , l.GetY1()  , l.GetY2() )

        if vmax is None :
            vmax = points.ymax() 
            for l in limits  : vmax = max ( vmax , l.GetY1()  , l.GetY2() )

            
        histo = ROOT.TH1F ( hID() , '' , 10 , 0 , len ( points ) )

        histo.GetXaxis().SetNdivisions(0)
        histo.SetMinimum ( vmin  )
        histo.SetMaximum ( vmax  )

    else :

        if vmin is None :
            vmin = points.xmin() 
            for l in limits  : vmin = min ( vmin , l.GetX1()  , l.GetX2() )
            
        if vmax is None :
            vmax = points.xmax() 
            for l in limits  : vmax = max ( vmax , l.GetX1()  , l.GetX2() )
            
        histo = ROOT.TH1F ( hID() , '', 10, vmin , vmax )

        histo.GetYaxis().SetNdivisions(0)
        histo.SetMinimum ( 0  )
        histo.SetMaximum ( len ( points )  )
        
    histo.draw()
    for box in boxes : box.draw()  
    points.draw( 'pe1')
    for l in limits : l.draw()


    return histo, points, limits , boxes   
    

#  ============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger import getLogger
    logger = getLogger( 'ostap.plotting.graph_summary' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    conf  = { 'marker_style' : 20 , 'marker_size' : 2 , 'marker_color' : 1 , 'line_color' : 1 }
    aconf = { 'marker_size' : 4 }
    
    data = [ Record (3 ,0.5,0.5 ,0.5, name = 'LHCb' , **conf ) ,
             Record ( VE(2.5,0.3**3)  , (-0.8,0.1) ,0.5, name = 'LHCb' , **conf ) ,
             Limit  ( 1.4 , 1.e-6 , **aconf )  
             ]

    ave = Average ( VE(2.0, 0.2**2 ) , (-0.3, 0.8) )  
    result1 = draw_summary (
        data + [ ave ] , average = ave  )

    result2 = draw_summary (
        data + [ ave ] , average = ave  , transpose = True )

    del result1
    del result2
    
# =============================================================================
##                                                                     The END
# =============================================================================

        
                
            
            
        
        
