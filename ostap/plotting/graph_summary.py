#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/plotting/graph_summary.py
# Prepare "summary" plot
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# result = draw_summary ( data  , vmax = 5 )
# @endcode
# Also one can add colored bands for "average"
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
# result = draw_summary ( data  , average  = average , vmax = 5 )
# @endcode
# "Average" can be also added into list of data points:
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
# result = draw_summary ( data + [ average ]  , average  = average , vmax = 5 )
# @endcode
# =============================================================================
"""Prepare ``summary'' plot
>>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
...          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
...          Limit  ( 2.5                    , label = 'BESIII'  ) ]
>>> result = draw_summary ( data  , vmax = 5 )

Also one can add colored bands for ``average'':

>>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
...          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
...          Limit  ( 2.5                    , label = 'BESIII'  ) ]
>>> average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
>>> result = draw_summary ( data  , average  = average , vmax = 5 )

``Average'' data  can be also added into list of data points:

>>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
... Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
... Limit  ( 2.5                    , label = 'BESIII'  ) ]
>>> average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
>>> result = draw_summary ( data + [ average ]  , average  = average , vmax = 5 )
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-04-11"
# =============================================================================
__all__     = (
    'Limit'         , ## graphical representation of the upper/lower limit (arrow) 
    'Record'        , ## graphical representation of data poins with several unceratities
    'Point'         , ## graphical representation of a single data point
    'Interval'      , ## graphical representation of an interval
    'Average'       , ## graphical representation of average  (box)
    'Separator'     , ## separator line for data blocks  
    'Summary'       , ## the  final summary graph 
    'make_summary'  , ## prepare summary graph 
    'draw_summary'  , ## draw summary graph
    'error_band'    , ## prepare drawing for the error band(s)
    )
# =============================================================================
from   ostap.core.ostap_types import num_types, string_types  
from   ostap.core.core        import VE , hID 
from   ostap.utils.cidict     import cidict , cidict_fun 
from   ostap.utils.utils      import vrange  
from   ostap.math.base        import pos_infinity, neg_infinity
import ostap.histos.graphs  
import ROOT
# =============================================================================
## Helper function to decode/pack/unpack/transform  errors and value
#  From sequence of values, get value, and sequece of positive
#  and sequence of negative errors
#  @code
#  v , ep , en = value_errors ( 0.1            , 0.2 , (-0.1,0.4) , (-0.2,0.1) , 0.3 ) 
#  v , ep , en = value_errors ( VE(0.1,0.2**2) , 0.2 , (-0.1,0.4) , (-0.2,0.1) , 0.3 )
#  @endcode 
def value_errors ( value , *errors ) :
    """Helper function to decode/pack/unpack/transform  errors and value
    
    >>> v , ep , en = value_errors ( 0.1            , 0.2 , (-0.1,0.4) , (-0.2,0.1) , 0.3 ) 
    >>> v , ep , en = value_errors ( VE(0.1,0.2**2) , 0.2 , (-0.1,0.4) , (-0.2,0.1) , 0.3 )
    """
    
    _covp  = 0.0
    _covn  = 0.0    
    _errsp = []
    _errsn = []
    
    if   isinstance ( value , num_types  ) :
        
        value = 1.0 * value
        
    elif isinstance ( value , VE ) and 0 <= value.cov2() :
        
        _covp += value.cov2()
        _covn += value.cov2()
        
        _errsp.append ( _covp **0.5 ) 
        _errsn.append ( _covn **0.5 )
        
        value = value.value()
        
    else :
        raise TypeError( 'Invalid value %s/%s ' % ( value , type ( value ) ) ) 
    
    for i, e in enumerate ( errors ) :
        
        if   isinstance ( e , num_types ) and 0 <= e    : ep , en =  e   , -e  
        elif 2 == len ( e ) and e[0] <= 0 and e[1] >= 0 : ep , en = e[1] , e[0] 
        elif 2 == len ( e ) and e[1] <= 0 and e[0] >= 0 : ep , en = e[0] , e[1] 
        else :
            raise TypeError( 'Invalid errors[%d]=%s ' % ( i , str ( e ) ) ) 
        
        _covp += ep * ep
        _covn += en * en
        
        _errsp.append ( _covp ** 0.5 ) 
        _errsn.append ( _covn ** 0.5 ) 

    return value , tuple ( _errsp ) , tuple (_errsn )

# =============================================================================
## Helper function to prepare drawing the error band
#  @code
#  objects = error_band2 ( value , positive_erorors , negative_errors , min_value = 0 , max_value = 10 , transpose = False )
#  for o in objects : o.draw()
#  @endcode
def error_band2 ( value , epos , eneg , min_value , max_value , **kwargs ) :
    """Helper function to prepare drawing the error band(s)
    >>> objects = error_band ( 1 , 0.2 , (-0.1,0.4) , (0.1,-0.3) , min_value = 0 , max_value = 10 , transpose = False )
    >>> for o in objects  : o.draw()
    """
    config = cidict ( transform = cidict_fun )
    config.update ( kwargs )

    transpose   = config.get ( 'transpose'   , False )
    transparent = config.get ( 'transparent' , -1    )

    if max_value < min_value :
        min_value , max_value = max_value , min_value
        
    boxes = []
        
    ## fill color 
    fcolor = config.get ( 'fill_color' , ROOT.kOrange )


    delta = max_value - min_value
    mnv   = min_value + 0.001 * delta
    mxv   = max_value - 0.001 * delta
    
    from  itertools import count
    for fc , ep , en in zip ( count ( fcolor , -1 ) , reversed ( epos ) ,  reversed ( eneg ) ) :
        
        if not transpose :
            box1 = ROOT.TBox ( value - en , mnv  , value + ep , mxv  )
            box2 = ROOT.TBox ( value - en , mnv  , value + ep , mxv  )
        else         :
            box1 = ROOT.TBox ( mnv  , value - en , mxv  , value + ep )
            box2 = ROOT.TBox ( mnv  , value - en , mxv  , value + ep )
            
        box1.set_fill_attributes ( **config ) 
        box2.set_fill_attributes ( **config ) 
        
        ## adjust fill color for transparency
        
        if 0 <= transparent <= 1 : box1.SetFillColorAlpha ( fc , transparent ) 
        else                     : box1.SetFillColor      ( fc ) 
        
        box2.SetLineColor ( fc - 1 )
        box2.SetFillStyle ( 0      )
        
        boxes.append ( box1 )
        boxes.append ( box2 )
                
    ##  mean value
    if not transpose : line  = ROOT.TLine ( value     , min_value , value     , max_value )
    else             : line  = ROOT.TLine ( min_value , value     , max_value , value     )
    
    line.set_line_attributes ( **config )
    line.SetLineWidth ( config.get ( 'average_width' , 3 ) )

    boxes.append ( line ) 
    
    return tuple ( boxes ) 


# =============================================================================
## Helper function to prepare drawing the error band
#  @code
#  objects = error_band ( 1 , 0.2 , (-0.1,0.4) , (0.1,-0.3) , min_value = 0 , max_value = 10 , transpose = False )
#  for o in objects : o.draw()
#  objects = error_band ( VE(1,0.5**2) , (0.1,-0.3) , min_value = 0 , max_value = 10 , transpose = False )
#  @endcode
def error_band ( value , *errors , **kwargs ) :
    """Helper function to prepare drawing the error band(s)
    >>> objects = error_band ( 1 , 0.2 , (-0.1,0.4) , (0.1,-0.3) , min_value = 0 , max_value = 10 , transpose = False )
    >>> for o in objects  : o.draw()
    """

    ## decode/pack/unpack/transform errors 
    v , epos , eneg = value_errors ( value , *errors )

    ## make colored bands 
    return error_band2 ( v , epos , eneg , **kwargs )

    
# =============================================================================
##  @class DrawConfig
#   Helper base class to keep the configuration of drawing attributes 
class DrawConfig(object) :
    """Helper base class to keep the configuration of drawingattributes 
    """    
    def __init__ ( self , **config ) :
        self.__config = cidict ( transform = cidict_fun )
        self.__config.update   ( config )
         
    @property
    def config ( self ) :
        """Configuration: line/fill/marker color/style/size/..."""
        return self.__config

# =============================================================================
## @class Label
#  Helper base class to keep/create  the label
class Label(DrawConfig) :
    """ Helper base class to keep/create  the label
    """
    def __init__  ( self , **config ) :
        super(Label,self).__init__  ( **config )

        if 'label' in self.config :
            if not 'text_font'  in self.config : self.config ['text_font' ] = 132
            if not 'text_align' in self.config : self.config ['text_align'] =  12
            
    # =========================================================================
    ## Create the label  <code>TLatex</code>
    #  @code
    #  point = ...
    #  label = point.label ( 5 )
    #  @endcode 
    def label ( self , vpos ) :
        """Create the label  <code>TLatex</code>
        >>> point = ...
        >>> label = point.label ( 5 )
        - see ROOT.TLatex
        """
        
        text     = self.config.get ( 'label' , ''         )
        if not text : return None
        
        xpos = self.config.get ( 'label_position' , 0 )
        
        label = ROOT.TLatex ( xpos , vpos , text )
        
        label.set_text_attributes ( **self.config )
        
        return label 
        

# =============================================================================
## @class  Limit
#  A graphical represenation of upper/lower limit as
#  arrow from <code>limit</endcode> to <code>zero</endcode>
#  @code
#  l = Limit ( 0.25 , line_color = 3  , arrow_size = 0.15 , ,arrow_style = '->' )  
#  @endcode 
class Limit(Label) :
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

    # =========================================================================
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
        name = self.config.get('name', self.config.get ( 'label', '' ))
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
class Record(Label) :
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
            
    ## construct a graph object for this data point/measurement  at givel level 
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
        
        name = self.config.get('name', self.config.get( 'label' , '') )
        if name : graph.SetName ( name ) 

        return graph

# ==============================================================================
## @class Point
#  Graphical representation of a single data point/measurement
#  @code
#  p1 = Point ( 15 , marker_style = 20 , marker_size = 4 )
#  @endcode
class Point(Record) :
    """Graphical representation of a single data point/measurement
    >>> p1 = Point ( 15 , marker_style = 20 , marker_size = 4 )
    """
    def __init__ ( self , value , **kwargs ) :

        ## initialize the base class 
        Record.__init__ ( self, float ( value ) , **kwargs ) 

# ==============================================================================
## @class Interval
#  Graphical representation of an interval point/measurement
#  @code
#  p1 = Interval ( 10 , 20 , line_color = 20 )
#  @endcode
class Interval(Record) :
    """Graphical representation of an interval point/measurement
    >>> p1 = Interval ( 10 , 20 , line_color = 20 )
    """
    def __init__ ( self , low , high , **kwargs ) :
        value = 0.5 *     ( high + low )
        error = 0.5 * abs ( high - low )
         
        Record.__init__ ( self, value , error , **kwargs ) 

        ## no marker 
        self.config [ 'marker_style' ] = 1

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
        
        if not 'fill_style'    in self.config : self.config [ 'fill_style'    ] = 1001 
        if not 'fill_color'    in self.config : self.config [ 'fill_color'    ] = ROOT.kOrange   
        if not 'line_color'    in self.config : self.config [ 'line_color'    ] = ROOT.kOrange + len ( self.positive_errors )    
        if not 'average_width' in self.config : self.config [ 'average_width' ] = 3 
        if not 'transparent'   in self.config : self.config [ 'transparent'   ] = 0.35 
    
    ## construct colored bands for the ``average with errors"
    def bands ( self , level ) :
        """Construct a graph object for this ``average-with-errors'' at givel level
        """
        return error_band2 ( self.value           ,
                             self.positive_errors ,
                             self.negative_errors ,
                             min_value    = 0     ,
                             max_value    = level ,
                             **self.config        )  


# =============================================================================
## @class Separator
#  the simplesy object - line, likely the separator 
class Separator(Label) :
    ""
    def __init__  ( self , low = None , high = None , **config ) :
        
        super(Label,self).__init__  ( **config )
        
        self.__low   = low 
        self.__high  = high 
        
        if not 'marker_style' in self.config : self.config [ 'marker_style' ] = 1
        if not 'line_style'   in self.config : self.config [ 'line_style'   ] = 2
        if not 'line_swidth'  in self.config : self.config [ 'line_width'   ] = 1
        
        if 'color' in self.config :
            if not 'line_color'   in self.config : self.config [ 'line_color'   ] = self.config [ 'color' ]

    # =========================================================================
    ## create an arrow for the upper limit at certain level 
    def line ( self , vpos , low = None , high = None ) :
        """ Create the arrow at certain level """
        
        style = self.config.get ( 'line_style' , 1 )
        width = self.config.get ( 'line_width' , 1 )


        l = low  if isinstance ( low  , num_types ) else self.__low
        h = high if isinstance ( high , num_types ) else self.__high 

        assert isinstance ( l , num_types ) , 'Invalid low-edge!'
        assert isinstance ( h , num_types ) , 'Invalid high-edge!'
        
        the_line = ROOT.TLine ( l , vpos , h , vpos )  
        
        the_line.set_line_attributes ( **self.config )

        return the_line 

    @property
    def low  ( self ) :
        "``low'': low edge"
        return self.__low
    @low.setter
    def low  ( self, value ) :
        assert value is None or isinatance ( value  , num_types ) , \
               "Invalid 'low' type!"
        self.__low = value

    @property
    def high ( self ) :
        "``high'': high  edge"
        return self.__high
    @high.setter
    def high ( self, value ) :
        assert value is None or isinatance ( value  , num_types ) , \
               "Invalid 'high' type!"
        self.__high = value

        
# ==============================================================================================
## @class Summary
#   The final graph for the ``summary plot''. It is a container of
#  - the histogram  *useful to define the range)
#  - the graph points
#  - the limit arrows
#  - the colored bands for ``average''
#  - the text(TLatex) labels
class Summary(object) :
    """The final graph for the ``summary plot''
    It is a container of
    - the histogram  *useful to define the range)
    - the graph points
    - the limit arrows
    - the colored bands for ``average''
    - the text(TLatex) labels
    """
    def __init__ ( self          ,
                   histo  = None ,
                   points = ()   ,
                   limits = ()   ,
                   bands  = ()   ,
                   labels = ()   ,
                   lines  = [] ) :

        self.__histo  = histo
        self.__points = tuple ( p for p in points ) 
        self.__limits = tuple ( l for l in limits )
        self.__bands  = tuple ( b for b in bands  )
        self.__labels = tuple ( l for l in labels )
        self.__lines  = tuple ( l for l in lines  ) 
        
    # ======================================================================================
    ##  draw the summary plot
    #   @code
    #   summary = ...
    #   summary .draw()
    #   @endcode
    def draw ( self  ) :
        """Draw the summary plot
        >>> summary = ...
        >>> summary .draw()
        """
        
        if   self.histo    :
            
            self.histo.draw ()            
            for b in self.bands  : b.draw()
            for p in self.points : p.draw('pe1')
            
        else :

            if self.points :
                self.points[0].draw('ape1')                
                for p in self.points[1:] :
                    p.draw('pe1')
                
            for b in self.bands  : b.draw()
        
        for l in self.lines      : l.draw()
        for l in self.limits     : l.draw()
        for l in self.labels     : l.draw()

    @property
    def histo ( self  ) :
        """``histo'' : histogram associated with graph"""
        return self.__histo
    @histo.setter
    def histo ( self , h ) :
        assert h is None or ( h and isinstance ( h , ROOT.TH1 ) ) ,\
               'Invalid histo type %s/%s' % ( h , type ( h ) )
        self.__histo = h
        
    @property
    def points ( self  ) :
        """``points'' : graph with points/measurements"""
        return self.__points
    
    @property
    def limits ( self  ) :
        """``limits'' : upper/lower limits"""
        return self.__limits

    @property
    def lines ( self  ) :
        """``lines'' : lines/separators"""
        return self.__lines

    @property
    def bands ( self  ) :
        """``bands'' : color bands associated with mean/averages"""
        return self.__bands

    @property
    def labels ( self  ) :
        """``labels'' : (la)text labels"""
        return self.__labels

    # =========================================================================
    ## number of data points in the summary plot
    #  @code
    #  summary = ...
    #  n = len ( summary )
    #  @endcode
    def __len__  ( self ) :
        """Number of data points in the summary plot
        >>> summary = ...
        >>> n = len ( summary )
        """
        return len ( self.points ) + len ( self.lines ) 
    # =========================================================================
    ## get xmin/xmax for this summary graph
    #  @code
    #  summary = ...
    #  xmin, xmax = summary.xminmax()
    #  @endcode  
    def xminmax ( self ) :
        """Get xmin/xmax for this summary graph 
        >>> summary = ...
        >>> xmin, xmax = summary.xminmax()
        """
        
        xmin = pos_infinity
        xmax = neg_infinity
        
        for p in self.points  :
            xmin = min ( xmin , p.xmin  ()  , p.xmax  () )
            xmax = max ( xmax , p.xmin  ()  , p.xmax  () )
            
        for l in self.limits  :
            xmin = min ( xmin , l.GetX1 ()  , l.GetX2 () )
            xmax = max ( xmax , l.GetX1 ()  , l.GetX2 () )
            
        for b in self.bands   :
            xmin = min ( xmin , b.GetX1 ()  , b.GetX2 () )
            xmax = max ( xmax , b.GetX1 ()  , b.GetX2 () )

        for b in self.lines   :
            xmin = min ( xmin , b.GetX1 ()  , b.GetX2 () )
            xmax = max ( xmax , b.GetX1 ()  , b.GetX2 () )
            
        return xmin, xmax

    # =========================================================================
    ## get ymin/ymax for this summary graph
    #  @code
    #  summary  = ...
    #  ymin, ymax = summary.yminmax()
    #  @endcode  
    def yminmax ( self ) :
        """Get ymin/ymax for this summary graph 
        >>> summary = ...
        >>> ymin, ymax = summary.yminmax()
        """
        ymin = pos_infinity
        ymax = neg_infinity
        
        for p in self.points  :
            ymin = min ( ymin , p.ymin  ()  , p.ymax  () )
            ymax = max ( ymax , p.ymin  ()  , p.ymax  () )
            
        for l in self.limits  :
            ymin = min ( ymin , l.GetY1 ()  , l.GetY2 () )
            ymax = max ( ymax , l.GetY1 ()  , l.GetY2 () )
            
        for b in self.bands   :
            ymin = min ( ymin , b.GetY1 ()  , b.GetY2 () )
            ymax = max ( ymax , b.GetY1 ()  , b.GetY2 () )

        for b in self.lines   :
            ymin = min ( ymin , b.GetY1 ()  , b.GetY2 () )
            ymax = max ( ymax , b.GetY1 ()  , b.GetY2 () )
            
        return ymin, ymax 

# =============================================================================================
## make summary plot  
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# result = make_summary ( data  , vmax = 5 )
# @endcode
# Also one can add colored bands for "average"
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
# result = make_summary ( data  , average  = average , vmax = 5 )
# @endcode
# "Average" can be alsoadded into list of data points:
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
# result = make_summary ( data + [ average ]  , average  = average , vmax = 5 )
# @endcode
def make_summary ( data               ,
                   average   = None   ,
                   transpose = False  ,
                   offset    = 0.5    ,
                   vmin      = None   ,
                   vmax      = None   ) :  
    """Prepare ``summary'' plot
    >>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
    ...          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
    ...          Limit  ( 2.5                    , label = 'BESIII'  ) ]
    >>> result = make_summary ( data )
    
    Also one can add colored bands for ``average'':
    
    >>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
    ...          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
    ...          Limit  ( 2.5                    , label = 'BESIII'  ) ]
    >>> average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
    >>> result = make_summary ( data  , average  = average )
    
    ``Average'' data  can be also added into list of data points:
    
    >>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
    ... Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
    ... Limit  ( 2.5                    , label = 'BESIII'  ) ]
    >>> average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
    >>> result = make_summary ( data + [ average ]  , average  = average  )
    """
    
    np      = len ( data )
    graphs  = []
    limits  = []
    labels  = []
    points  = []
    lines   = [] 


    data = reversed ( data ) if not transpose else data 


    vv = offset + 0.0 
    for i , record in enumerate ( data ) :

        vv = i + offset 
        
        if   isinstance ( record  , Record ) :
            
            point  = record.point ( vv )
            
            points.append ( point )
            
            label = record.label (  vv )
            if label : labels.append ( label )

        elif isinstance ( record , Limit ) :

            arrow , point = record. arrow ( vv ) 
            
            points.append ( point )
        
            limits.append ( arrow  )
            
            label = record.label (  vv )
            if label : labels.append ( label )
            
        elif isinstance ( record , Separator ) :

            line  = record.line ( vv , vmin , vmax )
            lines.append ( line )
            
        elif record and isinstance ( record , string_rypes ) and record.lower() in ( 'line' , 'sep' , 'separator' ) :

            line = Separator()
            line = line.line  ( vv , vmin , vmax )
            
            lines.append ( line )
            
        else :
            
            raise TypeError ( 'Invalid record #%d %s/%s' % ( i , str ( record ) , type ( record ) ) )

                
    assert not average or isinstance ( average , ( VE , Average ) ) ,\
           'Invalid average %s/%s'% ( average , type(average) )   

    bands = () 
    if isinstance( average , VE ) : average  = Average ( average  )
    if average  : 
        bands = average.bands ( np - 1.0 + 2 * offset ) 

    labels = tuple ( labels )  
    
    g = Summary ( points = points ,
                  limits = limits ,
                  bands  = bands  ,
                  labels = labels ,
                  lines  = lines  )  
    
    xmin, xmax  = g.xminmax()
    ymin, ymax  = g.yminmax()

    if transpose  :
        
        limits = tuple ( [ l.T()  for l in g.limits  ] )
        bands  = tuple ( [ b.T()  for b in g.bands   ] )
        points = tuple ( [ p.T()  for p in g.points  ] ) 
        labels = tuple ( [ l.T()  for l in g.labels  ] ) 

        g = Summary( points = points ,
                     limits = limits ,
                     bands  = bands  ,
                     labels = labels )  
        
        for p in g.points : 
            p.GetXaxis().SetNdivisions(0)
            p.SetMinimum ( xmin )
            p.SetMaximum ( xmax )
            
    else :
                
        
        for p in points : 
            p.GetYaxis().SetNdivisions(0)
            p.SetMinimum ( 0  )
            p.SetMaximum ( np )        
            

    return g

# =============================================================================================
##  Prepare and draw "summary" plot
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# result = draw_summary ( data  , vmax = 5 )
# @endcode
# Also one can add colored bands for "average"
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
# result = draw_summary ( data  , average  = average , vmax = 5 )
# @endcode
# "Average" can be alsoadded into list of data points:
# @code
# data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
#          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
#          Limit  ( 2.5                    , label = 'BESIII'  ) ]
# average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
# result = draw_summary ( data + [ average ]  , average  = average , vmax = 5 )
# @endcode
def draw_summary ( data      = []     ,
                   transpose = False  ,
                   average   = None   , 
                   vmin      = None   ,
                   vmax      = None   ,
                   offset    = 0.5    ) : 
    """Prepare and draw the ``summary'' plot
    >>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
    ...          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
    ...          Limit  ( 2.5                    , label = 'BESIII'  ) ]
    >>> result = draw_summary ( data  , vmax = 5 )
    
    Also one can add colored bands for ``average'':
    
    >>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
    ...          Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
    ...          Limit  ( 2.5                    , label = 'BESIII'  ) ]
    >>> average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
    >>> result = draw_summary ( data  , average  = average , vmax = 5 )
    
    ``Averag'' data  can be also added into list of data points:
    
    >>> data = [ Record ( 1.0 , 0.1 ,(-0.2, 0.5 ), label = 'LHCb'  , color = 4 ) ,
    ... Record ( 2.0 , 0.5 ,0.5         , label = 'Belle' , color = 3 , marker_style = 23 ) ,
    ... Limit  ( 2.5                    , label = 'BESIII'  ) ]
    >>> average = Average ( 2.2 , 0.3 , Label = 'PDG' ) 
    >>> result = draw_summary ( data + [ average ]  , average  = average , vmax = 5 )
    """
    
    summary = make_summary ( data                  ,
                             average   = average   ,
                             transpose = transpose ,
                             offset    = offset    ,
                             vmin      = vmin      ,
                             vmax      = vmax      ) 

    xmin, xmax = summary.xminmax()
    ymin, ymax = summary.xminmax()
    
    if transpose :
        
        if vmin is None : vmin = ymin
        if vmax is None : vmax = ymax

        histo = ROOT.TH1F ( hID() , '' , 10 , 0 , len ( summary ) - 1 + 2 * offset )

        histo.GetXaxis().SetNdivisions(0)
        histo.SetMinimum ( vmin  )
        histo.SetMaximum ( vmax  )

    else :
        
        if vmin is None : vmin = xmin
        if vmax is None : vmax = xmax
            
        histo = ROOT.TH1F ( hID() , '', 10, vmin , vmax )

        histo.GetYaxis().SetNdivisions(0)
        histo.SetMinimum ( 0  )
        histo.SetMaximum ( len ( summary ) - 1 + 2 * offset  )

    summary.histo = histo
    summary.draw()

    return summary
    

#  ============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger import getLogger
    logger = getLogger( 'ostap.plotting.graph_summary' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    
# =============================================================================
##                                                                     The END
# =============================================================================

        
                
            
            
        
        
