#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# =============================================================================
# @file test_plotting_sumamry_graph.py
# Test module for <code>graph_summary</code>
# =============================================================================
""" Test module for `graph_summary`
Prepare tests/examples for the  summary graphs
"""
# =============================================================================
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# =============================================================================
from   ostap.core.pyrouts           import VE 
from   ostap.plotting.graph_summary import  ( Average , Record   , draw_summary ,
                                              Point   , Interval , Limit  ) 
import ROOT, time
# =============================================================================
# logging
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ :
    logger = getLogger ( 'test_plotting_summary_graph' )
else :
    logger = getLogger ( __name__ )
# =============================================================================


# =============================================================================
def test_summary1 () :
    
    conf = { 'label_position' : 3861 , 'markersize' : 1.2 , 'text_size' : 0.05 , } ##  'line_width' : 2 } 
    
    ave = Average ( VE(3871.64,0.060**2) , label = 'New average' , **conf ) 
    
    data = [
        Record ( 3871.8  , 3.1         , 3.0   , label = 'D0'             , **conf ) ,
        Record ( 3873.0  , (+1.8,-1.6) , 1.3   , label = 'BaBar'          , **conf ) ,
        Record ( 3868.7  , 1.5         , 0.4   , label = 'BaBar-K^{0}'    , **conf ) ,
        Record ( 3871.4  , 0.6         , 0.1   , label = 'BaBar-K^{+}'    , **conf ) ,
        Record ( 3871.9  , 0.7         , 0.2   , label = 'BESIII'         , **conf ) ,
        Record ( 3871.95 , 0.48        , 0.12  , label = "LHCb'2010"      , **conf ) ,
        Record ( 3871.85 , 0.27        , 0.19  , label = 'Belle'          , **conf ) ,
        Record ( 3871.61 , 0.16        , 0.19  , label = 'CDF'            , **conf ) ,
        Record ( 3871.69 , 0.17                , label = "PDG\'2018"      , **conf ) ,
        Record ( 3871.70 , 0.067       , 0.068 , label = "LHCb'2020"      , color = ROOT.kRed , **conf ) ,
        Record ( 3871.59 , 0.060       , 0.030 , label = "LHCb'2020"      , color = ROOT.kRed , **conf ) ,
        ave   , 
        Record ( 3871.70 , 0.11                , label = 'm_{D^{0}}#plusm_{D^{*0}}' , **conf  )
        ]
    
    if ROOT.gStyle :
        ROOT.gStyle.SetEndErrorSize (5    )
        ROOT.gStyle.SetTickLength   (0.008)

    result = draw_summary ( data  , average  = ave , vmin = 3860 , vmax = 3877 , offset = 1.0 )  

    if ROOT.gPad :
        ROOT.gPad.RedrawAxis()
        
    time.sleep (3)

# =============================================================================
def test_summary2 ( ) :

    conf = { 'label_position' : 5.6 , 'markersize' : 1.2 , 'text_size' : 0.07 , } ##  'line_width' : 2 }
    
    data = [
        Limit  ( 2.1  , label  = 'Belle'                                     , **conf ) ,
        Record ( VE(3.4 ,1.4**2 ) ,        label = 'BaBar'                   , **conf ) ,
        Record ( VE(2.46,0.64**2) , 0.29 , label = 'LHCb'        , color = 2 , **conf ) ,
        Record ( 2.18 , ( 0.81, -0.77 )  , label = "Global fit"              , **conf ) ,
        Limit  ( 0.59 ,                    label = 'BESIII'                  , **conf ) ,
        ]
    
    if ROOT.gStyle :
        ROOT.gStyle.SetEndErrorSize (5    )
        ROOT.gStyle.SetTickLength   (0.008)
        
    result  = draw_summary ( data , vmin = 0 , vmax = 8 )
        
    if ROOT.gPad :
        ROOT.gPad.RedrawAxis()
        
    time.sleep (3)

# =============================================================================
def test_summary3 ( ) :

    conf = { 'label_position' : 5.6 , 'markersize' : 1.2 , 'text_size' : 0.07 , } ##  'line_width' : 2 }

    
    data = [
        Point    (    0          , **conf ) ,
        Point    (  +19          , **conf ) ,
        Interval (   -1  , +19   , **conf ) ,                
        Limit    (    0  , -40   , arrow_size = 0.01 , arrow_style = '|-|>', **conf ) ,                
        Point    (  +35          , **conf ) ,
        Point    (  +52          , **conf ) ,
        Interval (   -3  ,  -1   , **conf ) ,                
        Point    (  -15          , **conf ) ,
        Point    (  +91          , **conf ) ,
        Point    (  +60          , **conf ) ,
        Point    (  -79          , **conf ) ,
        Point    (  +88          , **conf ) ,
        Point    ( -215          , **conf ) ,
        Interval (  -70  , +124  , **conf ) ,
        Point    ( +102          , **conf ) ,
        Record   ( VE(7,12**2)   , **conf ) ,
        Limit    (    0  , -40   , arrow_size = 0.01 , arrow_style = '|-|>', **conf ) ,
        Point    ( +100          , **conf ) ,        
        Record   ( VE(25,90**2)  , **conf ) ,
        Point    ( -150          , **conf ) ,        
        Point    (    0          , **conf ) ,
        Point    (  -23          , **conf ) ,
        Point    (  +98          , **conf ) ,        
        Point    ( -149          , **conf ) ,
        Record   ( -3 , (+4,-15) , **conf ) ,
        Point    (  +53          , **conf ) ,
        Point    ( +166          , **conf ) ,
        Point    ( +260          , **conf ) ,
        Point    ( -182          , **conf ) ,
        Interval ( -250  , +2    , **conf ) ,
        Point    ( +13           , **conf ) ,
        Point    ( +164          , **conf ) ,
        ]
    
    
    if ROOT.gStyle :
        ROOT.gStyle.SetEndErrorSize (5    )
        ROOT.gStyle.SetTickLength   (0.008)
        
    result  = draw_summary ( data , vmin = -300 , vmax = 300 )
        
    if ROOT.gPad :
        ROOT.gPad.RedrawAxis()
        
    time.sleep (3)

# ==============================================================================
def test_smmary4() :

    conf  = { 'marker_style' : 20 , 'marker_size' : 2 , 'marker_color' : 1 , 'line_color' : 1 , 'label_position' : 4.0 }
    aconf = {}
    aconf.update ( conf )
    aconf.update ( { 'marker_size' : 4 , 'marker_style' : 1 } )
    
    
    data = [ Record (3 ,0.5,0.5 ,0.5, label = 'LHCb' , **conf ) ,
             Record ( VE(2.5,0.3**3)  , (-0.8,0.1) ,0.5, label = 'LHCb-2' , **conf ) ,
             Limit  ( 1.4 , 1.e-6 , label = 'BESIII' , **aconf )  
             ]

    ave = Average ( VE(2.0, 0.2**2 ) , (-0.3, 0.8) , label = 'PDG' , **conf )  

    data.append ( ave  )
    
    result1 = draw_summary ( data , average = ave  , vmax = 6 )

    if ROOT.gPad :
        ROOT.gPad.RedrawAxis()

    time.sleep (3)

    result2 = draw_summary ( data , average = ave  , transpose = True , vmax = 6 )
    
    if ROOT.gPad :
        ROOT.gPad.RedrawAxis()

    time.sleep (3)

      

# =============================================================================
if '__main__' == __name__ :

    test_summary1 ()
    test_summary2 ()
    test_summary3 ()
    test_summary4 ()
    
# =============================================================================
##                                                                      The END  
# =============================================================================

    
