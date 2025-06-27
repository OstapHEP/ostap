#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/frames/frame2histo.py
#  Helper utilities for frame -> histo reduction 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-06-16
# =============================================================================
""" Helper utilities for frame -> histo reduction """ 
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'DF_P2Model' , 
    'DF_P1Model' ,
    'DF_H3Model' ,
    'DF_H2Model' ,
    'DF_H1Model'
) 
# =============================================================================
import ostap.histos.histos
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.frames.frame2histos' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for frame -> histofram reductions' )
# =============================================================================

# =============================================================================
# Frame -> histogram prjections 
# =============================================================================
    
DF_H1Model = ROOT.ROOT.RDF.TH1DModel 
DF_H1Type  = ROOT.TH1D

DF_H2Model = ROOT.ROOT.RDF.TH2DModel 
DF_H2Type  = ROOT.TH2D

DF_H3Model = ROOT.ROOT.RDF.TH3DModel 
DF_H3Type  = ROOT.TH3D

DF_P1Model = ROOT.ROOT.RDF.TProfile1DModel 
DF_P1Type  = ROOT.TProfile

DF_P2Model = ROOT.ROOT.RDF.TProfile2DModel
DF_P2Type  = ROOT.TProfile2D

# =============================================================================
## convert 1D-histogram to "model" for usage with DataFrames 
def _h1_model_ ( histo ) :
    """ Convert 1D-histogram to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_H1Type ) :
        model = DF_H1Type()
        histo.Copy ( model )
    return DF_H1Model ( model ) 
# =============================================================================
## convert 2D-histogram to "model" for usage with DataFrames 
def _h2_model_ ( histo ) :
    """ Convert 2D-histogram to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_H2Type ) :
        model = DF_H2Type()
        histo.Copy ( model )
    return DF_H2Model ( model ) 
# =============================================================================
## convert 3D-histogram to "model" for usage with DataFrames 
def _h3_model_ ( histo ) :
    """ Convert 3D-histogram to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_H3Type ) :
        model = DF_H3Type()
        histo.Copy ( model )
    return DF_H3Model ( model ) 
# =============================================================================
## convert 1D-profile to "model" for usage with DataFrames 
def _p1_model_ ( histo ) :
    """ Convert 1D-profile to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_P1Type ) :
        model = DF_P1Type()
        histo.Copy ( model )
    return DF_P1Model ( model ) 
# =============================================================================
## convert 2D-profile to "model" for usage with DataFrames 
def _p2_model_ ( histo ) :
    """ Convert 2D-profile to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_P2Type ) :
        model = DF_P2Type()
        histo.Copy ( model )
    return P2Model ( model ) 
# ==============================================================================
ROOT.TH1.model        = _h1_model_
ROOT.TH2.model        = _h2_model_
ROOT.TH3.model        = _h3_model_
ROOT.TProfile  .model = _p1_model_
ROOT.TProfile2D.model = _p2_model_

# =============================================================================
_decorated_classes_ = (
    ROOT.TH1        , 
    ROOT.TH2        , 
    ROOT.TH3        , 
    ROOT.TProfile   , 
    ROOT.TProfile2D , 
)

_new_methods_       = (
    ROOT.TH1.model        , 
    ROOT.TH2.model        , 
    ROOT.TH3.model        , 
    ROOT.TProfile  .model , 
    ROOT.TProfile2D.model , 
)

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
