#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file pyrouts.py
#  Module with decoration of many ROOT objects for efficient use in python
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
#
# =============================================================================
"""Decoration of some ROOT objects for efficient use in python

Many native  root classes are equipped with new useful methods and operators,
in particular TH1(x) , TH2(x) , TAxis, TGraph(Errors), etc...
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
# =============================================================================
__all__     = (
    #
    'cpp'             , ## global C++ namespace
    'Ostap'           , ## C++ namespace Ostap 
    'rootID'          , ## construct the (global) unique ROOT identifier
    'funcID'          , ## construct the (global) unique ROOT identifier
    'funID'           , ## construct the (global) unique ROOT identifier
    'hID'             , ## construct the (global) unique ROOT identifier
    'histoID'         , ## construct the (global) unique ROOT identifier
    'dsID'            , ## construct the (global) unique ROOT identifier
    #
    'VE'              , ## Gaudi::Math::ValueWithError
    'histoGuess'      , ## guess the simple histo parameters
    #'useLL'           , ## use LL for histogram fit?
    #'allInts'         , ## natural histogram with natural entries?
    'SE'              , ## StatEntity
    'WSE'             , ## StatEntity with weights
    'iszero'          , ## is almost zero  ?
    'isequal'         , ## is almost equal ?
    #
    'binomEff'        , ## calculate binomial efficiency
    'binomEff2'       , ## calculate binomial efficiency
    'zechEff'         , ## calculete binomial efficiency using Zech's          prescription
    'wilsonEff'       , ## calculete binomial efficiency using Wilson's        prescription
    'agrestiCoullEff' , ## calculete binomial efficiency using Agresti-Coull's prescription
    #
    'binomEff_h1'     , ## calculate binomial efficiency for 1D-histos
    'binomEff_h2'     , ## calculate binomial efficiency for 2D-ihstos
    'binomEff_h3'     , ## calculate binomial efficiency for 3D-ihstos
    #
    'makeGraph'       , ## make ROOT Graph from input data
    'hToGraph'        , ## convert historgam to graph
    'hToGraph2'       , ## convert historgam to graph
    'hToGraph3'       , ## convert historgam to graph
    'lw_graph'        , ## make Lafferty-Wyatt's graph from the histo 
    'h1_axis'         , ## book 1D-histogram from axis 
    'h2_axes'         , ## book 2D-histogram from axes
    'h3_axes'         , ## book 3D-histogram from axes
    'axis_bins'       , ## convert list of bin edges to axis
    've_adjust'       , ## adjust the efficiency to be in physical range
    #
    )
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
logger = getLogger( __name__ )
if '__main__' == __name__ : logger = getLogger ( 'ostap.core.pyrouts' )
else                      : logger = getLogger ( __name__             )
# =============================================================================
# fixes
# =============================================================================
import ostap.fixes.fixes
from ostap.logger.utils import ROOTIgnore

# =============================================================================
logger.info ( 'Zillions of decorations for ROOT/RooFit objects')
# =============================================================================
with ROOTIgnore ( 2001 ) : 
    from ostap.core.core import ( cpp      , Ostap     , 
                                  ROOTCWD  , rootID    , 
                                  funcID   , funID     , fID             ,
                                  histoID  , hID       , dsID            ,
                                  cwd      , pwd       ,
                                  VE       , SE        , WSE             ,
                                  binomEff , binomEff2 ,
                                  zechEff  , wilsonEff , agrestiCoullEff ,
                                  iszero   , isequal   ,
                                  isint    , islong    , natural_entry   ) 
    
    
## ## silently load RooFit library trick...
## from ostap.logger.utils import mute
## with mute() : _tmp = ROOT.RooRealVar
## del mute


# =============================================================================
## decorate histograms 
# =============================================================================    
with ROOTIgnore( 2001 ) : 
    from ostap.histos.histos import ( binomEff_h1 , binomEff_h2 , binomEff_h3 ,
                                      h1_axis     , h2_axes     , h3_axes     ,
                                      axis_bins   , ve_adjust   , histoGuess  )
    
    
# =============================================================================
# Other decorations 
# =============================================================================
with ROOTIgnore( 2001 ) : 
    
    import ostap.trees.trees
    import ostap.trees.cuts
    
    import ostap.histos.param
    import ostap.histos.compare
    
    import ostap.io.root_file
    
    import ostap.math.polynomials
    import ostap.math.models
    import ostap.utils.hepdata 
    import ostap.utils.pdg_format 
    
    import ostap.plotting.canvas
    
    import ostap.fitting.minuit 
    import ostap.fitting.roofit 

    import ostap.plotting.draw_attributes 
    
# =============================================================================
## graphs 
# =============================================================================
with ROOTIgnore( 2001 ) :     
    from ostap.histos.graphs import ( makeGraph , lw_graph  , 
                                      hToGraph  , hToGraph2 , hToGraph3 )

    

del ROOTIgnore 
# =============================================================================
if '__main__' == __name__ :
            
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
