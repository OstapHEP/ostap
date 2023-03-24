// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/HistoInterpolation.h"
#include "Ostap/HistoInterpolators.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for 
 *   - class Ostap:::::Math::Histo1D
 *   - class Ostap:::::Math::Histo2D
 *   - class Ostap:::::Math::Histo3D
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-11-16
 */
// ============================================================================
/*  constructor with full  specification 
 *  @see Ostap::Math::HistoInterpolation
 *  @see Ostap::Math::HistoInterpolation::interpolate_1D
 *  @see Ostap::Math::HistoInterpolation::interpolate_2D
 *  @see Ostap::Math::HistoInterpolation::interpolate_3D
 */
// ============================================================================
Ostap::Math::HistoInterpolator::HistoInterpolator
( const bool edges       , 
  const bool extrapolate , 
  const bool density     ) 
  : m_edges       ( edges       ) 
  , m_extrapolate ( extrapolate )
  , m_density     ( density     ) 
{}
// ============================================================================
/*  constructor with full specification 
 *  @see Ostap::Math::HistoInterpolation
 *  @see Ostap::Math::HistoInterpolation::interpolate_1D
 */
// ===========================================================================
Ostap::Math::Histo1D::Histo1D 
( const TH1&                                  histo       ,
  const Ostap::Math::HistoInterpolation::Type t           ,
  const bool                                  edges       , 
  const bool                                  extrapolate ,
  const bool                                  density     ) 
  : Ostap::Math::HistoInterpolator ( edges , extrapolate , density ) 
  , m_h (   ) 
  , m_t ( t ) 
{
  const TH1* h = &histo ;
  Ostap::Assert  ( h && 1 == h->GetDimension()  , 
                   "Invalid type of ROOT::TH1"  , 
                   "Ostap::Math::Histo1D"       ) ;
  histo.Copy ( m_h ) ;
  m_h.SetDirectory ( nullptr ) ;
}
// ===========================================================================
// constructor from the histogram and predefined configuration 
// ===========================================================================
Ostap::Math::Histo1D::Histo1D 
( const TH1&                  histo ,
  const Ostap::Math::Histo1D& conf  ) 
  : Histo1D ( histo               , 
              conf.t           () ,
              conf.edges       () , 
              conf.extrapolate () , 
              conf.density     () )
{}
// ============================================================================
/*  constructor with full specification 
 *  @see Ostap::Math::HistoInterpolation
 *  @see Ostap::Math::HistoInterpolation::interpolate_2D
 */
// ===========================================================================
Ostap::Math::Histo2D::Histo2D 
( const TH2&                                  histo       ,
  const Ostap::Math::HistoInterpolation::Type tx          ,
  const Ostap::Math::HistoInterpolation::Type ty         ,
  const bool                                  edges       , 
  const bool                                  extrapolate ,
  const bool                                  density     ) 
  : Ostap::Math::HistoInterpolator ( edges , extrapolate , density ) 
  , m_h  (    ) 
  , m_tx ( tx ) 
  , m_ty ( ty ) 
{
  const TH2* h = &histo ;
  Ostap::Assert  ( h && 2 == h->GetDimension()  ,  
                   "Invalid type of ROOT::TH2"  , 
                   "Ostap::Math::Histo2D"       ) ;
  histo.Copy ( m_h ) ;
  m_h.SetDirectory ( nullptr ) ;
}
// ============================================================================
/*  constructor with full specification 
 *  @see Ostap::Math::HistoInterpolation
 *  @see Ostap::Math::HistoInterpolation::interpolate_3D
 */
// ===========================================================================
Ostap::Math::Histo3D::Histo3D 
( const TH3&                                  histo       ,
  const Ostap::Math::HistoInterpolation::Type tx          ,
  const Ostap::Math::HistoInterpolation::Type ty          ,
  const Ostap::Math::HistoInterpolation::Type tz          ,
  const bool                                  edges       , 
  const bool                                  extrapolate ,
  const bool                                  density     ) 
  : Ostap::Math::HistoInterpolator ( edges , extrapolate , density ) 
  , m_h  (    ) 
  , m_tx ( tx ) 
  , m_ty ( ty ) 
  , m_tz ( tz ) 
{
  Ostap::Assert  ( 3 == histo.GetDimension()   ,  
                   "Invalid type of ROOT::TH3" , 
                   "Ostap::Math::Histo3D"      ) ;
  histo.Copy ( m_h ) ;
  m_h.SetDirectory ( nullptr ) ;
}
// ============================================================================
Ostap::Math::Histo1D::Histo1D () 
  : Ostap::Math::HistoInterpolator () 
  , m_h ( ) 
  , m_t ( Ostap::Math::HistoInterpolation::Default )
{ m_h.SetDirectory ( nullptr ) ; }
// ============================================================================
Ostap::Math::Histo2D::Histo2D () 
  : Ostap::Math::HistoInterpolator () 
  , m_h  ( ) 
  , m_tx ( Ostap::Math::HistoInterpolation::Default )
  , m_ty ( Ostap::Math::HistoInterpolation::Default )
{ m_h.SetDirectory ( nullptr ) ; }
// ============================================================================
Ostap::Math::Histo3D::Histo3D () 
  : Ostap::Math::HistoInterpolator () 
  , m_h  ( ) 
  , m_tx ( Ostap::Math::HistoInterpolation::Default )
  , m_ty ( Ostap::Math::HistoInterpolation::Default )
  , m_tz ( Ostap::Math::HistoInterpolation::Default )
{ m_h.SetDirectory ( nullptr ) ; }
// ============================================================================
//                                                                      The END 
// ============================================================================
