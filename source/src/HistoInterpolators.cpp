// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/HistoInterpolation.h"
#include "Ostap/HistoInterpolators.h"
#include "Ostap/HistoHash.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "TH1.h"
#include "TAxis.h"
// ============================================================================
// Local
// ============================================================================
#include "Integrator1D.h"
#include "Integrator2D.h"
#include "Integrator3D.h"
#include "Exception.h"
#include "local_math.h"
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
namespace 
{
  // ==========================================================================
  /// split histogram if number of bins too large
  const unsigned short s_MAX_BINS = 25   ;
  const double         s_TINY_BIN = 0.01 ;
  // ==========================================================================
  /// has it very tiny bins between xmin and xmax ? 
  inline bool tiny_bin
  ( const TAxis& axis , 
    const double xmin , 
    const double xmax ) 
  {
    const double bw    = s_TINY_BIN * std::abs ( xmax - xmin ) ;
    const double range = axis.GetXmax() - axis.GetXmin() ; 
    const int    nbins = axis.GetNbins () ;
    //
    const TArrayD* a = axis.GetXbins() ;
    if ( !a || 0 == a->GetSize() )
    {
      // here we have uniform bins 
      const double binw = range / nbins ;
      return binw < bw ;                                        // RETURN 
    }
    //
    const int n1   = axis.FindFixBin ( xmin ) ;
    const int n2   = axis.FindFixBin ( xmax ) ;
    if ( n1 == n2 ) { return false ; }                          // RETURN
    //
    const int nmin = std::max ( std::min ( n1 , n2 ) , 1     ) ;
    if ( nbins < nmin ) { return false ; }                      // RETURN
    //
    const int nmax = std::min ( std::max ( n1 , n2 ) , nbins ) ;
    if ( nmax < 1     ) { return false ; }                      // RETURN
    //
    for ( int ibin = nmin ; ibin <= nmax ; ++ibin ) 
    { if ( axis.GetBinWidth ( ibin ) < bw ) { return true ; } } // RETURN 
    //
    return false ;
  }
  // ==========================================================================
}
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
  Ostap::Assert  ( 1 == histo.GetDimension()  , 
                   "Invalid type of ROOT::TH1"  , 
                   "Ostap::Math::Histo1D"       ) ;
  histo.Copy ( m_h ) ;
  m_h.SetDirectory ( nullptr ) ;
  //
  m_tag = Ostap::Utils::hash_histo ( *this ) ;
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
  Ostap::Assert  ( 2 == histo.GetDimension()  ,  
                   "Invalid type of ROOT::TH2"  , 
                   "Ostap::Math::Histo2D"       ) ;
  histo.Copy ( m_h ) ;
  m_h.SetDirectory ( nullptr ) ;
  //
  m_tag = Ostap::Utils::hash_histo ( *this ) ;
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
  //
  m_tag = Ostap::Utils::hash_histo ( *this ) ;
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

// ============================================================================
//  Integration 
// ============================================================================

// ============================================================================
// integral over whole histogram range 
// ============================================================================
double Ostap::Math::Histo1D::integral () const 
{
  const TAxis* xa = m_h.GetXaxis()  ;
  return integral ( xa->GetXmin() , xa->GetXmax() ) ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo1D::integral 
( const double low  , 
  const double high ) const
{
  //
  if      ( s_equal ( low  , high ) ) { return 0 ; }
  else if (           high < low    ) { return - integral ( high , low ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  //
  const double xmin = xa->GetXmin ()  ;
  if ( high <= xmin ) { return 0 ; }
  //
  const double xmax = xa->GetXmax ()  ;  
  if ( low  >= xmax ) { return 0 ; }
  //
  const double x_min = std::max ( low  , xmin ) ;
  const double x_max = std::min ( high , xmax ) ;
  //
  // split it if too large
  const unsigned int bin_min = xa->FindFixBin ( x_min ) ;
  const unsigned int bin_max = xa->FindFixBin ( x_max ) ;
  if ( bin_max > s_MAX_BINS + bin_min || tiny_bin ( *xa , x_min , x_max ) )
  {
    const int    bin_mid = ( bin_max + bin_min ) / 2 ;
    const double x_mid   = xa->GetBinCenter ( bin_mid ) ;
    return integral ( x_min , x_mid ) + integral ( x_mid , x_max ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Histo1D> s_integrator {} ;
  static char s_message[] = "Integral(Histo1D)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  //
  std::tie ( ierror , result , error ) = s_integrator.romberg_integrate
    ( tag () ,  
      &F     ,  
      x_min  , x_max ,                      // low & high edges
      workspace_romberg ( m_workspace ) ,   // workspace
      s_APRECISION_ROMBERG              ,   // absolute precision
      s_RPRECISION_ROMBERG              ,   // relative precision
      s_message                         , 
      __FILE__                          , 
      __LINE__                          ) ;
  //
  return result ;
}
// ============================================================================
// integral over whole histogram range 
// ============================================================================
double Ostap::Math::Histo2D::integral () const 
{
  const TAxis* xa = m_h.GetXaxis()  ;
  const TAxis* ya = m_h.GetYaxis()  ;
  return integral ( xa->GetXmin() , xa->GetXmax() , 
                    ya->GetXmin() , ya->GetXmax() ) ;
}
// ============================================================================
double Ostap::Math::Histo2D::integrateX
( const double y    ) const 
{
  const TAxis* xa = m_h.GetXaxis()  ;
  return integrateX ( y , xa->GetXmin() , xa->GetXmax() ) ;
}
// ============================================================================
double Ostap::Math::Histo2D::integrateY
( const double x    ) const 
{
  const TAxis* ya = m_h.GetYaxis()  ;
  return integrateY ( x , ya->GetXmin() , ya->GetXmax() ) ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo2D::integrateX  
( const double y     , 
  const double xmin  , 
  const double xmax  ) const
{
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( y <= ya->GetXmin () ) { return 0 ; }
  else if ( y >= ya->GetXmax () ) { return 0 ; }
  //
  if      ( s_equal ( xmin , xmax ) ) { return 0 ; }
  else if (           xmax < xmin   ) { return - integrateX ( y , xmax , xmin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( xmax <= xa->GetXmin () ) { return 0 ; }
  else if ( xmin >= xa->GetXmax () ) { return 0 ; }
  //
  const double x_min = std::max ( xmin , xa->GetXmin () ) ;
  const double x_max = std::min ( xmax , xa->GetXmax () ) ;
  //
  // split it if too large
  const unsigned int bin_min = xa->FindFixBin ( x_min ) ;
  const unsigned int bin_max = xa->FindFixBin ( x_max ) ;
  if ( bin_max > s_MAX_BINS + bin_min || tiny_bin ( *xa , x_min , x_max ) ) 
  {
    const int    bin_mid = ( bin_max + bin_min ) / 2 ;
    const double x_mid   = xa->GetBinCenter ( bin_mid ) ;
    return integrateX ( y , x_min , x_mid ) + integrateX ( y , x_mid , x_max ) ;
  }
  //
  // use GSL to evaluate the integral
  typedef Ostap::Math::IntegrateX2<Histo2D> IX ;
  static const Ostap::Math::GSL::Integrator1D<IX> s_integrator ;
  static const char s_message[] = "IntegrateX2(Histo2D)" ;
  //
  const IX fx { this , y } ;
  const auto F = s_integrator.make_function ( &fx ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  //
  std::tie ( ierror , result , error ) = s_integrator.romberg_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'Y' , y )                   , 
      &F     ,  
      x_min  , x_max ,                      // low & high edges
      workspace_romberg ( m_workspace ) ,   // workspace
      s_APRECISION_ROMBERG              ,   // absolute precision
      s_RPRECISION_ROMBERG              ,   // relative precision
      s_message                         , 
      __FILE__                          , 
      __LINE__                          ) ;
  //
  return result ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo2D::integrateY  
( const double x     , 
  const double ymin  , 
  const double ymax  ) const
{
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( x <= xa->GetXmin () ) { return 0 ; }
  else if ( x >= xa->GetXmax () ) { return 0 ; }
  //
  if      ( s_equal ( ymin , ymax ) ) { return 0 ; }
  else if (           ymax < ymin   ) { return - integrateY ( x , ymax , ymin ) ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( ymax <= ya->GetXmin () ) { return 0 ; }
  else if ( ymin >= ya->GetXmax () ) { return 0 ; }
  //
  const double y_min = std::max ( ymin , ya->GetXmin () ) ;
  const double y_max = std::min ( ymax , ya->GetXmax () ) ;
  //
  // split it if interval too large 
  const unsigned int bin_min = ya->FindFixBin ( y_min ) ;
  const unsigned int bin_max = ya->FindFixBin ( y_max ) ;
  if ( bin_max > s_MAX_BINS + bin_min || tiny_bin ( *ya , y_min , y_max ) ) 
  {
    const int    bin_mid = ( bin_max + bin_min ) / 2 ;
    const double y_mid   = ya->GetBinCenter ( bin_mid ) ;
    return integrateY ( x , y_min , y_mid ) + integrateY ( x , y_mid , y_max ) ;
  }
  //
  // use GSL to evaluate the integral
  typedef Ostap::Math::IntegrateY2<Histo2D> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char s_message[] = "IntegrateY2(Histo2D)" ;
  //
  const IY fx { this , x } ;
  const auto F = s_integrator.make_function ( &fx ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  //
  std::tie ( ierror , result , error ) = s_integrator.romberg_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , x )                   , 
      &F     ,  
      y_min  , y_max ,                      // low & high edges
      workspace_romberg ( m_workspace ) ,   // workspace
      s_APRECISION_ROMBERG              ,   // absolute precision
      s_RPRECISION_ROMBERG              ,   // relative precision
      s_message                         , 
      __FILE__                          , 
      __LINE__                          ) ;
  //
  return result ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo2D::integral
( const double xmin  , 
  const double xmax  , 
  const double ymin  , 
  const double ymax  ) const
{
  //
  if      ( s_equal ( xmin , xmax ) ) { return 0 ; }
  else if ( s_equal ( ymin , ymax ) ) { return 0 ; }
  else if ( xmax < xmin ) { return - integral ( xmax , xmin , ymin , ymax ) ; }
  else if ( ymax < ymin ) { return - integral ( xmin , xmax , ymax , ymin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( xmax <= xa->GetXmin () ) { return 0 ; }
  else if ( xmin >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( ymax <= ya->GetXmin () ) { return 0 ; }
  else if ( ymin >= ya->GetXmax () ) { return 0 ; }
  //
  const double x_min = std::max ( xmin , xa->GetXmin () ) ;
  const double x_max = std::min ( xmax , xa->GetXmax () ) ;
  const double y_min = std::max ( ymin , ya->GetXmin () ) ;
  const double y_max = std::min ( ymax , ya->GetXmax () ) ;
  //
  // use cubature   
  static const Ostap::Math::GSL::Integrator2D<Histo2D> s_cubature{} ;
  static const char s_message[] = "Integral(Histo2D)" ;
  const auto F = s_cubature.make_function ( this , x_min , x_max , y_min , y_max ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag ()              , 
      &F                  , 
      20000               ,
      s_APRECISION_CUBE2D , 
      s_RPRECISION_CUBE2D , 
      s_message           , 
      __FILE__            , 
      __LINE__            ) ;
  //
  return result ;
}


// ============================================================================
// integral over whole histogram range 
// ============================================================================
double Ostap::Math::Histo3D::integral () const 
{
  const TAxis* xa = m_h.GetXaxis ()  ;
  const TAxis* ya = m_h.GetYaxis ()  ;
  const TAxis* za = m_h.GetZaxis ()  ;
  return integral ( xa -> GetXmin () , xa -> GetXmax () , 
                    ya -> GetXmin () , ya -> GetXmax () ,
                    za -> GetXmin () , za -> GetXmax () ) ;
}
// ============================================================================
// integral over range 
// ============================================================================
double Ostap::Math::Histo3D::integrateXY 
( const double z     ) const 
{
  const TAxis* xa = m_h.GetXaxis ()  ;
  const TAxis* ya = m_h.GetYaxis ()  ;
  return integrateXY ( z , 
                       xa -> GetXmin () , xa -> GetXmax () , 
                       ya -> GetXmin () , ya -> GetXmax () ) ;
}
// ============================================================================
// integral over range 
// ============================================================================
double Ostap::Math::Histo3D::integrateYZ 
( const double x     ) const 
{
  const TAxis* ya = m_h.GetYaxis ()  ;
  const TAxis* za = m_h.GetZaxis ()  ;
  return integrateYZ ( x , 
                       ya -> GetXmin () , ya -> GetXmax () , 
                       za -> GetXmin () , za -> GetXmax () ) ;
}
// ============================================================================
// integral over range 
// ============================================================================
double Ostap::Math::Histo3D::integrateXZ 
( const double y     ) const 
{
  const TAxis* xa = m_h.GetXaxis ()  ;
  const TAxis* za = m_h.GetZaxis ()  ;
  return integrateXZ ( y , 
                       xa -> GetXmin () , xa -> GetXmax () , 
                       za -> GetXmin () , za -> GetXmax () ) ;
}
// ============================================================================
// integral over range 
// ============================================================================
double Ostap::Math::Histo3D::integrateX
( const double y     , 
  const double z     ) const 
{
  const TAxis* xa = m_h.GetXaxis ()  ;
  return integrateX ( y , z , xa -> GetXmin () , xa -> GetXmax () ) ;
}
// ============================================================================
// integral over range 
// ============================================================================
double Ostap::Math::Histo3D::integrateY
( const double x     , 
  const double z     ) const 
{
  const TAxis* ya = m_h.GetYaxis ()  ;
  return integrateX ( x , z , ya -> GetXmin () , ya -> GetXmax () ) ;
}
// ============================================================================
// integral over range 
// ============================================================================
double Ostap::Math::Histo3D::integrateZ
( const double x     , 
  const double y     ) const 
{
  const TAxis* za = m_h.GetZaxis ()  ;
  return integrateX ( x , y , za -> GetXmin () , za -> GetXmax () ) ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo3D::integral
( const double xmin  , 
  const double xmax  , 
  const double ymin  , 
  const double ymax  , 
  const double zmin  , 
  const double zmax  ) const
{
  //
  if      ( s_equal ( xmin , xmax ) ) { return 0 ; }
  else if ( s_equal ( ymin , ymax ) ) { return 0 ; }
  else if ( s_equal ( zmin , zmax ) ) { return 0 ; }
  else if ( xmax < xmin ) { return - integral ( xmax , xmin , ymin , ymax , zmin , zmax ) ; }
  else if ( ymax < ymin ) { return - integral ( xmin , xmax , ymax , ymin , zmin , zmax ) ; }
  else if ( zmax < zmin ) { return - integral ( xmin , xmax , ymin , ymax , zmax , zmin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( xmax <= xa->GetXmin () ) { return 0 ; }
  else if ( xmin >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( ymax <= ya->GetXmin () ) { return 0 ; }
  else if ( ymin >= ya->GetXmax () ) { return 0 ; }
  //
  const TAxis* za = m_h.GetZaxis ()  ;
  if      ( zmax <= za->GetXmin () ) { return 0 ; }
  else if ( zmin >= za->GetXmax () ) { return 0 ; }

  //
  const double x_min = std::max ( xmin , xa->GetXmin () ) ;
  const double x_max = std::min ( xmax , xa->GetXmax () ) ;
  const double y_min = std::max ( ymin , ya->GetXmin () ) ;
  const double y_max = std::min ( ymax , ya->GetXmax () ) ;
  const double z_min = std::max ( zmin , za->GetXmin () ) ;
  const double z_max = std::min ( zmax , za->GetXmax () ) ;
  //
  // use cubature   
  static const Ostap::Math::GSL::Integrator3D<Histo3D> s_cubature{} ;
  static const char s_message[] = "Integral(Histo3D)" ;
  const auto F = s_cubature.make_function ( this , 
                                            x_min , x_max , 
                                            y_min , y_max ,
                                            z_min , z_max ) ;
  //
  int     ierror =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag ()              , 
      &F                  , 
      50000               ,
      s_APRECISION_CUBE3D , 
      s_RPRECISION_CUBE3D , 
      s_message           , 
      __FILE__            , 
      __LINE__            ) ;
  //
  return result ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo3D::integrateXY 
( const double z    ,                          
  const double xmin , const double xmax ,
  const double ymin , const double ymax  ) const 
{
  //
  if      ( s_equal ( xmin , xmax ) ) { return 0 ; }
  else if ( s_equal ( ymin , ymax ) ) { return 0 ; }
  else if ( xmax < xmin ) { return - integrateXY ( z , xmax , xmin , ymin , ymax ) ; }
  else if ( ymax < ymin ) { return - integrateXY ( z , xmin , xmax , ymax , ymin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( xmax <= xa->GetXmin () ) { return 0 ; }
  else if ( xmin >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( ymax <= ya->GetXmin () ) { return 0 ; }
  else if ( ymin >= ya->GetXmax () ) { return 0 ; }
  //
  const TAxis* za = m_h.GetZaxis ()  ;
  if      ( z    <= za->GetXmin () ) { return 0 ; }
  else if ( z    >= za->GetXmax () ) { return 0 ; }
  //
  const double x_min = std::max ( xmin , xa->GetXmin () ) ;
  const double x_max = std::min ( xmax , xa->GetXmax () ) ;
  const double y_min = std::max ( ymin , ya->GetXmin () ) ;
  const double y_max = std::min ( ymax , ya->GetXmax () ) ;
  //
  // use 2D-cubature   
  //
  typedef Ostap::Math::IntegrateXY<Histo3D> FXY ;
  const FXY fxy ( this , z ) ;
  //
  static const Ostap::Math::GSL::Integrator2D<FXY> s_cubature{} ;
  static const char s_message[] = "IntegralXY(Histo3D)" ;
  const auto F = s_cubature.make_function ( &fxy  , 
                                            x_min , x_max , 
                                            y_min , y_max ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( Ostap::Utils::hash_combiner ( tag () , 'Z' , z ) , 
      &F , 
      50000 , 
      s_APRECISION_CUBE3D  , // absolute precision
      s_RPRECISION_CUBE3D  , // relative precision
      s_message            ,
      __FILE__             ,
      __LINE__             ) ;
  //
  return  result ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo3D::integrateYZ 
( const double x    ,                          
  const double ymin , const double ymax ,
  const double zmin , const double zmax  ) const 
{
  //
  if      ( s_equal ( zmin , zmax ) ) { return 0 ; }
  else if ( s_equal ( ymin , ymax ) ) { return 0 ; }
  else if ( ymax < ymin ) { return - integrateYZ ( x , ymax , ymin , zmin , zmax ) ; }
  else if ( zmax < zmin ) { return - integrateYZ ( x , ymin , ymax , zmax , zmin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( x    <= xa->GetXmin () ) { return 0 ; }
  else if ( x    >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( ymax <= ya->GetXmin () ) { return 0 ; }
  else if ( ymin >= ya->GetXmax () ) { return 0 ; }
  //
  const TAxis* za = m_h.GetZaxis ()  ;
  if      ( zmax <= za->GetXmin () ) { return 0 ; }
  else if ( zmin >= za->GetXmax () ) { return 0 ; }
  //
  const double y_min = std::max ( ymin , ya->GetXmin () ) ;
  const double y_max = std::min ( ymax , ya->GetXmax () ) ;
  const double z_min = std::max ( zmin , za->GetXmin () ) ;
  const double z_max = std::min ( zmax , za->GetXmax () ) ;
  //
  // use 2D-cubature   
  //
  typedef Ostap::Math::IntegrateYZ<Histo3D> FYZ ;
  const FYZ fyz ( this , x ) ;
  //
  static const Ostap::Math::GSL::Integrator2D<FYZ> s_cubature{} ;
  static const char s_message[] = "IntegralYZ(Histo3D)" ;
  const auto F = s_cubature.make_function ( &fyz  , 
                                            y_min , y_max , 
                                            z_min , z_max ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( Ostap::Utils::hash_combiner ( tag () , 'X' , x ) , 
      &F , 
      50000 , 
      s_APRECISION_CUBE3D  , // absolute precision
      s_RPRECISION_CUBE3D  , // relative precision
      s_message            ,
      __FILE__             ,
      __LINE__             ) ;
  //
  return  result ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo3D::integrateX 
( const double y    ,                          
  const double z    ,                          
  const double xmin , const double xmax  ) const 
{
  //
  if      ( s_equal ( xmin , xmax ) ) { return 0 ; }
  else if ( xmax < xmin ) { return - integrateX ( y , z , xmax , xmin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( xmax <= xa->GetXmin () ) { return 0 ; }
  else if ( xmin >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( y    <= ya->GetXmin () ) { return 0 ; }
  else if ( y    >= ya->GetXmax () ) { return 0 ; }
  //
  const TAxis* za = m_h.GetZaxis ()  ;
  if      ( z    <= za->GetXmin () ) { return 0 ; }
  else if ( z    >= za->GetXmax () ) { return 0 ; }
  //
  const double x_min = std::max ( xmin , xa->GetXmin () ) ;
  const double x_max = std::min ( xmax , xa->GetXmax () ) ;
  //
  // split it if too large
  const unsigned int bin_min = xa->FindFixBin ( x_min ) ;
  const unsigned int bin_max = xa->FindFixBin ( x_max ) ;
  if ( bin_max > s_MAX_BINS + bin_min || tiny_bin ( *xa , x_min , x_max ) ) 
  {
    const int    bin_mid = ( bin_max + bin_min ) / 2 ;
    const double x_mid   = xa->GetBinCenter ( bin_mid ) ;
    return integrateX ( y , z , x_min , x_mid ) + integrateX ( y , z , x_mid , x_max ) ;
  }
  //
  // use 1D integration
  //
  typedef Ostap::Math::IntegrateX3<Histo3D> FX ;
  const FX fx ( this , y , z ) ;
  //
  static const Ostap::Math::GSL::Integrator1D<FX> s_integrator ;
  static const char s_message[] = "IntegrateX(Histo3D)" ;
  const auto F = s_integrator.make_function ( &fx ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  //
  std::tie ( ierror , result , error ) = s_integrator.romberg_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'Y' , 'Z' , y , z  ) , 
      &F                                ,   // the function
      x_min    , x_max                  ,   // low & high edges
      workspace_romberg ( m_workspace ) ,   // workspace
      s_APRECISION_ROMBERG              ,   // absolute precision
      s_RPRECISION_ROMBERG              ,   // relative precision
      s_message                         , 
      __FILE__                          , 
      __LINE__                          ) ;
  //
  return  result ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo3D::integrateY 
( const double x    ,                          
  const double z    ,                          
  const double ymin , const double ymax  ) const 
{
  //
  if      ( s_equal ( ymin , ymax ) ) { return 0 ; }
  else if ( ymax < ymin ) { return - integrateY ( x , z , ymax , ymin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( x    <= xa->GetXmin () ) { return 0 ; }
  else if ( x    >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( ymax <= ya->GetXmin () ) { return 0 ; }
  else if ( ymin >= ya->GetXmax () ) { return 0 ; }
  //
  const TAxis* za = m_h.GetZaxis ()  ;
  if      ( z    <= za->GetXmin () ) { return 0 ; }
  else if ( z    >= za->GetXmax () ) { return 0 ; }
  //
  const double y_min = std::max ( ymin , ya->GetXmin () ) ;
  const double y_max = std::min ( ymax , ya->GetXmax () ) ;
  //
  // split it if too large
  const unsigned int bin_min = ya->FindFixBin ( y_min ) ;
  const unsigned int bin_max = ya->FindFixBin ( y_max ) ;
  if ( bin_max > s_MAX_BINS + bin_min || tiny_bin ( *ya , y_min , y_max ) ) 
  {
    const int    bin_mid = ( bin_max + bin_min ) / 2 ;
    const double y_mid   = ya->GetBinCenter ( bin_mid ) ;
    return integrateY ( x , z , y_min , y_mid ) + integrateY ( x , z , y_mid , y_max ) ;
  }
  //
  // use 1D integration
  //
  typedef Ostap::Math::IntegrateY3<Histo3D> FY ;
  const FY fy ( this , x , z ) ;
  //
  static const Ostap::Math::GSL::Integrator1D<FY> s_integrator ;
  static const char s_message[] = "IntegrateY(Histo3D)" ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  // 
  std::tie ( ierror , result , error ) = s_integrator.romberg_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , 'Z' , x , z  ) , 
      &F                                ,   // the function
      y_min    , y_max                  ,   // low & high edges
      workspace_romberg ( m_workspace ) ,   // workspace
      s_APRECISION_ROMBERG              ,   // absolute precision
      s_RPRECISION_ROMBERG              ,   // relative precision
      s_message                         , 
      __FILE__                          , 
      __LINE__                          ) ;
  //
  return  result ;
}
// ============================================================================
// integral between low and high
// ============================================================================
double Ostap::Math::Histo3D::integrateZ 
( const double x    ,                          
  const double y    ,                          
  const double zmin , const double zmax  ) const 
{
  //
  if      ( s_equal ( zmin , zmax ) ) { return 0 ; }
  else if ( zmax < zmin ) { return - integrateZ ( x , y , zmax , zmin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( x    <= xa->GetXmin () ) { return 0 ; }
  else if ( x    >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( y    <= ya->GetXmin () ) { return 0 ; }
  else if ( y    >= ya->GetXmax () ) { return 0 ; }
  //
  const TAxis* za = m_h.GetZaxis ()  ;
  if      ( zmax <= za->GetXmin () ) { return 0 ; }
  else if ( zmin >= za->GetXmax () ) { return 0 ; }
  //
  const double z_min = std::max ( zmin , za->GetXmin () ) ;
  const double z_max = std::min ( zmax , za->GetXmax () ) ;
  //
  // split it if too large
  const unsigned int bin_min = za->FindFixBin ( z_min ) ;
  const unsigned int bin_max = za->FindFixBin ( z_max ) ;
  if ( bin_max > s_MAX_BINS + bin_min || tiny_bin ( *za , z_min , z_max ) ) 
  {
    const int    bin_mid = ( bin_max + bin_min ) / 2 ;
    const double z_mid   = za->GetBinCenter ( bin_mid ) ;
    return integrateZ ( x , y , z_min , z_mid ) + integrateZ ( x , y , z_mid , z_max ) ;
  }
  //
  // use 1D integration
  //
  typedef Ostap::Math::IntegrateZ3<Histo3D> FZ ;
  const FZ fz ( this , x , y ) ;
  //
  static const Ostap::Math::GSL::Integrator1D<FZ> s_integrator ;
  static const char s_message[] = "IntegrateZ(Histo3D)" ;
  const auto F = s_integrator.make_function ( &fz ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  //
  std::tie ( ierror , result , error ) = s_integrator.romberg_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , 'Y' , x , y  ) , 
      &F                                ,   // the function
      z_min    , z_max                  ,   // low & high edges
      workspace_romberg ( m_workspace ) ,   // workspace
      s_APRECISION_ROMBERG              ,   // absolute precision
      s_RPRECISION_ROMBERG              ,   // relative precision
      s_message                         , 
      __FILE__                          , 
      __LINE__                          ) ;
  //
  return  result ;
}
// ============================================================================



// ============================================================================
//                                                                      The END 
// ============================================================================
