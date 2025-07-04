// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TH1.h"
#include "TAxis.h"
#include "TRandom.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/HistoInterpolation.h"
#include "Ostap/HistoInterpolators.h"
#include "Ostap/HistoHash.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local
// ============================================================================
#include "Integrator1D.h"
#include "Integrator2D.h"
#include "Integrator3D.h"
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

double Ostap::Math::Histo1D::xmin () const { return m_h.GetXaxis()->GetXmin() ; }
double Ostap::Math::Histo1D::xmax () const { return m_h.GetXaxis()->GetXmax() ; }
double Ostap::Math::Histo2D::xmin () const { return m_h.GetXaxis()->GetXmin() ; }
double Ostap::Math::Histo2D::xmax () const { return m_h.GetXaxis()->GetXmax() ; }
double Ostap::Math::Histo2D::ymin () const { return m_h.GetYaxis()->GetXmin() ; }
double Ostap::Math::Histo2D::ymax () const { return m_h.GetYaxis()->GetXmax() ; }
double Ostap::Math::Histo3D::xmin () const { return m_h.GetXaxis()->GetXmin() ; }
double Ostap::Math::Histo3D::xmax () const { return m_h.GetXaxis()->GetXmax() ; }
double Ostap::Math::Histo3D::ymin () const { return m_h.GetYaxis()->GetXmin() ; }
double Ostap::Math::Histo3D::ymax () const { return m_h.GetYaxis()->GetXmax() ; }
double Ostap::Math::Histo3D::zmin () const { return m_h.GetZaxis()->GetXmin() ; }
double Ostap::Math::Histo3D::zmax () const { return m_h.GetZaxis()->GetXmax() ; }

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
  if ( Ostap::Math::HistoInterpolation::Nearest == m_t && low <= xmin && xmax <= high ) 
  {
    // regular sum 
    double result = 0 ;
    const int nbins = xa->GetNbins() ;
    for ( int ibin = 1 ; ibin <= nbins ; ++ibin ) 
    { result += xa->GetBinWidth ( ibin ) * m_h.GetBinContent ( ibin ) ; }
    return result ;
  }
  //
  const double x_min = std::max ( low  , xmin ) ;
  const double x_max = std::min ( high , xmax ) ;
  //
  const unsigned int bin_min = std::max ( xa -> FindFixBin ( x_min ) , 1              ) ;
  const unsigned int bin_max = std::min ( xa -> FindFixBin ( x_max ) , xa->GetNbins() ) ;
  //
  if ( Ostap::Math::HistoInterpolation::Nearest == m_t ) 
  {
    double result = 0 ;
    if ( bin_min == bin_max  ) { return m_h.GetBinContent ( bin_min ) * ( x_max - x_min ) ; }  
    else 
    {
      // regular sum 
      for ( int ibin = bin_min + 1 ; ibin < bin_max ; ++ibin ) 
      { result += xa->GetBinWidth ( ibin ) * m_h.GetBinContent ( ibin ) ; }
      // first bin
      result   += ( xa -> GetBinUpEdge ( bin_min ) - x_min ) * m_h.GetBinContent ( bin_min ) ; 
      /// last bin 
      result   += ( x_max - xa->GetBinLowEdge ( bin_max )  ) * m_h.GetBinContent ( bin_max ) ; 
    }
    return result ;
  }
  //
  // split it if too large
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
  //
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
  const unsigned int bin_min = std::max ( xa->FindFixBin ( x_min ) , 1              ) ;
  const unsigned int bin_max = std::min ( xa->FindFixBin ( x_max ) , xa->GetNbins() ) ;
  //
  if ( Ostap::Math::HistoInterpolation::Nearest == m_tx ) 
  {
    double result = 0 ;
    if ( bin_min == bin_max  ) { result = (*this) ( xa->GetBinCenter ( bin_min ) , y )  * ( x_max - x_min ) ; }  
    else 
    {
      // regular sum 
      for ( int ibin = bin_min + 1 ; ibin < bin_max ; ++ibin ) 
      { result += (*this) ( xa -> GetBinCenter ( ibin    ) , y ) *   xa -> GetBinWidth  ( ibin )  ; }
      // first bin
      result   += (*this) ( xa -> GetBinCenter ( bin_min ) , y ) * ( xa -> GetBinUpEdge ( bin_min ) - x_min ) ;
      /// last bin 
      result   += (*this) ( xa -> GetBinCenter ( bin_max ) , y ) * ( x_max - xa->GetBinLowEdge ( bin_max )  ) ;
    }
    return result ;
  }
  //
  // split it if too large
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
  const unsigned int bin_min = std::max ( ya->FindFixBin ( y_min ) , 1              ) ;
  const unsigned int bin_max = std::min ( ya->FindFixBin ( y_max ) , ya->GetNbins() ) ;
  //
  if ( Ostap::Math::HistoInterpolation::Nearest == m_ty ) 
  {
    double result = 0 ;
    if ( bin_min == bin_max  ) { result = (*this) ( x , ya->GetBinCenter ( bin_min ) )  * ( y_max - y_min ) ; }  
    else 
    {
      // regular sum 
      for ( int ibin = bin_min + 1 ; ibin < bin_max ; ++ibin ) 
      { result += (*this) ( x , ya -> GetBinCenter ( ibin    ) ) *   ya -> GetBinWidth  ( ibin )  ; }
      // first bin
      result   += (*this) ( x , ya -> GetBinCenter ( bin_min ) ) * ( ya -> GetBinUpEdge ( bin_min ) - y_min ) ;
      /// last bin 
      result   += (*this) ( x , ya -> GetBinCenter ( bin_max ) ) * ( y_max - ya->GetBinLowEdge ( bin_max )  ) ;
    }
    return result ;
  }
  //
  // split it if interval too large 
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
  if ( Ostap::Math::HistoInterpolation::Nearest == m_tx && 
       Ostap::Math::HistoInterpolation::Nearest == m_ty && 
       xmin <= xa->GetXmin() && xa->GetXmax() <= xmax   && 
       ymin <= ya->GetXmin() && ya->GetXmax() <= ymax   ) 
  {
    // regular sum 
    double result = 0 ;
    const int nbx = xa->GetNbins() ;
    const int nby = ya->GetNbins() ;
    for ( int ix = 1 ; ix <= nbx ; ++ix ) 
    {
      const double bwx = xa -> GetBinWidth ( ix ) ;     
      for ( int iy = 1 ; iy <= nby ; ++iy ) 
      {
        const double bwy = ya -> GetBinWidth ( iy ) ;     
        result += m_h.GetBinContent ( ix , iy ) * bwx * bwy ;
      }
    }
    return result ;
  }
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
  //
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
  if ( Ostap::Math::HistoInterpolation::Nearest == m_tx && 
       Ostap::Math::HistoInterpolation::Nearest == m_ty && 
       Ostap::Math::HistoInterpolation::Nearest == m_tz && 
       xmin <= xa->GetXmin() && xa->GetXmax() <= xmax   && 
       ymin <= ya->GetXmin() && ya->GetXmax() <= ymax   &&
       zmin <= za->GetXmin() && za->GetXmax() <= zmax   ) 
  {
    // regular sum 
    double result = 0 ;
    const int nbx = xa -> GetNbins () ;
    const int nby = ya -> GetNbins () ;
    const int nbz = za -> GetNbins () ;
    for ( int ix = 1 ; ix <= nbx ; ++ix ) 
    {
      const double bwx = xa -> GetBinWidth ( ix ) ;     
      for ( int iy = 1 ; iy <= nby ; ++iy ) 
      {
        const double bwy = ya -> GetBinWidth ( iy ) ;     
        for ( int iz = 1 ; iz <= nbz ; ++iz ) 
        {
          const double bwz = za -> GetBinWidth ( iz ) ;     
          result += m_h.GetBinContent ( ix , iy , iz ) * bwx * bwy * bwz ;
        }
      }
    }
    return result ;
  }
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
double Ostap::Math::Histo3D::integrateXZ 
( const double y    ,                          
  const double xmin , const double xmax ,
  const double zmin , const double zmax  ) const 
{
  //
  if      ( s_equal ( zmin , xmax ) ) { return 0 ; }
  else if ( s_equal ( zmin , zmax ) ) { return 0 ; }
  else if ( xmax < xmin ) { return - integrateXZ ( y , xmax , xmin , zmin , zmax ) ; }
  else if ( zmax < zmin ) { return - integrateXZ ( y , xmin , xmax , zmax , zmin ) ; }
  //
  const TAxis* xa = m_h.GetXaxis ()  ;
  if      ( xmax <= xa->GetXmin () ) { return 0 ; }
  else if ( xmin >= xa->GetXmax () ) { return 0 ; }
  //
  const TAxis* za = m_h.GetZaxis ()  ;
  if      ( zmax <= za->GetXmin () ) { return 0 ; }
  else if ( zmin >= za->GetXmax () ) { return 0 ; }
  //
  const TAxis* ya = m_h.GetYaxis ()  ;
  if      ( y    <= ya->GetXmin () ) { return 0 ; }
  else if ( y    >= ya->GetXmax () ) { return 0 ; }
  //
  const double x_min = std::max ( xmin , xa->GetXmin () ) ;
  const double x_max = std::min ( xmax , xa->GetXmax () ) ;
  const double z_min = std::max ( zmin , za->GetXmin () ) ;
  const double z_max = std::min ( zmax , za->GetXmax () ) ;
  //
  // use 2D-cubature   
  //
  typedef Ostap::Math::IntegrateXY<Histo3D> FXZ ;
  const FXZ fxz ( this , y ) ;
  //
  static const Ostap::Math::GSL::Integrator2D<FXZ> s_cubature{} ;
  static const char s_message[] = "IntegralXZ(Histo3D)" ;
  const auto F = s_cubature.make_function ( &fxz  , 
                                            x_min , x_max , 
                                            z_min , z_max ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( Ostap::Utils::hash_combiner ( tag () , 'Y' , y ) , 
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
  const unsigned int bin_min = std::max ( xa->FindFixBin ( x_min ) , 1              ) ;
  const unsigned int bin_max = std::min ( xa->FindFixBin ( x_max ) , xa->GetNbins() ) ;
  //
  if ( Ostap::Math::HistoInterpolation::Nearest == m_tx ) 
  {
    double result = 0 ;
    if ( bin_min == bin_max  ) { result = (*this) ( xa->GetBinCenter ( bin_min ) , y , z )  * ( x_max - x_min ) ; }  
    else 
    {
      // regular sum 
      for ( int ibin = bin_min + 1 ; ibin < bin_max ; ++ibin ) 
      { result += (*this) ( xa -> GetBinCenter ( ibin    ) , y , z ) *   xa -> GetBinWidth  ( ibin )  ; }
      // first bin
      result   += (*this) ( xa -> GetBinCenter ( bin_min ) , y , z ) * ( xa -> GetBinUpEdge ( bin_min ) - x_min ) ;
      /// last bin 
      result   += (*this) ( xa -> GetBinCenter ( bin_max ) , y , z ) * ( x_max - xa->GetBinLowEdge ( bin_max )  ) ;
    }
    return result ;
  }
  //
  // split it if too large
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
  const unsigned int bin_min = std::max ( ya->FindFixBin ( y_min ) , 1              ) ;
  const unsigned int bin_max = std::min ( ya->FindFixBin ( y_max ) , ya->GetNbins() ) ;
  //
  if ( Ostap::Math::HistoInterpolation::Nearest == m_ty ) 
  {
    double result = 0 ;
    if ( bin_min == bin_max  ) { result = (*this) ( x , ya->GetBinCenter ( bin_min ) , z )  * ( y_max - y_min ) ; }  
    else 
    {
      // regular sum 
      for ( int ibin = bin_min + 1 ; ibin < bin_max ; ++ibin ) 
      { result += (*this) ( x , ya -> GetBinCenter ( ibin    ) , z ) *   ya -> GetBinWidth  ( ibin )  ; }
      // first bin
      result   += (*this) ( x , ya -> GetBinCenter ( bin_min ) , z ) * ( ya -> GetBinUpEdge ( bin_min ) - y_min ) ;
      /// last bin 
      result   += (*this) ( x , ya -> GetBinCenter ( bin_max ) , z ) * ( y_max - ya->GetBinLowEdge ( bin_max )  ) ;
    }
    return result ;
  }
  //
  // split it if too large
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
  const unsigned int bin_min = std::max ( za->FindFixBin ( z_min ) , 1              ) ;
  const unsigned int bin_max = std::min ( za->FindFixBin ( z_max ) , za->GetNbins() ) ;
  //
  if ( Ostap::Math::HistoInterpolation::Nearest == m_ty ) 
  {
    double result = 0 ;
    if ( bin_min == bin_max  ) { result = (*this) ( x , y , za->GetBinCenter ( bin_min ) )  * ( z_max - z_min ) ; }  
    else 
    {
      // regular sum 
      for ( int ibin = bin_min + 1 ; ibin < bin_max ; ++ibin ) 
      { result += (*this) ( x , y , za -> GetBinCenter ( ibin    ) ) *   za -> GetBinWidth  ( ibin )  ; }
      // first bin
      result   += (*this) ( x , y , za -> GetBinCenter ( bin_min ) ) * ( za -> GetBinUpEdge ( bin_min ) - z_min ) ;
      /// last bin 
      result   += (*this) ( x , y , za -> GetBinCenter ( bin_max ) ) * ( z_max - za->GetBinLowEdge ( bin_max )  ) ;
    }
    return result ;
  }
  //
  // split it if too large
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
// Random numbers 
// ============================================================================
/*  get a random number from this distribution
 *  - if maximum value is non-positive: generate uniform distribution
 *  - negative content is interpreted as zero
 *  @param rng  random generator
 *  @see TH1::GetRandom 
 *  @attention It can be rather inefficieint (e.g. for historgams 
 *  with large  number of empty vins 
 */
// ============================================================================
double
Ostap::Math::Histo1D::random ( TRandom* rng ) const
{
  double result = 0 ;
  random ( result , rng ) ;
  return result ;
}
// ============================================================================
// get random number
// ============================================================================
std::size_t Ostap::Math::Histo1D::random
( double&  result ,
  TRandom* rng    ) const 
{
  const double vmax   = m_h.GetBinContent ( m_h.GetMaximumBin() ) ;
  const double r      = rng ? rng->Rndm() : gRandom->Rndm();
  //
  const double xmin   = m_h.GetXaxis()->GetXmin()        ;
  const double xdelta = m_h.GetXaxis()->GetXmax() - xmin ;
  //
  double x = xmin + xdelta * r ;
  //
  std::size_t num = 1 ;
  // if empty or everything non-positive: uniform distribution 
  if ( ( vmax <= 0 ) || s_zero ( vmax ) )
    {
      result = x ;
      return num ;                              // RETURN
    }
  //
  double v = (*this) ( x ) ;
  double c = vmax * ( rng ? rng->Rndm() : gRandom->Rndm() ) ;
  num += 1 ; 
  while ( v < c )
    {
      x    = xmin + xdelta * ( rng ? rng->Rndm() : gRandom->Rndm() ) ; 
      v    = (*this) ( x ) ;
      c    = vmax          * ( rng ? rng->Rndm() : gRandom->Rndm() ) ;
      num += 2 ;
    }
  //
  result =  x ;
  return num ;
}

// ============================================================================
/* get a random number from this distribution
 *  - if maximum value is non-positive: generate uniform distribution
 *  - negative content is interpreted as zero
 *  @param rng  random generator
 *  @see TH1::GetRandom 
 *  @attention It can be rather inefficieint (e.g. for historgams 
 *  with large  number of empty vins 
 */
// ============================================================================
std::array<double,2>
Ostap::Math::Histo2D::random ( TRandom* rng ) const
{
  std::array<double,2> result { 0 , 0 } ;
  random ( result [ 0 ] , result [ 1 ] , rng ) ;
  return result ;
}
// ============================================================================
// get random number
// ============================================================================
std::size_t Ostap::Math::Histo2D::random
( double&  xr  ,
  double&  yr  ,
  TRandom* rng ) const 
{ 
  const double vmax   = m_h.GetBinContent ( m_h.GetMaximumBin() ) ;
  const double rx     = rng ? rng->Rndm() : gRandom->Rndm();
  const double ry     = rng ? rng->Rndm() : gRandom->Rndm();
  //
  const double xmin   = m_h.GetXaxis()->GetXmin()        ;
  const double xdelta = m_h.GetXaxis()->GetXmax() - xmin ;
  //
  const double ymin   = m_h.GetYaxis()->GetXmin()        ;
  const double ydelta = m_h.GetYaxis()->GetXmax() - xmin ;
  //
  double x = xmin + xdelta * rx ;
  double y = ymin + ydelta * ry ;
  //
  std::size_t num = 2 ;
  // if empty or everything non-positive: uniform distribution 
  if ( ( vmax <= 0 ) || s_zero ( vmax ) )
    {
      xr = x ;
      yr = y ;
      return num ;
    }
  //
  double v = (*this) ( x , y ) ;
  double c = vmax * ( rng ? rng->Rndm() : gRandom->Rndm() ) ;
  num += 1 ;
  while ( v < c )
    {
      x = xmin + xdelta * ( rng ? rng->Rndm() : gRandom->Rndm() ) ; 
      y = ymin + ydelta * ( rng ? rng->Rndm() : gRandom->Rndm() ) ; 
      v = (*this) ( x , y ) ;
      c = vmax          * ( rng ? rng->Rndm() : gRandom->Rndm() ) ;
      num += 3 ;
    }
  //
  xr = x ;
  yr = y ;
  return num ;
}
// ============================================================================
/* get a random number from this distribution
 *  - if maximum value is non-positive: generate uniform distribution
 *  - negative content is interpreted as zero
 *  @param rng  random generator
 *  @see TH1::GetRandom 
 *  @attention It can be rather inefficieint (e.g. for historgams 
 *  with large  number of empty vins 
 */
// ============================================================================
std::array<double,3>
Ostap::Math::Histo3D::random ( TRandom* rng ) const
{
  std::array<double,3> result { 0 , 0 , 0 } ;
  random ( result [ 0 ] , result [ 1 ] , result [ 2 ] , rng ) ;
  return result ;
}
// ============================================================================
// get random number
// ============================================================================
std::size_t Ostap::Math::Histo3D::random
( double&  xr  ,
  double&  yr  ,
  double&  zr  ,
  TRandom* rng ) const 
{ 
  //
  const double vmax   = m_h.GetBinContent ( m_h.GetMaximumBin() ) ;
  const double rx     = rng ? rng->Rndm() : gRandom->Rndm();
  const double ry     = rng ? rng->Rndm() : gRandom->Rndm();
  const double rz     = rng ? rng->Rndm() : gRandom->Rndm();
  //
  const double xmin   = m_h.GetXaxis()->GetXmin()        ;
  const double xdelta = m_h.GetXaxis()->GetXmax() - xmin ;
  //
  const double ymin   = m_h.GetYaxis()->GetXmin()        ;
  const double ydelta = m_h.GetYaxis()->GetXmax() - ymin ;
  //
  const double zmin   = m_h.GetZaxis()->GetXmin()        ;
  const double zdelta = m_h.GetZaxis()->GetXmax() - zmin ;
  //
  double x  = xmin + xdelta * rx ;
  double y  = ymin + ydelta * ry ;
  double z  = zmin + zdelta * rz ;
  //
  std::size_t num = 3 ;
  // if empty or everything non-positive: uniform distribution 
  if ( ( vmax <= 0 ) || s_zero ( vmax ) )
    {
      xr = x ;
      yr = y ;
      zr = z ;
      return num ;
    }
  //
  double v = (*this) ( x , y , z ) ;
  double c = vmax * ( rng ? rng->Rndm() : gRandom->Rndm() ) ;
  num     += 1 ;
  while ( v < c )
    {
      x    = xmin + xdelta * ( rng ? rng->Rndm() : gRandom->Rndm() ) ; 
      y    = ymin + ydelta * ( rng ? rng->Rndm() : gRandom->Rndm() ) ; 
      z    = zmin + zdelta * ( rng ? rng->Rndm() : gRandom->Rndm() ) ; 
      v    = (*this) ( x , y , z ) ;
      c    = vmax          * ( rng ? rng->Rndm() : gRandom->Rndm() ) ;
      num += 4 ;
    }
  //
  xr = x ;
  yr = y ;
  zr = z ;
  return num ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
