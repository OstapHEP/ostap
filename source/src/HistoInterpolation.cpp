// ============================================================================
// Include files 
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <array>
// ============================================================================
// ROOT 
// ============================================================================
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/HistoInterpolation.h"
// ============================================================================
/** @file 
 *  Implementation file for class : Gaudi::Math::HistoInterpolation
 *  @see  Gaudi::Math::HistoInterpolation
 *  Originally developed for Ostap in python
 *  Translated to C++ to get some gain in CPU 
 *
 *  @author Vanya Belyaev
 *  @date   2015-10-12
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double> s_equal{} ; // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>     s_zero{}  ; // zero for doubles
  // ==========================================================================
  // get bin content for 1D-histogram 
  inline Ostap::Math::ValueWithError _bin_ 
  ( const TH1&         h1              , 
    const unsigned int i               , 
    const bool         density = false ) 
  {    
    double v = h1.GetBinContent   ( i ) ;
    double e = h1.GetBinError     ( i ) ;
    //
    if ( density ) 
    {
      const double ibw = 1/h1.GetXaxis()->GetBinWidth( i ) ;
      v *= ibw ;
      e *= ibw ;
    }
    //
    return  Ostap::Math::ValueWithError ( v , e * e ) ;
  }
  // get bin content for 2D-histogram 
  inline Ostap::Math::ValueWithError _bin_ 
  ( const TH2&         h2              , 
    const unsigned int ix              , 
    const unsigned int iy              , 
    const bool         density = false ) 
  {    
    double v = h2.GetBinContent   ( ix , iy    ) ;
    double e = h2.GetBinError     ( ix , iy    ) ;
    //
    if ( density ) 
    {
      const double ibw = 1/( h2.GetXaxis()->GetBinWidth ( ix ) * 
                             h2.GetYaxis()->GetBinWidth ( iy ) ) ;
      v *= ibw ;
      e *= ibw ;
    }
    //
    return  Ostap::Math::ValueWithError ( v  , e * e ) ;
  }
  // get bin content for 3D-histogram 
  inline Ostap::Math::ValueWithError _bin_ 
  ( const TH3&         h3              , 
    const unsigned int ix              , 
    const unsigned int iy              ,
    const unsigned int iz              , 
    const bool         density = false ) 
  {    
    double v = h3.GetBinContent   ( ix , iy , iz ) ;
    double e = h3.GetBinError     ( ix , iy , iz ) ;
    //
    if ( density ) 
    {
      const double ibw = 1/( h3.GetXaxis()->GetBinWidth ( ix ) * 
                             h3.GetYaxis()->GetBinWidth ( iy ) * 
                             h3.GetZaxis()->GetBinWidth ( iz ) ) ;
      v *= ibw ;
      e *= ibw ;
    }
    //
    return  Ostap::Math::ValueWithError ( v  , e * e ) ;
  }
  // ==========================================================================
  // linear interpolation 
  inline Ostap::Math::ValueWithError _linear_ 
  ( const double                       x  , 
    const double                       x1 , 
    const double                       x2 , 
    const Ostap::Math::ValueWithError& v1 , 
    const Ostap::Math::ValueWithError& v2 )
  {
    const double dx = 1 / ( x1 - x2 )      ;
    const double c1 =     ( x1 - x  ) * dx ;
    const double c2 =     ( x  - x2 ) * dx ;
    //
    const double vv = v2.value() * c1      + v1.value() * c2      ;
    const double e2 = 
      ( 0 < v1.cov2() ? v1.cov2 () * c2 * c2 : 0.0 ) + 
      ( 0 < v2.cov2() ? v2.cov2 () * c1 * c1 : 0.0 ) ;
    //
    return Ostap::Math::ValueWithError ( vv , ( 0 >= e2 || s_zero ( e2 ) ) ? 0.0 : e2 ) ;
  }
  // ==========================================================================
  // quadratic interpolation 
  inline 
  Ostap::Math::ValueWithError _quadratic_ 
  ( const double                       x  , 
    const double                       x0 , 
    const double                       x1 , 
    const double                       x2 , 
    const Ostap::Math::ValueWithError& v0 , 
    const Ostap::Math::ValueWithError& v1 , 
    const Ostap::Math::ValueWithError& v2 )
  {
    const double dx0  = x - x0 ;
    const double dx1  = x - x1 ;
    const double dx2  = x - x2 ;
    //
    const double dx01 = x0 - x1 ;
    const double dx02 = x0 - x2 ;
    const double dx12 = x1 - x2 ;
    //
    const double c0   = dx1 * dx2 / ( dx01 * dx02 ) ;
    const double c1   = dx0 * dx2 / ( dx01 * dx12 ) ; 
    const double c2   = dx0 * dx1 / ( dx02 * dx12 ) ; 
    //
    const double vv = v0.value() * c0 + v2.value() * c2 - v1.value() * c1 ;
    //
    const double e2 = 
      ( 0 < v0.cov2() ? v0.cov2 () * c0 * c0 : 0.0 ) + 
      ( 0 < v1.cov2() ? v1.cov2 () * c1 * c1 : 0.0 ) +
      ( 0 < v2.cov2() ? v2.cov2 () * c2 * c2 : 0.0 ) ;
    //
    return Ostap::Math::ValueWithError ( vv , s_zero ( e2 ) ? 0.0 : e2 ) ;
  }
  inline 
  Ostap::Math::ValueWithError _quadratic2_ 
  ( const double                       x  , 
    const double                       x0 , 
    const double                       x1 , 
    const double                       x2 , 
    const double                       x3 , 
    const Ostap::Math::ValueWithError& v0 , 
    const Ostap::Math::ValueWithError& v1 , 
    const Ostap::Math::ValueWithError& v2 ,
    const Ostap::Math::ValueWithError& v3 )
  {
    return 
      x < x1 ? _quadratic_ ( x , x0 , x1 , x2 , v0 , v1 , v2 ) :
      x > x2 ? _quadratic_ ( x , x1 , x2 , x3 , v1 , v2 , v3 ) :
    0.5 *    ( _quadratic_ ( x , x0 , x1 , x2 , v0 , v1 , v2 ) + 
               _quadratic_ ( x , x1 , x2 , x3 , v1 , v2 , v3 ) ) ;  
  }
  // ==========================================================================
  // qubic interpolation 
  inline 
  Ostap::Math::ValueWithError _cubic_ 
  ( const double                       x  , 
    const double                       x0 , 
    const double                       x1 , 
    const double                       x2 , 
    const double                       x3 , 
    const Ostap::Math::ValueWithError& v0 , 
    const Ostap::Math::ValueWithError& v1 , 
    const Ostap::Math::ValueWithError& v2 ,
    const Ostap::Math::ValueWithError& v3 )
  {
    const double dx0  = x  - x0 ;
    const double dx1  = x  - x1 ;
    const double dx2  = x  - x2 ;
    const double dx3  = x  - x3 ;
    //
    const double dx01 = x0 - x1 ;
    const double dx02 = x0 - x2 ;
    const double dx03 = x0 - x3 ;
    const double dx12 = x1 - x2 ;
    const double dx13 = x1 - x3 ;
    const double dx23 = x2 - x3 ;
    //
    const double c0 =         dx1 * dx2 * dx3 / ( dx01 * dx02 * dx03 ) ;
    const double c1 = - dx0 *       dx2 * dx3 / ( dx01 * dx12 * dx13 ) ;
    const double c2 =   dx0 * dx1 *       dx3 / ( dx02 * dx12 * dx23 ) ;
    const double c3 = - dx0 * dx1 * dx2       / ( dx03 * dx13 * dx23 ) ;
    //
    const double vv = 
      v0.value() * c0 + 
      v1.value() * c1 + 
      v2.value() * c2 + 
      v3.value() * c3 ;
    //
    const double e2 = 
      ( 0 < v0.cov2() ? v0.cov2 () * c0 * c0 : 0.0 ) + 
      ( 0 < v1.cov2() ? v1.cov2 () * c1 * c1 : 0.0 ) +
      ( 0 < v2.cov2() ? v2.cov2 () * c2 * c2 : 0.0 ) +
      ( 0 < v3.cov2() ? v3.cov2 () * c3 * c3 : 0.0 ) ;
    //
    return Ostap::Math::ValueWithError ( vv , e2 ) ; 
  }
  // ==========================================================================
  // bilinear interpolation 
  inline Ostap::Math::ValueWithError _bilinear_ 
  ( const double                       x   , 
    const double                       y   , 
    const double                       x0  , 
    const double                       x1  , 
    const double                       y0  , 
    const double                       y1  , 
    //
    const Ostap::Math::ValueWithError& v00 , 
    const Ostap::Math::ValueWithError& v10 , 
    //
    const Ostap::Math::ValueWithError& v01 ,
    const Ostap::Math::ValueWithError& v11 )
  {
    return _linear_ 
      ( x  , 
        x0 , 
        x1 ,
        _linear_ ( y , y0 , y1 , v00 , v01 ) , 
        _linear_ ( y , y0 , y1 , v10 , v11 ) ) ;
  }
  // ==========================================================================
  // biquadratic interpolation 
  inline 
  Ostap::Math::ValueWithError _biquadratic_ 
  ( const double                       x   , 
    const double                       y   , 
    const double                       x0  , 
    const double                       x1  , 
    const double                       x2  , 
    const double                       y0  , 
    const double                       y1  , 
    const double                       y2  , 
    //
    const Ostap::Math::ValueWithError& v00 , 
    const Ostap::Math::ValueWithError& v10 , 
    const Ostap::Math::ValueWithError& v20 , 
    //
    const Ostap::Math::ValueWithError& v01 , 
    const Ostap::Math::ValueWithError& v11 , 
    const Ostap::Math::ValueWithError& v21 , 
    //
    const Ostap::Math::ValueWithError& v02 , 
    const Ostap::Math::ValueWithError& v12 , 
    const Ostap::Math::ValueWithError& v22 )
  {
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _quadratic_ ( y , y0 , y1 , y2 , v00 , v01 , v02 ) , 
        _quadratic_ ( y , y0 , y1 , y2 , v10 , v11 , v12 ) , 
        _quadratic_ ( y , y0 , y1 , y2 , v20 , v21 , v22 ) ) ;
  }
  // ==========================================================================
  // linear
  inline std::array<unsigned int,2> _linear_indices_ 
  ( const unsigned int ib , 
    const unsigned int nb , 
    const double       x  , 
    const double       xc ) 
  {
    const unsigned int ix0 = 
      1  >= ib ? 1      : 
      nb <= ib ? nb - 1 : 
      x  <  xc ? ib - 1 : ib ;
    return { ix0 , ix0 + 1 } ;
  }
  // quadratic 
  inline std::array<unsigned int,3> _quadratic_indices_ 
  ( const unsigned int    ib    , 
    const unsigned int    nb    , 
    const double       /* x  */ , 
    const double       /* xc */ ) 
  {
    const unsigned int ix0 = 
      1  >= ib     ? 1      : 
      nb <= ib + 1 ? nb - 2 : ib - 1 ;
    return { ix0 , ix0 + 1 , ix0 + 2 } ;
  }
  // cubic 
  inline std::array<unsigned int,4> _cubic_indices_   
  ( const unsigned int    ib    , 
    const unsigned int    nb    , 
    const double          x     , 
    const double          xc    ) 
  {
    const unsigned int ix0 = 
      2  >= ib     ? 1      : 
      nb <= ib + 1 ? nb - 3 :
      x < xc       ? ib - 2 : ib - 1 ;
    return { ix0 , ix0 + 1 , ix0 + 2 , ix0 + 3 } ;
  }
  // bi-quadratic: the same as  cubic 
  inline std::array<unsigned int,4> _quadratic2_indices_   
  ( const unsigned int    ib    , 
    const unsigned int    nb    , 
    const double          x     , 
    const double          xc    ) 
  { return _cubic_indices_ ( ib , nb , x , xc ) ; } 
  // ==========================================================================
}
// ============================================================================
// 1D-interpolation 
// ============================================================================
/** linear interpolation between two points 
 *  @param x  (INPUT) the x-value
 *  @param x0 (INPUT) x-coordinate of the first  point 
 *  @param x1 (INPUT) x-coordinate of the second point 
 *  @param y0 (INPUT) y-coordinate of the first  point 
 *  @param y1 (INPUT) y-coordinate of the second point 
 *  @return result of linear interpolation
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate
( const double                       x  , 
  const double                       x0 , 
  const double                       x1 ,
  const Ostap::Math::ValueWithError& y0 , 
  const Ostap::Math::ValueWithError& y1 ) 
{ return _linear_ ( x , x0, x1 , y0 , y1 ) ; }
// ============================================================================
/* quadratic (parabolic)  interpolation between three points 
 *  @param x  (INPUT) the x-value
 *  @param x0 (INPUT) x-coordinate of the first  point 
 *  @param x1 (INPUT) x-coordinate of the second point 
 *  @param x2 (INPUT) x-coordinate of the third  point 
 *  @param y0 (INPUT) y-coordinate of the first  point 
 *  @param y1 (INPUT) y-coordinate of the second point 
 *  @param x2 (INPUT) x-coordinate of the third  point 
 *  @return result of  quadratic (parabolic) interpolation
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate
( const double                       x  , 
  const double                       x0 , 
  const double                       x1 ,
  const double                       x2 ,
  const Ostap::Math::ValueWithError& y0 , 
  const Ostap::Math::ValueWithError& y1 , 
  const Ostap::Math::ValueWithError& y2 ) 
{ return _quadratic_( x , x0, x1 , x2 , y0 , y1 , y2 ) ; }
// ======================================================================
/*  qubic interpolation between four points 
 *  @param x  (INPUT) the x-value
 *  @param x0 (INPUT) x-coordinate of the first  point 
 *  @param x1 (INPUT) x-coordinate of the second point 
 *  @param x2 (INPUT) x-coordinate of the third  point 
 *  @param x3 (INPUT) x-coordinate of the fourth point 
 *  @param y0 (INPUT) y-coordinate of the first  point \f$ y(x_0) \f$  
 *  @param y1 (INPUT) y-coordinate of the second point \f$ y(x_1) \f$ 
 *  @param y2 (INPUT) x-coordinate of the third  point \f$ y(x_2) \f$ 
 *  @param y3 (INPUT) x-coordinate of the third  point \f$ y(x_3) \f$ 
 *  @return result of  quadratic (parabolic) interpolation  \f$ y(x) \f$
 */
// ======================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate
( const double                       x  , 
  const double                       x0 , 
  const double                       x1 ,
  const double                       x2 ,
  const double                       x3 ,
  const Ostap::Math::ValueWithError& y0 , 
  const Ostap::Math::ValueWithError& y1 , 
  const Ostap::Math::ValueWithError& y2 ,
  const Ostap::Math::ValueWithError& y3 ) 
{ return _cubic_( x , x0, x1 , x2 , x3 , y0 , y1 , y2 , y3 ) ; }
// ============================================================================
// 2D-interpolation 
// ============================================================================
/* bi-linear interpolation on grid 
 *  @param x   (INPUT) the x-value
 *  @param y   (INPUT) the x-value
 *  @param x0  (INPUT) the first  x-coordinate on the grid 
 *  @param x1  (INPUT) the second x-coordinate on the grid 
 *  @param y0  (INPUT) the first  y-coordinate on the grid 
 *  @param y1  (INPUT) the second y-coordinate on the grid 
 *  @param f00 (INPUT) function value for (x0,y0)
 *  @param f10 (INPUT) function value for (x1,y0)
 *  @param f01 (INPUT) function value for (x0,y1)
 *  @param f11 (INPUT) function value for (x1,y1)
 *  @return result of bi-linear interpolation
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate
( const double                       x   , 
  const double                       y   , 
  const double                       x0  , 
  const double                       x1  ,
  const double                       y0  , 
  const double                       y1  ,
  const Ostap::Math::ValueWithError& f00 , 
  const Ostap::Math::ValueWithError& f10 ,
  const Ostap::Math::ValueWithError& f01 , 
  const Ostap::Math::ValueWithError& f11 ) 
{
  return _bilinear_ ( x   , y   ,  
                      x0  , x1  , 
                      y0  , y1  , 
                      f00 , f10 , 
                      f01 , f11 ) ;
}
// ============================================================================
/* bi-quadratic interpolation on grid 
 *  @param x    (INPUT) the x-value
 *  @param y    (INPUT) the x-value
 *  @param x0   (INPUT) the first  x-coordinate on the grid 
 *  @param x1   (INPUT) the second x-coordinate on the grid 
 *  @param x2   (INPUT) the third  x-coordinate on the grid 
 *  @param y0   (INPUT) the first  y-coordinate on the grid 
 *  @param y1   (INPUT) the second y-coordinate on the grid 
 *  @param y2   (INPUT) the third  y-coordinate on the grid 
 *  @param f00  (INPUT) function value for (x0,y0)
 *  @param f10  (INPUT) function value for (x1,y0)
 *  @param f20  (INPUT) function value for (x2,y0)
 *  @param f01  (INPUT) function value for (x0,y1)
 *  @param f11  (INPUT) function value for (x1,y1)
 *  @param f21  (INPUT) function value for (x2,y1)
 *  @param f02  (INPUT) function value for (x0,y2)
 *  @param f12  (INPUT) function value for (x1,y2)
 *  @param f22  (INPUT) function value for (x2,y2)
 *  @return result of bi-quadratic interpolation
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate
( const double                       x   , 
  const double                       y   , 
  const double                       x0  , 
  const double                       x1  ,
  const double                       x2  ,
  const double                       y0  , 
  const double                       y1  ,
  const double                       y2  ,
  const Ostap::Math::ValueWithError& f00 , 
  const Ostap::Math::ValueWithError& f10 ,
  const Ostap::Math::ValueWithError& f20 ,
  const Ostap::Math::ValueWithError& f01 , 
  const Ostap::Math::ValueWithError& f11 ,
  const Ostap::Math::ValueWithError& f21 ,
  const Ostap::Math::ValueWithError& f02 , 
  const Ostap::Math::ValueWithError& f12 ,
  const Ostap::Math::ValueWithError& f22 ) 
{
  return _biquadratic_ 
    ( x   , y   , 
      x0  , x1  , x2  , 
      y0  , y1  , y2  , 
      f00 , f10 , f20 , 
      f01 , f11 , f21 , 
      f02 , f12 , f22 ) ;      
}
// ============================================================================
/* interpolate 1D historgam 
 *  @param h1          (INPUT) input histogram 
 *  @param x           (INPUT) the x-value 
 *  @param type        (INPUT) interpolation type 
 *  @param edges       (INPUT) use the special treatment of edges ? 
 *  @param extrapolate (INPUT) use extrapolation ? 
 *  @param density     (INPUT) rescale to density? 
 *  If "density" flag is activated, actually   the value of 
 *  density function, defined as a ratio of bin content over bin volume 
 *  is interpolated 
 *  @return value of interpolated function
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate_1D
( const TH1&                                  h1          ,   
  const double                                x           ,
  const Ostap::Math::HistoInterpolation::Type t           , 
  const bool                                  edges       , 
  const bool                                  extrapolate , 
  const bool                                  density     )  
{
  const TAxis* ax   = h1.GetXaxis()  ;
  if ( 0 == ax )                   { return ValueWithError() ; } // RETURN 
  //
  // get main parameters and treat the special "easy" cases:
  //
  const double       xmin  = ax->GetXmin  () ;
  if ( !extrapolate && xmin > x ) { return ValueWithError() ; } // RETURN 
  //
  const double       xmax  = ax->GetXmax  () ;
  if ( !extrapolate && xmax < x ) { return ValueWithError() ; } // RETURN 
  //
  if ( edges && !extrapolate && s_equal ( x , xmin ) ) 
  { return _bin_ ( h1 , 1 , density ) ; }                        // RETURN 
  //
  const unsigned int nbins = ax->GetNbins () ;
  //
  if ( edges && !extrapolate && s_equal ( x , xmax ) ) 
  { return _bin_ ( h1 , nbins , density ) ; }                     // RETURN 
  // 
  // adjust the interpolation type 
  //
  const Type itype = 
    ( t <= Nearest                  ) ? Nearest    : 
    ( 1 >= nbins                    ) ? Nearest    : 
    ( 2 == nbins  && t >= Linear    ) ? Linear     :
    ( 3 == nbins  && t >= Quadratic ) ? Quadratic  :
    ( 4 == nbins  && t >= Cubic     ) ? Cubic      :
    (                t >= Cubic     ) ? Cubic      : t ;
  //
  // the regular case 
  //
  unsigned int ib = ax->FindFixBin ( x ) ;
  //
  if      ( extrapolate &&         0 == ib ) { ib =     1 ; }
  else if ( extrapolate && nbins + 1 == ib ) { ib = nbins ; }
  //  
  if ( 0 == ib || nbins + 1 == ib ) { return ValueWithError() ;  } // RETURN
  //
  if ( Nearest == itype ) { return _bin_ ( h1 , ib , density ) ; } // RETURN
  //
  const double xc = ax->GetBinCenter ( ib ) ;
  //
  // use Nearest, if we at the bin center 
  if ( s_equal ( xc , x ) ) 
  { return _bin_ ( h1 , ib , density ) ; }                         // RETURN 
  //
  // treat left of the first and right half of the last bin separately:
  //
  if ( edges && !extrapolate && ( ( 1 == ib && x <= xc ) || ( nbins == ib && xc <= x ) ) )
  { return _bin_ ( h1 , ib , density ) ; }                        // RETURN  
  //
  // linear interpolation 
  //
  if ( Linear == itype )
  {
    // CORRECT LINEAR INTERPOLATION IN X 
    const std::array<unsigned int,2> ix = _linear_indices_ ( ib , nbins , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _bin_ ( h1 , ix[0] , density ) , 
        _bin_ ( h1 , ix[1] , density ) ) ;
  }
  //
  if ( Quadratic == itype && 3 == nbins ) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ib , nbins , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _bin_ ( h1 , ix[0] , density ) , 
        _bin_ ( h1 , ix[1] , density ) , 
        _bin_ ( h1 , ix[2] , density ) ) ;
  }
  if ( Quadratic == itype && 4 <= nbins ) 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ib , nbins , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _bin_ ( h1 , ix[0] , density ) , 
        _bin_ ( h1 , ix[1] , density ) , 
        _bin_ ( h1 , ix[2] , density ) ,
        _bin_ ( h1 , ix[3] , density ) ) ; 
  }
  //
  // CORRECT CUBIC INTERPOLATION IN X 
  const std::array<unsigned int,4> ix = _cubic_indices_ ( ib , nbins , x , xc ) ;
  //
  const double x0 = ax->GetBinCenter ( ix[0] ) ;
  const double x1 = ax->GetBinCenter ( ix[1] ) ;
  const double x2 = ax->GetBinCenter ( ix[2] ) ;
  const double x3 = ax->GetBinCenter ( ix[3] ) ;
  //
  return _cubic_ 
    ( x , x0 , x1 , x2 , x3 , 
      _bin_ ( h1 , ix[0] , density ) , 
      _bin_ ( h1 , ix[1] , density ) , 
      _bin_ ( h1 , ix[2] , density ) ,
      _bin_ ( h1 , ix[3] , density ) ) ;
}
// ============================================================================
/*  interpolate 2D historgam 
 *  @param  h2         (INPUT) input histogram 
 *  @param  x          (INPUT) the value 
 *  @param  y          (INPUT) the value 
 *  @param tx          (INPUT) interpolation type in x-direction
 *  @param ty          (INPUT) interpolation type in y-direction
 *  @param edges       (INPUT) use the special treatment of edges ? 
 *  @param extrapolate (INPUT) use extrapolation ? 
 *  @param density     (INPUT) rescale to density? 
 *  If "density" flag is activated, actually   the value of 
 *  density function, defined as a ratio of bin content over bin volume 
 *  is interpolated 
 *  @return valeu of interpolated function
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate_2D
( const TH2&                                  h2          , 
  const double                                x           ,
  const double                                y           ,
  const Ostap::Math::HistoInterpolation::Type tx          , 
  const Ostap::Math::HistoInterpolation::Type ty          , 
  const bool                                  edges       , 
  const bool                                  extrapolate , 
  const bool                                  density     )
{
  const TAxis* ax   = h2.GetXaxis()  ;
  if ( 0 == ax )                   { return ValueWithError() ; } // RETURN 
  const TAxis* ay   = h2.GetYaxis()  ;
  if ( 0 == ay )                   { return ValueWithError() ; } // RETURN 
  //
  // get main parameters and treat the special "easy" cases:
  //
  const double       xmin  = ax->GetXmin  () ;
  if ( !extrapolate && xmin > x )  { return ValueWithError() ; } // RETURN 
  //
  const double       xmax  = ax->GetXmax  () ;
  if ( !extrapolate && xmax < x )  { return ValueWithError() ; } // RETURN 
  //
  const double       ymin  = ay->GetXmin  () ;
  if ( !extrapolate && ymin > y )  { return ValueWithError() ; } // RETURN 
  //
  const double       ymax  = ay->GetXmax  () ;
  if ( !extrapolate && ymax < y )  { return ValueWithError() ; } // RETURN 
  //  
  const unsigned int nbx = ax->GetNbins () ;
  const unsigned int nby = ay->GetNbins () ;
  //
  // adjust the interpolation type 
  //
  Type itypex = 
    ( tx <= Nearest                ) ? Nearest    : 
    ( 1  >= nbx                    ) ? Nearest    : 
    ( 2  == nbx && tx >= Linear    ) ? Linear     :
    ( 3  == nbx && tx >= Quadratic ) ? Quadratic  :
    ( 4  == nbx && tx >= Cubic     ) ? Cubic      :
    (              tx >= Cubic     ) ? Cubic      : tx ;  
  //
  Type itypey = 
    ( ty <= Nearest                ) ? Nearest    : 
    ( 1  >= nby                    ) ? Nearest    : 
    ( 2  == nby && ty >= Linear    ) ? Linear     :
    ( 3  == nby && ty >= Quadratic ) ? Quadratic  :
    ( 4  == nby && ty >= Cubic     ) ? Cubic      :
    (              ty >= Cubic     ) ? Cubic      : ty ;
  //
  // find the bin 
  //
  unsigned int ibx = ax -> FindFixBin ( x ) ;
  unsigned int iby = ay -> FindFixBin ( y ) ;
  //
  if      ( 0       == ibx && s_equal ( x , xmin ) ) { ibx += 1 ; }
  else if ( nbx + 1 == ibx && s_equal ( x , xmax ) ) { ibx -= 1 ; }
  if      ( 0       == iby && s_equal ( y , ymin ) ) { iby += 1 ; }
  else if ( nby + 1 == iby && s_equal ( y , ymax ) ) { iby -= 1 ; }
  //
  if      ( extrapolate &&       0 == ibx ) { ibx =   1 ; }
  else if ( extrapolate && nbx + 1 == ibx ) { ibx = nbx ; }
  if      ( extrapolate &&       0 == iby ) { iby =   1 ; }
  else if ( extrapolate && nby + 1 == iby ) { iby = nby ; }
  //
  if ( 0 == ibx || nbx < ibx ) { return ValueWithError () ; }      // RETURN 
  if ( 0 == iby || nby < iby ) { return ValueWithError () ; }      // RETURN 
  //
  if ( Nearest == itypex && Nearest == itypey ) 
  { return _bin_ ( h2 , ibx , iby , density ) ; }                  // RETURN 
  //
  // get bin centres 
  //
  const double xc = ax->GetBinCenter ( ibx ) ;
  const double yc = ay->GetBinCenter ( iby ) ;
  //
  // special treatment of edges 
  if ( edges && !extrapolate && ( ( 1 == ibx && x <= xc ) || ( nbx == ibx && xc <= x ) ) ) { itypex = Nearest ; }
  if ( edges && !extrapolate && ( ( 1 == iby && y <= yc ) || ( nby == iby && yc <= y ) ) ) { itypey = Nearest ; }
  //
  // adjust the interpolation rules if needed 
  if ( Nearest != itypex && s_equal ( xc , x ) ) { itypex = Nearest ; }
  if ( Nearest != itypey && s_equal ( yc , y ) ) { itypey = Nearest ; }
  //
  if ( Nearest == itypex && Nearest == itypey )     // (1) 
  { return _bin_ ( h2 , ibx , iby , density ) ; }                 // RETURN
  //
  //
  if      ( Nearest == itypex && Nearest == itypey )     // (1) 
  { return _bin_ ( h2 , ibx , iby , density ) ; }
  //
  else if ( Nearest == itypex && Linear  == itypey )     // (2) 
  {
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _linear_ 
      ( y , y0 , y1 , 
        _bin_ ( h2 , ibx , iy[0] , density ) ,
        _bin_ ( h2 , ibx , iy[1] , density ) ) ;
  }
  // 
  else if ( Nearest == itypex && Quadratic == itypey && 3 == nby)  // (3) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _quadratic_ 
      ( y , y0 , y1 , y2 , 
        _bin_ ( h2 , ibx , iy[0] , density ) , 
        _bin_ ( h2 , ibx , iy[1] , density ) , 
        _bin_ ( h2 , ibx , iy[2] , density ) ) ;
  }
  else if ( Nearest == itypex && Quadratic == itypey )              // (3') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic2_ 
      ( y , y0 , y1 , y2 , y3 , 
        _bin_ ( h2 , ibx , iy[0] , density ) , 
        _bin_ ( h2 , ibx , iy[1] , density ) , 
        _bin_ ( h2 , ibx , iy[2] , density ) ,
        _bin_ ( h2 , ibx , iy[3] , density ) ) ;
  }
  //
  else if ( Nearest == itypex && Cubic == itypey )   //  (4) 
  {
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _cubic_ 
      ( y , y0 , y1 , y2 , y3 , 
        _bin_ ( h2 , ibx , iy[0] , density ) , 
        _bin_ ( h2 , ibx , iy[1] , density ) , 
        _bin_ ( h2 , ibx , iy[2] , density ) ,
        _bin_ ( h2 , ibx , iy[3] , density ) ) ;
  }
  //
  else if ( Linear == itypex && Nearest == itypey )    // (5)
  {
    // CORRECT LINEAR INTERPOLATION IN X 
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _bin_ ( h2 , ix[0] , iby , density ) ,
        _bin_ ( h2 , ix[1] , iby , density ) ) ;
  }
  //
  else if ( Linear == itypex && Linear == itypey )    // (6) 
  {
    // CORRECT LINEAR INTERPOLATION IN X 
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _bilinear_ 
      ( x    ,   y   , 
        x0 , x1 , 
        y0 , y1 ,
        _bin_ ( h2 , ix[0] , iy[0] , density ) , 
        _bin_ ( h2 , ix[1] , iy[0] , density ) , 
        _bin_ ( h2 , ix[0] , iy[1] , density ) , 
        _bin_ ( h2 , ix[1] , iy[1] , density ) ) ;
  }
  //
  else if ( Linear == itypex && Quadratic == itypey &&  3 == nby )   // (7)  
  {   
    // CORRECT LINEAR INTERPOLATION IN X 
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _linear_ 
      (  x , x0 , x1 , 
         _quadratic_ ( y , y0 , y1 , y2 , 
                       _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                       _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                       _bin_ ( h2 , ix[0] , iy[2] , density ) ) , 
         _quadratic_ ( y , y0 , y1 , y2 , 
                       _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                       _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                       _bin_ ( h2 , ix[1] , iy[2] , density ) ) ) ;
  }
  //
  else if ( Linear == itypex && Quadratic == itypey    )          // (7')  
  {   
    // CORRECT LINEAR INTERPOLATION IN X 
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _linear_ 
      (  x , x0 , x1 , 
         _quadratic2_ ( y , y0 , y1 , y2 , y3 , 
                        _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                        _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                        _bin_ ( h2 , ix[0] , iy[2] , density ) ,
                        _bin_ ( h2 , ix[0] , iy[3] , density ) ) ,
         _quadratic2_ ( y , y0 , y1 , y2 , y3 , 
                        _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                        _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                        _bin_ ( h2 , ix[1] , iy[2] , density ) ,
                        _bin_ ( h2 , ix[1] , iy[3] , density ) ) ) ;
  }
  //
  else if ( Linear == itypex && Cubic == itypey )   // (8)  
  {
    // CORRECT LINEAR INTERPOLATION IN X 
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[0] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[1] , iy[3] , density ) ) ) ;
  }
  else if ( Quadratic == itypex && Nearest == itypey &&  3 == nbx )  // (9) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double        x0 = ax->GetBinCenter ( ix[0] ) ;
    const double        x1 = ax->GetBinCenter ( ix[1] ) ;
    const double        x2 = ax->GetBinCenter ( ix[2] ) ;
    //    
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _bin_ ( h2 , ix[0] , iby , density ) , 
        _bin_ ( h2 , ix[1] , iby , density ) , 
        _bin_ ( h2 , ix[2] , iby , density ) ) ;
  }
  //
  else if ( Quadratic == itypex && Nearest == itypey )  // (9') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double        x0 = ax->GetBinCenter ( ix[0] ) ;
    const double        x1 = ax->GetBinCenter ( ix[1] ) ;
    const double        x2 = ax->GetBinCenter ( ix[2] ) ;
    const double        x3 = ax->GetBinCenter ( ix[3] ) ;
    //    
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _bin_ ( h2 , ix[0] , iby , density ) , 
        _bin_ ( h2 , ix[1] , iby , density ) , 
        _bin_ ( h2 , ix[2] , iby , density ) ,
        _bin_ ( h2 , ix[3] , iby , density ) ) ;
  }
  //
  else if ( Quadratic == itypex &&  Linear == itypey && 3 == nbx )  // (10)
  {
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h2 , ix[0] , iy[0] , density ) ,
                   _bin_ ( h2 , ix[0] , iy[1] , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                   _bin_ ( h2 , ix[1] , iy[1] , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                   _bin_ ( h2 , ix[2] , iy[1] , density ) ) ) ;
  }
  //
  else if ( Quadratic == itypex &&  Linear == itypey )  // (10')
  {
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h2 , ix[0] , iy[0] , density ) ,
                   _bin_ ( h2 , ix[0] , iy[1] , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                   _bin_ ( h2 , ix[1] , iy[1] , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                     _bin_ ( h2 , ix[2] , iy[1] , density ) ) ,
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h2 , ix[3] , iy[0] , density ) , 
                   _bin_ ( h2 , ix[3] , iy[1] , density ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && 3 == nbx && 3 == nby )  // (11) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _biquadratic_ 
      ( x    , y  , 
        x0   , x1 , x2 , 
        y0   , y1 , y2 , 
        _bin_ ( h2 , ix[0] , iy[0] , density ) , 
        _bin_ ( h2 , ix[1] , iy[0] , density ) , 
        _bin_ ( h2 , ix[2] , iy[0] , density ) , 
        _bin_ ( h2 , ix[0] , iy[1] , density ) , 
        _bin_ ( h2 , ix[1] , iy[1] , density ) , 
        _bin_ ( h2 , ix[2] , iy[1] , density ) , 
        _bin_ ( h2 , ix[0] , iy[2] , density ) , 
        _bin_ ( h2 , ix[1] , iy[2] , density ) , 
        _bin_ ( h2 , ix[2] , iy[2] , density ) ) ;
  }
  else if ( Quadratic == itypex && Quadratic == itypey && 3 == nbx )  // (11') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,
                       _bin_ ( h2 , ix[0] , iy[0] , density ) ,
                       _bin_ ( h2 , ix[0] , iy[1] , density ) ,
                       _bin_ ( h2 , ix[0] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[0] , iy[3] , density ) ) ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,
                       _bin_ ( h2 , ix[1] , iy[0] , density ) ,
                       _bin_ ( h2 , ix[1] , iy[1] , density ) ,
                       _bin_ ( h2 , ix[1] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[1] , iy[3] , density ) ) ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,
                       _bin_ ( h2 , ix[2] , iy[0] , density ) ,
                       _bin_ ( h2 , ix[2] , iy[1] , density ) ,
                       _bin_ ( h2 , ix[2] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[2] , iy[3] , density ) ) ) ;
  }
  else if ( Quadratic == itypex && Quadratic == itypey && 3 == nby )  // (11'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic_ ( y , y0 , y1 , y2 ,
                      _bin_ ( h2 , ix[0] , iy[0] , density ) ,
                      _bin_ ( h2 , ix[0] , iy[1] , density ) ,
                      _bin_ ( h2 , ix[0] , iy[2] , density ) ) ,
        _quadratic_ ( y , y0 , y1 , y2 , 
                      _bin_ ( h2 , ix[1] , iy[0] , density ) ,
                      _bin_ ( h2 , ix[1] , iy[1] , density ) ,
                      _bin_ ( h2 , ix[1] , iy[2] , density ) ) ,
        _quadratic_ ( y , y0 , y1 , y2 ,
                      _bin_ ( h2 , ix[2] , iy[0] , density ) ,
                      _bin_ ( h2 , ix[2] , iy[1] , density ) ,
                      _bin_ ( h2 , ix[2] , iy[2] , density ) ) ,
        _quadratic_ ( y , y0 , y1 , y2 ,
                      _bin_ ( h2 , ix[3] , iy[0] , density ) ,
                      _bin_ ( h2 , ix[3] , iy[1] , density ) ,
                      _bin_ ( h2 , ix[3] , iy[2] , density ) ) ) ;
  }
  else if ( Quadratic == itypex && Quadratic == itypey )  // (11''') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,
                       _bin_ ( h2 , ix[0] , iy[0] , density ) ,
                       _bin_ ( h2 , ix[0] , iy[1] , density ) ,
                       _bin_ ( h2 , ix[0] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[0] , iy[3] , density ) ) ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,
                       _bin_ ( h2 , ix[1] , iy[0] , density ) ,
                       _bin_ ( h2 , ix[1] , iy[1] , density ) ,
                       _bin_ ( h2 , ix[1] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[1] , iy[3] , density ) ) ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,
                       _bin_ ( h2 , ix[2] , iy[0] , density ) ,
                       _bin_ ( h2 , ix[2] , iy[1] , density ) ,
                       _bin_ ( h2 , ix[2] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[2] , iy[3] , density ) ) ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,
                       _bin_ ( h2 , ix[3] , iy[0] , density ) ,
                       _bin_ ( h2 , ix[3] , iy[1] , density ) ,
                       _bin_ ( h2 , ix[3] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[3] , iy[3] , density ) ) );  
  }
  //
  else if ( Quadratic == itypex && Cubic == itypey && 3 == nbx )  // (12) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic_ 
      ( x  , x0 , x1 , x2 , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[0] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[1] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[2] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[2] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[2] , iy[3] , density ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Cubic == itypey )  // (12') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic2_ 
      ( x  , x0 , x1 , x2 , x3 , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[0] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[1] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[2] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[2] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[2] , iy[3] , density ) ) ,
        _cubic_ ( y , y0 , y1 , y2 , y3 , 
                  _bin_ ( h2 , ix[3] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[3] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[3] , iy[2] , density ) ,
                  _bin_ ( h2 , ix[3] , iy[3] , density ) ) ) ;
  }
  else if ( Cubic == itypex && Nearest == itypey )  // (13) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _bin_ ( h2 , ix[0] , iby , density ) , 
        _bin_ ( h2 , ix[1] , iby , density ) , 
        _bin_ ( h2 , ix[2] , iby , density ) , 
        _bin_ ( h2 , ix[3] , iby , density ) ) ;
  }
  //
  else if ( Cubic == itypex && Linear == itypey )  // (14) 
  {
    // CORRECT CUBIC INTERPOLATION IN X
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,
        _linear_ ( y, y0 , y1 , 
                   _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                   _bin_ ( h2 , ix[0] , iy[1] , density ) ) , 
        _linear_ ( y, y0 , y1 , 
                   _bin_ ( h2 , ix[1] , iy[0] , density ) ,
                   _bin_ ( h2 , ix[1] , iy[1] , density ) ) , 
        _linear_ ( y, y0 , y1 , 
                   _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                   _bin_ ( h2 , ix[2] , iy[1] , density ) ) , 
        _linear_ ( y, y0 , y1 , 
                   _bin_ ( h2 , ix[3] , iy[0] , density ) , 
                   _bin_ ( h2 , ix[3] , iy[1] , density ) ) ) ;
  }
  else if ( Cubic == itypex && Quadratic == itypey && 3 == nby)  // (15) 
  {
    // CORRECT CUBIC INTERPOLATION IN X
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,
        _quadratic_ ( y, y0 , y1 , y2 , 
                      _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                      _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                      _bin_ ( h2 , ix[0] , iy[2] , density ) ) , 
        _quadratic_ ( y, y0 , y1 , y2 , 
                      _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                      _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                      _bin_ ( h2 , ix[1] , iy[2] , density ) ) , 
        _quadratic_ ( y, y0 , y1 , y2 , 
                      _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                      _bin_ ( h2 , ix[2] , iy[1] , density ) , 
                      _bin_ ( h2 , ix[2] , iy[2] , density ) ) , 
        _quadratic_ ( y, y0 , y1 , y2 , 
                      _bin_ ( h2 , ix[3] , iy[0] , density ) , 
                      _bin_ ( h2 , ix[3] , iy[1] , density ) , 
                      _bin_ ( h2 , ix[3] , iy[2] , density ) ) ) ;
  }
  //
  else if ( Cubic == itypex && Quadratic == itypey )  //        (15') 
  {
    // CORRECT CUBIC INTERPOLATION IN X
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,
        _quadratic2_ ( y, y0 , y1 , y2 , y3 , 
                       _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                       _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                       _bin_ ( h2 , ix[0] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[0] , iy[3] , density ) ) , 
        _quadratic2_ ( y, y0 , y1 , y2 , y3 , 
                       _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                       _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                       _bin_ ( h2 , ix[1] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[1] , iy[3] , density ) ) , 
        _quadratic2_ ( y, y0 , y1 , y2 , y3 , 
                       _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                       _bin_ ( h2 , ix[2] , iy[1] , density ) , 
                       _bin_ ( h2 , ix[2] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[2] , iy[3] , density ) ) ,
        _quadratic2_ ( y, y0 , y1 , y2 , y3  , 
                       _bin_ ( h2 , ix[3] , iy[0] , density ) , 
                       _bin_ ( h2 , ix[3] , iy[1] , density ) , 
                       _bin_ ( h2 , ix[3] , iy[2] , density ) ,
                       _bin_ ( h2 , ix[3] , iy[3] , density ) ) ) ;
  }
  else if ( Cubic == itypex && Cubic == itypey )  // (16) 
  {
    // CORRECT CUBIC INTERPOLATION IN X
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,
        _cubic_ ( y , y0 , y1 , y2 , y3 ,  
                  _bin_ ( h2 , ix[0] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[2] , density ) , 
                  _bin_ ( h2 , ix[0] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 ,  
                  _bin_ ( h2 , ix[1] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[2] , density ) , 
                  _bin_ ( h2 , ix[1] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 ,  
                  _bin_ ( h2 , ix[2] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[2] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[2] , iy[2] , density ) , 
                  _bin_ ( h2 , ix[2] , iy[3] , density ) ) , 
        _cubic_ ( y , y0 , y1 , y2 , y3 ,  
                  _bin_ ( h2 , ix[3] , iy[0] , density ) , 
                  _bin_ ( h2 , ix[3] , iy[1] , density ) , 
                  _bin_ ( h2 , ix[3] , iy[2] , density ) , 
                  _bin_ ( h2 , ix[3] , iy[3] , density ) ) ) ;
  }
  //
  return _bin_ ( h2 , ibx , iby , density ) ;
}
// ============================================================================
/*  interpolate 3D historgam 
 *  @param  h3         (INPUT) input histogram 
 *  @param  x          (INPUT) the x-value 
 *  @param  y          (INPUT) the y-value 
 *  @param  y          (INPUT) the z-value 
 *  @para   type       (INPUT) interpolation type 
 *  @param edges       (INPUT) use the special treatment of edges ? 
 *  @param extrapolate (INPUT) use extrapolation ? 
 *  @param density     (INPUT) rescale to density? 
 *  If "density" flag is activated, actually   the value of 
 *  density function, defined as a ratio of bin content over bin volume 
 *  is interpolated 
 *  @return value of interpolated function
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::HistoInterpolation::interpolate_3D
( const TH3&                                  h3          , 
  const double                                x           ,
  const double                                y           ,
  const double                                z           ,
  const Ostap::Math::HistoInterpolation::Type tx          , 
  const Ostap::Math::HistoInterpolation::Type ty          , 
  const Ostap::Math::HistoInterpolation::Type tz          , 
  const bool                                  edges       , 
  const bool                                  extrapolate , 
  const bool                                  density     )
{
  const TAxis* ax   = h3.GetXaxis()  ;
  if ( 0 == ax )                   { return ValueWithError() ; } // RETURN 
  const TAxis* ay   = h3.GetYaxis()  ;
  if ( 0 == ay )                   { return ValueWithError() ; } // RETURN 
  const TAxis* az   = h3.GetZaxis()  ;
  if ( 0 == az )                   { return ValueWithError() ; } // RETURN 
  //
  // get main parameters and treat the special "easy" cases:
  const double       xmin  = ax->GetXmin  () ;
  if ( !extrapolate && xmin > x )  { return ValueWithError() ; } // RETURN 
  //
  const double       xmax  = ax->GetXmax  () ;
  if ( !extrapolate && xmax < x )  { return ValueWithError() ; } // RETURN 
  //
  const double       ymin  = ay->GetXmin  () ;
  if ( !extrapolate && ymin > y )  { return ValueWithError() ; } // RETURN 
  //
  const double       ymax  = ay->GetXmax  () ;
  if ( !extrapolate && ymax < y )  { return ValueWithError() ; } // RETURN 
  //  
  const double       zmin  = az->GetXmin  () ;
  if ( !extrapolate && zmin > z )  { return ValueWithError() ; } // RETURN 
  //
  const double       zmax  = az->GetXmax  () ;
  if ( !extrapolate && zmax < z )  { return ValueWithError() ; } // RETURN 
  //  
  const unsigned int nbx = ax->GetNbins () ;
  const unsigned int nby = ay->GetNbins () ;
  const unsigned int nbz = az->GetNbins () ;
  //
  // adjust the interpolation type 
  //
  Type itypex = 
    ( tx <= Nearest                ) ? Nearest    : 
    ( 1  >= nbx                    ) ? Nearest    : 
    ( 2  == nbx && tx >= Linear    ) ? Linear     :
    ( 3  == nbx && tx >= Quadratic ) ? Quadratic  :
    ( 4  == nbx && tx >= Cubic     ) ? Cubic      :
    (              tx >= Cubic     ) ? Cubic      : tx ;  
  //
  Type itypey = 
    ( ty <= Nearest                ) ? Nearest    : 
    ( 1  >= nby                    ) ? Nearest    : 
    ( 2  == nby && ty >= Linear    ) ? Linear     :
    ( 3  == nby && ty >= Quadratic ) ? Quadratic  :
    ( 4  == nby && ty >= Cubic     ) ? Cubic      :
    (              ty >= Cubic     ) ? Cubic      : ty ;
  //
  Type itypez = 
    ( tz <= Nearest                ) ? Nearest    : 
    ( 1  >= nbz                    ) ? Nearest    : 
    ( 2  == nbz && tz >= Linear    ) ? Linear     :
    ( 3  == nbz && tz >= Quadratic ) ? Quadratic  :
    ( 4  == nbz && tz >= Cubic     ) ? Cubic      :
    (              tz >= Cubic     ) ? Cubic      : tz ;
  //
  // find the bin 
  //
  unsigned int ibx = ax -> FindFixBin ( x ) ;
  unsigned int iby = ay -> FindFixBin ( y ) ;
  unsigned int ibz = az -> FindFixBin ( z ) ;
  //
  if      ( 0       == ibx && s_equal ( x , xmin ) ) { ibx += 1 ; }
  else if ( nbx + 1 == ibx && s_equal ( x , xmax ) ) { ibx -= 1 ; }
  if      ( 0       == iby && s_equal ( y , ymin ) ) { iby += 1 ; }
  else if ( nby + 1 == iby && s_equal ( y , ymax ) ) { iby -= 1 ; }
  if      ( 0       == ibz && s_equal ( z , zmin ) ) { ibz += 1 ; }
  else if ( nbz + 1 == ibz && s_equal ( z , zmax ) ) { ibz -= 1 ; }
  //
  if      ( extrapolate &&       0 == ibx ) { ibx =   1 ; }
  else if ( extrapolate && nbx + 1 == ibx ) { ibx = nbx ; }
  if      ( extrapolate &&       0 == iby ) { iby =   1 ; }
  else if ( extrapolate && nby + 1 == iby ) { iby = nby ; }
  if      ( extrapolate &&       0 == ibz ) { ibz =   1 ; }
  else if ( extrapolate && nbz + 1 == ibz ) { ibz = nbz ; }
  //
  if ( 0 == ibx || nbx < ibx ) { return ValueWithError () ; }      // RETURN 
  if ( 0 == iby || nby < iby ) { return ValueWithError () ; }      // RETURN 
  if ( 0 == ibz || nbz < ibz ) { return ValueWithError () ; }      // RETURN 
  //
  if ( Nearest == itypex && Nearest == itypey && Nearest == itypez ) 
  { return _bin_ ( h3 , ibx , iby , ibz , density ) ; }                      // RETURN 
  //
  // get bin centres 
  const double xc = ax->GetBinCenter ( ibx ) ;
  const double yc = ay->GetBinCenter ( iby ) ;
  const double zc = az->GetBinCenter ( ibz ) ;
  //
  // special treatment of edges 
  if ( edges && !extrapolate && ( ( 1 == ibx && x <= xc ) || ( nbx == ibx && xc <= x ) ) ) { itypex = Nearest ; }
  if ( edges && !extrapolate && ( ( 1 == iby && y <= yc ) || ( nby == iby && yc <= y ) ) ) { itypey = Nearest ; }
  if ( edges && !extrapolate && ( ( 1 == ibz && z <= zc ) || ( nbz == ibz && zc <= z ) ) ) { itypez = Nearest ; }
  //
  // redefine the interpolation rules  if we at the bin centres 
  if ( Nearest != itypex && s_equal ( xc , x ) ) { itypex = Nearest ; }
  if ( Nearest != itypey && s_equal ( yc , y ) ) { itypey = Nearest ; }
  if ( Nearest != itypez && s_equal ( zc , z ) ) { itypez = Nearest ; }
  //
  if      ( Nearest == itypex && Nearest == itypey && Nearest == itypez ) 
  { return _bin_ ( h3 , ibx , iby , ibz , density ) ; }                    // (1)   RETURN
  //
  else if ( Linear  == itypex && Nearest == itypey && Nearest == itypez )  // (2) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _bin_ ( h3 , ix[0] , iby , ibz , density ) , 
        _bin_ ( h3 , ix[1] , iby , ibz , density ) ) ;  
  }
  else if ( Quadratic == itypex && Nearest == itypey && Nearest == itypez && 3 == nbx )  // (3) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _bin_ ( h3 , ix[0] , iby , ibz , density ) , 
        _bin_ ( h3 , ix[1] , iby , ibz , density ) ,
        _bin_ ( h3 , ix[2] , iby , ibz , density ) ) ;  
  }
  else if ( Quadratic == itypex && Nearest == itypey && Nearest == itypez )  // (3') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _bin_ ( h3 , ix[0] , iby , ibz , density ) , 
        _bin_ ( h3 , ix[1] , iby , ibz , density ) ,
        _bin_ ( h3 , ix[2] , iby , ibz , density ) ,
        _bin_ ( h3 , ix[3] , iby , ibz , density ) ) ;  
  }
  else if ( Cubic    == itypex && Nearest == itypey && Nearest == itypez )  // (4) 
  {
     // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _bin_ ( h3 , ix[0] , iby , ibz , density ) , 
        _bin_ ( h3 , ix[1] , iby , ibz , density ) ,
        _bin_ ( h3 , ix[2] , iby , ibz , density ) ,
        _bin_ ( h3 , ix[3] , iby , ibz , density ) ) ;  
  }
  else if ( Nearest == itypex && Linear == itypey && Nearest == itypez )   // (5)
  {
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _linear_ 
      ( y , y0 , y1 , 
        _bin_ ( h3 , ibx , iy[0] , ibz , density ) , 
        _bin_ ( h3 , ibx , iy[1] , ibz , density ) ) ; 
  }
  else if ( Linear == itypex && Linear == itypey && Nearest == itypez )   // (6)
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 ,
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) ) ) ;
  }
  else if ( Quadratic == itypex && Linear == itypey && Nearest == itypez && 3 == nbx )   // (7)
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) ) ) ;
  }
  else if ( Quadratic == itypex && Linear == itypey && Nearest == itypez )   // (7')
  {
    // CORRECT  BI- QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) ) ,
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) ) ) ;
  }
  else if ( Cubic == itypex && Linear == itypey && Nearest == itypez )   // (8)
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) ) , 
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) ) ,
        _linear_ ( y , y0 , y1 , 
                   _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
                   _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) ) ) ;
  }
  else if ( Nearest == itypex && Quadratic == itypey && Nearest == itypez && 3 == nby )   // (9)
  {
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _quadratic_ 
      ( y , y0 , y1 , y2 ,  
        _bin_ ( h3 , ibx , iy[0] , ibz , density ) , 
        _bin_ ( h3 , ibx , iy[1] , ibz , density ) , 
        _bin_ ( h3 , ibx , iy[2] , ibz , density ) ) ;
  }
  else if ( Nearest == itypex && Quadratic == itypey && Nearest == itypez )   // (9')
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic2_ 
      ( y , y0 , y1 , y2 ,  y3 , 
        _bin_ ( h3 , ibx , iy[0] , ibz , density ) , 
        _bin_ ( h3 , ibx , iy[1] , ibz , density ) , 
        _bin_ ( h3 , ibx , iy[2] , ibz , density ) ,
        _bin_ ( h3 , ibx , iy[3] , ibz , density ) ) ;
  }
  else if ( Linear == itypex && Quadratic == itypey && Nearest == itypez && 3 == nby)   // (10)
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 ,
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ) , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ) ) ; 
  }
  else if ( Linear == itypex && Quadratic == itypey && Nearest == itypez )   // (10')
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 , 
                       _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) , 
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 , 
                       _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) ) ; 
  }
  else if ( Quadratic == itypex && Quadratic == itypey && Nearest == itypez && 3 == nbx && 3 == nby )   // (11)
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _quadratic_
      ( x , x0 , x1 , x2 , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ) , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ) , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ) ) ; 
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Nearest == itypez && 3 == nbx )   // (11')
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic_
      ( x , x0 , x1 , x2 , 
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 ,
                       _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) , 
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 , 
                       _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) , 
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 , 
                       _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[2] , iy[3] , ibz , density ) ) ) ; 
  }
  else if ( Quadratic == itypex && Quadratic == itypey && Nearest == itypez && 3 == nby )   // (11'')
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _quadratic2_
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ) , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ) , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ) ,
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[3] , iy[2] , ibz , density ) ) ) ; 
  }
  else if ( Quadratic == itypex && Quadratic == itypey && Nearest == itypez )   // (11''')
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic2_
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,  
                       _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) , 
        _quadratic2_ ( y , y0 , y1 , y2 , y3 , 
                       _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) , 
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,  
                       _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[2] , iy[3] , ibz , density ) ) ,
        _quadratic2_ ( y , y0 , y1 , y2 , y3 ,  
                       _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[3] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[3] , iy[3] , ibz , density ) ) ) ; 
  }
  //
  else if ( Cubic == itypex && Quadratic == itypey && Nearest == itypez && 3 == nby )   // (12)
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    return _cubic_
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ) , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ) , 
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ) ,
        _quadratic_ ( y , y0 , y1 , y2 ,  
                      _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
                      _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) , 
                      _bin_ ( h3 , ix[3] , iy[2] , ibz , density ) ) ) ; 
  }
  //
  else if ( Cubic == itypex && Quadratic == itypey && Nearest == itypez )   // (12')
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _cubic_
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 , 
                       _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) , 
                       _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) , 
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 , 
                       _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) , 
                       _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) , 
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 , 
                       _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ,
                       _bin_ ( h3 , ix[2] , iy[3] , ibz , density ) ) ,
        _quadratic2_ ( y , y0 , y1 , y2 ,  y3 , 
                       _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
                       _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) , 
                       _bin_ ( h3 , ix[3] , iy[2] , ibz , density ) , 
                       _bin_ ( h3 , ix[3] , iy[3] , ibz , density ) ) ) ; 
  }
  else if ( Nearest == itypex && Cubic == itypey && Nearest == itypez )   // (13)
  {
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _cubic_ 
      ( y , y0 , y1 , y2 , y3 ,  
        _bin_ ( h3 , ibx , iy[0] , ibz , density ) , 
        _bin_ ( h3 , ibx , iy[1] , ibz , density ) , 
        _bin_ ( h3 , ibx , iy[2] , ibz , density ) ,
        _bin_ ( h3 , ibx , iy[3] , ibz , density ) ) ;
  }
  else if ( Linear == itypex && Cubic == itypey && Nearest == itypez )   // (14)
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //    
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 ,
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) ,   
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) ) ;  
  }  
  else if ( Quadratic == itypex && Cubic == itypey && Nearest == itypez && 3 == nbx )   // (15)
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) , 
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) ,
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[2] , iy[3] , ibz , density ) ) ) ;  
  }  
  else if ( Quadratic == itypex && Cubic == itypey && Nearest == itypez )   // (15')
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 ,  x3 ,  
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) , 
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) ,
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[2] , iy[3] , ibz , density ) ) ,
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[3] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[3] , iy[3] , ibz , density ) ) ) ;  
  }  
  else if ( Cubic == itypex && Cubic == itypey && Nearest == itypez )   // (16)
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[0] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[0] , iy[3] , ibz , density ) ) , 
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[1] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[1] , iy[3] , ibz , density ) ) ,
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[2] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[2] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[2] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[2] , iy[3] , ibz , density ) ) , 
        _cubic_ 
        ( y , y0 , y1 , y2 , y3 ,  
          _bin_ ( h3 , ix[3] , iy[0] , ibz , density ) , 
          _bin_ ( h3 , ix[3] , iy[1] , ibz , density ) , 
          _bin_ ( h3 , ix[3] , iy[2] , ibz , density ) ,
          _bin_ ( h3 , ix[3] , iy[3] , ibz , density ) ) ) ;  
  }  
  else if ( Nearest  == itypex && Nearest == itypey && Linear == itypez )  // (17) 
  {
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _linear_ 
      ( z , z0 , z1 , 
        _bin_ ( h3 , ibx , iby , iz[0] , density ) , 
        _bin_ ( h3 , ibx , iby , iz[1] , density ) ) ;  
  }
  else if ( Linear  == itypex && Nearest == itypey && Linear == itypez )  // (18) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _linear_
      ( x , x0 , x1 ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) ) ) ;
  }
  else if ( Quadratic  == itypex && Nearest == itypey && Linear == itypez &&  3 == nbx )  // (19) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic_
      ( x , x0 , x1 , x2 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) ) ) ;
  }
  else if ( Quadratic  == itypex && Nearest == itypey && Linear == itypez )  // (19') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic2_
      ( x , x0 , x1 , x2 , x3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) ) ) ;
  }
  else if ( Cubic  == itypex && Nearest == itypey && Linear == itypez )  // (20) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _cubic_
      ( x , x0 , x1 , x2 , x3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) ) ) ;
  }
  else if ( Nearest  == itypex && Linear == itypey && Linear == itypez )  // (21) 
  {
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _linear_
      ( y , y0 , y1 ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) ) ) ;
  }
  else if ( Linear  == itypex && Linear == itypey && Linear == itypez )  // (22) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 ,
        _linear_
        ( y , y0 , y1 ,  
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 ,  
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) ) ) ;
  }
  else if ( Quadratic  == itypex && Linear == itypey && Linear == itypez && 3 == nbx )  // (23) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _linear_
        ( y , y0 , y1 ,  
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) ) ) ; 
  }
  else if ( Quadratic  == itypex && Linear == itypey && Linear == itypez )  // (23') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _linear_
        ( y , y0 , y1 ,  
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) ) ) ) ; 
  }
  else if ( Cubic  == itypex && Linear == itypey && Linear == itypez )  // (24) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,
        _linear_
        ( y , y0 , y1 ,  
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) ) ) ) ; 
  }
  else if ( Nearest  == itypex && Quadratic == itypey && Linear == itypez && 3 == nby )  // (25) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic_
      ( y , y0 , y1 , y2 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) ) ) ;
  }

  else if ( Nearest  == itypex && Quadratic == itypey && Linear == itypez )  // (25') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic2_
      ( y , y0 , y1 , y2 , y3 ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) ) ) ;
  }
  // 
  else if ( Linear == itypex && Quadratic == itypey && Linear == itypez && 3 == nby )  // (26) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 ,
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ) ) ;  
  }
  // 
  else if ( Linear == itypex && Quadratic == itypey && Linear == itypez )  // (26') 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) ) ;  
  }
  // 
  else if ( Quadratic == itypex && Quadratic == itypey && Linear == itypez && 3 == nbx && 3 == nby )  // (27) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
   //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ) ) ;
  }
  // 
  else if ( Quadratic == itypex && Quadratic == itypey && Linear == itypez && 3 == nby )  // (27') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
   //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Linear == itypez && 3 == nbx )  // (27'') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
   //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Linear == itypez )  // (27''') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) ) ) ) ;
  }
  // 
  else if ( Cubic == itypex && Quadratic == itypey && Linear == itypez &&  3 == nby )  // (28) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ; 
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) ) ) ) ;
  }
  // 
  else if ( Cubic == itypex && Quadratic == itypey && Linear == itypez )  // (28') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ; 
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) ) ) ) ;  
  }
  // 
  else if ( Nearest  == itypex && Cubic == itypey && Linear == itypez )  // (29) 
  {
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
   //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) ) ) ;
  }
  else if ( Linear  == itypex && Cubic == itypey && Linear == itypez )  // (30) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
   //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) , 
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) ) ;
  }
  // 
  else if ( Quadratic == itypex && Cubic == itypey && Linear == itypez &&  3 == nbx )  // (31) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
   //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) , 
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) , 
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) ) ) ) ;
  }

  else if ( Quadratic == itypex && Cubic == itypey && Linear == itypez )  // (31') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
   //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) ) , 
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) ) ,
          _linear_ 
          ( z , z0 , z1 , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) ) ) ) ;
  }
  //
  else if ( Cubic == itypex && Cubic == itypey && Linear == itypez )  // (32) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
   //
    // CORRECT LINEAR INTERPOLATION IN Z
    const std::array<unsigned int,2> iz = _linear_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,  
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) ) ) , 
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) ) ) , 
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) ) ) ,
      _cubic_
      ( y , y0 , y1 , y2 , y3 , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) ) , 
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) ) ,
        _linear_ 
        ( z , z0 , z1 , 
          _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) ) ) ) ;
  }
  else if ( Nearest  == itypex && Nearest == itypey && Quadratic == itypez &&3 ==  nbz )  // (33) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic_ 
      ( z , z0 , z1 , z2 ,  
        _bin_ ( h3 , ibx , iby , iz[0] , density ) , 
        _bin_ ( h3 , ibx , iby , iz[1] , density ) , 
        _bin_ ( h3 , ibx , iby , iz[2] , density ) ) ;
  }
  else if ( Nearest  == itypex && Nearest == itypey && Quadratic == itypez )  // (33') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( z , z0 , z1 , z2 , z3 ,  
        _bin_ ( h3 , ibx , iby , iz[0] , density ) , 
        _bin_ ( h3 , ibx , iby , iz[1] , density ) , 
        _bin_ ( h3 , ibx , iby , iz[2] , density ) ,
        _bin_ ( h3 , ibx , iby , iz[3] , density ) ) ;
  }
  else if ( Linear  == itypex && Nearest == itypey && Quadratic == itypez && 3 == nbz )  // (34) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ) ) ;    
  }
  //
  else if ( Linear  == itypex && Nearest == itypey && Quadratic == itypez )  // (34') 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,   
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,    
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) ) ;    
  }
  // 
  else if ( Quadratic == itypex && Nearest == itypey && Quadratic == itypez && 3 == nbx && 3 == nbz )  // (35) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ) ) ;    
  }  
  else if ( Quadratic == itypex && Nearest == itypey && Quadratic == itypez && 3 ==  nbx )  // (35') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 ,  z3 , 
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ) ) ;    
  }
  else if ( Quadratic == itypex && Nearest == itypey && Quadratic == itypez &&  3 == nbz )  // (35'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 ,  x3 ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[2] , density ) ) ) ;    
  }  
  //
  else if ( Quadratic == itypex && Nearest == itypey && Quadratic == itypez )  // (35'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ,    
          _bin_ ( h3 , ix[2] , iby , iz[3] , density ) ) ,    
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 , 
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[2] , density ) ,   
          _bin_ ( h3 , ix[3] , iby , iz[3] , density ) ) ) ;    
  }
  // 
  else if ( Cubic == itypex && Nearest == itypey && Quadratic == itypez && 3 == nbz )  // (36) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[2] , density ) ) ) ;    
  }
  // 
  else if ( Cubic == itypex && Nearest == itypey && Quadratic == itypez )  // (36') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[2] , iby , iz[3] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 , 
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[2] , density ) ,    
          _bin_ ( h3 , ix[3] , iby , iz[3] , density ) ) ) ;    
  }
  //
  else if ( Nearest == itypex && Linear == itypey && Quadratic == itypez && 3 == nbz )  // (37) 
  {
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _linear_
      ( y , y0 , y1 , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ) ) ;
  }
  //
  else if ( Nearest == itypex && Linear == itypey && Quadratic == itypez )  // (37') 
  {
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_
      ( y , y0 , y1 , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 , 
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) ) ;
  }
  // 
  else if ( Linear == itypex && Linear == itypey && Quadratic == itypez && 3 == nbz )  // (38)
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ) ) ;
  }
  // 
  else if ( Linear == itypex && Linear == itypey && Quadratic == itypez )  // (38')
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) ) ;
  }
  //  
  else if ( Quadratic == itypex && Linear == itypey && Quadratic == itypez && 3 ==  nbx && 3 == nbz )  // (39) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Linear == itypey && Quadratic == itypez && 3 == nbx )  // (39') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,   
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ) ) ;
  }
  //  
  else if ( Quadratic == itypex && Linear == itypey && Quadratic == itypez && nbz == 3 )  // (39'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ) ) ;
  }
  //  
  else if ( Quadratic == itypex && Linear == itypey && Quadratic == itypez )  // (39''') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ) ,
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ) ) ;
  }
  //
  else if ( Cubic == itypex && Linear == itypey && Quadratic == itypez && 3 ==  nbz )  // (40) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ) ) ;
  }
  //
  else if ( Cubic == itypex && Linear == itypey && Quadratic == itypez )  // (40') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ) , 
        _linear_
        ( y , y0 , y1 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 ,  z3 , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Nearest == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nby && 3 == nbz )  // (41) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic_
      ( y , y0 , y1 , y2 ,  
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ) ) ;
  }
  //
  else if ( Nearest == itypex && Quadratic == itypey && Quadratic == itypez &&  3 == nby )  // (41') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_
      ( y , y0 , y1 , y2 ,  
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 , 
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[2] , iz[3] , density ) ) ) ;
  }
  //
  else if ( Nearest == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nbz )  // (41'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic2_
      ( y , y0 , y1 , y2 , y3 ,   
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[2] , density ) ) ) ;
  }
  //
  else if ( Nearest == itypex && Quadratic == itypey && Quadratic == itypez )  // (41''') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_
      ( y , y0 , y1 , y2 , y3 , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,   
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,   
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[3] , iz[2] , density ) ) ) ;
  }
  //  
  else if ( Linear == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nby && 3 == nbz )  // (42) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Linear == itypex && Quadratic == itypey && Quadratic == itypez && 3 ==  nby )  // (42') 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) ) ;  
  }
  //
  else if ( Linear == itypex && Quadratic == itypey && Quadratic == itypez &&  3 ==  nbz )  // (42'') 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Linear == itypex && Quadratic == itypey && Quadratic == itypez )  // (42''') 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,   
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,   
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,   
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,   
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  // 
  // SKIP IT HERE move it to the end 
  // else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez )  // (43) 
  //
  else if ( Cubic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nby && 3 == nbz )  // (44) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Cubic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nby )  // (44') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Cubic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nbz )  // (44'') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ) ) , 
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Cubic == itypex && Quadratic == itypey && Quadratic == itypez )  // (44''') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,   
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Nearest == itypex && Cubic == itypey && Quadratic == itypez && 3 == nbz )  // (45) 
  {
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _cubic_
      ( y , y0 , y1 , y2 , y3 ,  
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ) , 
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ) ,
        _quadratic_ 
        ( z , z0 , z1 , z2 ,  
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[2] , density ) ) ) ;
  }
  // 
  else if ( Nearest == itypex && Cubic == itypey && Quadratic == itypez )  // (45') 
  {
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_
      ( y , y0 , y1 , y2 , y3 ,  
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) , 
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ) ,
        _quadratic2_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[3] , iz[3] , density ) ) ) ;
  }
  // 
  else if ( Linear == itypex && Cubic == itypey && Quadratic == itypez && 3 == nbz )  // (46) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Linear == itypex && Cubic == itypey && Quadratic == itypez )  // (46') 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,   
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 ,  z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) ) ;  
  }
  // 
  else if ( Quadratic == itypex && Cubic == itypey && Quadratic == itypez && 3 == nbx && 3 == nbz )  // (47) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
   //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ) ) , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Quadratic == itypex && Cubic == itypey && Quadratic == itypez && 3 == nbx )  // (47') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
   //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,   
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ) ;  
  }
  // 
  else if ( Quadratic == itypex && Cubic == itypey && Quadratic == itypez && 3 == nbz )  // (47'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
   //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ) ) , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Quadratic == itypex && Cubic == itypey && Quadratic == itypez )  // (47''') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
   //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 ,z3 ,   
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) ,  
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Cubic == itypex && Cubic == itypey && Quadratic == itypez && 3 == nbz )  // (48) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ) ) , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Cubic == itypex && Cubic == itypey && Quadratic == itypez )  // (48') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ,
        _cubic_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  // 
  else if ( Nearest  == itypex && Nearest == itypey && Cubic == itypez )  // (49) 
  {
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( z , z0 , z1 , z2 , z3 ,  
        _bin_ ( h3 , ibx , iby , iz[0] , density ) , 
        _bin_ ( h3 , ibx , iby , iz[1] , density ) , 
        _bin_ ( h3 , ibx , iby , iz[2] , density ) ,
        _bin_ ( h3 , ibx , iby , iz[3] , density ) ) ;
  }
  //
  else if ( Linear  == itypex && Nearest == itypey && Cubic == itypez )  // (50) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //   
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) ) ;    
  }
  //
  else if ( Quadratic  == itypex && Nearest == itypey && Cubic == itypez && 3 == nbx )  // (51) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //   
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[2] , iby , iz[3] , density ) ) ) ;
  }
  //
  else if ( Quadratic  == itypex && Nearest == itypey && Cubic == itypez )  // (51') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //   
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[2] , iby , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[3] , iby , iz[3] , density ) ) ) ;
  }
  //
  else if ( Cubic  == itypex && Nearest == itypey && Cubic == itypez )  // (52) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //   
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iby , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iby , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[2] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[2] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[2] , iby , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[3] , iby , iz[0] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[1] , density ) , 
          _bin_ ( h3 , ix[3] , iby , iz[2] , density ) ,
          _bin_ ( h3 , ix[3] , iby , iz[3] , density ) ) ) ;
  }
  else if ( Nearest  == itypex && Linear == itypey && Cubic == itypez )  // (53) 
  {
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( y , y0, y1 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) ) ;    
  }
  //
  else if ( Linear == itypex && Linear == itypey && Cubic == itypez )  // (54) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) ,
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) ) ;    
  }
  //
  else if ( Quadratic == itypex && Linear == itypey && Cubic == itypez && 3 == nbx )  // (55) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) ,
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) , 
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Linear == itypey && Cubic == itypez )  // (55') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) ,
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) , 
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ) ,
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) ) ) ;
  }
  //
  else if ( Cubic == itypex && Linear == itypey && Cubic == itypez )  // (56) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT LINEAR INTERPOLATION IN Y
    const std::array<unsigned int,2> iy = _linear_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ) ,
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ) , 
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ) , 
        _linear_ 
        ( y , y0, y1 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Nearest  == itypex && Quadratic == itypey && Cubic == itypez && 3 == nby )  // (57) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( y , y0, y1 , y2 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[2] , iz[3] , density ) ) ) ;
  }
  // 
  else if ( Nearest  == itypex && Quadratic == itypey && Cubic == itypez )  // (57') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( y , y0, y1 , y2 , y3 ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[2] , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[3] , iz[3] , density ) ) ) ;
  }
  // 
  else if ( Linear  == itypex && Quadratic == itypey && Cubic == itypez &&  3 == nby )  // (58) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ) ) ;
  }
  //
  else if ( Linear  == itypex && Quadratic == itypey && Cubic == itypez )  // (58') 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_ 
      ( x , x0 , x1 , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ,
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 ,  
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Quadratic  == itypex && Quadratic == itypey && Cubic == itypez && 3 ==  nbx && 3 == nby )  // (59) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Quadratic  == itypex && Quadratic == itypey && Cubic == itypez && 3 == nbx )  // (59') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Quadratic  == itypex && Quadratic == itypey && Cubic == itypez && 3 == nby )  // (59'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3  ,  
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ) ,
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Quadratic  == itypex && Quadratic == itypey && Cubic == itypez )  // (59''') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 ,  
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ,
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 ,  
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Cubic  == itypex && Quadratic == itypey && Cubic == itypez && 3 ==  nby )  // (60) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ) ,
        _quadratic_ 
        ( y , y0, y1 , y2 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Cubic  == itypex && Quadratic == itypey && Cubic == itypez )  // (60') 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( x , x0 , x1 , x2 ,  x3 , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ,
        _quadratic2_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  // 
  else if ( Nearest  == itypex && Cubic == itypey && Cubic == itypez )  // (61) 
  {
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_ 
      ( y , y0, y1 , y2 , y3 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[0] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[0] , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[1] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[2] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ibx , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ibx , iy[3] , iz[2] , density ) ,
          _bin_ ( h3 , ibx , iy[3] , iz[3] , density ) ) ) ;
  }
  //
  else if ( Linear == itypex && Cubic == itypey && Cubic == itypez )  // (62) 
  {
    // CORRECT LINEAR INTERPOLATION IN X
    const std::array<unsigned int,2> ix = _linear_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _linear_
      ( x , x0 , x1 ,        
      _cubic_ 
      ( y , y0, y1 , y2 , y3 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
          _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) , 
      _cubic_ 
      ( y , y0, y1 , y2 , y3 , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
        _cubic_ 
        ( z , z0 , z1 , z2 , z3 ,  
          _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
          _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,
          _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Cubic == itypey && Cubic == itypez && 3 == nbx )  // (63) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_
      ( x , x0 , x1 , x2 ,         
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Cubic == itypey && Cubic == itypez )  // (63') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_
      ( x , x0 , x1 , x2 , x3 ,          
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) ,
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  //
  else if ( Cubic == itypex && Cubic == itypey && Cubic == itypez )  // (64) 
  {
    // CORRECT CUBIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _cubic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _cubic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT CUBIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _cubic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _cubic_
      ( x , x0 , x1 , x2 , x3 ,         
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[3] , iz[3] , density ) ) ) , 
        _cubic_ 
        ( y , y0, y1 , y2 , y3 , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) ,
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[2] , iz[3] , density ) ) , 
          _cubic_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[3] , iz[3] , density ) ) ) ) ;
  }
  // ==================================================================================
  // MOVED HERE 
  // ==================================================================================
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nbx && 3 ==  nby && 3 == nbz )  // (43) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nbx && 3 == nby )  // (43') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ) ) ;  
  }
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez && 3 ==  nby && 3 == nbz )  // (43'') 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nbx && 3 == nbz )  // (43''') 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nbx )  // (43^4) 
  {
    // CORRECT QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,3> ix = _quadratic_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic_ 
      ( x , x0 , x1 , x2 ,  
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez && 3 ==  nby )  // (43^5) 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,3> iy = _quadratic_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,   
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) ) , 
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ) ,
        _quadratic_
        ( y , y0 , y1 , y2 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[3] , density ) ) ) ) ;  
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez && 3 == nbz )  // (43^6) 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,3> iz = _quadratic_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _quadratic2_
        ( y , y0 , y1 , y2 , y3,   
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) ) ) , 
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) ) , 
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ,
          _quadratic_ 
          ( z , z0 , z1 , z2 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;
  }
  //
  else if ( Quadratic == itypex && Quadratic == itypey && Quadratic == itypez )  // (43^7) 
  {
    // CORRECT BI-QUADRATIC INTERPOLATION IN X 
    const std::array<unsigned int,4> ix = _quadratic2_indices_ ( ibx , nbx , x , xc ) ;
    //
    const double x0 = ax->GetBinCenter ( ix[0] ) ;
    const double x1 = ax->GetBinCenter ( ix[1] ) ;
    const double x2 = ax->GetBinCenter ( ix[2] ) ;
    const double x3 = ax->GetBinCenter ( ix[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Y 
    const std::array<unsigned int,4> iy = _quadratic2_indices_ ( iby , nby , y , yc ) ;
    //
    const double y0 = ay->GetBinCenter ( iy[0] ) ;
    const double y1 = ay->GetBinCenter ( iy[1] ) ;
    const double y2 = ay->GetBinCenter ( iy[2] ) ;
    const double y3 = ay->GetBinCenter ( iy[3] ) ;
    //
    // CORRECT BI-QUADRATIC INTERPOLATION IN Z
    const std::array<unsigned int,4> iz = _quadratic2_indices_ ( ibz , nbz , z , zc ) ;
    //
    const double z0 = az->GetBinCenter ( iz[0] ) ;
    const double z1 = az->GetBinCenter ( iz[1] ) ;
    const double z2 = az->GetBinCenter ( iz[2] ) ;
    const double z3 = az->GetBinCenter ( iz[3] ) ;
    //
    return _quadratic2_ 
      ( x , x0 , x1 , x2 , x3 ,  
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[0] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[0] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[2] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[2] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[0] , iy[3] , iz[2] , density ) ,
            _bin_ ( h3 , ix[0] , iy[3] , iz[3] , density ) ) ) ,
        _quadratic2_
        ( y , y0 , y1 , y2 , y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[1] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[2] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[1] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[2] , density ) , 
            _bin_ ( h3 , ix[1] , iy[3] , iz[3] , density ) ) ) , 
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[2] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[2] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[2] , iy[2] , iz[2] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[2] , iy[3] , iz[2] , density ) ) ) , 
        _quadratic2_
        ( y , y0 , y1 , y2 ,  y3 ,  
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[0] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[2] , density ) , 
            _bin_ ( h3 , ix[3] , iy[0] , iz[3] , density ) ) , 
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[1] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[1] , iz[2] , density ) ,
            _bin_ ( h3 , ix[3] , iy[1] , iz[3] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[3] , iy[2] , iz[2] , density ) ) ,
          _quadratic2_ 
          ( z , z0 , z1 , z2 , z3 ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[0] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[1] , density ) , 
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ,  
            _bin_ ( h3 , ix[3] , iy[3] , iz[2] , density ) ) ) ) ;  
  }
  //

  //
  // should never go here.
  // 
  return _bin_( h3 , ibx , iby , ibz , density ) ;  // RETURN 
}

// ============================================================================
// The END 
// ============================================================================
