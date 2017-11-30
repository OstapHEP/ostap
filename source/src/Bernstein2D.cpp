// ============================================================================
// Include files 
// ============================================================================
// STD& STL
// ============================================================================
#include <cmath>
#include <array>
#include <climits>
#include <cassert>
#include <numeric>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Bernstein2D.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
// ============================================================================
/** @file 
 *  Implementation file for functions, related to Bernstein's polynomnials 
 *
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein2D::Bernstein2D
( const unsigned short      nX   ,
  const unsigned short      nY   ,
  const double              xmin ,
  const double              xmax ,
  const double              ymin ,
  const double              ymax )
  : m_nx   ( nX ) 
  , m_ny   ( nY )
//
  , m_pars ( ( nX + 1 ) * ( nY + 1 ) , 0.0 )
//
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
  , m_ymin ( std::min ( ymin , ymax ) )
  , m_ymax ( std::max ( ymin , ymax ) )
    //
  , m_bx   () 
  , m_by   ()
{
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= nX ; ++ix ) 
  { m_bx.push_back ( Bernstein ( BB ( ix , nX ) , xmin , xmax ) ) ; }
  //
  for ( unsigned short iy = 0 ; iy <= nY ; ++iy ) 
  { m_by.push_back ( Bernstein ( BB ( iy , nY ) , ymin , ymax ) ) ; }
  //
}
// ============================================================================
// from symmetric variant
// ============================================================================
Ostap::Math::Bernstein2D::Bernstein2D
( const Ostap::Math::Bernstein2DSym& right ) 
  : m_nx ( right.nX() ) 
  , m_ny ( right.nY() ) 
  , m_pars ( ( right.nX () + 1 ) * ( right.nY () + 1 ) , 0.0 ) 
  , m_xmin ( right.xmin() ) 
  , m_xmax ( right.xmax() ) 
  , m_ymin ( right.xmin() ) 
  , m_ymax ( right.xmax() ) 
  , m_bx   () 
  , m_by   () 
{
  //
  m_bx.reserve ( m_nx ) ;
  m_by.reserve ( m_ny ) ;
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix ) 
  { m_bx.push_back ( Bernstein ( BB ( ix , m_nx ) , m_xmin , m_xmax ) ) ; }
  //
  for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy ) 
  { m_by.push_back ( Bernstein ( BB ( iy , m_ny ) , m_ymin , m_ymax ) ) ; }
  //
  for ( unsigned short ix = 0 ; ix <= m_nx  ; ++ix ) 
  {
    setPar ( ix , ix , right.par ( ix , ix ) ) ;
    for ( unsigned short iy = 0 ; iy < ix  ; ++iy ) 
    {
      const double p = right.par ( ix ,iy ) ;
      setPar ( ix , iy , p ) ;
      setPar ( iy , ix , p ) ;
    }
  }
}
// ============================================================================
// Move constructor 
// ============================================================================
Ostap::Math::Bernstein2D::Bernstein2D
(       Ostap::Math::Bernstein2D&& right ) 
  : m_nx   ( std::move ( right.m_nx   ) ) 
  , m_ny   ( std::move ( right.m_ny   ) ) 
  , m_pars ( std::move ( right.m_pars ) ) 
  , m_xmin ( std::move ( right.m_xmin ) ) 
  , m_xmax ( std::move ( right.m_xmax ) ) 
  , m_ymin ( std::move ( right.m_ymin ) ) 
  , m_ymax ( std::move ( right.m_ymax ) ) 
  , m_bx   ( std::move ( right.m_bx   ) ) 
  , m_by   ( std::move ( right.m_by   ) ) 
{}
// ============================================================================
// swap  two 2D-polynomials 
// ============================================================================
void Ostap::Math::Bernstein2D::swap ( Ostap::Math::Bernstein2D& right ) 
{
  std::swap ( m_nx   , right.m_nx    ) ;
  std::swap ( m_ny   , right.m_ny    ) ;
  std::swap ( m_pars , right.m_pars  ) ;
  std::swap ( m_xmin , right.m_xmin  ) ;
  std::swap ( m_xmax , right.m_xmax  ) ;
  std::swap ( m_ymin , right.m_ymin  ) ;
  std::swap ( m_ymax , right.m_ymax  ) ;
  std::swap ( m_bx   , right.m_bx    ) ;
  std::swap ( m_by   , right.m_by    ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::Bernstein2D::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
    { result += par ( ix , iy ) * fx[ix] * fy[iy] ; } 
  }
  //
  const double scalex = ( m_nx + 1 ) / ( xmax () - xmin () ) ;
  const double scaley = ( m_ny + 1 ) / ( ymax () - ymin () ) ;
  //
  return result * scalex * scaley ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein2D::evaluate 
( const double x ,
  const double y ) const
{
  /// the trivial cases
  if ( x < m_xmin || x > m_xmax ) { return 0.0        ; }
  if ( y < m_ymin || y > m_ymax ) { return 0.0        ; }
  //
  //
  if      ( 0 == npars ()       ) { return 0.0        ; }
  else if ( 1 == npars ()       ) 
  { 
    const double scalex = ( m_nx + 1 ) / ( xmax() - xmin() ) ;
    const double scaley = ( m_ny + 1 ) / ( ymax() - ymin() ) ;
    return m_pars [0] * scalex * scaley ; 
  }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i )
  { fx[i] = m_bx[i] ( x )  ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i )
  { fy[i] = m_by[i] ( y )  ; }
  //
  return calculate ( fx  ,  fy ) ;
}
// ============================================================================
/** get the integral over 2D-region 
 *  \f[  x_min < x < x_max, y_min< y< y_max\f] 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integral() const 
{ return std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ; }
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( s_equal ( xlow  , m_xmin ) &&
            s_equal ( xhigh , m_xmax ) &&
            s_equal ( ylow  , m_ymin ) &&
            s_equal ( yhigh , m_ymax )  )  { return integral () ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  else if ( xhigh <  xmin () || xlow >  xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow >  ymax() ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i )
  { fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i )
  { fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  return calculate (  fx , fy ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param x     variable 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateX 
( const double y    , 
  const double xlow , const double xhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh  ) { return -1*integrateX ( y , xhigh , xlow ) ; }
  else if ( xhigh <= xmin () || xlow >= xmax() ) { return 0 ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( s_equal ( xlow  , m_xmin ) &&
            s_equal ( xhigh , m_xmax )         ) { return integrateX ( y ) ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i )
  { fy[i] = m_by[i] ( y ) ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i )
  { fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  return calculate (  fx , fy ) ;
}
// ============================================================================
/*  integral over y-dimension 
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param x     variable 
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  >  yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( yhigh <= ymin () || ylow >= ymax() ) { return 0 ; }
  else if ( s_equal ( ylow  , m_ymin ) &&
            s_equal ( yhigh , m_ymax )         ) { return integrateY ( x ) ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i )
  { fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i )
  { fx[i] = m_bx[i] ( x ) ; }
  //
  return calculate (  fx , fy ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param x     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateX ( const double y ) const 
{
  if ( y < ymin () || y > ymax() ) { return 0 ; }
  //
  const std::vector<double> fx ( m_nx + 1 , (xmax()  -  xmin() ) / ( m_nx  + 1 ) ) ;
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i )
  { fy[i] = m_by[i] ( y ) ; }
  //
  return calculate (  fx , fy ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param y     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateY 
( const double x ) const 
{
 if ( x < xmin () || x > xmax() ) { return 0 ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i )
  { fx[i] = m_bx[i] ( x ) ; }
  //
  const std::vector<double> fy ( m_ny + 1 , (ymax()  -  ymin() ) / ( m_ny  + 1 ) ) ;
  //
  return calculate ( fx ,  fy ) ;
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Bernstein2D::setPar
( const unsigned int   k     , 
  const double         value )
{
  if ( k >= npars() )                     { return false ; }
  if ( s_equal ( m_pars [ k ] , value ) ) { return false ; }
  m_pars [ k ] = value ;
  return true ;
}
// ============================================================================
// Operators 
// ============================================================================
Ostap::Math::Bernstein2D&
Ostap::Math::Bernstein2D::operator+=( const double a )
{
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2D&
Ostap::Math::Bernstein2D::operator*=( const double a )
{
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2D&
Ostap::Math::Bernstein2D::operator-=( const double a )
{
  if ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2D&
Ostap::Math::Bernstein2D::operator/=( const double a )
{
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::operator-() const
{
  Bernstein2D b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::__add__   ( const double value ) const
{ return (*this) + value ; }
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::__radd__  ( const double value ) const
{ return value + (*this) ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::__mul__   ( const double value ) const
{ return (*this) * value ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::__rmul__  ( const double value ) const
{ return value * (*this) ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::__sub__  ( const double value ) const
{ return (*this) - value ; }
// ============================================================================
// Subtract Bernstein polynomial from a constant
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::__rsub__ ( const double value ) const
{ return value - (*this) ; }
// ============================================================================
// Divide Benrstein polynomial by a constant
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D:: __div__   ( const double value ) const
{ return (*this) / value ; }
// ============================================================================
// Negate Bernstein polynomial
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein2D::__neg__ ()  const
{ return -(*this); }
// ============================================================================





  

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein2DSym::Bernstein2DSym
( const unsigned short      n    ,
  const double              xmin ,
  const double              xmax )
  : m_n    ( n ) 
//
  , m_pars ( ( n + 1 ) * ( n + 2 ) / 2 , 0.0 )
//
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
//
  , m_b    () 
{
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short i = 0 ; i <= n ; ++i ) 
  { m_b.push_back ( Bernstein ( BB ( i , n ) , xmin , xmax ) ) ; }
  //
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::Bernstein2DSym::Bernstein2DSym
(    Ostap::Math::Bernstein2DSym&& right )
  : m_n    ( std::move ( right.m_n    )  )
  , m_pars ( std::move ( right.m_pars )  )
  , m_xmin ( std::move ( right.m_xmin )  )
  , m_xmax ( std::move ( right.m_xmax )  )
  , m_b    ( std::move ( right.m_b    )  )
{}
// ============================================================================
// swap
// ============================================================================
void Ostap::Math::Bernstein2DSym::swap( Ostap::Math::Bernstein2DSym& right )
{
  std::swap ( m_n    ,  right.m_n    ) ;
  std::swap ( m_pars ,  right.m_pars ) ;
  std::swap ( m_xmin ,  right.m_xmin ) ;
  std::swap ( m_xmax ,  right.m_xmax ) ;
  std::swap ( m_b    ,  right.m_b    ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::Bernstein2DSym::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_n ; ++ix )
  {
    result   += par ( ix , ix ) * fx[ix] * fy[ix] ;
    for  ( unsigned short iy = 0 ; iy < ix ; ++iy )
    { result += par ( ix , iy ) * ( fx[ix] * fy[iy] + fx[iy] * fy[ix] ) ; } 
  }
  //
  const double scalex = ( m_n  + 1 ) / ( xmax () - xmin () ) ;
  const double scaley = scalex ;
  //
  return result * scalex * scaley ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein2DSym::evaluate
( const double x ,
  const double y ) const
{
  /// the trivial cases
  if ( x < xmin () || x > xmax () ) { return 0.0        ; }
  if ( y < ymin () || y > ymax () ) { return 0.0        ; }
  //
  if      ( 0 == npars ()       ) { return 0.0 ; }
  else if ( 1 == npars ()       ) 
  {
    const double scale = ( m_n + 1 ) / ( xmax() - xmin() ) ;
    return m_pars [0] * ( scale * scale ) ;
  }
  ///
  std::vector<double> fy ( m_n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_n ; ++i )
  { fy[i] = m_b[i] ( y ) ; }
  //
  std::vector<double> fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i )
  { fx[i] = m_b[i] ( x ) ; }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}} 
 *  \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integral 
( const double xlow , const double xhigh ,
  const double ylow , const double yhigh ) const 
{
  //
  if      ( xlow  > xhigh   ) { return -integral ( xhigh , xlow    , ylow   , yhigh  ) ; }
  else if ( ylow  > yhigh   ) { return -integral ( xlow  , xhigh   , yhigh  , ylow   ) ; }
  //
  else if ( xlow  < xmin () ) { return  integral ( xmin() , xhigh  , ylow   , yhigh  ) ; }
  else if ( xhigh > xmax () ) { return  integral ( xlow   , xmax() , ylow   , yhigh  ) ; }
  else if ( ylow  < ymin () ) { return  integral ( xlow   , xhigh  , ymin() , yhigh  ) ; }
  else if ( yhigh > ymax () ) { return  integral ( xlow   , xhigh  , ylow   , ymax() ) ; }
  //
  else if ( s_equal ( xlow , xhigh ) || s_equal ( ylow , yhigh ) ) { return 0 ; }
  //
  std::vector<double> fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i )
  { fx[i] = m_b[i].integral ( xlow , xhigh ) ; }
  //
  std::vector<double> fy ( m_n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_n ; ++i )
  { fy[i] = m_b[i].integral ( ylow , yhigh ) ; }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param x     variable 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateX 
( const double y    ,  
  const double xlow , const double xhigh ) const 
{ return integrateY ( y , xlow , xhigh ) ; }
// ============================================================================
/*  integral over y-dimension 
 *  \f[ \int_{y_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param y     variable 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateY
( const double x    ,
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -integrateY ( x , yhigh , ylow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow >  ymax() ) { return 0 ; }
  else if ( s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () )  ) { return integrateY ( x ) ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_n ; ++i )
  { fy[i] = m_b[i].integral ( y_low , y_high ) ; }
  //
  std::vector<double> fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i )
  { fx[i] = m_b[i] ( x ) ; }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}} 
 *  \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 */
// ============================================================================
double  Ostap::Math::Bernstein2DSym::integral   () const 
{
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_n ; ++ix )
  {
    result   += par ( ix , ix ) ;
    for  ( unsigned short iy = 0 ; iy < ix ; ++iy )
    { result += 2 * par ( ix , iy ) ; } 
  }
  //
  return result ;
}
// ============================================================================
/* integral over x-dimension 
 *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param x     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateX ( const double y ) const 
{ return integrateY ( y ) ; }
// ============================================================================
/*  integral over y-dimension 
 *  \f[ \int_{y_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param y     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateY ( const double x ) const 
{
  //
  if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  //
  std::vector<double>       fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i ) { fx[i] = m_b[i] ( x ) ; }
  const std::vector<double> fy ( m_n + 1 , ( ymax() - ymin () ) / ( m_n + 1 ) ) ;
  //
  return  calculate ( fx , fy ) ;
}
// ============================================================================
// set (k)-parameter
// ============================================================================
bool Ostap::Math::Bernstein2DSym::setPar
( const unsigned int   k     , 
  const double         value )
{
  //
  if ( k >= npars() )                     { return false ; }
  if ( s_equal ( m_pars [ k ] , value ) ) { return false ; }
  m_pars [ k ] = value ;
  //
  return true ;
}
// ============================================================================
Ostap::Math::Bernstein2DSym&
Ostap::Math::Bernstein2DSym::operator+=( const double a )
{
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2DSym&
Ostap::Math::Bernstein2DSym::operator*=( const double a )
{
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2DSym&
Ostap::Math::Bernstein2DSym::operator-=( const double a )
{
  if ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2DSym&
Ostap::Math::Bernstein2DSym::operator/=( const double a )
{
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::operator-() const
{
  Bernstein2DSym b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::__add__   ( const double value ) const
{ return (*this) + value ; }
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::__radd__  ( const double value ) const
{ return value + (*this) ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::__mul__   ( const double value ) const
{ return (*this) * value ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::__rmul__  ( const double value ) const
{ return value * (*this) ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::__sub__  ( const double value ) const
{ return (*this) - value ; }
// ============================================================================
// Subtract Bernstein polynomial from a constant
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::__rsub__ ( const double value ) const
{ return value - (*this) ; }
// ============================================================================
// Divide Benrstein polynomial by a constant
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym:: __div__   ( const double value ) const
{ return (*this) / value ; }
// ============================================================================
// Negate Bernstein polynomial
// ============================================================================
Ostap::Math::Bernstein2DSym
Ostap::Math::Bernstein2DSym::__neg__ ()  const
{ return -(*this); }
// ============================================================================


// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive2D::Positive2D
( const unsigned short      nX   ,
  const unsigned short      nY   ,
  const double              xmin ,
  const double              xmax ,
  const double              ymin ,
  const double              ymax )
  : m_bernstein (   nX , nY , xmin , xmax , ymin , ymax ) 
  , m_sphere    ( ( nX + 1 ) * ( nY + 1 ) - 1 )
{
  updateBernstein () ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::Positive2D::Positive2D
(       Ostap::Math::Positive2D&& right ) 
  : m_bernstein ( std::move ( right.m_bernstein ) ) 
  , m_sphere    ( std::move ( right.m_sphere    ) ) 
{}
// ============================================================================
// swap 
// ============================================================================
void Ostap::Math::Positive2D::swap ( Ostap::Math::Positive2D& right ) 
{
  Ostap::Math::swap ( m_bernstein ,  right.m_bernstein ) ;
  Ostap::Math::swap ( m_sphere    ,  right.m_sphere    ) ;  
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive2D::setPar 
( const unsigned int k     , 
  const double       value )
{
  //
  const bool update = m_sphere.setPhase ( k , value ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive2D::updateBernstein ()
{
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , m_sphere.x2 ( ix ) ) ;
    update = updated || update ;  
  }
  //
  return update ;
}
// ============================================================================
// get the parameter value
// ============================================================================
double Ostap::Math::Positive2D::par ( const unsigned int k ) const 
{ return m_sphere.phase ( k ) ; }
// ============================================================================
/*  get the integral over 2D-region           
 *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}} 
 *        \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 */
// ============================================================================
double  Ostap::Math::Positive2D::integral   () const { return 1 ; }
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Positive2D::integral   
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{ 
  return
    s_equal ( xlow  , xmin() ) && 
    s_equal ( xhigh , xmax() ) && 
    s_equal ( ylow  , ymin() ) && 
    s_equal ( yhigh , ymax() )  ?  1.0 :
    m_bernstein.integral ( xlow , xhigh , ylow , yhigh ) ; 
}
// ============================================================================




// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive2DSym::Positive2DSym
( const unsigned short      N    ,
  const double              xmin ,
  const double              xmax )
  : m_bernstein (   N , xmin , xmax ) 
  , m_sphere    ( ( N + 1 ) * ( N + 2 ) / 2 - 1  )
{
  updateBernstein () ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::Positive2DSym::Positive2DSym
(       Ostap::Math::Positive2DSym&& right ) 
  : m_bernstein ( std::move ( right.m_bernstein ) )
  , m_sphere    ( std::move ( right.m_sphere    ) ) 
{}
// ============================================================================
// swap 
// ============================================================================
void Ostap::Math::Positive2DSym::swap ( Ostap::Math::Positive2DSym& right ) 
{
  Ostap::Math::swap ( m_bernstein ,  right.m_bernstein ) ;
  Ostap::Math::swap ( m_sphere    ,  right.m_sphere    ) ;  
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive2DSym::setPar 
( const unsigned int k     , 
  const double       value )
{
  //
  const bool update = m_sphere.setPhase ( k , value ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive2DSym::updateBernstein ()
{
  //
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , m_sphere.x2 ( ix ) ) ; 
    update = updated || update ; 
  }
  //
  return update ;
}
// ============================================================================
// get the parameter value
// ============================================================================
double Ostap::Math::Positive2DSym::par ( const unsigned int  k ) const 
{ return m_sphere.phase ( k ) ; }
// ============================================================================
/*  get the integral over 2D-region 
 *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}} 
 *   \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y \f] 
 */
// ============================================================================
double  Ostap::Math::Positive2DSym::integral   () const { return 1 ; }
// ============================================================================
/*  get the integral over 2D-region 
 *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} 
 *   \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y \f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Positive2DSym::integral
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const
{
  return 
    s_equal ( xlow  , xmin () ) &&
    s_equal ( xhigh , xmax () ) &&
    s_equal ( ylow  , ymin () ) &&
    s_equal ( yhigh , ymax () ) ?  1.0 :
    m_bernstein.integral ( xlow , xhigh , ylow , yhigh ) ; 
}
// ============================================================================
/* integral over x-dimension 
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param x     variable 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ======================================================================
double  Ostap::Math::Positive2DSym::integrateX 
( const double y    , 
  const double xlow , const double xhigh ) const 
{ return m_bernstein.integrateX ( y , xlow , xhigh ) ; }
// ======================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param y     variable 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 */
// ======================================================================
double Ostap::Math::Positive2DSym::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{ return m_bernstein.integrateY ( x , ylow , yhigh ) ; }
// ======================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param x     variable 
 */
// ======================================================================
double Ostap::Math::Positive2DSym::integrateX ( const double y ) const 
{ return m_bernstein.integrateX ( y ) ; }
// ======================================================================
/*  integral over y-dimension 
 *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param y     variable 
 */
// ======================================================================
double Ostap::Math::Positive2DSym::integrateY ( const double x ) const 
{ return m_bernstein.integrateY ( x ) ; }
// ======================================================================


// ============================================================================
// The END 
// ============================================================================
