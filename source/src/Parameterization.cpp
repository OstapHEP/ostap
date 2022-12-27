// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Polynomials.h"
#include "Ostap/Parameterization.h"
// ============================================================================
// Local
// ============================================================================
#include "local_math.h"
// ============================================================================
/** @file 
 *  Implemenation of functions and classes
 *  from the file Ostap/Parameterization.h
 *  @author Vanya Belyaev  Ivan.Belyaev@itep.ru
 *  @date 2019-06-30
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  //  precompute values of Legendre polynomials P_i(x) 
  inline void _legendre_values
  ( std::vector<double>& values , const long double x ) 
  { Ostap::Math::legendre_values ( values.begin() , values.end () , x ) ; }
  // ==========================================================================
  //  precompute values of Legendre integrals P_i(x) 
  void _legendre_integrals
  ( std::vector<double>& values , 
    const long double    xlow   , 
    const long double    xhigh  ) 
  {
    const unsigned long N = values.size() ;
    if ( 0 == N ) { return ; }
    values [0] =         xhigh - xlow ;
    if ( 1 == N ) { return ; }
    values [1] = 0.5 * ( xhigh - xlow ) * ( xhigh + xlow ) ;
    if ( 2 == N ) { return ; }
    //
    long double p0_h = 1     ;
    long double p1_h = xhigh ;
    long double pi_h = 0     ;
    //
    long double p0_l = 1     ;
    long double p1_l = xlow  ;
    long double pi_l = 0     ;
    //
    for ( unsigned short i = 2 ; i < N ; ++i )
    {
      pi_h = ( ( 2 * i - 1 ) * xhigh * p1_h - ( i - 1 ) * p0_h ) / i ;
      pi_l = ( ( 2 * i - 1 ) * xlow  * p1_l - ( i - 1 ) * p0_l ) / i ;
      //
      const long double ii_h =  xhigh * pi_h - p1_h ;
      const long double ii_l =  xlow  * pi_l - p1_l ;
      //
      values [ i ] = ( ii_h - ii_l ) / ( i + 1 ) ;
      //
      p0_h = p1_h ;
      p1_h = pi_h ;
      p0_l = p1_l ;
      p1_l = pi_l ; 
    }
  }
  // ==========================================================================
}
// ============================================================================
// Negation operators 
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum2::operator-() const
{ LegendreSum2 s(*this) ; Ostap::Math::negate ( s.m_pars ) ; return s ; }
// ============================================================================
Ostap::Math::LegendreSum3 
Ostap::Math::LegendreSum3::operator-() const
{ LegendreSum3 s(*this) ; Ostap::Math::negate ( s.m_pars ) ; return s ; }
// ============================================================================
Ostap::Math::LegendreSum4 
Ostap::Math::LegendreSum4::operator-() const
{ LegendreSum4 s(*this) ; Ostap::Math::negate ( s.m_pars ) ; return s ; }
// ============================================================================


// ============================================================================
// 2- Legendre parameterization 
// ============================================================================
// constructor 
// ============================================================================
Ostap::Math::LegendreSum2::LegendreSum2 
( const unsigned short NX , 
  const unsigned short NY , 
  const double   xmin     , 
  const double   xmax     , 
  const double   ymin     , 
  const double   ymax     ) 
  : Parameters ( ( NX + 1 ) * ( NY + 1 ) ) 
    //
  , m_NX   ( NX ) 
  , m_NY   ( NY ) 
    //
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
    //
  , m_ymin ( std::min ( ymin , ymax ) )
  , m_ymax ( std::max ( ymin , ymax ) )
    //
  , m_cache_x ( ( m_NX + 1 ) , 0.0 )
  , m_cache_y ( ( m_NY + 1 ) , 0.0 )
{}
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::LegendreSum2::LegendreSum2 
( const std::vector<double>& pars ,
  const unsigned short       NX   , 
  const unsigned short       NY   , 
  const double               xmin , 
  const double               xmax , 
  const double               ymin , 
  const double               ymax )
  : Ostap::Math::LegendreSum2 ( NX , NY , xmin , xmax , ymin , ymax )
{
  setPars ( pars  ) ;
}
// ============================================================================
/* constructor orm the product of two Legendre sums
 *  \f$ S(x,y) = S_x(x)\times S_y(y) \f$ 
 *  @param sx (INPUT) the first  Legendre sum 
 *  @param sy (INPUT) the second Legendre sum 
 */
// ============================================================================
Ostap::Math::LegendreSum2::LegendreSum2
( const Ostap::Math::LegendreSum& sx , 
  const Ostap::Math::LegendreSum& sy )
  : Parameters ( ( sx.degree() + 1 ) * ( sy.degree() + 1 ) ) 
    //
  , m_NX   ( sx.degree() ) 
  , m_NY   ( sy.degree() ) 
    //
  , m_xmin ( sx.xmin () )
  , m_xmax ( sx.xmax () )
    //
  , m_ymin ( sy.xmin () )
  , m_ymax ( sy.xmax () )
    //
  , m_cache_x ( ( m_NX + 1 ) , 0.0 )
  , m_cache_y ( ( m_NY + 1 ) , 0.0 )
{
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  {
    for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    {
      const unsigned int k =  ix * ( m_NY + 1 ) + iy ;
      m_pars [ k ] = sx.par ( ix ) * sy.par ( iy ) ;
    }
  }
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::LegendreSum2::evaluate
( const double x , 
  const double y ) const 
{
  //
  long double xx     = tx ( x ) ;
  long double yy     = ty ( y ) ;
  //
  _legendre_values ( m_cache_x , xx ) ;
  _legendre_values ( m_cache_y , yy ) ;
  //
  return calculate () ;
}
// ============================================================================
/*  update  the Legendre expansion by addition of one "event" with 
 *  the given weight
 *  @code
 *  LegendreSum2 sum = ... ;
 *  for ( auto x : .... ) { sum.fill( x , y ) ; }
 *  @endcode
 */
// ============================================================================
bool Ostap::Math::LegendreSum2::fill 
( const double x      , 
  const double y      , 
  const double weight ) 
{
  // no update 
  if      ( x < m_xmin || x > m_xmax ) { return false ; }
  else if ( y < m_ymin || y > m_ymax ) { return false ; }
  else if ( s_zero ( weight )        ) { return true  ; }
  //
  const long double w   = weight * 4.0L / 
    ( ( m_ymax - m_ymin ) * ( m_xmax - m_xmin ) ) ;
  //
  const double xx  =  tx ( x ) ;
  const double yy  =  ty ( y ) ;
  //
  _legendre_values ( m_cache_x , xx ) ;
  _legendre_values ( m_cache_y , yy ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { const unsigned int k = index ( ix , iy ) ;
      m_pars[k] += w * m_cache_x[ix] * m_cache_y[iy] * ( ix + 0.5L ) * ( iy + 0.5L ) ; } }
  //
  return true ;
} 
// ============================================================================
// Integrals and projections 
// ============================================================================
/*  get the integral 
 *  \f$ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}} f(x,y){\mathrm{d}} x {\mathrm{d}} y \f$
 */
// ============================================================================
double Ostap::Math::LegendreSum2::integral
( const double xlow  ,   
  const double xhigh , 
  const double ylow  ,   
  const double yhigh ) const 
{
  //
  if ( s_equal ( xlow , xhigh ) || s_equal ( ylow , yhigh ) ) { return 0 ; }
  else if ( xlow > xhigh ) { return -integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow > yhigh ) { return -integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  // 
  const double xl = std::max ( xlow  , m_xmin ) ;
  const double xh = std::min ( xhigh , m_xmax ) ;
  const double yl = std::max ( ylow  , m_ymin ) ;
  const double yh = std::min ( yhigh , m_ymax ) ;
  //
  if      ( xh <= m_xmin || xl >= m_xmax ) { return 0 ; }
  else if ( yh <= m_ymin || yl >= m_ymax ) { return 0 ; }
  //
  if  ( s_equal  ( xl , m_xmin ) && s_equal  ( xh , m_xmax ) && 
        s_equal  ( yl , m_ymin ) && s_equal  ( yh , m_ymax ) ) { return integral() ; }
  //
  const double txl =  tx ( xl ) ;
  const double txh =  tx ( xh ) ;
  //
  const double tyl =  ty ( yl ) ;
  const double tyh =  ty ( yh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_x , txl , txh ) ;
  _legendre_integrals ( m_cache_y , tyl , tyh ) ;
  //
  return calculate()  * ( m_xmax - m_xmin ) * ( m_ymax - m_ymin ) * 1./4 ;
}
// ============================================================================
/*  get the integral 
 *  \f$ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f(x,y){\mathrm{d}} x {\mathrm{d}} y \f$
 */
// ============================================================================
double Ostap::Math::LegendreSum2::integral   () const 
{ return m_pars[0] * ( m_xmax - m_xmin ) * ( m_ymax - m_ymin ); }
// ============================================================================
/*  integrate over x dimension 
 *  \f$ f(y) =  \int_{x_{min}}^{x_{max}} F(x,y) {\mathrm{d}} x \f$
 */
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum2::integralX  () const 
{
  // vector of coefficients 
  std::vector<double> pars ( m_NY + 1 , 0.0 ) ;
  //
  // get the integral 
  const unsigned short ix = 0 ;
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
  { pars [ iy ] = m_pars [ index ( ix , iy )  ] ; }
  //
  Ostap::Math::scale ( pars , m_xmax - m_xmin ) ;
  return LegendreSum ( pars , m_ymin , m_ymax ) ;
}
// ============================================================================
/*  integrate over y dimension 
 *  \f$ f(x) =  \int_{y_{min}}^{y_{max}} F(x,y) {\mathrm{d}} y \f$
 */
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum2::integralY  () const 
{
  // vector of coefficients 
  std::vector<double> pars ( m_NX + 1 , 0.0 ) ;
  //
  // get the integral 
  const unsigned short iy = 0 ;
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { pars [ ix ] = m_pars [ index ( ix , iy )  ] ; }
  //
  Ostap::Math::scale ( pars , m_ymax - m_ymin ) ;
  return LegendreSum ( pars , m_xmin , m_xmax ) ;
}
// ============================================================================
/* integrate over x dimension 
 *  \f$ f_x(y) =  \int_{x_{low}}^{x_{high}} f(x,y) {\mathrm{d}} x \f$
 */
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum2::integralX
( const double xlow  , 
  const double xhigh ) const 
{
  //
  if ( s_equal ( xlow , xhigh   ) ) { return LegendreSum ( 0 , m_ymin  , m_ymax ) ; }
  if (           xlow > xhigh     ) { return -integralX  ( xhigh , xlow ) ; }
  //
  const double xl  = std::max ( xlow  , m_xmin ) ;
  const double xh  = std::min ( xhigh , m_xmax  );
  if  ( s_equal ( xl , m_xmin ) && s_equal ( xh , m_xmax ) ) { return integralX () ; }
  //  
  const double txl = tx ( xl ) ;
  const double txh = tx ( xh ) ;
  //
  // prepare cache 
  _legendre_integrals ( m_cache_x , txl ,  txh ) ;
  //
  // vector of coefficients 
  std::vector<double> pars ( m_NY + 1  , 0.0 ) ;
  //
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
  { for ( unsigned short ix = 0 ;  ix <= m_NX ; ++ix ) 
  { pars[iy] += 0.5 * m_pars [ index ( ix , iy )  ] * m_cache_x[ix] ; } }
  //
  Ostap::Math::scale ( pars , m_xmax - m_xmin ) ;
  return LegendreSum ( pars , m_ymin , m_ymax ) ; 
}
// ============================================================================
/* integrate over y dimension 
 *  \f$ f(x) =  \int_{y_{low}}^{y_{high}} f(x,y) {\mathrm{d}} y \f$
 */
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum2::integralY
( const double ylow  , 
  const double yhigh ) const 
{
  //
  if ( s_equal ( ylow , yhigh   ) ) { return LegendreSum ( 0 , m_xmin  , m_xmax ) ; }
  if (           ylow > yhigh     ) { return -integralY  ( yhigh , ylow ) ; }
  //
  const double yl = std::max ( ylow  , m_ymin ) ;
  const double yh = std::min ( yhigh , m_ymax  );
  if  ( s_equal ( yl , m_ymin ) && s_equal ( yh , m_ymax ) ) { return integralY () ; }
  //  
  const double tyl = ty ( yl ) ;
  const double tyh = ty ( yh ) ;
  //
  // prepare cache 
  _legendre_integrals ( m_cache_y , tyl ,  tyh ) ;
  //
  // vector of coefficients 
  std::vector<double> pars ( m_NX + 1  , 0.0 ) ;
  //
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
  { for ( unsigned short ix = 0 ;  ix <= m_NX ; ++ix ) 
      { pars[ix] += 0.5 * m_pars [ index ( ix , iy )  ] * m_cache_y[iy] ; } } // <<---- !!!
  //
  Ostap::Math::scale ( pars , m_ymax - m_ymin ) ;
  return LegendreSum ( pars , m_xmin , m_xmax ) ; 
}
// ============================================================================
// "transpose": \f$ S(x,y) \equiv F(y,x) \f$
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum2::transpose() const 
{
  LegendreSum2 t ( m_NY , m_NX  , m_ymin , m_ymax , m_xmin , m_xmax ) ;
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
  { for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
    { t.m_pars [ t.index ( iy , ix ) ] = par ( ix , iy ) ; } }
  //
  return t ; 
}
// ============================================================================
// 3D
// ============================================================================
Ostap::Math::LegendreSum3::LegendreSum3 
( const unsigned short NX , 
  const unsigned short NY , 
  const unsigned short NZ , 
  const double   xmin     , 
  const double   xmax     , 
  const double   ymin     , 
  const double   ymax     ,
  const double   zmin     , 
  const double   zmax     ) 
  : Parameters ( ( NX + 1 ) * (  NY + 1 ) * ( NZ + 1 ) )
    //
  , m_NX   ( NX ) 
  , m_NY   ( NY ) 
  , m_NZ   ( NZ ) 
    //
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
    //
  , m_ymin ( std::min ( ymin , ymax ) )
  , m_ymax ( std::max ( ymin , ymax ) )
    //
  , m_zmin ( std::min ( zmin , zmax ) )
  , m_zmax ( std::max ( zmin , zmax ) )
    //
  , m_cache_x ( NX + 1 , 0.0 ) 
  , m_cache_y ( NY + 1 , 0.0 ) 
  , m_cache_z ( NZ + 1 , 0.0 )
{}
// ============================================================================
Ostap::Math::LegendreSum3::LegendreSum3 
( const std::vector<double>& pars ,
  const unsigned short       NX   , 
  const unsigned short       NY   , 
  const unsigned short       NZ   , 
  const double               xmin , 
  const double               xmax , 
  const double               ymin , 
  const double               ymax ,
  const double               zmin , 
  const double               zmax )
  : Ostap::Math::LegendreSum3 ( NX   , NY   , NZ ,
                                xmin , xmax ,
                                ymin , ymax ,
                                zmin , zmax )
{
  setPars ( pars ) ;
  
}
// ============================================================================
/*  constructor orm the product of two Legendre sums
 *  \f$ S(x,y,z) = S_x(x)\times S_y(y) \times S_z(z) \f$ 
 *  @param sx (INPUT) the first  Legendre sum 
 *  @param sy (INPUT) the second Legendre sum 
 *  @param sz (INPUT) the third Legendre sum 
 */
// ============================================================================
Ostap::Math::LegendreSum3::LegendreSum3
( const Ostap::Math::LegendreSum& sx , 
  const Ostap::Math::LegendreSum& sy ,
  const Ostap::Math::LegendreSum& sz )
  : Parameters ( ( sx.degree() + 1 ) * ( sy.degree() + 1 ) * ( sz.degree() + 1 ) )
    //
  , m_NX   ( sx.degree () ) 
  , m_NY   ( sy.degree () ) 
  , m_NZ   ( sz.degree () ) 
    //
  , m_xmin ( sx.xmin   () )
  , m_xmax ( sx.xmax   () )
    //
  , m_ymin ( sy.xmin   () )
  , m_ymax ( sy.xmax   () )
    //
  , m_zmin ( sz.xmin   () )
  , m_zmax ( sz.xmax   () )
    //
  , m_cache_x ( ( m_NX + 1 ) , 0.0 )
  , m_cache_y ( ( m_NY + 1 ) , 0.0 )
  , m_cache_z ( ( m_NZ + 1 ) , 0.0 )
{
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { m_pars [ index ( ix , iy , iz ) ] = sx.par ( ix ) * sy.par ( iy ) *  sz.par ( iz ) ; } } }
}
// ============================================================================
/*  constructor form the product of two Legendre sums
 *  \f$ S(x,y,z) = S_{xy}(x,y) \times S_z(z) \f$ 
 *  @param sxy (INPUT) the first  Legendre sum 
 *  @param sz  (INPUT) the second Legendre sum 
 */
// ============================================================================
Ostap::Math::LegendreSum3::LegendreSum3
( const Ostap::Math::LegendreSum2& sxy ,
  const Ostap::Math::LegendreSum&  sz  ) 
  : Parameters ( ( sxy.nx() + 1 ) * ( sxy.ny() + 1 ) * ( sz.degree() + 1 ) )
    //
  , m_NX   ( sxy.nx ()     ) 
  , m_NY   ( sxy.ny ()     ) 
  , m_NZ   ( sz .degree () ) 
    //
  , m_xmin ( sxy.xmin() )
  , m_xmax ( sxy.xmax() )
    //
  , m_ymin ( sxy.ymin() )
  , m_ymax ( sxy.ymax() )
    //
  , m_zmin ( sz .xmin() )
  , m_zmax ( sz .xmax() )
    //
  , m_cache_x ( ( m_NX + 1 ) , 0.0 )
  , m_cache_y ( ( m_NY + 1 ) , 0.0 )
  , m_cache_z ( ( m_NZ + 1 ) , 0.0 )
{
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { m_pars [ index ( ix , iy , iz ) ] = sxy.par ( ix , iy ) * sz.par ( iz ) ; } } }
}
// ============================================================================
/*  constructor form the product of two Legendre sums
 *  \f$ S(x,y,z) = S_x(x) \times S_{yz}(y,z)\f$ 
 *  @param sx   (INPUT) the first  Legendre sum 
 *  @param syz  (INPUT) the second Legendre sum 
 */
// ============================================================================
Ostap::Math::LegendreSum3::LegendreSum3
( const Ostap::Math::LegendreSum&  sx  ,
  const Ostap::Math::LegendreSum2& syz ) 
  : Parameters ( ( sx.degree() + 1 ) * ( syz.nx() + 1 ) * ( syz.ny() + 1 ) )
    //
  , m_NX   ( sx  .degree () )
  , m_NY   ( syz .nx     () )
  , m_NZ   ( syz .ny     () )
    //
  , m_xmin ( sx  .xmin   () )
  , m_xmax ( sx  .xmax   () )
    //
  , m_ymin ( syz .xmin   () )
  , m_ymax ( syz .xmax   () )
    //
  , m_zmin ( syz .ymin   () )
  , m_zmax ( syz .ymax   () )
    //
  , m_cache_x ( ( m_NX + 1 ) , 0.0 )
  , m_cache_y ( ( m_NY + 1 ) , 0.0 )
  , m_cache_z ( ( m_NZ + 1 ) , 0.0 )
{
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { m_pars [ index ( ix , iy , iz ) ] = sx.par ( ix ) * syz.par ( iy , iz ) ; } } }
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::LegendreSum3::evaluate    
( const double x , 
  const double y , 
  const double z ) const 
{
  //
  long double xx     = tx ( x ) ;
  long double yy     = ty ( y ) ;
  long double zz     = tz ( z ) ;
  //
  _legendre_values ( m_cache_x , xx ) ;
  _legendre_values ( m_cache_y , yy ) ;
  _legendre_values ( m_cache_z , zz ) ;
  //
  return calculate () ;
}
// ============================================================================
/*  update  the Legendre expansion by addition of one "event" with 
 *  the given weight
 *  @code
 *  LegendreSum3 sum = ... ;
 *  for ( auto x : .... ) { sum.fill ( x , y , z ) ; }
 *  @endcode
 *  This is a useful function to make an unbinned parameterization 
 *  of certain distribution and/or efficiency 
 */
// ============================================================================
bool Ostap::Math::LegendreSum3::fill
( const double x , 
  const double y , 
  const double z , 
  const double weight ) 
{
  // no update 
  if      ( x < m_xmin || x > m_xmax ) { return false ; }
  else if ( y < m_ymin || y > m_ymax ) { return false ; }
  else if ( z < m_zmin || z > m_zmax ) { return false ; }
  else if ( s_zero ( weight )        ) { return true  ; }
  //
  const long double w = weight * 8.0L / 
    ( ( m_zmax - m_zmin ) * ( m_ymax - m_ymin ) * ( m_xmax - m_xmin ) ) ;
  //
  const double xx  =  tx ( x ) ;
  const double yy  =  ty ( y ) ;
  const double zz  =  tz ( z ) ;
  //
  _legendre_values ( m_cache_x , xx ) ;
  _legendre_values ( m_cache_y , yy ) ;
  _legendre_values ( m_cache_z , zz ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { const unsigned int k = index ( ix , iy , iz ) ;
        m_pars[k] += w * m_cache_x[ix] * m_cache_y[iy] * m_cache_z[iz] 
          * ( ix + 0.5L ) * ( iy + 0.5L ) * ( iz + 0.5L ) ; } }}
  //
  return true ;
}
// ============================================================================
/*  integrate over x dimension 
 *  \f$ f_x(y,z) =  \int_{x_{min}}^{x_{max}} f(x,y,z) {\mathrm{d}} x \f$
 */
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum3::integralX () const 
{
  LegendreSum2 r ( m_NY , m_NZ , m_ymin , m_ymax , m_zmin , m_zmax );
  //
  const unsigned short ix = 0 ;
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { r.setPar ( iy , iz , m_pars [ index ( ix, iy, iz ) ] ) ; } }
  //
  r *= ( m_xmax - m_xmin ) ;
  //
  return r ;
}
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum3::integralY () const 
{
  LegendreSum2 r ( m_NX , m_NZ , m_xmin , m_xmax , m_zmin , m_zmax );
  //
  const unsigned short iy = 0 ;
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { r.setPar ( ix , iz , m_pars [ index ( ix, iy, iz ) ] ) ; } }
  //
  r *= ( m_ymax - m_ymin ) ;
  //
  return r ;
}
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum3::integralZ () const 
{
  LegendreSum2 r ( m_NX , m_NY , m_xmin , m_xmax , m_ymin , m_ymax );
  //
  const unsigned short iz = 0 ;
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { r.setPar ( ix , iy , m_pars [ index ( ix, iy, iz ) ] ) ; } }
  //
  r *= ( m_zmax - m_zmin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over x dimension 
 *  \f$ f(y,z) =  \int_{x_{low}}^{x_{high}} F(x,y,z) {\mathrm{d}} x \f$
 */
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum3::integralX 
( const double xlow  , 
  const double xhigh ) const 
{
  LegendreSum2 r ( m_NY , m_NZ , m_ymin , m_ymax , m_zmin , m_zmax );
  //
  if ( s_equal (  xlow , xhigh ) ) { return  r ; }
  if (            xlow > xhigh   ) { return -integralX ( xhigh , xlow ) ; }
  //
  const  double xl = std::max ( xlow  , m_xmin ) ;
  const  double xh = std::min ( xhigh , m_xmax ) ;
  //
  if  ( xh <= m_xmin || xl >= m_xmax ) { return r ; }
  //
  if  ( s_equal ( xl , m_xmin ) && s_equal ( xh , m_xmax ) ) { return integralX() ; }
  //
  const double txl =  tx ( xl ) ;
  const double txh =  tx ( xh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_x , txl , txh ) ;
  //
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy )
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { 
      double value = 0 ;
      for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
      { value += 0.5 * m_pars [ index ( ix, iy, iz ) ] * m_cache_x [ ix ] ; }
      r.setPar ( iy , iz , value ) ; 
    }
  }
  //
  r *= ( m_xmax - m_xmin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over y dimension 
 *  \f$ f(x,z) =  \int_{y_{low}}^{y_{high}} F(x,y,z) {\mathrm{d}} y \f$
 */
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum3::integralY 
( const double ylow  , 
  const double yhigh ) const 
{
  LegendreSum2 r ( m_NX , m_NZ , m_xmin , m_xmax , m_zmin , m_zmax );
  //
  if ( s_equal ( ylow , yhigh ) ) { return  r ; }
  if (           ylow > yhigh   ) { return -integralY ( yhigh , ylow ) ; }
  //
  const  double yl = std::max ( ylow  , m_ymin ) ;
  const  double yh = std::min ( yhigh , m_ymax ) ;
  //
  if  ( yh <= m_ymin || yl >= m_ymax ) { return r  ; }
  //
  if  ( s_equal ( yl , m_ymin ) && s_equal ( yh , m_ymax ) ) { return integralY () ; }
  //
  const double tyl =  ty ( yl ) ;
  const double tyh =  ty ( yh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_y , tyl , tyh ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix )
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { 
      double value = 0 ;
      for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
      { value += 0.5 * m_pars [ index ( ix, iy, iz ) ] * m_cache_y [ iy ] ; }
      r.setPar ( ix , iz , value ) ; 
    }
  }
  //
  r *= ( m_ymax - m_ymin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over z dimension 
 *  \f$ f(x,y) =  \int_{z_{low}}^{z_{high}} f(x,y,z) {\mathrm{d}} z \f$
 */
// ============================================================================
Ostap::Math::LegendreSum2 
Ostap::Math::LegendreSum3::integralZ 
( const double zlow  , 
  const double zhigh ) const 
{
  LegendreSum2 r ( m_NX , m_NY , m_xmin , m_xmax , m_ymin , m_ymax );
  //
  if ( s_equal ( zlow , zhigh ) ) { return  r ; }
  if (           zlow > zhigh   ) { return -integralZ ( zhigh , zlow ) ; }
  //
  const  double zl = std::max ( zlow  , m_zmin ) ;
  const  double zh = std::min ( zhigh , m_zmax ) ;
  //
  if  ( zh <= m_zmin || zl >= m_zmax ) { return r  ; }
  //
  if  ( s_equal ( zl , m_zmin ) && s_equal ( zh , m_zmax ) ) { return integralZ () ; }
  //
  const double tzl =  tz ( zl ) ;
  const double tzh =  tz ( zh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_z , tzl , tzh ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix )
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { 
      double value = 0 ;
      for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { value += 0.5 * m_pars [ index ( ix, iy, iz ) ] * m_cache_z [ iz ] ; }
      r.setPar ( ix , iy , value ) ; 
    }
  }
  //
  r *= ( m_zmax - m_zmin ) ;
  //
  return r ;
}
// ============================================================================
/*  integral 
 *  \f$ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}}
 *    \int_{z_{low}}^{z_{high}} f(x,y,z) {\mathrm{d}} x {\mathrm{d}} y {\mathrm{d}} z \f$ 
 */
// ============================================================================
double Ostap::Math::LegendreSum3::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh , 
  const double zlow , const double zhigh ) const 
{
  if ( s_equal ( xlow , xhigh ) || 
       s_equal ( ylow , yhigh ) ||  
       s_equal ( zlow , zhigh ) ) { return 0 ; }
  else if ( xlow > xhigh ) { return -integral ( xhigh , xlow  , ylow  , yhigh , zlow  , zhigh ) ; }
  else if ( ylow > yhigh ) { return -integral ( xlow  , xhigh , yhigh , ylow  , zlow  , zhigh ) ; }
  else if ( zlow > zhigh ) { return -integral ( xlow  , xhigh , ylow  , yhigh , zhigh , zlow  ) ; }
  // 
  const double xl = std::max ( xlow  , m_xmin ) ;
  const double xh = std::min ( xhigh , m_xmax ) ;
  const double yl = std::max ( ylow  , m_ymin ) ;
  const double yh = std::min ( yhigh , m_ymax ) ;
  const double zl = std::max ( zlow  , m_zmin ) ;
  const double zh = std::min ( zhigh , m_zmax ) ;
  //
  if      ( xh <= m_xmin || xl >= m_xmax ) { return 0 ; }
  else if ( yh <= m_ymin || yh >= m_ymax ) { return 0 ; }
  else if ( zh <= m_zmin || zh >= m_zmax ) { return 0 ; }
  //
  if  ( s_equal  ( xl , m_xmin ) && s_equal  ( xh , m_xmax ) && 
        s_equal  ( yl , m_ymin ) && s_equal  ( yh , m_ymax ) && 
        s_equal  ( zl , m_zmin ) && s_equal  ( zh , m_zmax ) ) { return integral() ; }
  //
  const double txl =  tx ( xl ) ;
  const double txh =  tx ( xh ) ;
  //
  const double tyl =  ty ( yl ) ;
  const double tyh =  ty ( yh ) ;
  //
  const double tzl =  tz ( zl ) ;
  const double tzh =  tz ( zh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_x , txl , txh ) ;
  _legendre_integrals ( m_cache_y , tyl , tyh ) ;
  _legendre_integrals ( m_cache_z , tzl , tzh ) ;
  //
  return calculate()  * ( m_xmax - m_xmin ) * ( m_ymax - m_ymin ) * ( m_zmax - m_zmin ) * 1./8 ;
}
// ============================================================================
/*  integral 
 *  \f$ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}
 *    \int_{z_{min}}^{z_{max}} f(x,y,z) {\mathrm{d}} x {\mathrm{d}} y {\mathrm{d}} z \f$ 
 */
// ============================================================================
double Ostap::Math::LegendreSum3::integral () const 
{ return m_pars[0] * ( m_xmax - m_xmin ) * ( m_ymax - m_ymin ) * ( m_zmax - m_zmin ); }

// ============================================================================
// 4D
// ============================================================================
Ostap::Math::LegendreSum4::LegendreSum4 
( const unsigned short NX , 
  const unsigned short NY , 
  const unsigned short NZ , 
  const unsigned short NU , 
  const double   xmin     , 
  const double   xmax     , 
  const double   ymin     , 
  const double   ymax     ,
  const double   zmin     , 
  const double   zmax     ,
  const double   umin     , 
  const double   umax     ) 
  : Parameters ( ( NX + 1 ) * (  NY + 1 ) * ( NZ + 1 ) * ( NU + 1 ) ) 
    //
  , m_NX   ( NX ) 
  , m_NY   ( NY ) 
  , m_NZ   ( NZ ) 
  , m_NU   ( NU ) 
    //
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
    //
  , m_ymin ( std::min ( ymin , ymax ) )
  , m_ymax ( std::max ( ymin , ymax ) )
    //
  , m_zmin ( std::min ( zmin , zmax ) )
  , m_zmax ( std::max ( zmin , zmax ) )
    //
  , m_umin ( std::min ( umin , umax ) )
  , m_umax ( std::max ( umin , umax ) )
    //
  , m_cache_x ( m_NX + 1 , 0.0 ) 
  , m_cache_y ( m_NY + 1 , 0.0 ) 
  , m_cache_z ( m_NZ + 1 , 0.0 )
  , m_cache_u ( m_NU + 1 , 0.0 )
    //
{}
// ============================================================================
Ostap::Math::LegendreSum4::LegendreSum4 
( const std::vector<double>& pars ,
  const unsigned short       NX   , 
  const unsigned short       NY   , 
  const unsigned short       NZ   , 
  const unsigned short       NU   , 
  const double               xmin , 
  const double               xmax , 
  const double               ymin , 
  const double               ymax ,
  const double               zmin , 
  const double               zmax ,
  const double               umin , 
  const double               umax ) 
  : Ostap::Math::LegendreSum4 ( NX   , NY   ,
                                NZ   , NU   ,
                                xmin , xmax ,
                                ymin , ymax ,
                                zmin , zmax ,
                                umin , umax )
{
  setPars  ( pars ) ;
}
// ============================================================================
/*  constructor orm the product of two Legendre sums
 *  \f$ S(x,y,z) = S_x(x)\times S_y(y) \times S_z(z) \f$ 
 *  @param sx (INPUT) the first  Legendre sum 
 *  @param sy (INPUT) the second Legendre sum 
 *  @param sz (INPUT) the third  Legendre sum 
 *  @param su (INPUT) the fourth Legendre sum 
 */
// ============================================================================
Ostap::Math::LegendreSum4::LegendreSum4 
( const LegendreSum&  sx , 
  const LegendreSum&  sy ,
  const LegendreSum&  sz ,
  const LegendreSum&  su ) 
  : Parameters ( ( sx.degree() + 1 ) * 
                 ( sy.degree() + 1 ) * 
                 ( sz.degree() + 1 ) * 
                 ( su.degree() + 1 ) )
    //
  , m_NX   ( sx.degree () ) 
  , m_NY   ( sy.degree () ) 
  , m_NZ   ( sz.degree () ) 
  , m_NU   ( su.degree () ) 
    //
  , m_xmin ( sx.xmin   () )
  , m_xmax ( sx.xmax   () )
    //
  , m_ymin ( sy.xmin   () )
  , m_ymax ( sy.xmax   () )
    //
  , m_zmin ( sz.xmin   () )
  , m_zmax ( sz.xmax   () )
    //
  , m_umin ( su.xmin   () )
  , m_umax ( su.xmax   () )
    //
  , m_cache_x ( ( m_NX + 1 ) , 0.0 )
  , m_cache_y ( ( m_NY + 1 ) , 0.0 )
  , m_cache_z ( ( m_NZ + 1 ) , 0.0 )
  , m_cache_u ( ( m_NU + 1 ) , 0.0 )
{
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
        { m_pars [ index ( ix , iy , iz , iu ) ] = sx.par ( ix ) * sy.par ( iy ) * sz.par ( iz ) * su.par ( iu ) ; } } 
    }
  } 
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::LegendreSum4::evaluate    
( const double x , 
  const double y , 
  const double z , 
  const double u ) const 
{
  //
  const long double xx = tx ( x ) ;
  const long double yy = ty ( y ) ;
  const long double zz = tz ( z ) ;
  const long double uu = tu ( u ) ;
  //
  _legendre_values ( m_cache_x , xx ) ;
  _legendre_values ( m_cache_y , yy ) ;
  _legendre_values ( m_cache_z , zz ) ;
  _legendre_values ( m_cache_u , uu ) ;
  //
  return calculate () ;
}
// ===========================================================================
/** update  the Legendre expansion by addition of one "event" with 
 *  the given weight
 *  @code
 *  LegendreSum3 sum = ... ;
 *  for ( auto x : .... ) { sum.fill ( x , y , z , u ) ; }
 *  @endcode
 *  This is a useful function to make an unbinned parameterization 
 *  of certain distribution and/or efficiency 
 */
// ===========================================================================
bool Ostap::Math::LegendreSum4::fill 
( const double x      , 
  const double y      , 
  const double z      , 
  const double u      , 
  const double weight ) 
{
  // no update 
  if      ( x < m_xmin || x > m_xmax ) { return false ; }
  else if ( y < m_ymin || y > m_ymax ) { return false ; }
  else if ( z < m_zmin || z > m_zmax ) { return false ; }
  else if ( u < m_umin || u > m_umax ) { return false ; }
  else if ( 0 == weight              ) { return true  ; }
  //
  const long double w   = weight * 16.0L / 
    ( ( m_umax - m_umin ) * ( m_zmax - m_zmin ) *
      ( m_ymax - m_ymin ) * ( m_xmax - m_xmin ) ) ;
  //
  const double xx  =  tx ( x ) ;
  const double yy  =  ty ( y ) ;
  const double zz  =  tz ( z ) ;
  const double uu  =  tu ( u ) ;
  //
  _legendre_values ( m_cache_x , xx ) ;
  _legendre_values ( m_cache_y , yy ) ;
  _legendre_values ( m_cache_z , zz ) ;
  _legendre_values ( m_cache_u , uu ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
        { const unsigned int k = index ( ix , iy , iz , iu ) ;
          m_pars[k] += w 
            * m_cache_x[ix] * m_cache_y[iy] 
            * m_cache_z[iz] * m_cache_u[iu] 
            * ( ix + 0.5L ) * ( iy + 0.5L ) 
            * ( iz + 0.5L ) * ( iu + 0.5L ) ; } } } }
  //
  return true ;
}
// ============================================================================
/*  integrate over x dimension 
 *  \f$ f_x(y,z,u) =  \int_{x_{min}}^{x_{max}} f(x,y,z,u) {\mathrm{d}} x \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3  
Ostap::Math::LegendreSum4::integralX  () const 
{
  LegendreSum3 r ( m_NY   , m_NZ   , m_NU , 
                   m_ymin , m_ymax , 
                   m_zmin , m_zmax , 
                   m_umin , m_umax );
  //
  const unsigned short ix = 0 ;
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
      { r.setPar ( iy , iz , iu ,  m_pars [ index ( ix, iy, iz , iu ) ] ) ; } } }
  //
  r *= ( m_xmax - m_xmin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over y dimension 
 *  \f$ f_x(x,z,u) =  \int_{y_{min}}^{y_{max}} f(x,y,z,u) {\mathrm{d}} y \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3  
Ostap::Math::LegendreSum4::integralY  () const 
{
  LegendreSum3 r ( m_NX   , m_NZ   , m_NU , 
                   m_xmin , m_xmax , 
                   m_zmin , m_zmax , 
                   m_umin , m_umax );
  //
  const unsigned short iy = 0 ;
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
      { r.setPar ( ix , iz , iu , m_pars [ index ( ix, iy, iz , iu ) ] ) ; } } }
  //
  r *= ( m_ymax - m_ymin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over z dimension 
 *  \f$ f_x(x,y,u) =  \int_{z_{min}}^{z_{max}} f(x,y,z,u) {\mathrm{d}} z \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3  
Ostap::Math::LegendreSum4::integralZ  () const 
{
  LegendreSum3 r ( m_NX   , m_NY   , m_NU , 
                   m_xmin , m_xmax , 
                   m_ymin , m_ymax , 
                   m_umin , m_umax );
  //
  const unsigned short iz = 0 ;
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
      { r.setPar ( ix , iy , iu , m_pars [ index ( ix , iy, iz , iu ) ] ) ; } } }
  //
  r *= ( m_zmax - m_zmin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over u dimension 
 *  \f$ f_x(x,y,z) =  \int_{u_{min}}^{u_{max}} f(x,y,z,u) {\mathrm{d}} u \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3  
Ostap::Math::LegendreSum4::integralU  () const 
{
  LegendreSum3 r ( m_NX   , m_NY   , m_NZ , 
                   m_xmin , m_xmax , 
                   m_ymin , m_ymax , 
                   m_zmin , m_zmax );
  //
  const unsigned short iu = 0 ;
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { r.setPar ( ix , iy , iz , m_pars [ index ( ix, iy, iz , iu ) ] ) ; } } }
  //
  r *= ( m_umax - m_umin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over x dimension 
 *  \f$ f(y,z,u) =  \int_{x_{min}}^{x_{max}} F(x,y,z,u) {\mathrm{d}} x \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3 
Ostap::Math::LegendreSum4::integralX 
( const double xlow  , 
  const double xhigh ) const 
{
  LegendreSum3 r ( m_NY   , m_NZ   , m_NU , 
                   m_ymin , m_ymax , 
                   m_zmin , m_zmax ,
                   m_umin , m_umax );
  //
  if ( s_equal (  xlow , xhigh ) ) { return  r ; }
  if (            xlow > xhigh   ) { return -integralX ( xhigh , xlow ) ; }
  //
  const  double xl = std::max ( xlow  , m_xmin ) ;
  const  double xh = std::min ( xhigh , m_xmax ) ;
  //
  if  ( xh <= m_xmin || xl >= m_xmax ) { return r ; }
  //
  if  ( s_equal ( xl , m_xmin ) && s_equal ( xh , m_xmax ) ) { return integralX() ; }
  //
  const double txl =  tx ( xl ) ;
  const double txh =  tx ( xh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_x , txl , txh ) ;
  //
  for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy )
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
      { 
        double value = 0 ;
        for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
        { value += 0.5 * m_pars [ index ( ix, iy, iz , iu ) ] * m_cache_x [ ix ] ; }
        r.setPar ( iy , iz , iu , value ) ; } } }
  //
  r *= ( m_xmax - m_xmin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over y dimension 
 *  \f$ f(x,z,u) =  \int_{y_{low}}^{y_{high}} F(x,y,z,u) {\mathrm{d}} y \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3 
Ostap::Math::LegendreSum4::integralY 
( const double ylow  , 
  const double yhigh ) const 
{
  LegendreSum3 r ( m_NX   , m_NZ , m_NU , 
                   m_xmin , m_xmax , 
                   m_zmin , m_zmax , 
                   m_umin , m_umax );
  //
  if ( s_equal ( ylow , yhigh ) ) { return  r ; }
  if (           ylow > yhigh   ) { return -integralY ( yhigh , ylow ) ; }
  //
  const  double yl = std::max ( ylow  , m_ymin ) ;
  const  double yh = std::min ( yhigh , m_ymax ) ;
  //
  if  ( yh <= m_ymin || yl >= m_ymax ) { return r  ; }
  //
  if  ( s_equal ( yl , m_ymin ) && s_equal ( yh , m_ymax ) ) { return integralY () ; }
  //
  const double tyl =  ty ( yl ) ;
  const double tyh =  ty ( yh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_y , tyl , tyh ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix )
  { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
    { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
      { 
        double value = 0 ;
        for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
        { value += 0.5 * m_pars [ index ( ix, iy, iz , iu ) ] * m_cache_y [ iy ] ; }
        r.setPar ( ix , iz , iu , value ) ; } } }
  //
  r *= ( m_ymax - m_ymin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over z dimension 
 *  \f$ f(x,y,u) =  \int_{z_{low}}^{z_{high}} F(x,y,z,u) {\mathrm{d}} z \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3 
Ostap::Math::LegendreSum4::integralZ 
( const double zlow  , 
  const double zhigh ) const 
{
  LegendreSum3 r ( m_NX , m_NY , m_NU ,
                   m_xmin , m_xmax ,
                   m_ymin , m_ymax , 
                   m_umin , m_umax );
  //
  if ( s_equal ( zlow , zhigh ) ) { return  r ; }
  if (           zlow > zhigh   ) { return -integralZ ( zhigh , zlow ) ; }
  //
  const  double zl = std::max ( zlow  , m_zmin ) ;
  const  double zh = std::min ( zhigh , m_zmax ) ;
  //
  if  ( zh <= m_zmin || zl >= m_zmax ) { return r  ; }
  //
  if  ( s_equal ( zl , m_zmin ) && s_equal ( zh , m_zmax ) ) { return integralY () ; }
  //
  const double tzl =  tz ( zl ) ;
  const double tzh =  tz ( zh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_z , tzl , tzh ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix )
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
      { 
        double value = 0 ;
        for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
        { value += 0.5 * m_pars [ index ( ix, iy, iz , iu ) ] * m_cache_z [ iz ] ; }
        r.setPar ( ix , iy , iu , value ) ; } } }
  //
  r *= ( m_zmax - m_zmin ) ;
  //
  return r ;
}
// ============================================================================
/*  integrate over u dimension 
 *  \f$ f(x,y,z) =  \int_{u_{low}}^{u_{high}} F(x,y,z,u) {\mathrm{d}} u \f$
 */
// ============================================================================
Ostap::Math::LegendreSum3 
Ostap::Math::LegendreSum4::integralU 
( const double ulow  , 
  const double uhigh ) const 
{
  LegendreSum3 r ( m_NX , m_NY , m_NZ ,
                   m_xmin , m_xmax ,
                   m_ymin , m_ymax , 
                   m_zmin , m_zmax );
  //
  if ( s_equal ( ulow , uhigh ) ) { return  r ; }
  if (           ulow > uhigh   ) { return -integralU ( uhigh , ulow ) ; }
  //
  const  double ul = std::max ( ulow  , m_umin ) ;
  const  double uh = std::min ( uhigh , m_umax ) ;
  //
  if  ( uh <= m_umin || ul >= m_umax ) { return r  ; }
  //
  if  ( s_equal ( ul , m_umin ) && s_equal ( uh , m_umax ) ) { return integralU () ; }
  //
  const double tul =  tu ( ul ) ;
  const double tuh =  tu ( uh ) ;
  //
  // prepare cache
  //
  _legendre_integrals ( m_cache_u , tul , tuh ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix )
  { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
    { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
      { 
        double value = 0 ;
        for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
        { value += 0.5 * m_pars [ index ( ix, iy, iz , iu ) ] * m_cache_u [ iu ] ; }
        r.setPar ( ix , iy , iz , value ) ; } } }
  //
  r *= ( m_umax - m_umin ) ;
  //
  return r ;
}
 
// ============================================================================
// The END 
// ============================================================================
