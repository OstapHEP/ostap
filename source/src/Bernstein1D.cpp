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
#include "Ostap/Choose.h"
#include "Ostap/Bernstein1D.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "bernstein_utils.h"
// ============================================================================
/** @file 
 *  Implementation file for functions, related to Bernstein's polynomnials 
 *
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
// Bernstein EVEN 
// ============================================================================
/*  constructor
 *  the actual degree of polynomial will be 2*N
 *  @param N  parameter that defiend the order of polynomial (2*N)
 *  @param xmin low edge 
 *  @param xmax high edge 
 */
// ============================================================================
Ostap::Math::BernsteinEven::BernsteinEven 
( const unsigned short N    , 
  const double         xmin ,
  const double         xmax ) 
  : m_N         (   N                   )   
  , m_bernstein ( 2*N + 1 , xmin , xmax )
{}
// ============================================================================
/*  constructor from list of coefficients 
 *  @param xmin low edge 
 *  @param xmax high edge 
 */
// ============================================================================
Ostap::Math::BernsteinEven::BernsteinEven 
( const std::vector<double>& pars , 
  const double               xmin ,
  const double               xmax ) 
  : m_N         (   pars.size()                   )   
  , m_bernstein ( 2*pars.size() + 1 , xmin , xmax )
{
  for ( unsigned short i = 0 ; i < pars.size() ; ++i ) { setPar ( i , pars[i] ) ; }
}
// ============================================================================
/* set k-parameter
 *  @param k index
 *  @param value new value 
 *  @return true if parameter is actually changed 
 */
// ============================================================================
bool Ostap::Math::BernsteinEven::setPar
( const unsigned short k , const double value ) 
{
  if ( npars() <= k ) { return false ; }
  const bool b1 = m_bernstein.setPar (             k , value ) ;
  const bool b2 = m_bernstein.setPar ( 2*m_N + 1 - k , value ) ;
  return b1 || b2 ;
}
// ============================================================================
// get all parameters (by value!!! COPY!!)
// ============================================================================
std::vector<double>
Ostap::Math::BernsteinEven::pars () const 
{
  return std::vector<double>( m_bernstein.pars().begin()           , 
                              m_bernstein.pars().begin() + m_N + 1 ) ;                            
}
// ============================================================================
//  Sum of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__add__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp += value ; return tmp ; }
// ============================================================================
//  Sum of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__radd__  ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp += value ; return tmp ; }
// ============================================================================
//  Product of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__mul__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp *= value ; return tmp ; }
// ============================================================================
//  Product of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__rmul__  ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp *= value ; return tmp ; }
// ============================================================================
//  Subtract a constant from Bernstein polynomial
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__sub__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp -= value ; return tmp ; }
// ============================================================================
//  Subtract (right) a constant and Bernstein polynomial
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__rsub__  ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp *= -1 ; tmp += value ; return tmp ; }
// ============================================================================
//  Division of Bernstein polynomial and constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__truediv__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp /= value ; return tmp ; }
// ============================================================================



// ============================================================================
// POSITIVE 
// ============================================================================

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive::Positive
( const unsigned short      N     ,
  const double              xmin  ,
  const double              xmax  )
  : m_bernstein ( N , xmin , xmax )
  , m_sphereA   ( 1                      ) 
  , m_sphereR   ( std::max ( 1 , N - 1 ) )
  , m_rs        ( std::max ( 1 , N - 1 ) , 0.0  ) 
  , m_v1        ( N + 1 , 0.0 ) 
  , m_v2        ( N + 1 , 0.0 ) 
  , m_aux       ( N + 1 , 0.0 ) 
{
  updateBernstein () ;
}
// // ============================================================================
// // constructor from the list of phases 
// // ============================================================================
// Ostap::Math::Positive::Positive
// ( const std::vector<double>& pars ,
//   const double               xmin ,
//   const double               xmax )
//   : m_bernstein ( pars.size() , xmin , xmax )
//   , m_sphere    ( pars , 3 ) 
// {
//   updateBernstein () ;
// }
// // ============================================================================
// // constructor from the sphere with coefficients  
// // ============================================================================
// Ostap::Math::Positive::Positive
// ( const Ostap::Math::NSphere& sphere , 
//   const double                xmin   , 
//   const double                xmax   )
//   : m_bernstein ( sphere.dim() , xmin , xmax )
//   , m_sphere    ( sphere ) 
// {
//   updateBernstein () ;
// }
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive::updateBernstein ()
{
  //
  bool update = false ;
  //
  // degree 
  const unsigned short o    = degree() ;
  //
  const long double    norm = m_bernstein.npars() / 
    ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  // few simple cases 
  //
  if       ( 0 == o ) { return m_bernstein.setPar ( 0 , norm ) ; }
  else if  ( 1 == o )  
  {
    const bool updated0 = m_bernstein.setPar ( 0 , m_sphereA.x2 ( 0 ) * norm ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , m_sphereA.x2 ( 1 ) * norm ) ;
    update              = updated1 || update ;
    //
    return updated0 || updated1 ;
  }
  else if  ( 2 == o )
  {
    // get the root from R-sphere 
    const long double r  = m_sphereR.x2 ( 0 ) ;
    const long double ri = 1 - r ;
    /// get alpha and beta from A-sphere 
    long double alpha     = m_sphereA . x2 ( 0 ) ;
    long double beta      = m_sphereA . x2 ( 1 ) ;
    ///
    const long double c0 =   r  * r  * alpha ;
    const long double c1 = - r  * ri * alpha ;
    const long double c2 =   ri * ri * alpha ;
    //
    const bool updated1 = m_bernstein.setPar ( 0 , c0            ) ;
    const bool updated2 = m_bernstein.setPar ( 1 , c1 + 2 * beta ) ;
    const bool updated3 = m_bernstein.setPar ( 2 , c2            ) ;
    //
    m_bernstein *= ( norm / ( c0 + c1 + c2 + 2 * beta ) ) ;
    //
    return updated1 || updated2 || updated3 ;
  }
  // ==========================================================================
  // generic case 
  // ==========================================================================
  /// get alpha and beta from A-sphere 
  const long double alpha     = m_sphereA . x2 ( 0 ) ;
  const long double beta      = m_sphereA . x2 ( 1 ) ;
  ///
  /// get root-parameters from R-sphere and integrate them to get the roots 
  const unsigned nR = m_rs.size() ;
  for ( unsigned short iR = 0 ; iR < nR ; ++iR ) { m_rs [ iR ] = m_sphereR.x2 ( iR ) ; }
  std::partial_sum ( m_rs.begin() , m_rs.end () , m_rs.begin() ) ;
  ///
  const bool even = ( 0 == o % 2 );
  //
  std::array<long double,3> br ;   // helper second order polynomial 
  /// "half-order"
  const unsigned short m = even ? o / 2 :  ( o - 1 ) / 2 ;
  //
  unsigned short n1 = 2 ;
  unsigned short n2 = 2 ;
  //
  if ( even )
  {
    m_v1 [ 0 ] = alpha ;                                           n1 = 1 ;    
    m_v2 [ 0 ] = 0     ; m_v2 [ 1 ] = 0.5 *beta ; m_v2 [ 2 ] = 0 ; n2 = 3 ;
  }
  else 
  {
    m_v1 [ 0 ] = beta  ; m_v1 [ 1 ] = 0     ; n1 = 2 ;
    m_v2 [ 0 ] = 0     ; m_v2 [ 1 ] = alpha ; n2 = 2 ;
  }
  //
  for ( unsigned short iR = 0 ; iR < nR ; ++iR ) 
  {
    const long double r  = m_rs [ iR ] ;
    const long double ri = 1 - r ;
    //
    br[ 0 ] =   r  * r ;
    br[ 1 ] =  -r  * ri ;
    br[ 2 ] =   ri * ri ;
    //
    if ( 0 == iR % 2 ) 
    {
      Ostap::Math::Utils::b_multiply 
        ( m_v1.begin () , m_v1.begin () + n1 , br.begin () , br.end () , m_aux.begin () );
      n1  += 2 ;
      std::swap  ( m_v1 , m_aux ) ;
    }
    else 
    {
      Ostap::Math::Utils::b_multiply 
        ( m_v2.begin () , m_v2.begin () + n2 , br.begin () , br.end () , m_aux.begin () );
      n2  += 2 ;
      std::swap  ( m_v2 , m_aux ) ;      
    }
    //
  }
  //
  const unsigned short nP = m_bernstein.npars() ;
  for ( unsigned short iP = 0 ; iP < nP ; ++iP )
  { update |= m_bernstein.setPar ( iP , m_v1 [ iP ] + m_v2 [ iP ] ) ; }
  //
  if ( update ) 
  {
    const long double s1 = std::accumulate ( m_bernstein.pars() . begin () , 
                                             m_bernstein.pars() . end   () , 0.0L ) ;
    m_bernstein *= norm / s1 ;
  }
  //
  return update ;
}
// =============================================================================
// get the integral between low and high 
// =============================================================================
double Ostap::Math::Positive::integral
( const double low , const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_bernstein.integral ( low , high )  ; 
}
// ============================================================================
// get all parameters (phases on sphere)
// ============================================================================
std::vector<double> Ostap::Math::Positive::pars  () const
{
  std::vector<double> r ( npars() , 0.0 ) ;
  r [ 0 ] = m_sphereA.par( 0 ) ;
  const std::vector<double>& p = m_sphereR.pars () ;
  std::copy  ( p.begin() , p.end() , r.begin() + 1 ) ;
  return r  ;  
}
   





// ============================================================================
// POSITIVE EVEN
// ============================================================================


// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::PositiveEven::PositiveEven
( const unsigned short      N    ,
  const double              xmin ,
  const double              xmax )
  : m_even   ( N , xmin , xmax )
  , m_sphere ( N , 3 ) 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the list of phases 
// ============================================================================
Ostap::Math::PositiveEven::PositiveEven
( const std::vector<double>& pars ,
  const double               xmin ,
  const double               xmax )
  : m_even      ( pars.size() , xmin , xmax )
  , m_sphere    ( pars , 3 ) 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the sphere with coefficients  
// ============================================================================
Ostap::Math::PositiveEven::PositiveEven
( const Ostap::Math::NSphere& sphere , 
  const double                xmin   , 
  const double                xmax   )
  : m_even    ( sphere.dim() , xmin , xmax )
  , m_sphere  ( sphere ) 
{
  updateBernstein () ;
}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::PositiveEven::updateBernstein ()
{
  ///
  bool         update = false ;
  /// degree 
  const unsigned short o = m_even.degree() ;
  //
  const double   norm    = m_even.npars() / 
    ( m_even.xmax() -  m_even.xmin () ) ;
  //
  // few simple cases 
  //
  if       ( 0 == o ) { return m_even.setPar( 0 , norm ) ; }
  //
  // get the parameters of "global" non-negative symmetric parabola 
  //
  const double a0 = m_sphere.x2(0)      ;
  const double a1 = m_sphere.x2(1) - a0 ;
  const double a2 = a0                  ;
  //
  // "elevate to degree of bernstein"
  const unsigned short N =  m_even.bernstein().degree()  ;
  std::vector<long double> v ( N + 1 ) ;
  v[0] = a0 ;
  v[1] = a1 ;
  v[2] = a2 ;
  std::fill ( v.begin() + 3 , v.end() , a2 ) ;
  // repeate the elevation cycles: 
  for ( unsigned short   n = 2  ; n < N ; ++n ) 
  {
    // "current" degree 
    for ( unsigned short k = n ;  1<= k ; --k ) 
    {
      v[k]  = ( n + 1 - k ) * v[k] + k * v[k-1] ;
      v[k] /=   n + 1  ;
    } 
  }
  //
  // now  we have a non-negative symmetric parabola coded.
  //   - add a proper positive polynomial to it.
  const unsigned short nV = v.size() ;
  const unsigned short nX = m_sphere.nX() ;
  for ( unsigned short ix = 2 ; ix < nX ; ++ix ) 
  {
    const double x = m_sphere.x2 ( ix ) ;
    v[      ix - 2 ] += x ; 
    v[ nV - ix + 1 ] += x ; // keep symmetry 
  }
  //
  const double isum = norm / std::accumulate ( v.begin() , v.end() , 0.0L ) ;
  //
  const unsigned short nE = m_even.npars() ;
  //
  for ( unsigned short ix = 0 ; ix < nE ; ++ix ) 
  {
    const bool updated = m_even.setPar ( ix , 2 * v[ix] * isum ) ;
    update = updated || update ;
  }
  //
  return update ;
}
// =============================================================================
// get the integral between low and high 
// =============================================================================
double Ostap::Math::PositiveEven::integral
( const double low , const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_even.integral ( low , high )  ; 
}
// ============================================================================


// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Monotonic::Monotonic 
( const unsigned short      N          ,
  const double              xmin       ,
  const double              xmax       , 
  const bool                increasing ) 
  : m_bernstein  ( N , xmin , xmax ) 
  , m_sphere     ( N , 3           )
  , m_increasing ( increasing      )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Monotonic::Monotonic 
( const std::vector<double>& pars       ,
  const double               xmin       ,
  const double               xmax       ,
  const bool                 increasing ) 
  : m_bernstein  ( pars.size () , xmin , xmax ) 
  , m_sphere     ( pars , 3   ) 
  , m_increasing ( increasing )  
{
  updateBernstein () ;
}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Monotonic::updateBernstein ()
{
  //
  bool   update = false ;
  //
  // get sphere coefficients 
  std::vector<double> v ( m_sphere.nX() ) ;
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { v[ix] = m_sphere.x2 ( ix ) * ( ix + 1 ) ; }
  //
  // integrate them and to get new coefficients
  if   ( m_increasing ) { std::partial_sum ( v. begin() , v. end() ,  v. begin() ) ; }
  else                  { std::partial_sum ( v.rbegin() , v.rend() ,  v.rbegin() ) ; }
  //
  const double isum = m_bernstein.npars() 
    / std::accumulate ( v.begin() , v.end() , 0.0 ) 
    / ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// get the minimal value of function 
// ============================================================================
double Ostap::Math::Monotonic::fun_min () const
{
  const std::vector<double>& ps = m_bernstein.pars() ;
  return  std::min ( ps.front() , ps.back() ) ;
}
// ============================================================================
// get the maximal value of function 
// ============================================================================
double Ostap::Math::Monotonic::fun_max () const
{
  const std::vector<double>& ps = m_bernstein.pars() ;
  return  std::max( ps.front() , ps.back() ) ;
}
// ============================================================================
// get the integral between low and high 
// =============================================================================
double Ostap::Math::Monotonic::integral
( const double low , const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_bernstein.integral ( low , high )  ; 
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Convex::Convex
( const unsigned short      N          ,
  const double              xmin       ,
  const double              xmax       , 
  const bool                increasing ,
  const bool                convex     ) 
  : m_bernstein  ( N , xmin , xmax ) 
  , m_sphere     ( N , 3      ) 
  , m_increasing ( increasing )
  , m_convex     ( convex     )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Convex::Convex
( const std::vector<double>& pars       ,
  const double               xmin       ,
  const double               xmax       ,
  const bool                 increasing ,
  const bool                 convex     ) 
  : m_bernstein  ( pars.size() , xmin, xmax ) 
  , m_sphere     ( pars , 3   ) 
  , m_increasing ( increasing )
  , m_convex     ( convex     )  
{
  updateBernstein () ;
}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Convex::updateBernstein ()
{
  //
  bool   update = false ;
  //  
  /// degree 
  const unsigned short o = degree() ;
  //
  const double   norm    = m_bernstein.npars() / 
    ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  // few simple cases 
  //
  if       ( 0 == o ) { return m_bernstein.setPar ( 0 , norm ) ; }
  else if  ( 1 == o )  
  {
    const double a  =  m_sphere.x2 ( 0 ) ;
    const double b  =  m_sphere.x2 ( 1 ) ;
    const double x0 =  m_increasing ? a     : a + b ;
    const double x1 =  m_increasing ? a + b : a     ;
    //
    const bool updated0 = m_bernstein.setPar ( 0 , x0 * norm / ( 2 * a + b ) ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , x1 * norm / ( 2 * a + b ) ) ;
    update              = updated1 || update ;
    //
    return update ;
  }
  //
  // get sphere coefficients 
  //
  std::vector<double>  v ( m_sphere.nX() ) ;
  const unsigned short vs = v.size()    ;
  //
  const std::array<double,2> a = { { m_sphere.x2(0) , m_sphere.x2(1) } };
  for ( unsigned short ix = 2 ; ix < vs ; ++ix ) 
  { v[ix] = m_sphere.x2 ( ix ) ; }
  //
  // integrate them twice and to get new coefficients
  std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
  std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
  //
  if ( !m_convex ) 
  {
    const  double last = v.back() ;
    for ( unsigned short k = 0 ; k < vs; ++k) 
    { 
      v[k] = last  - v[k] ; 
      if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
    }
  }
  //
  if ( m_increasing != m_convex )
  { std::reverse ( v.begin() , v.end() ) ; }
  //
  // add a positive linear function 
  //
  const unsigned short d = degree() ;
  for ( unsigned short k = 0 ; k < vs ; ++k ) 
  {
    const double r1 =  double(k) / d ;
    //
    v[k] +=  
      m_increasing ?  
      a[0] +       r1   * a[1] :
      a[0] + ( 1 - r1 ) * a[1] ;
    //
    if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
  }
  //
  const double isum = m_bernstein.npars() 
    / std::accumulate ( v.begin() , v.end() , 0.0 ) 
    / ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  for ( unsigned short ix = 0 ; ix < vs ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;  
}
// ============================================================================
// get the minimal value of function 
// ============================================================================
double Ostap::Math::Convex::fun_min () const
{
  const std::vector<double>& ps = m_bernstein.pars() ;
  return  std::min ( ps.front() , ps.back() ) ;
}
// ============================================================================
// get the maximal value of function 
// ============================================================================
double Ostap::Math::Convex::fun_max () const
{
  const std::vector<double>& ps = m_bernstein.pars() ;
  return  std::max( ps.front() , ps.back() ) ;
}
// ============================================================================
// get the integral between low and high 
// =============================================================================
double Ostap::Math::Convex::integral
( const double low , const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_bernstein.integral ( low , high )  ; 
}




// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ConvexOnly::ConvexOnly 
( const unsigned short      N     ,
  const double              xmin   ,
  const double              xmax   , 
  const bool                convex ) 
  : m_bernstein ( N , xmin , xmax ) 
  , m_sphere    ( N , 3  )
  , m_convex    ( convex )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ConvexOnly::ConvexOnly 
( const std::vector<double>& pars    ,
  const double               xmin    ,
  const double               xmax    ,
  const bool                 convex  ) 
  : m_bernstein ( pars.size()  , xmin , xmax ) 
  , m_sphere    ( pars , 3  )
  , m_convex    ( convex )  
{
  updateBernstein () ;
}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::ConvexOnly::updateBernstein ()
{
  // 
  bool   update = false ;
  // 
  const unsigned short o = degree() ;
  //
  const double   norm    = m_bernstein.npars() / 
    ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  if       ( 0 == o ) { return m_bernstein.setPar( 0 , norm ) ; }
  else if  ( 1 == o )  
  {
    const bool updated0 = m_bernstein.setPar ( 0 , m_sphere.x2 ( 0 ) * norm ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , m_sphere.x2 ( 1 ) * norm ) ;
    update              = updated1 || update ;
    //
    return update ;
  }
  //
  // get sphere coefficients 
  //
  std::vector<double>  v ( m_sphere.nX() ) ;
  const unsigned short vs = v.size()    ;
  //
  // get parameters from the sphere:
  //
  if ( !m_convex ) 
  {
    const std::array<double,2> a = { { m_sphere.x2(0) , m_sphere.x2(1) } };
    for ( unsigned short ix = 2 ; ix < vs ; ++ix ) 
    { v[ix] = m_sphere.x2 ( ix ) ; }
    //
    // integrate them twice and to get new coefficients
    std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
    std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
    //
    {
      const  double last = v.back() ;
      for ( unsigned short k = 0 ; k < vs; ++k) 
      { 
        v[k] = last  - v[k] ; 
        if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
      }
    }
    //
    // subtract the linear component and 
    // add positive linear function
    //
    const double v1 = a[0] - v.front() ;
    const double v2 = a[1] - v.back()  ;
    const unsigned int   d = degree() ;
    for ( unsigned short k = 0 ; k < vs ; ++k ) 
    {
      const double r1 =  double(k)  / d ;
      v[k] +=  ( 1 - r1 ) * v1  + r1 * v2 ;
      if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
    }
  }
  else 
  { 
    std::array<double,3> a = { { m_sphere.x2(0) , 
                                 m_sphere.x2(1) , 
                                 m_sphere.x2(2) } };
    for ( unsigned short ix = 3 ; ix < vs ; ++ix ) 
    { v[ix] = m_sphere.x2 ( ix ) ; }
    // integrate them twice and to get new coefficients
    std::partial_sum ( v.  begin() + 3 , v.  end()     ,  v.  begin() + 3 ) ; 
    std::partial_sum ( v.  begin() + 3 , v.  end()     ,  v.  begin() + 3 ) ; 
    //    
    const double a0 = a[0] ;
    const double a2 = a[2] ;
    const double a1_min = -1*std::sqrt ( a0 * a2 ) ;
    const double a1_max = 0.5 * ( a0 + a2 ) ;
    //
    const double a1 = a1_min + a[1] * ( a1_max - a1_min ) ;
    //
    const double c0 = a0         ;
    const double c1 = 2*(a1-a0)  ;
    const double c2 = a0+a2-2*a1 ; 
    //
    const unsigned int   d = degree() ;
    for ( unsigned short k = 0 ; k < vs ; ++k ) 
    {
      double vv = c0 ;
      const double r1 =  double(k) / d ;
      if ( 0 != k ) { vv += r1 * c1 ; }
      if ( 1 <  k ) { vv += r1 * ( k - 1 ) * c2 / ( d - 1 ) ; }
      v[k] +=  vv ;
      if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
    }
  }
  //
  const double isum = m_bernstein.npars() 
    / std::accumulate ( v.begin() , v.end() , 0.0 ) 
    / ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  for ( unsigned short ix = 0 ; ix < vs ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// get the integral between low and high 
// =============================================================================
double Ostap::Math::ConvexOnly::integral
( const double low , const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_bernstein.integral ( low , high )  ; 
}



// ============================================================================
//                                                                      The END 
// ============================================================================
