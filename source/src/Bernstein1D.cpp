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
#include <tuple>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Choose.h"
#include "Ostap/Bernstein1D.h"
#include "Ostap/Polynomials.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "bernstein_utils.h"
#include "syncedcache.h"
// ============================================================================
/** @file 
 *  Implementation file for functions, related to Bernstein's polynomnials 
 *
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /** @var s_PHASE2  
   *  useful delta-phase for 2D-positive polynomials 
   */
  const long double s_PHASE2 = std::asin ( std::sqrt ( 1.0L / 3.0L ) ) ;
  // ==========================================================================
  /** @var s_PHASE2  
   *  useful delta-phase for 2D-positive polynomials 
   */
  const std::vector<double> s_PHASES2 ( 1 , s_PHASE2) ;   
  // ==========================================================================
  
  // ==========================================================================
  typedef std::tuple<double,std::vector<double> >  PPROOTS ;
  const PPROOTS& deltas_for_pproots  ( const unsigned short N ) 
  {
    //
    typedef std::map<unsigned short, PPROOTS> MAP   ;
    typedef SyncedCache<MAP>                  CACHE ;
    // 
    // ========================================================================
    //
    static CACHE s_cache {} ;
    // ========================================================================
    { // look into the cache ==================================================
      CACHE::Lock lock { s_cache.mutex() } ;
      auto it = s_cache->find  ( N ) ;
      if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
      // ======================================================================
    } // ======================================================================
    // ========================================================================
    
    std::vector<double> pproots ( N + 1 ) ;
    double alpha = Ostap::Math::Utils::positive_pseudo_roots ( N , pproots ) ;
    pproots.push_back ( 1 ) ;
    //
    std::adjacent_difference ( pproots.begin() , pproots.end() , pproots.begin() );
    //
    std::transform ( pproots.begin() , 
                     pproots.end()   , 
                     pproots.begin() , 
                     [] (  const double x ) -> double 
                     { return std::sqrt ( std::max ( x , 0.0 ) ) ; } ) ;
    // ========================================================================
    // convert x to   phis 
    pproots = Ostap::Math::NSphere::phis ( pproots ) ;
    //
    alpha = std::acos ( std::sqrt ( alpha ) ) ;
    //    
    PPROOTS result = std::make_tuple ( alpha , pproots ) ;
    // ========================================================================
    { // update the cache =====================================================
      CACHE::Lock lock  { s_cache.mutex() } ;
      // update the cache
      s_cache->insert ( std::make_pair ( N , result ) ) ;
      auto it = s_cache->find  ( N ) ;
      return it->second ;
    } // ======================================================================
  }
  // ==========================================================================
}
// ============================================================================
// Bernstein EVEN 
// ============================================================================
/*  constructor
 *  @param N    the order of even polynomial 
 *  @param xmin low edge 
 *  @param xmax high edge 
 */
// ============================================================================
Ostap::Math::BernsteinEven::BernsteinEven 
( const unsigned short N    , 
  const double         xmin ,
  const double         xmax ) 
  : m_bernstein ( 0 == N % 2 ? N : N - 1 , xmin , xmax )
{}
// ============================================================================
/*  constructor from the list of coefficients 
 *  @param xmin low edge 
 *  @param xmax high edge 
 */
// ============================================================================
Ostap::Math::BernsteinEven::BernsteinEven 
( const std::vector<double>& pars , 
  const double               xmin ,
  const double               xmax ) 
  : m_bernstein ( pars.empty() ? 0 : 2 * pars.size() - 2 , xmin , xmax )
{
  setPars ( pars.begin() , pars.end() ) ;
}
// ============================================================================
/*  set k-parameter
 *  @param k index
 *  @param value new value 
 *  @return true if parameter is actually changed 
 */
// ============================================================================
bool Ostap::Math::BernsteinEven::setPar
( const unsigned short k     , 
  const double         value ) 
{
  //
  const unsigned short np  =  npars  () ;
  const unsigned short d   =  degree () ;
  //
  if      ( np <=     k ) { return false ; }
  else if ( d  == 2 * k ) { return m_bernstein.setPar ( k , value ) ; }
  //
  const bool b1 = m_bernstein.setPar (     k , value ) ;
  const bool b2 = m_bernstein.setPar ( d - k , value ) ;
  //
  return b1 || b2 ;
  //
}
// ============================================================================
// get all parameters (by value!!! COPY!!)
// ============================================================================
std::vector<double>
Ostap::Math::BernsteinEven::pars () const 
{
  return std::vector<double> ( m_bernstein.pars().begin() , 
                               m_bernstein.pars().begin() + npars() ) ;
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
double Ostap::Math::Utils::positive_pseudo_roots 
( const unsigned short N       , 
  std::vector<double>& pproots )
{
  //
  if      ( 2 >  N ) { pproots = {{     }} ; return 1.0   ; }
  else if ( 2 == N ) { pproots = {{ 0.5 }} ; return 1.0/3 ; }
  //
  pproots.resize ( N - 1 ) ;
  //
  const unsigned short K    =         N / 2   ;
  const bool           even =  ( 0 == N % 2 ) ;
  //
  if  ( even ) 
  {
    // roots here are K roots of T_{K} and (K-1) roots of U_{K-1}
    for ( unsigned short n = 0 ; n < N - 1 ; ++n ) 
    {
      const unsigned short k = n / 2 ;
      const double pp = 0 == n % 2 ? 
        - std::cos ( ( 2 * k + 1 ) * M_PIl / ( 2 * K ) ) :
        - std::cos ( (     k + 1 ) * M_PIl / (     K ) ) ; 
      pproots [ n ] = 0.5 * ( pp + 1 ) ;
    }
  }
  else 
  {
    // roots here are K roots of V_{K} and K roots of W_{K}
    for ( unsigned short n = 0 ; n < N - 1 ; ++n ) 
    {
      const unsigned short k = n / 2 ;
      const double pp = ( 0 == n % 2 ) ? 
        std::cos ( ( 2 * K - 2 * k     ) * M_PIl / ( 2 * K + 1 ) ) :
        std::cos ( ( 2 * K - 2 * k - 1 ) * M_PIl / ( 2 * K + 1 ) ) ;
      pproots [ n ] = 0.5 * ( pp + 1 ) ;
    }  
  }
  //
  if  ( even ) 
  {
    //
    const double mu    = 1 / ( 1.0 - N * N )     ;
    const double kappa = ( 1 + mu ) / ( 1 - mu ) ;
    const double alpha = kappa / ( 1 + kappa )   ;
    //
    return alpha ;
  }
  //
  return 0.5 ; 
}

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive::Positive
( const unsigned short       N    ,
  const double               xmin ,
  const double               xmax )
  : m_bernstein ( N , xmin , xmax )
  , m_sphereA   ( 0 == N ? 0 :     1 ) 
  , m_sphereR   ( N <  2 ? 0 : N - 1 )
  , m_rs        ( N <  2 ? 0 : N - 1 , 0.0 ) 
  , m_v1        ( N <  2 ? 0 : N + 1 , 0.0 ) 
  , m_v2        ( N <  2 ? 0 : N + 1 , 0.0 ) 
  , m_aux       ( N <  2 ? 0 : N + 1 , 0.0 ) 
{
  if ( 2 <= N ) 
  {
    const PPROOTS& pp = deltas_for_pproots ( N ) ;
    m_sphereA = NSphere ( "" , {{ std::get<0> ( pp ) }} ) ; 
    m_sphereR = NSphere ( "" ,    std::get<1> ( pp )    ) ;
  }
  updateBernstein () ;
}
// ============================================================================
// constructor from the list of parameters/phases 
// ============================================================================
 Ostap::Math::Positive::Positive
 ( const std::vector<double>& pars ,
   const double               xmin ,
   const double               xmax )
   : Positive ( pars.begin (), pars.end () , xmin , xmax ) 
 {}
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
  // few simple cases, treated explicitely  
  //
  if       ( 0 == o ) { return m_bernstein.setPar ( 0 , norm ) ; }
  else if  ( 1 == o )  
  {
    /// get alpha and beta from A-sphere 
    const long double alpha     = m_sphereA . x2 ( 0 ) ;
    const long double beta      = 1 - alpha ;
    ///
    const bool updated0 = m_bernstein.setPar ( 0 , alpha * norm ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , beta  * norm ) ;
    update              = updated1 || update ;
    //
    return updated0 || updated1 ;
  }
  else if  ( 2 == o )
  {
    /// get alpha and beta from A-sphere 
    const long double alpha     = m_sphereA . x2 ( 0 ) ;
    const long double beta      = m_sphereA . x2 ( 1 ) ;
    //
    // get the root from R-sphere 
    const long double r         = m_sphereR . x2 ( 0 ) ;
    ///
    std::array<long double,3> v2 = { { 0.0 , 0.0 , 0.0 } } ;
    Ostap::Math::Utils::bernstein2_from_roots ( r , r , v2 ) ;
    //
    const long double sv = v2 [ 0 ] + v2 [ 1 ] + v2 [ 2 ] ;
    const long double na = alpha * norm / sv ;
    //
    const long double c0 = v2 [ 0 ] * na                ;
    const long double c1 = v2 [ 1 ] * na  + beta * norm ;
    const long double c2 = v2 [ 2 ] * na                ;
    //
    const bool updated1 = m_bernstein.setPar ( 0 , c0 ) ;
    const bool updated2 = m_bernstein.setPar ( 1 , c1 ) ;
    const bool updated3 = m_bernstein.setPar ( 2 , c2 ) ;
    //
    return updated1 || updated2 || updated3 ;
    //
  }
  // ==========================================================================
  // generic case 
  // ==========================================================================
  /// get alpha and beta from A-sphere 
  const long double alpha     = m_sphereA . x2 ( 0 ) ;
  const long double beta      = 1 - alpha ;
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
    m_v1 [ 0 ] = alpha                                      ; n1 = 1 ;    
    m_v2 [ 0 ] = 0     ; m_v2 [ 1 ] = beta ; m_v2 [ 2 ] = 0 ; n2 = 3 ;
  }
  else 
  {
    m_v1 [ 0 ] = alpha ; m_v1 [ 1 ] = 0                         ; n1 = 2 ;
    m_v2 [ 0 ] = 0     ; m_v2 [ 1 ] = beta                      ; n2 = 2 ;
  }
  //
  for ( unsigned short iR = 0 ; iR < nR ; ++iR ) 
  {
    const long double r  = m_rs [ iR ] ;    
    Ostap::Math::Utils::bernstein2_from_roots ( r , r , br ) ;
    //
    if ( 0 == iR % 2 ) 
    {
      Ostap::Math::Utils::b_multiply ( m_v1.begin () , m_v1.begin () + n1 , br , m_aux.begin () );
      n1  += 2 ;
      std::swap  ( m_v1 , m_aux ) ;
    }
    else 
    {
      Ostap::Math::Utils::b_multiply ( m_v2.begin () , m_v2.begin () + n2 , br , m_aux.begin () );
      n2  += 2 ;
      std::swap  ( m_v2 , m_aux ) ;      
    }
  }
  //
  const long double s1 = norm * alpha / std::accumulate ( m_v1.begin () , m_v1.end() , 0.0L ) ;
  const long double s2 = norm * beta  / std::accumulate ( m_v2.begin () , m_v2.end() , 0.0L ) ;
  //
  const unsigned short nP = m_bernstein.npars() ;
  for ( unsigned short iP = 0 ; iP < nP ; ++iP )
  { 
    const bool updated = m_bernstein.setPar ( iP , m_v1 [ iP ] * s1 + m_v2 [ iP ] * s2 ) ; 
    update = updated || update ;
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
  std::vector<double>::iterator t = r.begin() ;
  //
  const std::vector<double>& pa = m_sphereA.pars() ;
  t = std::copy ( pa.begin () , pa.end() , t ) ;
  const std::vector<double>& pr = m_sphereR.pars() ;
  t = std::copy ( pr.begin () , pr.end() , t ) ;
  //
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
  : m_even     ( N     , xmin , xmax )
  , m_positive ( N / 2 , xmin , xmax )  
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
  : m_even      ( 2 * pars.size() , xmin , xmax )
  , m_positive  (     pars        , xmin , xmax ) 
{
  updateBernstein () ;
}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::PositiveEven::updateBernstein ()
{
  const std::vector<double>& p = m_positive.bernstein().pars() ;
  return m_even.setPars ( p.begin() , p.end() ) ;
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
  : m_bernstein  ( N                  , xmin , xmax ) 
  , m_positive   ( 1 <= N ? N - 1 : 0 , xmin , xmax ) 
  , m_sphere     ( 1 <= N ?     1 : 0 )
  , m_increasing ( increasing )  
  , m_aux        ( N + 1 ) 
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
  : Monotonic ( pars.begin() , pars.end() , xmin , xmax , increasing )
{}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Monotonic::updateBernstein ()
{
  //
  bool update = false ;
  //
  const unsigned short o = degree() ;
  const long double    norm = m_bernstein.npars() / 
    ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  if      ( 0 == o ) { return m_bernstein.setPar ( 0 , norm ) ; }
  else if ( 1 == o ) 
  {
    /// get alpha from sphere 
    const long double a     = m_sphere . x2 ( 0 ) ;
    const long double alpha = m_increasing ? a / ( 1 + a ) : 1 / ( 1 + a ) ;
    const long double beta  = 1 - alpha     ;
    //
    const bool updated0 = m_bernstein.setPar ( 0 , alpha * norm ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , beta  * norm ) ;
    update              = updated1 || update ;
    return updated0 || updated1 ;
  }
  // 
  // generic case  
  //
  // get alpha and beta from sphere 
  const long double alpha  = m_sphere . x2 ( 0 ) ;
  const long double beta   = 1 - alpha ;
  //
  const std::vector<double>& ppars = m_positive.bernstein().pars() ;
  if ( m_increasing ) { std::partial_sum ( ppars. begin () , ppars. end () , m_aux. begin () + 1 ) ; }
  else                { std::partial_sum ( ppars.rbegin () , ppars.rend () , m_aux.rbegin () + 1 ) ; }
  //
  const unsigned short nx = m_aux.size() ;
  const double s1 = alpha * norm / std::accumulate ( m_aux.begin() , m_aux.end() , 0.0 ) ;
  const double s2 = beta  * norm / nx    ;
  //
  for ( unsigned short i = 0 ; i < nx ; ++i ) 
  { 
    const bool updated = m_bernstein.setPar ( i , s2 + s1 * m_aux [ i ] ) ; 
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
// get all parameters (by value)
// ============================================================================
std::vector<double> Ostap::Math::Monotonic::pars  () const
{
  std::vector<double> r { m_positive.pars() } ; 
  const std::vector<double>& pa = m_sphere  .pars() ;
  r.insert ( r.begin() , pa.begin() , pa.end() ) ;
  return r  ;  
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
  : m_bernstein  ( N                  , xmin , xmax ) 
  , m_positive   ( 2 <= N ? N - 2 : 0 , xmin , xmax ) // needed for  2<= N 
  , m_sphereA    ( 2 <= N ? 1 : 0  )                  // needed for 2 <= N  
  , m_sphereI    ( 1 <= N ? 1 : 0  )                  // needed for 1 <= N 
  , m_increasing ( increasing )
  , m_convex     ( convex     )  
  , m_aux        ( N + 1      )
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
  : Convex ( pars.begin() , pars.end(),  xmin , xmax , increasing , convex ) 
{}
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
    /// get alpha from "linear" sphere 
    const long double a     = m_sphereI . x2 ( 0 ) ;
    const long double alpha = m_increasing ? a / ( 1 + a ) : 1 / ( 1 + a ) ;
    const long double beta  = 1 - alpha     ;
    //
    const bool updated0 = m_bernstein.setPar ( 0 , alpha * norm ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , beta  * norm ) ;
    update              = updated1 || update ;
    return updated0 || updated1 ;
  }
  //
  const std::vector<double>& b_pars = m_positive.bpars() ;
  //
  m_aux [ 0 ] = 0 ; m_aux [ 1 ] = 0 ;
  std::copy ( b_pars.begin() , b_pars.end() , m_aux.begin () + 2 ) ;
  //
  // integrate them twice and to get new coefficients
  std::partial_sum ( m_aux. begin () + 2 , m_aux. end () , m_aux. begin () + 2 ) ; 
  std::partial_sum ( m_aux. begin () + 2 , m_aux. end () , m_aux. begin () + 2 ) ; 
  //
  if ( !m_convex ) 
  { Ostap::Math::scale_and_shift ( m_aux.begin () , m_aux.end () , -1.0L , m_aux.back () ) ; }  
  //
  if ( m_increasing != m_convex ) { std::reverse ( m_aux.begin() , m_aux.end() ) ; }
  //
  /// get alpha & beta from A-sphere 
  const long double alpha = m_sphereA . x2 ( 0 ) ;
  const long double beta  = 1 - alpha ;
  //
  const long double is1 = alpha * norm / std::accumulate ( m_aux.begin() , m_aux.end() , 0.0L) ;
  
  // add a positive linear function 
  const long double a  =  m_sphereI.x2 ( 0 ) ;
  const long double a0 = m_increasing ? a / ( 1 + a ) : 1 / ( 1 + a ) ;
  const long double a1 = 1 - a0 ;
  
  const unsigned short np  = m_bernstein.npars() ;
  const long double    is2 = 2 * beta * norm / np ;
  for ( unsigned short k   = 0 ; k < np ; ++k ) 
  {
    const long double w = double ( k ) / o ;
    const long double f = a0 * ( 1 - w ) + a1 * w ;
    //
    const bool updated = m_bernstein.setPar ( k , is1 * m_aux [ k ] + is2 * f ) ;
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
// get all parameters (by value)
// ============================================================================
std::vector<double> Ostap::Math::Convex::pars  () const
{
  std::vector<double> r { m_positive.pars() } ; 
  const std::vector<double>& pi = m_sphereI .pars() ;
  r.insert ( r.begin() , pi.begin() , pi.end() ) ;
  const std::vector<double>& pa = m_sphereA .pars() ;
  r.insert ( r.begin() , pa.begin() , pa.end() ) ;
  return r  ;  
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
