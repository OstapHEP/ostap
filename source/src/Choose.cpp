// $Id$ 
// ============================================================================
// Include files
// ============================================================================
// STD/STL 
// ============================================================================
#include <cstdint>
#include <cmath>
#include <map>
#include <climits>
#include <algorithm>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Clenshaw.h"
#include "Ostap/Choose.h"
#include "Ostap/Math.h"
// ============================================================================
//  local 
// ============================================================================
#include "syncedcache.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Calculaate binbomial coefficients and related quantities 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<std::uintmax_t>::is_specialized,
		  "numeric_limits<std::uintmax_t> is not specialzaed!" ) ;
  static_assert ( std::numeric_limits<std::size_t>::is_specialized,
		  "numeric_limits<std::size_t> is not specialzaed!" ) ;
  static_assert ( std::numeric_limits<unsigned long long>::is_specialized,
		  "numeric_limits<usigned long long> is not specialzaed!" ) ;
  static_assert ( std::numeric_limits<unsigned long>::is_specialized,
		  "numeric_limits<usigned long> is not specialzaed!" ) ;
  static_assert ( std::numeric_limits<unsigned int>::is_specialized,
		  "numeric_limits<usigned int> is not specialzaed!" ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned int>::max ()       <= std::numeric_limits<unsigned long>::max ()      && 		  
		  std::numeric_limits<unsigned long>::max ()      <= std::numeric_limits<unsigned long long>::max () &&  
		  std::numeric_limits<unsigned long long>::max () <= std::numeric_limits<std::uintmax_t>::max ()     && 
		  std::numeric_limits<std::size_t>::max ()        <= std::numeric_limits<std::uintmax_t>::max () ,
		  "Invalid std::numeric_limits<*UINT*>::max () hierarchy!" ) ;
  // ==========================================================================================================
  const std::uintmax_t s_ULLMAX = std::numeric_limits<std::uintmax_t>::max () - 1 ;
  const long double    s_emax   = std::log ( 0.4L * ( s_ULLMAX / 2 ) ) ;
  // ==========================================================================
  /** calculate the binomial coefficient C(k,n) =  n!/((n-k)!*k!)
   *  In case of overflow std::numeric_limits<unsigned long long>::max is returned 
   */
  inline std::uintmax_t choose 
  ( unsigned short n , 
    unsigned short k ) 
  {
    //
    // =========================================================================
    //
    if      ( k > n                ) { return 0 ; }
    else if ( 0 == k || n == k     ) { return 1 ; }
    else if ( 1 == k || n == k + 1 ) { return n ; }
    else if ( 2 == k || n == k + 2 ) { return 1ull * n * ( n - 1 ) / 2 ; }
    //
    k = std::min ( k , (unsigned short) ( n - k ) ) ;
    //
    std::uintmax_t r = 1  ;
    const std::uintmax_t L = s_ULLMAX / n ; 
    for ( unsigned short d = 1 ; d <= k  ; ++d , --n ) 
    {
      // if ( r > s_ULLMAX / n * d ) { return s_ULLMAX ; }  //  RETURN
      if ( r > L * d ) { return s_ULLMAX ; }  //  RETURN
      r = ( r / d ) * n + ( r % d ) * n / d;
    }
    return r ;
  }
  // =====================================================================
  /** calculate the inverse binomial coefficient 
   */
  inline long double ichoose 
  ( const unsigned short n , 
    const unsigned short k ) 
  {
    //
    if      ( k > n                ) { return 0        ; }
    else if ( 0 == k || n == k     ) { return 1        ; }
    else if ( 1 == k || n == k + 1 ) { return 1.0L / n ; }
    else if ( 2 == k || n == k + 2 ) { return 2.0L / ( 1ull * n * ( n - 1 ) ) ; }
    else if ( n <= 67              ) { return 1.0L / Ostap::Math::choose ( n , k ) ; }
    //
    const unsigned short kk = std::min ( k , (unsigned short) ( n - k ) ) ;
    //
    if ( 120 <= kk )
    {
      const long double logr =
	std::lgamma ( 1.0L +      kk ) +
	std::lgamma ( 1.0L + n  - kk ) -
	std::lgamma ( 1.0L + n       ) ;
      return  std::exp ( logr ) ;
    }
    //
    long double r = 1  ;
    for ( unsigned short i = 1 ; i <= kk  ; ++i ) 
    {
      r *= i ;
      r /= ( n + 1 - i ) ;
      if ( !r ) { return r ; }
    }
    return r ;
  }
  // ======================================================================
  /*  calculate the logarithm of binomial coefficient
   *  \f$ \log C^n_k \f$
   */
  // ======================================================================
  inline long double log_choose
  ( const unsigned short n ,
    const unsigned short k ) 
  {
    if      ( k <= 1 || k >= n ) { return 0 ; } //
    else if ( n <= 67 )          { return std::log ( 1.0L * Ostap::Math::choose ( n , k ) ) ; }
    //
    return 
      std::lgamma ( 1.0L + n     ) - 
      std::lgamma ( 1.0L + k     ) - 
      std::lgamma ( 1.0L + n - k ) ;
  }
  // ======================================================================
  /// evaluate the binomial coefficient as long double C(k,n) = n!/((n-k)!*k!)
  inline long double choose_long_double 
  ( const unsigned short n ,
    const unsigned short k ) 
  {
    //
    // ====================================================================
    //
    if      ( k > n                ) { return 0 ; }
    else if ( 0 == k || n == k     ) { return 1 ; }
    else if ( 0 == k || n == k     ) { return 1 ; }
    else if ( 1 == k || n == k + 1 ) { return n ; }
    else if ( 2 == k || n == k + 2 ) { return 1.0L * n * ( n - 1 )                                     /   2 ; }
    else if ( 3 == k || n == k + 3 ) { return 1.0L * n * ( n - 1 ) * ( n - 2 )                         /   6 ; }
    else if ( 4 == k || n == k + 4 ) { return 1.0L * n * ( n - 1 ) * ( n - 2 ) * ( n - 3 )             /  24 ; }
    else if ( 5 == k || n == k + 5 ) { return 1.0L * n * ( n - 1 ) * ( n - 2 ) * ( n - 3 ) * ( n - 4 ) / 120 ; }
    //
    else if ( n <= Ostap::Math::N_CHOOSE_MAX ) { return 1.0L * Ostap::Math::choose ( n , k ) ; }
    //
    const unsigned  short k1 =  2 * k < n ? k : n - k ;
    //
    long double a = std::lgamma ( 1.0L + n     ) ;
    a            -= std::lgamma ( 1.0L + n - k ) ;
    a            -= std::lgamma ( 1.0L + k     ) ; 
    //
    return std::exp ( a ) ;
  }
  // ==========================================================================
  inline unsigned long long _stirling_ 
  ( const unsigned short n ,
    const unsigned short k ) 
  {
    return  
      0 == n && n == k ? 1 :
      0 == n || 0 == k ? 0 :
      _stirling_ ( n - 1 , k ) * ( n - 1 ) + _stirling_ ( n - 1 , k - 1 ) ;
  }
  // ==========================================================================
  inline long double _stirling_double_ 
  ( const unsigned short n ,
    const unsigned short k ) 
  {
    return  
      0 == n && n == k ? 1 :
      0 == n || 0 == k ? 0 :
      _stirling_double_ ( n - 1 , k ) * ( n - 1 ) + 
      _stirling_double_ ( n - 1 , k - 1 ) ;
  }
  // ==========================================================================
  /// zero for doubles  
  const Ostap::Math::Zero<double> s_zero{}  ;       // zero for doubles
  // ==========================================================================
}
// ============================================================================
/* calculate the binomial coefficient C(n,k) = n!/((n-k)!*k!)
 * the result is exact for all n,k<=67 
 * @warning In case of overflow std::numeric_limits<unsigned long long>::max is returned 
 * @author Vanya BELYAEV Ivan.Belyaev@irep.ru
 * @date 2015-03-08
 */
// ============================================================================
std::uintmax_t
Ostap::Math::choose 
( const unsigned short n , 
  const unsigned short k ) 
{
  //
  if      ( k > n                ) { return 0 ; }
  else if ( 0 == k || n == k     ) { return 1 ; }
  else if ( 1 == k || n == k + 1 ) { return n ; }
  else if ( 2 == k || n == k + 2 ) { return 1ull * n * ( n - 1 )                                     /  2 ; }
  else if ( 3 == k || n == k + 3 ) { return 1ull * n * ( n - 1 ) * ( n - 2 )                         /  6 ; }
  else if ( 4 == k || n == k + 4 ) { return 1ull * n * ( n - 1 ) * ( n - 2 ) * ( n - 3 )             / 24 ; }
  //
  Ostap::Assert ( n <= N_CHOOSE_MAX ,
		  "choose: `n' is too large to fit result into 64-bit integer, use `choose_double'" ,
		  "Ostap::Math::choose" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const unsigned short m  = n < 2 * k ?  ( n - k ) : k ;
  //
  typedef std::pair<unsigned short,unsigned short> KEY    ;
  typedef unsigned long long                       RESULT ;
  typedef std::map<KEY,RESULT>                     MAP    ;
  typedef SyncedCache<MAP>                         CACHE  ;
  /// the cache
  static CACHE                        s_CACHE     {} ; // the cache
  static const std::size_t            s_MAX_CACHE { 2500 } ;
  //
  // 
  const KEY key { n , m } ;
  //
  { // (1) check a value already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // (2) calculate it using Pascal's rule : 
  //
  const RESULT p1 = choose ( n - 1 , m - 1 ) ;
  if ( s_ULLMAX <= p1 ) { return p1 ; }
  //
  const RESULT p2 = choose ( n - 1 , m     ) ;
  if ( s_ULLMAX <= p2 ) { return p2 ; }
  //
  const RESULT result = p1 + p2 ;
  if ( s_ULLMAX <= result ) { return result ; }
  //
  {  // (3) add valid calculated value into the cache 
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear () ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;
  }  
  //
  return result ;
}
// ============================================================================
/*  calculate the inverse binomial coefficient 
 *  \f$ a = C(n,k)^{-1} = \frac{ (n-k)!k!}{n!}\f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2020-01-31
 */
// ============================================================================
double Ostap::Math::ichoose 
( const unsigned short n , 
  const unsigned short k ) 
{
  //
  if      ( k > n                ) { return 0        ; }
  else if ( 0 == k || n == k     ) { return 1        ; }
  else if ( 1 == k || n == k + 1 ) { return 1.0L / n ; }
  else if ( 2 == k || n == k + 2 ) { return 2.0L / ( 1ull * n * ( n - 1 ) ) ; }
  //
  const unsigned short m  = n < 2 * k ?  ( n - k ) : k ;
  //
  typedef std::pair<unsigned short,unsigned short> KEY    ;
  typedef double                                   RESULT ; 
  typedef std::map<KEY,RESULT>                     MAP    ;
  typedef SyncedCache<MAP>                         CACHE  ;
  /// the cache
  static CACHE                        s_CACHE     {} ; // the cache
  static const std::size_t            s_MAX_CACHE { 10000 } ;
  //
  const KEY key { n , m } ;
  //
  { // (1) check a value already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // (2) calculate it!
  const double result = ::ichoose ( n , m ) ;
  //
  { // (3) add calculated value into the cache
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;
  }  
  //
  return result ;
}
// ============================================================================
/** calculate the binomial coefficient C(k,n) = n!/((n-k)!*k!)
 *  @author Vanya BELYAEV Ivan.Belyaev@irep.ru
 *  @date 2015-03-08
 */
// ============================================================================
long double
Ostap::Math::choose_double 
( const unsigned short n , 
  const unsigned short k ) 
{
  //
  if      ( k > n                ) { return 0 ; }
  else if ( 0 == k || n == k     ) { return 1 ; }
  else if ( 1 == k || n == k + 1 ) { return n ; }
  else if ( 2 == k || n == k + 2 ) { return 1.0L * n * ( n - 1 )                                     /   2 ; }
  else if ( 3 == k || n == k + 3 ) { return 1.0L * n * ( n - 1 ) * ( n - 2 )                         /   6 ; }
  else if ( 4 == k || n == k + 4 ) { return 1.0L * n * ( n - 1 ) * ( n - 2 ) * ( n - 3 )             /  24 ; }
  else if ( 5 == k || n == k + 5 ) { return 1.0L * n * ( n - 1 ) * ( n - 2 ) * ( n - 3 ) * ( n - 4 ) / 120 ; }
  //
  else if ( n <= N_CHOOSE_MAX    ) { return choose ( n , k ) ; }
  //
  const unsigned short m  = n < 2 * k ?  ( n - k ) : k ;
  //
  typedef std::pair<unsigned short,unsigned short> KEY    ; 
  typedef long double                              RESULT ;
  typedef std::map<KEY,RESULT>                     MAP    ;
  typedef SyncedCache<MAP>                         CACHE  ;
  /// the cache
  static CACHE                        s_CACHE     {} ; // the cache
  static const std::size_t            s_MAX_CACHE { 10000 } ;
  // 
  const KEY key { n , m } ;
  //
  // 
  { // (1) check a value already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // (2) calculate it!
  const RESULT result = ::choose_long_double ( n , m ) ;
  //
  { // (3) add calculated value into the cache
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;
  }  
  //
  return result ;
}
// ============================================================================
/*  calculate the generalized binomial coefficient C(a,n) 
 *  \f$C(\alpha,k) = \frac{\alpha}{k}\frac{\alpha-1}{k-1}...\f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
double Ostap::Math::gen_choose 
( const double         a ,
  const unsigned short k )
{
  //
  if      ( 0 == k        ) { return 1 ; } 
  else if ( 1 == k        ) { return a ; }
  else if ( s_zero ( a )  ) { return 0 ; }
  //
  long double r = 1 ;
  long double b = a ;
  for ( unsigned short d = k ; 0 < d  ; --d ) 
  {
    // if ( s_zero( b ) ) { return 0 ; }  // RETURN  
    r *= b ;
    r /= d ; 
    b -= 1 ;
  }
  return r ;
}
// ============================================================================
/*  calculate the generalized binomial coefficient C(n/2,k) 
 *  \f$C(n,k) = \frac{n/2}{k}\frac{n/2-1}{k-1}...\f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
double Ostap::Math::choose_half 
( const int            n ,
  const unsigned short k )
{
  if      ( 0 == k              ) { return 1                    ; }
  else if ( 0 < n && 0 == n % 2 ) { return choose ( n * 2 , k ) ; }
  else if ( 1 == k              ) { return 0.5L * n             ; }  // attention! 
  else if ( 0 == n              ) { return 0                    ; } 
  //
  long double r = 1 ;
  int         N = n ;
  for ( unsigned short d = k ; 0 < d  ; --d ) 
  {
    r *= N ;
    r /= d ; 
    N -= 2 ;  // ATTENTION 
  }
  //
  return 
    k < 63 ?
    r / std::pow ( 2L   , k ) : 
    r / std::pow ( 2.0L , k ) ;
}
// ============================================================================
/*  calculate the logarithm of binomial coefficient
 *  \f$ \log C^n_k \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
long double Ostap::Math::log_choose 
( const unsigned short n ,
  const unsigned short k ) 
{
  if      ( k <= 1 || k >= n ) { return 0 ; } //
  //
  const unsigned short m  = n < 2 * k ?  ( n - k ) : k ;
  //
  typedef std::pair<unsigned short,unsigned short> KEY    ;
  typedef long double                              RESULT ;
  typedef std::map<KEY,RESULT>                     MAP    ;
  typedef SyncedCache<MAP>                         CACHE  ;
  /// the cache
  static CACHE                        s_CACHE     {} ; // the cache
  static const std::size_t            s_MAX_CACHE { 10000 } ;
  //
  const KEY key { n , m } ;
  //
  // 
  { // (1) check a value already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // (2) calculate it!
  const RESULT result = ::log_choose ( n , m ) ;
  //
  { // (3) add calculated value into the cache
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;
  }  
  //
  return result ;
}
// ============================================================================
/* calculate unsigned Stirling number of 1st kind 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
std::uintmax_t Ostap::Math::stirling1 
( const unsigned short n ,
  const unsigned short k ) 
{ return _stirling_ ( n , k ) ; }
// ============================================================================
/*  calculate unsigned Stirling number of 1st kind 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
double Ostap::Math::stirling1_double
( const unsigned short n ,
  const unsigned short k ) 
{ return _stirling_double_ ( n , k ) ; }  
// ============================================================================

// ============================================================================
/* Eulerian number A(n,k)
 *  @see https://en.wikipedia.org/wiki/Eulerian_number
 *  @param n   \f$ 0 \le n \f$
 *  @param k   \f$ 0 \le k \le n \f$ 
 *  @return euleria number A(n,k)  
 */
// ==========================================================================
std::uintmax_t
Ostap::Math::eulerian 
( const unsigned short n , 
  const unsigned short k ) 
{
  if      ( !n && !k               ) { return 1 ; }
  else if (  n <=  k               ) { return 0 ; }
  else if (  0 ==  k || n == k + 1 ) { return 1 ; }
  //
  Ostap::Assert ( n <= N_EULERIAN_MAX ,
		  "`n' is too larger to fit result into 64-bit integer, use `eurelian_double'" ,
		  "Ostap::Math::eurelian" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  typedef std::pair<unsigned short,unsigned short> KEY    ;
  typedef std::uintmax_t                           RESULT ;
  typedef std::map<KEY,RESULT>                     MAP    ;
  typedef SyncedCache<MAP>                         CACHE  ;
  /// the cache
  static CACHE                        s_CACHE     {} ; // the cache
  static const std::size_t            s_MAX_CACHE { 10000 } ;
  //
  const KEY key { n , k } ;
  // 
  { // (1) check a value already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // (2) calculate it!
  const RESULT A1     = eulerian ( n - 1 , k - 1 ) ;
  const RESULT A2     = eulerian ( n - 1 , k     ) ;
  const RESULT result = ( n - k ) * A1 + ( k + 1 ) * A2 ;
  //
  { // (3) add calculated value into the cache
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;
  }  
  //
  return result ;
}
// ============================================================================
/*  Eulerian number A(n,k)
 *  @see https://en.wikipedia.org/wiki/Eulerian_number
 *  @param n   \f$ 0 \le n \f$
 *  @param k   \f$ 0 \le k \le n \f$ 
 *  @return euleria number A(n,k)  
 */
// ==========================================================================
double
Ostap::Math::eulerian_double 
( const unsigned short n , 
  const unsigned short k )
{
  if      ( !n && !k               ) { return 1 ; }
  else if (  n <=  k               ) { return 0 ; }
  else if (  0 ==  k || n == k + 1 ) { return 1 ; }
  //
  // ATTENTION! 
  if ( n <= N_EULERIAN_MAX ) { return eulerian ( n , k ) ; }
  //
  typedef std::pair<unsigned short,unsigned short> KEY    ;
  typedef double                                   RESULT ;
  typedef std::map<KEY,RESULT>                     MAP    ;
  typedef SyncedCache<MAP>                         CACHE  ;
  /// the cache
  static CACHE                        s_CACHE     {} ; // the cache
  static const std::size_t            s_MAX_CACHE { 10000 } ;
  //
  const KEY key { n , k } ;
  // 
  { // (1) check a value already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // (2) calculate it!
  const RESULT A1     = eulerian_double ( n - 1 , k - 1 ) ;
  const RESULT A2     = eulerian_double ( n - 1 , k     ) ;
  const RESULT result = ( n - k ) * A1 + ( k + 1 ) * A2 ;
  //
  { // (3) add calculated value into the cache
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;
  }  
  //
  return result ;
}
// ============================================================================


// ============================================================================
/* Eulerian polynomials
 *  \f$  A_n(t) = \sum_k  A(k,k) t^k \f$ 
 *  @see https://en.wikipedia.org/wiki/Eulerian_number
 */
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::Eulerian::Eulerian
( const unsigned short N )
  : m_N { N }
{}
// ============================================================================
// evaluate the Eulerian polynomial 
// ============================================================================
double
Ostap::Math::Eulerian::evaluate
( const double x ) const
{
  //
  if      ( 2 > m_N            ) { return 1     ; }
  else if ( 2 == m_N           ) { return 1 + x ; }
  else if ( 1 < std::abs ( x ) ) { return std::pow ( x , m_N - 1 ) * evaluate ( 1 / x ) ; }
  //
  const std::vector<double>& row = eulerian ( m_N ) ;
  return Ostap::Math::Clenshaw::monomial_sum ( row.begin() , row.end() , x * 1.0L ).first ; 
}
// ============================================================================
// get the derivative
// ============================================================================
double
Ostap::Math::Eulerian::derivative 
( const double x ) const
{
  const std::vector<double>& row = eulerian ( m_N ) ;
  return Ostap::Math::Clenshaw::monomial_sum ( row.begin() , row.end() , x * 1.0L ).second ; 
}
// ============================================================================

// ============================================================================
template <unsigned short N>
using RW = Ostap::Math::EulerianRow_<N>;
// ============================================================================

// ============================================================================
/* Get a row of Eulerian numbers for given N
 *  @see https://en.wikipedia.org/wiki/Eulerian_number
 *  @param n   \f$ 0 \le n \f$
 *  @param k   \f$ 0 \le k \le n \f$ 
 *  @return euleria number A(n,k)  
 */
// ============================================================================
const std::vector<double>&
Ostap::Math::eulerian
( const unsigned short n )
{
  //
  typedef unsigned short         KEY    ;
  typedef std::vector<double>    RESULT ;
  typedef std::map<KEY,RESULT>   MAP    ;
  typedef SyncedCache<MAP>       CACHE  ;
  /// the cache
  static CACHE                   s_CACHE     {}       ; 
  static const std::size_t       s_MAX_CACHE { 1000 } ;
  //
  const KEY key { n } ;
  //
  { // (1) check if a value is already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_CACHE->empty () )
      {
      s_CACHE->insert ( std::make_pair (  0 , RESULT{ 1.0 } ) ) ;
      s_CACHE->insert ( std::make_pair (  1 , RESULT( RW<1>  ::row.begin() , RW<1>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  2 , RESULT( RW<2>  ::row.begin() , RW<2>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  3 , RESULT( RW<3>  ::row.begin() , RW<3>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  4 , RESULT( RW<4>  ::row.begin() , RW<4>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  5 , RESULT( RW<5>  ::row.begin() , RW<5>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  6 , RESULT( RW<6>  ::row.begin() , RW<6>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  7 , RESULT( RW<7>  ::row.begin() , RW<7>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  8 , RESULT( RW<8>  ::row.begin() , RW<8>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair (  9 , RESULT( RW<9>  ::row.begin() , RW<9>  ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 10 , RESULT( RW<10> ::row.begin() , RW<10> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 11 , RESULT( RW<11> ::row.begin() , RW<11> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 12 , RESULT( RW<12> ::row.begin() , RW<12> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 13 , RESULT( RW<13> ::row.begin() , RW<13> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 14 , RESULT( RW<14> ::row.begin() , RW<14> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 15 , RESULT( RW<15> ::row.begin() , RW<15> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 16 , RESULT( RW<16> ::row.begin() , RW<16> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 17 , RESULT( RW<17> ::row.begin() , RW<17> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 18 , RESULT( RW<18> ::row.begin() , RW<18> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 19 , RESULT( RW<19> ::row.begin() , RW<19> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 20 , RESULT( RW<20> ::row.begin() , RW<20> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 21 , RESULT( RW<21> ::row.begin() , RW<21> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 22 , RESULT( RW<22> ::row.begin() , RW<22> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 23 , RESULT( RW<23> ::row.begin() , RW<23> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 24 , RESULT( RW<24> ::row.begin() , RW<24> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 25 , RESULT( RW<25> ::row.begin() , RW<25> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 26 , RESULT( RW<21> ::row.begin() , RW<26> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 27 , RESULT( RW<22> ::row.begin() , RW<27> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 28 , RESULT( RW<23> ::row.begin() , RW<28> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 29 , RESULT( RW<24> ::row.begin() , RW<29> ::row.end() ) ) ) ;
      s_CACHE->insert ( std::make_pair ( 30 , RESULT( RW<25> ::row.begin() , RW<30> ::row.end() ) ) ) ;
    }
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // calculate the previous row:
  const RESULT& prev = eulerian ( n - 1 ) ;
  // current row
  RESULT result ( n , 0.0 ) ;
  // A(N,0) = 1 
  result [ 0    ] = 1 ; // A (N,0)=1 
  result [ n -1 ] = 1 ; // A (N,N-1) = 1 
  //
  for ( unsigned short k = 1 ; k + 1 < n ; ++k )
  { result [ k ] = ( n - k ) * prev [ k - 1 ] + ( k + 1 ) * prev [ k ] ; }
  //
  { // (3) add calculated value into the cache
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    //
    return s_CACHE->insert ( std::make_pair ( key , result ) ).first->second ;
  }  
}
// ============================================================================


// ============================================================================
/** Stieltjes constants 
 *  @see https://en.wikipedia.org/wiki/Stieltjes_constants
 *  @see http://www.plouffe.fr/simon/constants/stieltjesgamma.txt
 */
// ============================================================================
double Ostap::Math::stieltjes 
( const unsigned short N )
{
  static const std::array<long double,N_STIELTJES+1> s_STIELTJES {
.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495146314472498070824809605040144865428362241739976449235362535003337429373377376739427925952582470949160087352039481656708532331517766115286211995015079847937L     ,
-.7281584548367672486058637587490131913773633833433795259900655974140143357151148487808692824484401460407720727888674475946468021218957925934357889512942805222349515279549606111502095344684131064742443502359366559538469614380259381045738025493664557651677028e-1L ,
-.9690363192872318484530386035212529359065806101340749880701365451850755382280414171978197381374537319286223858587946823125346497950962571972526493648255902282233425225376030952454650839481623874159222094228123258327703499005518661324301612485127438367384026e-2L , 
.2053834420303345866160046542753384285715804445410618245481483336913834492112970053570557166228566702972851647965455643681856998286631518491198828868652927507879634329919223749694521937588516086128932900100887582794153232910617743621951263465535646918550222e-2L  , 
.2325370065467300057468170177526068000904469413784850990758040907124841005315521900301678059039306360827843551103225522196456664007682437918861215754578332362815081585973614343456254349981107648311254210559424139592687296840743156065574648306992970376230434e-2L  ,
.7933238173010627017533348774444448307315394045848870757342562698231482118017152023797200635876308162719168578623728950693993277932862610560981556563334528255145467565669705027727800328021356905705583506896973325036121197655863337023278457363824591927681467e-3L  , 
-.2387693454301996098724218419080042777837151563580786314764253073910675599929638714368611141285111024780673075909533038020656912269595124699717764644539646672994359920016312185594700697656216552224200467213493615215357142501593344546382399161688323344633459e-3L ,
-.5272895670577510460740975054788582819962534729698953310134042268856827324651411821440413807979996094255347137014747313651473070509002813384872173815912447148812297961891402918386808326741994669807197014678473723857395756624902062083941316454527809794259124e-3L ,
-.3521233538030395096020521650012087417291805337923503566573315073642817765060653010801409187200115151480243957720152737286032860943501971136758373996071542419096556138574571993689232705602040409831370056942984921502852878319887965827056475310881183605291025e-3L ,
-.3439477441808804817791462379822739062078953859444162975929190484315010334446152837095754389345718990812939708751106559999723650791217087757476384733602840643363955402889741500455608351851986205711990705391789569312430556033918346536209762876458392962080569e-4L ,
.2053328149090647946837222892370653029598537741667643038402087143530090240710691751984960510609028168654630716692437652938609455701636670314763888351107401248788037757774478758676467051586842465830033648294891261556376188509817851906503013130840649822003101e-3L  ,
.2701844395439035266729020820679556738278420586884025039737358031367999909642929802378444607202380543341892065433898331362430539345163022251660977905645827208257807545662673885178807828093358372485533868496985893326451504654464211037386959473080948618383043e-3L  ,
.1672729121051401933535015433411834466078066328055658280477909376512195970327407625539042834685454860827918008130473096673223305416815315435148622550674827969874566633280326548214334209784161729704490822545521988490380808808815410064684131411221199776802669e-3L  ,
-.2746380660376015886000760369335518152678533767039553609283308916757051860700887286766227018014141226063021263616344526017942940219226384963138784380546201092996086073554616887396220116191432159038003288891331098518783382994530877461169787303572031210085377e-4L ,
-.2092092620592999458371396973445849578315442115060695624342083257187577618413479238849643630363170068978952103732270941025164202162406217483524889032192155723532675386434812392932022134837802510909806951284079809204130775034571304740902331856694608730746506e-3L , 
-.2834686553202414466429344749971269770687029807176752539699432929676256905331671583502753305306099359983296641743658578448670204816395094978958650260714071371637925247680241849293162055953372734373236183678636434144439111527617109957199974180603721462395347e-3L , 
-.1996968583089697747077845632032403919157649740340612798596671625543805947413858376185052901045229228152050749963179893268016715236270037554630588742666500589733359678401153567212058999564688677646468080389369421514423787946034894284978610156860905576783237e-3L ,
.2627703710991833669946659763051012281607869292911406079711751835228318283659619892378138323561685615911623560201540792781969494575153922147852554290415989679216348097499078673743176921693516030034263511744512517839540600579140351515313919384582936117575383e-4L  ,
.3073684081492528265927547519486256455238112907314616910811036523148083902097289937656175556967012404059095414514102941951119674092289417995555858007852262090239169157743457758869727907667465462418873038971056792798849296000885600799035155113134156928909482e-3L  ,
.5036054530473556290555964377171600353212698076494978373237909270104380947646223189527024912899340799999802379685872863039429327259976022887863525558688062335605764854361061913564269294911630304700024202875309737842304685438543203315544401816345818487684388e-3L  ,
.4663435615115594494005948244335505251131434739256889976707266280985445821300329007010634781743483037676276776893978664788490376098061008113591441797256336823512267157959337408681128870085117167176652714745827612534802328769644422488731525593967185587090795e-3L  ,
.1044377697560001158107956743677204910444282507055467478343714867390804411994132201971857629970557612477080276816174297098442819046793799873530097156786216211170154344305780995660062733094223612700160596820798224374455841978660825820242682142793915007842212e-3L  ,
-.5415995822039977016551961731741055845438609287007488018391913163842120274727846608325326997754074546503753106668640607229437483993682169764132200837068069923525151709033207834551632712814364366899723899724930759778803264577711209436089086578419457274846364e-3L ,
-.1243962090408245779299741599537165809147028113964637716532971108378003227300447476688043928770660618906589245188584132576689392153326370316493290511316121734714766939840656684143340621090548028988852205331608344386554756230080950750377804966175863406283901e-2L ,
-.1588511278903561561906196611521115857318722822144129067478194125480950383379748723977182688294225715704915149922400750354649555457261240500959167389177194357862621241155775647357272837529119488611489255817396821745210322189035753496755771351316733341251102e-2L ,
-.1074591952738488824724291987353173089273979331453170361409902581781368852879725083510487626853780543878731137936481194598448964551855691857004868428077492623871354826190616571737826265422531415866113138027144599171360216617932498664012392456911554072571146e-2L ,
.6568035186371544315047730033562152488860650604775373760992825250891170965034776663910579342625757859715232062111062915856525206876176042963292134187616474914067354036499371011291177322406733034202470236023762935497158074335393463790043477652282574744043603e-3L  ,
.3477836913618538209007359574258811547662915663885919292269369089163601152516386090184507808325259934643247527620919916504590034421507603998744879095106816233978110202439753349023737964681719231706674211339160697022354031902450736661875203742977403043595923e-2L  ,
.6400068531700629458107228221945863666637198144588475220590330625837093213897134487336264083586702624780911771466881249158837079615155947036988364224391544873825936603455801551951739436482654980238076705978473654568449071615955207841466403251011797640636588e-2L  ,
.7371151770472239134412402423559402157841327488512840153180311593260094964080124902874271317941480450169928391633082168004655647838467136598965148773299440893911326321952210876313787678510143647647900141023350228313526634026921869524825211534867202695510451e-2L  ,
.3557728855573160947913537748908402610809650649522125076131381742604895375347768268411760486644490119698820111281432172758850410765212613111042978612558180006655129918249475202003030282427138985184474816070832944619701532581247303970119806166784524080913113e-2L  ,
-.7513325997815228933135160081576145616636587418058437829820653897726656851364608736093406575498735199953375557645227668649877454170501745904952874352569833828340361752923325058378630036421872568935701095205574194985993332448759010555982110176282553449891202e-2L ,
-.2570372910842040179348788378034991655408420213709272996040553488073373925092915194869402974674313960730940827885805986079503641578184306369509419528819255356778241434806326040117004842490643574687237090396478809065180424485158191278852413572924082385053692e-1L ,
-.4510673410808021990498284969956376278181476720118751914813317935089473799365486780079426205395715048231336142622172311472080346051129198009798869840229135326813527696281135231523051833221278728337033589841096529250915588639261001568781069882956884511058368e-1L ,
-.5112692802150846442507582003801215757213099693370947031283705080222204757165161391299384814318941989350085938186440763932854104808109240365643509465403766309112448633408990761899675968647705676410731088908908770878621830867075270142641783507042859178782092e-1L ,
-.2037304360386131270575189730254735763855514540127568404226904414244063724030988503350463213893919466606979495327877395077626203087646069193281918128083913334868246722812020110757918056153158263567517772557098298030827761720020741589113319214710546324868855e-1L ,
.7248215881681133373380044422036406713494156746516600301109791139709584630115609076184768490496514066910360733891536130570598692290165869195273978736246742908618723458316488977615459103966319968155863913816498268488587418654995797804155697655649958847357895e-1L  ,
.2360263822743015027209817621993795633447977318546636606795696524336922591383506319852601197035035699203623244114544400619644905611166062223529956334728610781278059048805132264984375254618970256448247285280843274150107341682480156517967257924219186788145330L     ,
.4289634463848091527368615465389604104951855560836298508881511467174782165406052598052578837703479513149605271888650666596964536012997065102084281647972845368239473903765489789230269010448415666019284057957391635619051037937314897060412487116366199455816836L     ,
.5179218426929237189788930575162080963032869116382440388334761465713068965780043079965859932817502694727319695029057167172255798781972489650526343539907081526696209591694159811807914277744490491038137744209287798148227736083094114236529502428204747874931539L     ,
.2487215593946154650844919104403834212137661679529456385350509654042314398657765133858509357834519851660462686780038493302932372007287500470374521580805939391035285372957104115822097373824399880947615656580184641803366839892091737835298438895206210595752375L     ,
-.7195748469013003506888739112196854239261748831904649437785449263509786356293101897746613034720604271142999671606071191146849003123224065620009891082439621461705375083652807042878644516387056361921714149226264455652409110895061126053434853106969343380563178L    ,
-2.638794927335734535788281675648811310351782987035211499393243485907022239259505477040826570395791839498801399715723709634295474856557592397489594583251684015695249594247301673922474234026168332480339091413459364363081675429428395332081292613817734336565494L    ,
-5.264930312355023828811032859580336606977042993956806905911825107752096822582823647571129789899774007449077183706697896011698161033442136588696911829176438099711663253191742572376731229116446238288895074259071273499782678262582441265329247305947655986144330L    ,
-7.188745889503527282342094824577664215932574859281338805833810538667337501527011967915372049590948156198532169549971970898863800579381929252662783608393442431305718388587696952733141648850458812574971021050868587199490079795392192132425737225940758618726619L    ,
-5.072344589916372492298940404798877527548898547920341433576114007970208835005861247222957450067358896625318542197749232212662512002598509927706415764052961176009738698115190727795940434449848894227264093298941761645790234028006626255339596450726680973432486L    ,
6.609915609096965813839975106591453976121624142676258962466161263230138247955950680501240606225730339947058738765837375315012758818853754069379249379465003451583201983190721503025985727108804064481134263041435677655221268593105524179459213774579064542886341L     ,
34.03977498215874824766115211222887543481585062040728845880766195415091850111097098102234038334063582255146136304735662018166909676792633649172593400903096852662977956340472884566589366552943189154254333330735508711818454888919225114589304435439875918866231L     ,
78.68247976324258495603848420938847587592859021606368068347311248131245045816477978354222279127347778484079585199133838207604437335391341578617869858249504960215493811609322620525071575891947562789820590429774578807866291149245067581168982559771813175765951L     ,
125.8443876319784690933640869642634548210830426469704980639576282947523289773055415013393636381818820365174280709651417996658374631556841875987387400798841785568945454815405004576022607269753629090377417076709892483708691350430892879564872582173161084734278L     ,
126.8236026513227165967252536486575555384835759448901701980956993469732544041271419075472839202481342000655527466907773263202970132939623809589591513912865283594558393685656496366403863160270130724104654611609335672863609966882450434670623623037390894982865L     ,
-19.19691187302785580049992278889461000488617138521125094955804692768632497077431894010336425512641087536957362345904530571260115521167316919248856770615789826449707900700147413646023179865419434659545981369663884566897807702752807188599160262261877584723048L    ,
-463.1889230267168108085343582392362931366116957180006901624675536855778259610180865548529674806333795445784720929015557429014573440317655361459554522755159806774751854061019646746026001962717761365205206395869989959129560170027580385380148222072396574437614L    ,
-1340.659144376892189727101963668908904813680425823094277360983924970123545905951320804164099765412339187206633285729948584751525184352621143592938185171889188832432976967119173814211291335944639480016298975541814145863037680692144453474530167024849384966632L    ,
-2572.454740404435516763574381743194537144349074129635190025656675261689968560135282533680795608528937664815193978544465137589608940353479663547865399362584761794796720486829475401161510062152263672635951091948051219547961625385885682988539181752335571182500L    ,
-3457.141208645389953809462307608548197267392571262945382531901198419835265449139114139752724233485244785535339119102051777751947298714580782741473753259866854277965873748702018770556207132668621431246868940695237459978570056725932042526772704517088887009050L    ,
-2055.275816231974304572825156345659341714378586681402878963865132625780739633561673788853189681295671720531311167392431020178456949604170241477352885604106843954614982380859102212646287075146135017827460000558232950054762819030973497931729190368748362000642L    , 
5372.282213203191284181833545959014329555896642905224954409370153194551639289251004905214862802153702614891093139187266989686294985772286939348850809250994703278188368030520223602114426496922514344077289177925193404535132017379204347643038237561684778349535L     ,
24019.38937760698818204861035985736409058929002953145915799769330693343585546963061640387520046045702663578653894052569085070940959579119305407994509817271201272546721845612222693895881739318046735749910581720490135833336177178720592740120790968916641082098L     ,
57424.31929696407550657104394311255386400824370256752620991647966953367902572223109625676082012930154382556432735615395354528456397119780281123421908197993868853546446676051766559695367292479004717927284373185020644824412329349572910288157946256744189407140L     , 
98543.25459014604211093211428941711531650104882341554197930622783608489036287602715386941804646663591106764763550763071863725791968111497783011529453607923260637137817777048517149527447460321197286999081445397775020004411932959031469803884925185616975390011L     ,
111670.9578149410793387893403211563974853622621861860070492558245106321152197978896356784986709659246398447028682342515620904866717282773920131267534144689106562607331771415956711044043588809723539995371146389025477672830672938535201055031258424853698622521L     ,
5333.665210500764343403713072261136718914593293565738489694853347359352589030809494678259291885891979613064318888829142263676324267800058016249399128627054877573090524750544029388489126257672223385886733854495868967399358045126339726416940001003776331780669L     ,
-390972.6873133963957775428164347301012865910198508746441638769290837847013047456217619412022961983123163124272712579711729696456704271694618158088390339217519866909714445427931569604142174394070839878737071843151424174567821468441597847504269831119565319909L    ,
-1303180.712532519808573578260908777678044948332570067129963421421848029593041641458697133568977980576700883270277981455463933985088650175081474094621803307108859599115906956541023262677524307005797980356671824629484862625887654008181910920206214383524051245L    ,
-2845076.552608612116523243100237370870167655675924216498192869458404187329953640177299560183025897558495297082006826021694564155486292072145777913453704712722953138646764731475253810351084080951541575328521843225182807547286142617638437012547115115046403931L    ,
-4540526.609737724100509220867003971251048957338640918478864644401431771260187110948708426225730696673484910608169992239039359208562183905400100070421973273676747789972100524122177059432652884762115926751661713705281886452057166399772740063853326419021621033L    ,
-4341905.139001516198415670264383693842260958150127609160422048895643211404947195718257471803919996372167852148973302363686896828779612377060187735478109745765995707365214973502468316183453023901759425272748017178186339655775494652737519189536753913543805680L    ,
2871566.945972460509905357306147525690106604257063853632662792335691123095953849148759187546180165016153744782976866579448776215707008569628936807475513691359949909887658409137518529362425919225421409999326669652492263291715447395177074780931961343784615131L     ,
26604908.54668677092425521912451936480795299150646426379923783074812038437400291966231913301672803147593025933983520666997287407005167498500146308540816725205740171309465570488989279619206587140609236965698637546585522770792996633444924580867051639968204913L     ,
79321663.11929906059211001875506986109078678605929299821733167638956228811094702496866006501239037925742409911928107256247873459490010385482553313165643858925233075421940693131005923399599328970140465301608533561129008200715549168186165980180247842021792257L     ,
166215134.0468254355839825083764775026922493156873793437567024004040720676074041053803255270535873330605351293137839271608605683647864515142266769019278900135634372572577287785416673156462232216791854730153546523662690892362166893419994060430067089918268485L     ,
255153258.3082389934909848078277768987902985514378076503823427021825601645614601862599207385980737560906754309716204100869275599148675590710569782697323107792348115172035957067240704301036234325777779946135648388560892047021038514968620844931091344338610286L     ,
212655631.6918540369512019248627158396680614022195978324871483336108798420688927107279671233050594954058347952103381180030465863386349837153315936212092445328817755413706057708197992558132570343461837062889786366607647313628680221334335862993688051677647064L     ,
-298767089.4311661838592139845456232588222876086992170490259044012663702238324435706623905293153630798969724334390118682919837048033136392454203199753628085167331437818516158612262299723479084451715092046370319754099231070567632457508936168937828037690257857L    ,
-1919487427.732802396229627757729327898023244830602085729304679210428118470905660169503320644717346154413202421403219588184732013410426970230132448744061230223078726672201628641965688174011380164991661246163863144544756484425183731247748949286364250669027383L    ,
-5515574258.129220269427426037306452066936170169171072411797916588895750116790116120042606241013411261739570410392314324043817608753827749236198030823299522969460857196981095390178299797003363459101202590130435740820533907989476579686558014505421730147475846L    , 
-11483450987.92625603025962308750354819358989636858236343418922095146003529652099454932564500092640974706164032396286399805090630030297913074096080965875780082646299473016444439508070872009844923178137154154440180991389231196962826618318021120909673249402227L    ,
-17570152277.77726386163608615618525785039204899986664599956286229502140896936467083284521308786084527317826035812378168370902350943468812664885992884478165346500022035578688139196146186122588937724007564250471442030547825810788316951289967818338033762342923L    
} ;


  Ostap::Assert ( N < s_STIELTJES.size ()  , 
                  "Index for Stielties's constant is too large" , 
                  "Ostap::Math::stieltjes" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;

  return s_STIELTJES [ N ] ; 
} 

// ============================================================================
// The END 
// ============================================================================
