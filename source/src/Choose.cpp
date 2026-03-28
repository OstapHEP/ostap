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
    const unsigned  short k1 = 2*k < n ? k : n - k ;
    if ( k1 * std::log2 ( M_E * n / k1 ) < 63 ) 
      { return std::log ( 1.0L * Ostap::Math::choose ( n , k ) ) ; }
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
    if      ( k > n            ) { return 0 ; }
    else if ( 0 == k || n == k ) { return 1 ; }
    else if ( n <= 67          ) { return 1.0L * Ostap::Math::choose ( n , k ) ; }
    //
    const unsigned  short k1 =  2 * k < n ? k : n - k ;
    if ( k1 * std::log2 ( M_E * n / k1 ) < 63 )
      { return Ostap::Math::choose ( n , k1 ) ; }
    //
    long double a = std::lgamma ( ( long double )     n + 1 )  ;
    if ( a < s_emax ) { return Ostap::Math::choose ( n , k1 ) ; }
    a            -= std::lgamma ( ( long double ) n - k + 1 ) ;
    if ( a < s_emax ) { return Ostap::Math::choose ( n , k1 ) ; }
    a            -= std::lgamma ( ( long double )     k + 1 ) ;
    if ( a < s_emax ) { return Ostap::Math::choose ( n , k1 ) ; }
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
  else if ( 2 == k || n == k + 2 ) { return 1ull * n * ( n - 1 ) / 2 ; }
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
  if      ( k > n            ) { return 0 ; }
  else if ( 0 == k || n == k ) { return 1 ; }
  else if ( n <= 67          ) { return Ostap::Math::choose ( n , k ) ; }
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
double Ostap::Math::log_choose 
( const unsigned short n ,
  const unsigned short k ) 
{
  if      ( k <= 1 || k >= n ) { return 0 ; } //
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
  else if (  1 ==  k || n == k + 1 ) { return 1 ; }
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
/* Eulerian number A(n,k)
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
  else if (  1 ==  k || n == k + 1 ) { return 1 ; }
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
/* Get a row of Eulerian numbers for given N
 *  @see https://en.wikipedia.org/wiki/Eulerian_number
 *  @param n   \f$ 0 \le n \f$
 *  @param k   \f$ 0 \le k \le n \f$ 
 *  @return euleria number A(n,k)  
 */
// ============================================================================
/**
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
  static CACHE                   s_CACHE = { { 0 , { 1.0 } } } ; // the cache
  static const std::size_t       s_MAX_CACHE { 1000 } ;
  //
  const KEY key { n } ;
  //
  { // (1) check a value already calculated
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  // calculate the row 
  RESULT result ; result.reserve ( n + 1 ) ;
  //
  //
  { // (3) add calculated value into the cache
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;    
  }  
  //
  return result ;
}
*/

// ============================================================================
// The END 
// ============================================================================
