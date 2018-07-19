// $Id$ 
// ============================================================================
// Include files
// ============================================================================
// STD/STL 
// ============================================================================
#include <cmath>
#include <climits>
#include <algorithm>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Clenshaw.h"
#include "Ostap/Choose.h"
#include "Ostap/Math.h"
// ============================================================================
/** @file 
 *  Cacualte binbomial coefficients and related quantities 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  typedef std::numeric_limits<unsigned long long> ULLTYPE ;
  static_assert( ULLTYPE::is_specialized, 
                 "numeric_limits<unsigned long long> is not specialzaed!" ) ;
  // ==========================================================================  
  const unsigned long long s_ullmax = ULLTYPE::max () ;
  const long double        s_emax   = std::log ( 0.2L * s_ullmax )  ;
  const unsigned short     s_digits = ULLTYPE::digits - 2         ;
  // ==========================================================================
  /** calculate the binomial coefficient C(k,n) = n!/((n-k)!*k!)
   *  In case of overflow std::numeric_limits<unsigned long long>::max is returned 
   */
  inline unsigned long long  
  _choose_ ( unsigned short n , unsigned short k ) 
  {
    //
    if      ( k > n            ) { return 0 ; }
    else if ( 0 == k || n == k ) { return 1 ; }
    //
    k = std::min ( k , (unsigned short) ( n - k ) ) ;
    unsigned long long r = 1  ;
    for ( unsigned short d = 1 ; d <= k  ; ++d , --n ) 
    {
      if ( r > s_ullmax / n * d ) { return s_ullmax ; }  //  RETURN
      // r *= n ;
      // r /= d ; 
      r = ( r / d ) * n + ( r % d ) * n / d;
    }
    return r ;
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
unsigned long long  
Ostap::Math::choose 
( const unsigned short n , 
  const unsigned short k ) { return _choose_ ( n , k ) ; }
// ============================================================================
/** calculate the binomial coefficient C(k,n) = n!/((n-k)!*k!)
 *  @author Vanya BELYAEV Ivan.Belyaev@irep.ru
 *  @date 2015-03-08
 */
// ============================================================================
double Ostap::Math::choose_double 
( const unsigned short n ,
  const unsigned short k ) 
{
  //
  if      ( k > n            ) { return 0 ; }
  else if ( 0 == k || n == k ) { return 1 ; }
  else if ( n < s_digits     ) { return _choose_ ( n , k ) ; }
  else if ( n <= 67          ) { return _choose_ ( n , k ) ; }
  //
  const unsigned  short k1 = 2*k < n ? k : n - k ;
  if ( k1 * std::log2 ( M_E * n / k1 ) < 63 ) 
  { return _choose_ ( n , k ) ; }
  // 
  long double a = std::lgamma ( (long double)     n + 1 )  ;
  if ( a < s_emax ) { return _choose_ ( n , k ) ; }
  a            -= std::lgamma ( (long double) n - k + 1 ) ;
  if ( a < s_emax ) { return _choose_ ( n , k ) ; }
  a            -= std::lgamma ( (long double)     k + 1 ) ;
  if ( a < s_emax ) { return _choose_ ( n , k ) ; }
  //
  return std::exp ( a ) ;
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
  else if ( 1 == k              ) { return 0.5 * n              ; }  // attention! 
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
  else if ( n <= 67 )          { return std::log( (long double) _choose_ ( n , k ) ) ; }
  //
  const unsigned  short k1 = 2*k < n ? k : n - k ;
  if ( k1 * std::log2 ( M_E * n / k1 ) < 63 ) 
  { return std::log( (long double) _choose_ ( n , k ) ) ; }
  // 
  return 
    std::lgamma ( (long double) ( n     + 1 ) ) - 
    std::lgamma ( (long double) ( k     + 1 ) ) - 
    std::lgamma ( (long double) ( n - k + 1 ) ) ;
}
// ============================================================================
/* calculate unsigned Stirling number of 1st kind 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
unsigned long long Ostap::Math::stirling1 ( const unsigned short n ,
                                            const unsigned short k ) 
{ return _stirling_ ( n , k ) ; }
// ============================================================================
/*  calculate unsigned Stirling number of 1st kind 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-08
 */
// ============================================================================
double Ostap::Math::stirling1_double ( const unsigned short n ,
                                       const unsigned short k ) 
{ return _stirling_double_ ( n , k ) ; }  
// ========================================================================



// ============================================================================
// The END 
// ============================================================================
