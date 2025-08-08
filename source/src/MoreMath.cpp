// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <array>
#include <climits>
#include <complex>
#include <cassert>
#include <algorithm>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_math.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_airy.h"
#include "gsl/gsl_sf_fermi_dirac.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_log.h"
#include "gsl/gsl_sf_psi.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_ellint.h"
#include "gsl/gsl_sf_lambert.h"
#include "gsl/gsl_sf_dilog.h"
#include "gsl/gsl_sf_zeta.h"
#include "gsl/gsl_sf_clausen.h"
// ============================================================================
// LHCbMath
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/MakeArray.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Clausen.h"
#include "Ostap/Choose.h"
#include "Ostap/Interpolants.h"
#include "Ostap/Bernstein.h"
#include "Ostap/ChebyshevApproximation.h"
// ============================================================================
// Local
// ============================================================================
#include "GSL_sentry.h"
#include "Faddeeva.hh"
#include "gauss.h"
#include "Ostap/Workspace.h"
#include "local_math.h"
#include "local_gsl.h"
#include "Integrator1D.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  implementation file for function from file LHCbMath/MoreFunctions.h
 *  @see LHCbMath/MoreFunctions.h
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<float>::is_specialized            , 
                  "mumeric_limits are not specialized for flota  "      ) ;
  static_assert ( std::numeric_limits<double>::is_specialized           , 
                  "mumeric_limits are not specialized for doubles"      ) ;
  static_assert ( std::numeric_limits<long double>::is_specialized      , 
                  "mumeric_limits are not specialized for long doubles" ) ;
  // ==========================================================================
  /// equality criteria for doubles
  /// const Ostap::Math::Equal_To<double> s_equal{} ;       // equality criteria for doubles
  /// zero for doubles  
  /// const Ostap::Math::Zero<double>     s_zero {} ;       // zero for doubles
  /// "almost infinity" 
  const long double s_infinity = 0.9 * std::numeric_limits<double>::max () ;
  /// "relatively large value" 
  const int         s_imax     = -32   ;
  const long double s_large    = std::pow ( 2.0 , -s_imax ) ;
  // epsilon
  const long double s_epsilon  =  
    std::min ( (long double) std::numeric_limits<double>     ::epsilon () , 
               10000       * std::numeric_limits<long double>::epsilon () ) ;
  // ==========================================================================
  /// get a factorial 
  inline long double 
  _factorial_d_ ( const unsigned short N ) 
  {
    return 
      0 == N ?  1 : 
      1 == N ?  1 : 
      2 == N ?  2 : 
      3 == N ?  6 : 
      4 == N ? 24 : N * _factorial_d_ ( N - 1 ) ;
  }
  // ==========================================================================
}
// ============================================================================
namespace
{
  // ==========================================================================
  /** sum of N-terms in the exponential expansion 
   *  \f$ f = \sum_{i=0}^{N} \frac{x^k}{k!}\f$
   *  @param x  INPUT the argument 
   *  @param N  INPUT N-terms to be used 
   *  @return partial exponential sum 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-03-26
   */
  inline long double _exp_N_ 
  ( const long double    x , 
    const unsigned short N ) 
  {
    long double r = 1 ;
    long double t = 1 ;
    for ( unsigned short n = 1 ; n <= N ; ++n ) 
    {
      t *= x ;
      t /= n ;
      r += t ;
      if ( r > s_infinity ) { return s_infinity ; }  // RETURN 
    }
    return r ;
  }
  // ==========================================================================
}
// ============================================================================
/*  sum of N-terms in the exponential expansion 
 *  \f$ f = \sum_{i=0}^{N} \frac{x^k}{k!}\f$
 *  @param x  INPUT the argument 
 *  @param N  INPUT N-terms to be used 
 *  @return partial expoenntial sum 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */
// ============================================================================
double Ostap::Math::exp_N ( const double x , const unsigned short N ) 
{
  return 
    0 == N       ? 1 :
    1 == N       ? 1 + x : 
    2 == N       ? 1 + x * ( 1 + x *   0.5 ) : 
    3 == N       ? 1 + x * ( 1 + x * ( 0.5 + x / 6.0 ) ) : 
    s_zero ( x ) ? 1. : _exp_N_ ( x , N ) ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  /** calculate "relative or reduced exponent"
   *  \f$f(x) = N! ( e^{x} - \sum_{k=0}^{N-1} \frac{x^k}{k!})/x^{N} \f$ 
   *  @param x  INPUT the argument 
   *  @param N  INPUT N-terms to be used 
   *  @return the value of "reduced exponent"
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-03-26
   */
  inline long double _exp_rel_N_  
  ( const long double  x , 
    const unsigned int N ) 
  {
    //
    // 1. calculate regular part 
    //
    long double hm1 = 1  ; 
    long double km1 = 0  ;
    long double h0  = 0  ;
    long double k0  = 1  ;
    long double hp1 = h0 ;
    long double kp1 = k0 ;
    for ( unsigned long n = 2 ; n <= 100000 ; ++n ) 
    {
      const long double   an = ( 0 == n % 2 ?  x * n / 2 :  -x*( N + ( n - 1 ) / 2 ) ) ;
      const unsigned long bn = n + N ;
      //
      hp1 = bn * h0 + an * hm1 ;
      kp1 = bn * k0 + an * km1 ;  
      //
      hm1 = h0  ;
      km1 = k0  ;
      h0  = hp1 ;
      k0  = kp1 ;
      //
      if ( std::abs ( hp1 ) > s_large || std::abs ( kp1 ) > s_large ) 
      {
        //
        h0  = std::ldexp ( h0  , s_imax ) ;
        k0  = std::ldexp ( k0  , s_imax ) ;
        hm1 = std::ldexp ( hm1 , s_imax ) ;
        km1 = std::ldexp ( km1 , s_imax ) ;
      }
      // time-to-time check the convergency 
      if ( 0 == n % 5 ) 
      {
        const long double delta  =  ( hm1 / km1 ) / ( h0 / k0 ) - 1 ;
        if ( std::abs ( delta ) <= 2 * s_epsilon ) { break ; }
      }
    }
    const long double result = h0/k0 ;
    //
    // add the "irregular part"
    //
    return 1 / ( 1 - x / ( N + 1 + result )  ) ;
  }
  // ==========================================================================
}
// ============================================================================
/* "relative or reduced exponent"
 *  \f$f(x) = N! ( e^{x} - \sum_{k=0}^{N-1} \frac{x^k}{k!})/x^{N} \f$ 
 *  @param x  INPUT the argument 
 *  @param N  INPUT N-terms to be used 
 *  @return the value of "reduced exponent"
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */
// ============================================================================
double Ostap::Math::exp_rel_N ( const double x , const unsigned short N ) 
{
  // todo: switch to GSL implementation??
  const long double y = x ;
  return 
    0 == N ? std::exp ( y ) :
    1 == N ? exprel   ( x ) :
    s_zero ( x ) ? 1 : 
    _exp_rel_N_ ( y , N ) ;  
}
// ============================================================================
/*  compute \f$ f(x) = \frac{e^x-1}{x}\f$
 *  @return the value of psi function 
 *  @see exp_rel_N 
 */        
// ============================================================================
double Ostap::Math::exprel ( const double x ) 
{
  //
  const long double y = x ;
  return
    x < GSL_LOG_DBL_MIN ? -1.0 / y             : // RETURN
    x > GSL_LOG_DBL_MAX ? s_infinity           : // RETURN 
    1 > std::abs ( x )  ? std::expm1 ( y ) / y : // RETURN 
    ( std::exp ( y ) - 1 ) / y ;
}
// ============================================================================
namespace
{
  // ==========================================================================
  /** regularized incomplete gamma function 
   *  \f$ \gamma^{\ast}(a,x) = \frac{x^{-a}}{\Gamma(a)\gamma(a,x) }\f$, 
   *  where 
   *  \f$ \gamma(a,x) = \Gamma(a) - \Gamma(a,x) = \int_0^x e^{-t}t^{a-1}dt\f$, 
   *  @param a INPUT a-parameter 
   *  @param x INPUT x-argument 
   *  @return the value of regularized incomplete gamma function 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-03-27
   */
  inline long double _gamma_star_ 
  ( const long double a , 
    const long double x )  
  {
    //
    long double t  = 1     ; 
    long double r  = t / a ;
    //
    for ( unsigned long n = 1 ; n <1000000; ++n ) 
    {
      //
      t *= -x ; 
      t /=  n ;
      //
      if ( 0 == a + n ) { break ; }   // BREAK  :-( 
      r += t / ( a + n ) ;
      //
      if      ( std::abs ( t ) <= 2 * s_epsilon ) { break ; } // BREAK 
    }
    //
    long double Ga = std::tgamma ( a ) ;
    r /= Ga ;
    //
    return 
      r < -s_infinity ? -s_infinity :
      r >  s_infinity ?  s_infinity : r ;
  }
  // ==========================================================================
  /** regularized incomplete gamma function 
   *  \f$ \gamma^{\ast}(a,x) = \frac{x^{-a}}{\Gamma(a)\gamma(a,x) }\f$, 
   *  where 
   *  \f$ \gamma(a,x) = \Gamma(a) - \Gamma(a,x) = \int_0^x e^{-t}t^{a-1}dt\f$, 
   *  @param a INPUT a-parameter 
   *  @param x INPUT x-argument 
   *  @return the value of regularized incomplete gamma function 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-03-27
   */
  inline long double _gamma_star_2_ 
  ( const long double a , 
    const long double x )  
  {
    //
    long double t = 1 / a ; 
    long double r = t     ;
    //
    for ( unsigned long n = 1 ; n <1000000; ++n ) 
    {
      //
      if ( 0 == a + n ) { break ; }   // BREAK :-( 
      t *=   x       ;
      t /= ( a + n ) ;
      //
      r += t ;
      if ( std::abs ( t ) <= 2 * s_epsilon ) { break ; } // BREAK 
    }
    //
    return r * std::exp ( -x ) / std::tgamma ( a ) ;
    //
  }
  // ==========================================================================
}
// ============================================================================
/** regularized incomplete gamma function 
 *  \f$ \gamma^{\ast}(a,x) = \frac{x^{-a}}{\Gamma(a)\gamma(a,x) }\f$, 
 *  where 
 *  \f$ \gamma(a,x) = \Gamma(a) - \Gamma(a,x) = \int_0^x e^{-t}t^{a-1}dt\f$, 
 *  @param a INPUT a-parameter 
 *  @param x INPUT x-argument 
 *  @return the value of regularized incomplete gamma function 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-27
 */
// ============================================================================
double Ostap::Math::gamma_star ( const double a , const double x ) 
{
  //
  if ( Ostap::Math::isint ( a ) || std::abs ( a - Ostap::Math::round ( a ) ) < 1.e-4 ) 
  {
    const int         n = Ostap::Math::round ( a ) ;
    const long double y = x ;
    if ( n <= 0 ) { return 1./ std::pow ( y , std::abs ( n ) ) ; }
  }
  //
  // const double r3 = 
  //   0 < x && 0 < a ? gsl_sf_gamma_inc_P ( a , x ) * std::pow ( x , -a ) : -1 ;
  //
  return
    1.1 < x ? _gamma_star_2_ ( a , x ) :  _gamma_star_   ( a , x ) ;    
}
// ============================================================================
/** regularized incomplete gamma function 
 *  \f$ \gamma^{\ast}(n,x) = \frac{x^{-n}}{\Gamma(n)\gamma(n,x) }\f$, 
 *  where 
 *  \f$ \gamma(n,x) = \Gamma(n) - \Gamma(n,x) = \int_0^x e^{-t}t^{n-1}dt\f$, 
 *  @param n INPUT n-parameter 
 *  @param x INPUT x-argument 
 *  @return the value of regularized incomplete gamma function 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-27
 */
// ============================================================================
double Ostap::Math::gamma_star ( const int n , const double x ) 
{
  const long double y = x ;
  return 
    n   <= 0 ? 1./ std::pow ( y , std::abs ( n ) ) : 
    1.1 <  x ? _gamma_star_2_ ( n ,  y ) : _gamma_star_  ( n , y ) ;  
}
// ============================================================================
/*  normalized incomplete gamma function 
 *  \f$ Q(a,x) = \frac { \Gamma ( a , x ) }{\Gamma(a) } \f$, 
 *  where \f$ \Gamma(a,x) =  \int_x^{+\infty} t^{a-1} e^{-t}dt \f$ 
 *  is an incomplete uppper Gamma function
 *  @return the value of normalized incomplete gamma function 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2022-07-12
 */
// ============================================================================
double Ostap::Math::gamma_inc_Q
( const double a , 
  const double x ) 
{
  if ( 0 < a && 0 <= x && s_zero ( x ) ) { return 1 ; }
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_gamma_inc_Q_e ( a  , x, &result) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_gamma_inc_Q_e function" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN() ; }
    //
  }
  return result.val ;
}
// ============================================================================
/* normalized incomplete gamma function 
 *  \f$ P(a,x) = 1 - Q  (a, x ) = 
 *    = \frac{ \int_0^{x} t^{a-1} e^{-t}dt } { \Gamma(a) } \f$ 
 *  @return the value of normalized incomplete gamma function 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2022-07-12
 */
// ============================================================================
double Ostap::Math::gamma_inc_P
( const double a , 
  const double x ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_gamma_inc_P_e ( a  , x, &result) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_gamma_inc_P_e function" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN() ; }
    //
  }
  return result.val ;
}
// ============================================================================
/*  alpha_n 
 *  \f$\alpha_n(x) = \int_1^\inf t^n e^{-tz}dt \f$
 *  @param n INPUT n-parameter 
 *  @param x INPUT x-argument 
 *  @return the function value 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-27
 */
// ============================================================================
double Ostap::Math::alpha_N 
( const unsigned short n , const double x ) 
{
  const long double z = x ;
  //
  long double result  = _factorial_d_ ( n ) ;
  result /= std::pow ( z , n + 1 ) ;
  //
  return result * std::exp ( -z ) * _exp_N_ ( z , n ) ;
}
// ============================================================================
/* alpha'_n 
 *  \f$\alpha^{\prime}_n(x) = \int_0^1 t^n e^{-tx}dt \f$
 *  @param n INPUT n-parameter 
 *  @param x INPUT x-argument 
 *  @return the function value 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-27
 */
// ============================================================================
double Ostap::Math::alpha_prime_N ( const unsigned short n , const double x ) 
{
  const long double z   = x     ;
  const long double np1 = n + 1 ;
  //
  return
    s_zero   (  x ) ?                         1 / np1 :
    std::exp ( -z ) * _exp_rel_N_ ( z , n + 1 ) / np1 ;
}
// ============================================================================
/*  beta_n 
 *  \f$\beta_n(x) = \int_{-1}^{+1} t^n e^{-tx}dt \f$
 *  @param n INPUT n-parameter 
 *  @param x INPUT x-argument 
 *  @return the function value 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-27
 */
// ============================================================================
double Ostap::Math::beta_N ( const unsigned short n , const double x ) 
{
  return 
    0 == n % 2 ? 
    alpha_prime_N ( n , x ) + alpha_prime_N ( n , -x ) : 
    alpha_prime_N ( n , x ) - alpha_prime_N ( n , -x ) ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  /* evaluate "simple" continued fraction:
   *  \f$f(x) = a_0 + \frac{1}{ a_1 + \frac{1}{ a_2 + ...} } \f$
   */
  template <class ITERATOR>
  inline long double 
  _simple_cf_ 
  ( ITERATOR begin , 
    ITERATOR end   ) 
  {
    if ( begin == end ) { return 0 ; }
    //
    long double hm1 = 0  ; 
    long double km1 = 1  ;
    long double h0  = 1  ;
    long double k0  = 0  ;
    long double hp1 = h0 ;
    long double kp1 = k0 ;
    for ( ; begin != end ; ++begin ) 
    {
      const double a = *begin ;
      hp1 = a * h0 + hm1 ;
      kp1 = a * k0 + km1 ;  
      //
      hm1 = h0  ;
      km1 = k0  ;
      h0  = hp1 ;
      k0  = kp1 ;
      //
      if ( std::abs ( hp1 ) > s_large || std::abs ( kp1 ) > s_large ) 
      {
        h0  = std::ldexp ( h0  , s_imax ) ;
        k0  = std::ldexp ( k0  , s_imax ) ;
        hm1 = std::ldexp ( hm1 , s_imax ) ;
        km1 = std::ldexp ( km1 , s_imax ) ;
      }
    }
    //
    return h0 / k0 ;
  }
  // ==========================================================================
  /* evaluate "simple" continued fraction:
   *  \f$f(x) = \frac{b_0}{ 1 + \frac{b_1}{ 1 +  ...} } \f$
   */
  template <class ITERATOR>
  inline long double 
  _simple_cf_b_ 
  ( ITERATOR begin , 
    ITERATOR end   ) 
  {
    if ( begin == end ) { return 0 ; }
    //
    long double hm1 = 1  ; 
    long double km1 = 0  ;
    long double h0  = 0  ;
    long double k0  = 1  ;
    long double hp1 = h0 ;
    long double kp1 = k0 ;
    for ( ; begin != end ; ++begin ) 
    {
      const double b = *begin ;
      hp1 = h0 + b * hm1 ;
      kp1 = k0 + b * km1 ;  
      //
      hm1 = h0  ;
      km1 = k0  ;
      h0  = hp1 ;
      k0  = kp1 ;
      //
      if ( std::abs ( hp1 ) > s_large || std::abs ( kp1 ) > s_large ) 
      {
        h0  = std::ldexp ( h0  , s_imax ) ;
        k0  = std::ldexp ( k0  , s_imax ) ;
        hm1 = std::ldexp ( hm1 , s_imax ) ;
        km1 = std::ldexp ( km1 , s_imax ) ;
      }
    }
    //
    return h0 / k0 ;
  }
  // ==========================================================================
  /* evaluate "simple" continued fraction:
   *  \f$f(x) = \frac{a_1}{ b_1 + \frac{a_2}{ b_2 +  ...} } \f$
   */
  template <class ITERATOR1, class ITERATOR2>
  inline long double 
  _simple_cf_ 
  ( ITERATOR1 ab , 
    ITERATOR1 ae , 
    ITERATOR2 bb , 
    ITERATOR2 be ) 
  {
    if ( ab == ae || bb == be ) { return 0 ; }
    //
    long double hm1 = 1  ; 
    long double km1 = 0  ;
    long double h0  = 0  ;
    long double k0  = 1  ;
    long double hp1 = h0 ;
    long double kp1 = k0 ;
    for ( ; ab != ae && bb != be ;  ++ab, ++bb  ) 
    {
      //
      const double a = *ab ;
      const double b = *bb ;
      //
      hp1 = b * h0 + a * hm1 ;
      kp1 = b * k0 + a * km1 ;  
      //
      hm1 = h0  ;
      km1 = k0  ;
      h0  = hp1 ;
      k0  = kp1 ;
      //
      if ( std::abs ( hp1 ) > s_large || std::abs ( kp1 ) > s_large ) 
      {
        h0  = std::ldexp ( h0  , s_imax ) ;
        k0  = std::ldexp ( k0  , s_imax ) ;
        hm1 = std::ldexp ( hm1 , s_imax ) ;
        km1 = std::ldexp ( km1 , s_imax ) ;
      }
    }
    //
    return h0 / k0 ;
  }
  // ==========================================================================
}
// ============================================================================
/* evaluate "simple" continued fraction 
 *  \f$f(x) = a_0 + \frac{1}{ a_1 + \frac{1}{ a_2 + ...} } \f$
 *  @param a  INPUT  coefficients  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::continued_fraction_simple 
( const std::vector<double>& a ) 
{ return a.empty() ? 0 : _simple_cf_ ( a.begin() , a.end() ) ; }
// ============================================================================
/* evaluate "simple" continued fraction 
 *  \f$f(x) = \frac{b_0}{ 1 + \frac{b_1}{ 1 + ...}} \f$
 *  @param b  INPUT  coefficients  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::continued_fraction_b 
( const std::vector<double>& b )
{ 
  return 
    b.empty()            ? 0 :
    s_zero ( b.front() ) ? 0 : _simple_cf_b_ ( b.begin() , b.end() ) ; 
}
// ============================================================================
/* evaluate the continued fraction 
 *  \f$f(x) =  [b_0+] \frac{a_1}{ b_1 + \frac{a_2}{ b_2 + ...}} \f$
 *  @param a  INPUT  coefficients (len = N) 
 *  @param b  INPUT  coefficients (len = N or N+1)
 *  @attention shortest sequence is extended ith zeroes 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::continued_fraction
( const std::vector<double>& a , 
  const std::vector<double>& b ) 
{
  //
  if      ( a.size()     == b.size() ) 
  { return             _simple_cf_ ( a.begin() , a.end() , b.begin()     , b.end() ) ; }
  else if ( a.size() + 1 == b.size() ) 
  { return b.front() + _simple_cf_ ( a.begin() , a.end() , b.begin() + 1 , b.end() ) ; }
  //
  return std::numeric_limits<double>::quiet_NaN();
}
// ============================================================================
/*  confluent hypergeometrical function  1F1 aka Kummer's function
 *  \f$ f(a,b,x) = \sum_i  \frac{(a,i)}{((b,i)}\frac{x^i}{i!}\$f 
 *  @param a INPUT a-parameter 
 *  @param b INPUT b-argument  (b>0)
 *  @param x argument
 *  @retutrn value of Kummer function
 */
// ============================================================================
double Ostap::Math::kummer 
( const unsigned short a ,
  const unsigned short b , 
  const double         x ) 
{
  //
  // simple cases 
  //
  if      ( 0 == a || s_zero ( x ) ) { return 1 ; }
  else if ( a == b ) 
  {
    const long double z = x ;
    return std::abs ( x ) < 0.3 ? std::expm1 ( z ) + 1.0L : std::exp ( z ) ;
  }
  else if ( 1 == a && a < b        ) { return exp_rel_N  ( x , b - 1 ) ; }
  else if ( a + 1 == b ) 
  {
    long double gs = gamma_star ( a , -x );
    return gs * _factorial_d_ ( a ) ;
  }
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror  = gsl_sf_hyperg_1F1_int_e ( a  , b , x, &result) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from hyperg_1F1_int_e function" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN() ; }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/*  scaled complementary error function 
 *  \f$ 1 -  erf (x) = e^{-x^2} erfcx(x)  \f$ 
 *  @param x  the argument 
 *  @return the value of the scaled complementary error function 
 *  @attention  overflow happens for x<-26.6
 *  The actual implementation is copied from http://ab-initio.mit.edu/Faddeeva
 *  @see http://ab-initio.mit.edu/Faddeeva
 *  @see https://en.wikipedia.org/wiki/Error_function
 *  @see https://en.wikipedia.org/wiki/Faddeeva_function
 */
// ============================================================================
double Ostap::Math::erfcx ( const double x ) { return Faddeeva::erfcx ( x ) ; }
// ============================================================================
/*  scaled complementary error function 
 *  \f$ 1 -  erf (x) = e^{-x^2} erfcx(x)  \f$ 
 *  @param x  the argument 
 *  @return the value of the scaled complementary error function 
 *  @attention  overflow happens for x<-26.6
 *  The actual implementation is copied from http://ab-initio.mit.edu/Faddeeva
 *  @see http://ab-initio.mit.edu/Faddeeva
 *  @see https://en.wikipedia.org/wiki/Error_function
 *  @see https://en.wikipedia.org/wiki/Faddeeva_function
 */
// ============================================================================
std::complex<double> 
Ostap::Math::erfcx ( const std::complex<double>& x )
{ return Faddeeva::erfcx ( x ) ; }
// ============================================================================
/*  compute Faddeeva "w" function:
 *  w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
 *  @return the value of the scaled complementary error function 
 *  The actual implementation is copied from http://ab-initio.mit.edu/Faddeeva
 *  @see http://ab-initio.mit.edu/Faddeeva
 *  @see https://en.wikipedia.org/wiki/Error_function
 *  @see https://en.wikipedia.org/wiki/Faddeeva_function
 */
// ============================================================================
std::complex<double> 
Ostap::Math::faddeeva_w ( const std::complex<double>& x ) 
{ return Faddeeva::w ( x ) ; }
// ============================================================================
/*  complex error function (the error function of complex arguments)
 *  @param x  the argument 
 *  @return the value of the complex error function 
 *  The actual implementation is copied from http://ab-initio.mit.edu/Faddeeva
 *  @see http://ab-initio.mit.edu/Faddeeva
 *  @see https://en.wikipedia.org/wiki/Error_function
 *  @see https://en.wikipedia.org/wiki/Faddeeva_function
 */
// ============================================================================
std::complex<double>
Ostap::Math::erf   ( const std::complex<double>& x ) 
{ return Faddeeva::erf ( x ) ; }
// ============================================================================
/*  complementary complex error function 
 *  \f$ 1 -  erf (x) = erfc(x)  \f$         
 *  @param x  the argument 
 *  @return the value of the complementary complex error function 
 *  The actual implementation is copied from http://ab-initio.mit.edu/Faddeeva
 *  @see http://ab-initio.mit.edu/Faddeeva
 *  @see https://en.wikipedia.org/wiki/Error_function
 *  @see https://en.wikipedia.org/wiki/Faddeeva_function
 */
// ============================================================================
std::complex<double>  
Ostap::Math::erfc  ( const std::complex<double>& x ) 
{ return Faddeeva::erfc ( x ) ; }
// ============================================================================
/* imaginary error function 
 *  \f$ erfi(x) = -i \mathrm{erf}(ix) = \frac{2}{\sqrt{\pi}} \int_0^x e^{t^2}dt\f$ 
 *  @param x the argument
 *  @return the value of the imaginary error function 
 *  The actual implementation is copied from http://ab-initio.mit.edu/Faddeeva
 *  @see http://ab-initio.mit.edu/Faddeeva
 *  @see https://en.wikipedia.org/wiki/Error_function
 */
// ============================================================================
double Ostap::Math::erfi ( const double x ) 
{ return Faddeeva::erfi ( x ) ; }  
// ============================================================================
/* imaginary error function 
 *  \f$ erfi(x) = -i \mathrm{erf}(ix) = \frac{2}{\sqrt{\pi}} \int_0^x e^{t^2}dt\f$ 
 *  @param x the argument
 *  @return the value of the imaginary error function 
 *  The actual implementation is copied from http://ab-initio.mit.edu/Faddeeva
 *  @see http://ab-initio.mit.edu/Faddeeva
 *  @see https://en.wikipedia.org/wiki/Error_function
 */
// ============================================================================
std::complex<double> Ostap::Math::erfi ( const std::complex<double>& x ) 
{ return Faddeeva::erfi ( x ) ; }  
// ============================================================================
/** Inverse scaled error function for \f$ 0 < x \f$ 
 *  @return value of the inverse scaled error function 
 *  @see Ostap::Math:::erfcx 
 */
// ============================================================================
double Ostap::Math::erfcxinv   ( const double  x  ) 
{
  if ( x <= 0 ) { std::numeric_limits<double>::quiet_NaN (); }
  //
  static const double s_C = 2.0 / std::sqrt ( M_PI ) ;
  //
  double v =  ( x <= 1 ) ? s_SQRTPIi / x : - std::sqrt ( std::log ( x ) ) ;
  //
  for ( unsigned short i = 0 ; i < 25 ; ++i ) 
  {
    const double ev = Ostap::Math::erfcx ( v ) ;
    const double dv = ( ev - x ) / ( 2 * v * ev - s_C ) ;
    if  ( 1 < i && s_equal ( v , v - dv ) ) { return v - dv ; }
    v -= dv ;
  }
  //
  return v ;
}
// ============================================================================
/*  compute sech fuction 
 *  \$f f(x) = \frac{1}{\cosh x} = \frac{2}{ e^{x}+e^{-x} }\f$
 *  @return the value of sech function 
 */
// ============================================================================
double Ostap::Math::sech ( const double x ) 
{ return 700 < std::abs ( x )  ? 0.0 : 2.0 / ( std::exp(x)+std::exp(-x) ) ; }
// ============================================================================
/*  compute sech function 
 *  \$f f(x) = \frac{1}{\cosh x} = \frac{2}{ e^{x}+e^{-x} }\f$
 *  @return the value of sech function 
 */
// ============================================================================
std::complex<double> Ostap::Math::sech 
( const std::complex<double>& x )
{ return 700 < std::abs ( x.real() ) ? 
    std::complex<double>(0,0) : 2.0 / ( std::exp(x)+std::exp(-x) ) ; }
// ============================================================================


// ============================================================================
// Gamma function and friends 
// ============================================================================
/*  Gamma function \f$ \Gamma ( x )\ f$ 
 *  @see Ostap::Math::gamma 
 */
// ============================================================================
double Ostap::Math::tgamma ( const double x ) { return std::tgamma ( x )  ; }
// ============================================================================
/*  Gamma function \f$ \Gamma ( x )\ f$ 
 *  @see Ostap::Math::gamma 
 */
// ============================================================================
double Ostap::Math::gamma  ( const double x ) { return std::tgamma ( x )  ; }
// ============================================================================
/* logarithm of gamma function
 *  \f$ \log \Gamma ( x ) \f$ 
 */
// ============================================================================
double Ostap::Math::lgamma ( const double x ) { return std::lgamma ( x )  ; }
// ============================================================================
/*  compute inverse Gamma function 
 *  \$f f(x) = \frac{1}{\Gamma(x)}\f$
 *  @return the value of inverse Gamma functions 
 */
// ============================================================================
double Ostap::Math::igamma ( const double x ) 
{
  if ( x > 170 || ( x<=0 && Ostap::Math::isint ( x ) ) ) { return 0 ; }  // RETURN 
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_gammainv_e ( x , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_gammainv_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/* Logarithm of gamma function for complex argument 
 *  \f$ \log \Gamma ( x ) $
 * Note that the imaginary part (arg) is not well-determined 
 *   when |z| is very large, due to inevitable roundoff in restricting 
 *   to (-\pi,\pi]. This will result in a GSL_ELOSS error when it occurs. 
 *   The absolute value part (lnr), however, never suffers from loss of precision.
 */
// ============================================================================
std::complex<double> 
Ostap::Math::lgamma ( const std::complex<double>& x ) 
{ 
  // simple case 
  if ( s_zero ( x.imag() ) && 0 < x.real() && !s_zero ( x.real() ) ) 
  { return std::lgamma ( x.real() ) ; }
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result r ;
  gsl_sf_result a ;
  //
  const int ierror = gsl_sf_lngamma_complex_e ( x.real() , x.imag() , &r , &a ) ;
  //
  if ( ierror && GSL_ELOSS != ierror ) 
  {
    // 
    gsl_error ( "Error from gsl_sf_gammainv_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return std::complex<double>( r.val , a.val ) ;  
}
// ===========================================================================
/*  Gamma function of complex argument 
 *  \f$ \Gamma ( x ) \f$ 
 */
// ===========================================================================
std::complex<double> 
Ostap::Math::gamma 
( const std::complex<double>& x ) 
{
  // simple case 
  if ( s_zero ( x.imag() ) && 0 < x.real() && !s_zero ( x.real() ) ) 
  { return std::tgamma ( x.real() ) ; }
  //
  return std::exp ( Ostap::Math::lgamma ( x ) ) ;
}

// ===========================================================================
/*  Gamma function of complex argument 
 *  \f$ \Gamma ( x ) \f$ 
 */
// ===========================================================================
std::complex<double> 
Ostap::Math::tgamma 
( const std::complex<double>& x ) 
{
  // simple case 
  if ( s_zero ( x.imag() ) && 0 < x.real() && !s_zero ( x.real() ) ) 
  { return std::tgamma ( x.real() ) ; }
  //
  return std::exp ( Ostap::Math::lgamma ( x ) ) ;
}

// ============================================================================
/*  compute psi function 
 *  \$f f(x) = \frac{d}{dx}\ln \Gamma(x)\f$
 *  @return the value of psi function 
 */
// ============================================================================G
double Ostap::Math::psi ( const double x ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ( false )  ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_psi_e ( x , &result ) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_psi_e" , __FILE__ , __LINE__ , ierror ) ;
      if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN(); }
    }
  //
  return result.val ;
}
// ============================================================================
/* compute polygamma function 
 * \f$ \psi^{(n)}(x) = \left(\frac{d}{dx}\right)^{(n)}\psi(x) = 
 * = \left(\frac{d}{dx}\right)^{(n)} \log \Gamma (x) \f$ 
 * @return value of polygamma function
 * @see Ostap::Math::polygamma
 * @see Ostap::Math::digamma
 * @see Ostap::Math::trigamma
 */
// ============================================================================
double Ostap::Math::psi
( const double         x , 
  const unsigned short n )
{
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ( false )  ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_psi_n_e ( n , x , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_psi_n_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ===========================================================================  
  
  

// ============================================================================
// Beta function
// ============================================================================
/*  beta function for 
 *  \f$ B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
 *  - \f$ 0<x\f$
 *  - \f$ 0<y\f$ 
 *  @return value of beta function 
 */
// ============================================================================
double Ostap::Math::beta
( const double x ,
  const double y ) 
{ 
  //
  if ( x <  0 && s_zero ( x ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  if ( y <  0 && s_zero ( y ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  if ( s_equal ( x , 1 )      ) { return 1 / y ; }
  if ( s_equal ( y , 1 )      ) { return 1 / x ; }
  //
  if ( 0.9 < x && x <= 50 && isushort ( x ) ) { return beta ( (unsigned short) round ( x ) , y ) ; }
  if ( 0.9 < y && y <= 50 && isushort ( y ) ) { return beta ( x , (unsigned short) round ( y ) ) ; }
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ( false )  ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_beta_e ( x , y , &result ) ;
  if ( ierror ) 
  {
    // 
    if ( GSL_EUNDRFLW == ierror ) { return 0 ;}
    // 
    gsl_error ( "Error from gsl_sf_beta_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/*  beta function for 
 *  \f$ B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
 *  - \f$ 0<x\f$
 *  - \f$ 0<y\f$ 
 *  @return value of beta function 
 */
// ============================================================================
double Ostap::Math::beta 
( const unsigned short x , 
  const unsigned short y )
{
  if ( x < 1 || y < 1 ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  const unsigned short i = std::min ( x , y ) ;
  const unsigned short j = std::max ( x , y ) ;
  //
  if ( 51 <= i ) { return beta ( 1.0 * x , 1.0 * y ) ; }
  //
  double result = 1.0 / j ;
  for ( unsigned short  k = 1 ; k < i ; ++k ) { result *= k * 1.0 / ( k + j ) ; }
  return result ;
}
// ============================================================================
/*  beta function for 
 *  \f$ B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
 *  - \f$ 0<x\f$
 *  - \f$ 0<y\f$ 
 *  @return value of beta function 
 */
// ============================================================================
double Ostap::Math::beta 
( const unsigned short x , 
  const double         y )
{
  if ( x < 1 || y <= 0  ) { return std::numeric_limits<double>::quiet_NaN () ; }
  if ( 51 <= x ) { return beta ( 1.0 * x , y ) ; }
  //
  double result = 1.0 / y ;
  for ( unsigned short k = 1 ; k < x ; ++k ) { result *= k * 1.0 / ( k + y ) ; }
  return result ;
}
// ============================================================================
/*  beta function for 
 *  \f$ B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
 *  - \f$ 0<x\f$
 *  - \f$ 0<y\f$ 
 *  @return value of beta function 
 */
// ============================================================================
double Ostap::Math::beta 
( const double         x , 
  const unsigned short y ) { return beta ( y , x ) ; }
// ============================================================================
/* natural logarithm of beta function 
 *  \f$ \log B(x,y) = \log \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
 *  - \f$ 0<x\f$
 *  - \f$ 0<y\f$ 
 *  @return value of logarith of beta function 
 */
// ============================================================================
double Ostap::Math::lnbeta ( const double x , const double y ) 
{ 
  //
  if  ( x <  0 && s_zero ( x ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  if  ( y <  0 && s_zero ( y ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  if  ( s_equal ( x , 1 )      ) { return - std::log ( y )  ; }
  if  ( s_equal ( y , 1 )      ) { return - std::log ( x )  ; }
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ( false )  ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_lnbeta_e ( x , y , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_lnbeta_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/*  get the gaussian integral
 *  \f[ f = \int_a^b \exp { -\alpha^2 x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param low   the low  integration limit
 *  @param high  the high integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ============================================================================
double Ostap::Math::gaussian_integral
( const double alpha ,
  const double beta  ,
  const double low   ,
  const double high  ) 
{
  // note the difference in the arguments! 
  return details::gaussian_int ( alpha * alpha , beta , low , high ) ;
}
// ============================================================================
/*  get the gaussian integral
 *  \f[ f = \int_{a}^{_\inf} \exp { -\alpha^2 x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param low   the low  integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ============================================================================
double Ostap::Math::gaussian_integral_right
( const double alpha ,
  const double beta  ,
  const double low   ) 
{
  // note the difference in the arguments! 
  return details::gaussian_int_R ( alpha * alpha , beta , low ) ;
}
// ============================================================================
/*  get the gaussian integral
 *  \f[ f = \int_{-\inf}^b \exp { -\alpha^2 x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param high  the high integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ============================================================================
double Ostap::Math::gaussian_integral_left
( const double alpha ,
  const double beta  ,
  const double high  ) 
{
  // note the difference in the arguments! 
  return details::gaussian_int_L ( alpha * alpha , beta , high ) ;
}
// ============================================================================
namespace
{
  // ==========================================================================
  // standard Gauss PDF 
  inline double _gauss_pdf_ ( const double x )
  {
    static const double s_norm  = 1.0 / std::sqrt ( 2.0 * M_PI ) ;
    const double arg = -0.5 * x * x ;
    if ( arg < s_EXP_UNDERFLOW ) { return 0 ; } // RETURN
    return s_norm * std::exp ( arg ) ;
  } // ========================================================================
  // ==========================================================================
  /// standard Gauss CDF
  inline double _gauss_cdf_ ( const double x ) 
  {
    static const double s_isqrt2    = 1.0 / std::sqrt ( 2.0 ) ;
    static const double s_high      =  6                * s_isqrt2 ;
    static const double s_underflow = -s_ERFC_UNDERFLOW * s_isqrt2 ;
    //
    if      ( x >= s_high      ) { return 1 ; }
    else if ( x <= s_underflow ) { return 0 ; }
    //  
    const double y = x * s_isqrt2 ;
    //
    return y < -2 ? 0.5 * std::erfc ( -y ) : 0.5 * ( 1 + std::erf ( y ) ) ;
  } // ========================================================================
  // ==========================================================================
  /// Gauss integral
  inline double _gauss_int_ ( const double a , const double b )
  {
    //
    static const double s_isqrt2    = 1.0 / std::sqrt ( 2.0 ) ;
    static const double s_high      = 6                * s_isqrt2 ;
    static const double s_underflow = s_ERFC_UNDERFLOW * s_isqrt2 ;
    //     
    const double xmin = std::min ( a , b ) ;
    const double xmax = std::max ( a , b ) ;
    //
    if      ( xmax <= -s_underflow || xmin   >= s_underflow ) { return 0 ; }
    else if ( xmin <= -s_high      && s_high <= xmax        ) { return 1 ; }
    //
    const double ya = a * s_isqrt2 ;
    const double yb = b * s_isqrt2 ;
    //
    if      ( xmax <= -2 ) { return 0.5 * ( std::erfc ( std::abs ( yb ) ) - std::erfc ( std::abs ( ya ) ) ) ; } 
    else if ( xmin >= +2 ) { return 0.5 * ( std::erfc (            ya   ) - std::erfc (            yb   ) ) ; }
    //
    return 0.5 * ( std::erf ( yb ) - std::erf ( ya ) ) ;
  } // =======================================================================
  // ==========================================================================
} //                                             The end of anonymous namespace 
// ============================================================================
/*  get the standard gaussian pdf 
 *  @see https://en.wikipedia.org/wiki/Normal_distribution
 *  @param x x-value  
 *  @param sigma sigma (width)  
 *  @param mu    mu (location)
 *  @return the value of gaussian pdf 
 */
// ============================================================================
double Ostap::Math::gauss_pdf
( const double x     ,
  const double mu    ,
  const double sigma )
{
  const double aisigma = 1.0 / std::abs ( sigma ) ;
  const double dx  = ( x  - mu ) * aisigma ;
  return _gauss_pdf_ ( dx ) * aisigma ;
}
// ============================================================================
/*  get the standard gaussian cdf 
 *  @see https://en.wikipedia.org/wiki/Normal_distribution
 *  \f$ f(x) = \frac{1}{2} \left( 1 + erf ( \frac{x} { \sqrt{2} } ) \right) \f$ 
 *  @param x x-value  
 *  @return the value of gaussian cdf 
 */
// ============================================================================
double Ostap::Math::gauss_cdf 
( const double x     ,
  const double mu    ,
  const double sigma )
{
  const double y = ( x - mu ) / std::abs ( sigma ) ;
  return _gauss_cdf_ ( y ) ;
}
// ============================================================================
/*  get the Gaussian integral 
 *  @see https://en.wikipedia.org/wiki/Normal_distribution
 *  \f$ f(x) = \frac{1}{2} \left( 1 + erf ( \frac{x} { \sqrt{2} } ) \right) \f$ 
 *  \f[ f(a,b;\mu,\sigma = \int_a^b \frac{1}{\sqrt{2\pi}\sigma}
 *     \mathrm{e}^{-\frac{1}{2} \left( \frac{x-\mu}{\sigma}\right)^2}dx \f]
 *  @param a low integration limit
 *  @param b high integration limit
 *  @param mu location of Gaussian
 *  @param sigma width of the Gaussian
 */
// ============================================================================
double Ostap::Math::gauss_int
( const double a     ,
  const double b     ,
  const double mu    ,
  const double sigma )
{
  if ( s_equal ( a , b ) ) { return 0 ; }
  //
  const double iasigma = 1.0 / std::abs ( sigma ) ;
  const double as      = ( a - mu ) * iasigma ;
  const double bs      = ( b - mu ) * iasigma ;
  //
  return s_equal ( as , bs ) ? 0 : _gauss_int_ ( as , bs ) ;
}
// ============================================================================
namespace
{
  // ==========================================================================
  /** trivial Cavalieri's integral 
   *  \f[ I = \int\limits_{x_{low}}^{x_{high}} x^n dx   \f]
   *  @see https://en.wikipedia.org/wiki/Cavalieri%27s_quadrature_formula
   *  - for n<0 , it is assumed that function is continuous between xlow and xhigh 
   */
  inline double _cavalieri_ 
  ( const int    n     ,
    const double xlow  ,
    const double xhigh )
  {
    // (1) trivial checks 
    if      ( !n                       ) { return xhigh - xlow ; }
    else if ( s_equal ( xlow , xhigh ) ) { return 0 ; }
    else if ( xhigh < xlow             ) { return -_cavalieri_ ( n , xhigh , xlow ) ; }
    //
    // (0) log-case
    if ( -1 == n ) { return std::log ( std::abs ( xhigh / xlow ) ) ; }
    //
    if ( 0 > n && ( s_zero ( xlow ) || s_zero ( xhigh ) ) ) 
      { return std::numeric_limits<double>::quiet_NaN() ; }   // return NaN 
    //
    // the basic/best case 
    if      ( 0 <= xlow  ) { return ( std::pow ( xhigh , n + 1 ) - std::pow ( xlow  , n + 1 ) ) / ( n + 1 ) ;}
    else if ( 0 >= xhigh )
      {     
        const double result = ( std::pow ( std::abs ( xhigh ) , n + 1 ) - 
                                std::pow ( std::abs ( xlow  ) , n + 1 ) ) / ( n + 1 ) ;
        return ( 0 == n % 2 ? +1 : -1 ) * result ;   
      }
    else if ( n < 0 ) { return std::numeric_limits<double>::quiet_NaN() ; }   // return NaN 
    //
    const double v1 = std::pow (            xhigh  , n + 1 ) ;
    const double v2 = std::pow ( std::abs ( xlow ) , n + 1 ) ;
    //
    return ( v1 + ( 0 == n % 2 ? +1 : -1 ) * v2 ) / ( n + 1 ) ;
  }
  // ==========================================================================
  /** trivial Cavalieri's integral 
   *  \f[ I = \int\limits_{x_{low}}^{x_{high}} x^n dx   \f]
   *  @see https://en.wikipedia.org/wiki/Cavalieri%27s_quadrature_formula
   *  - for n<0 , it is assumed that function is continuous between x_low and xhigh 
   */
  inline double _cavalieri_ 
  ( const double n     ,
    const double xlow  ,
    const double xhigh )
  {
    // (1) trivial checks 
    if      ( !n                       ) { return xhigh - xlow ; }
    else if ( s_equal ( xlow , xhigh ) ) { return 0 ; }
    else if ( xhigh < xlow             ) { return -_cavalieri_ ( n , xhigh , xlow ) ; }
    //
    // (2) integer argument: covers also the log-case
    if ( Ostap::Math::isint ( n ) )
      {
        const int nn = std::lround ( n ) ;
        return _cavalieri_ ( nn , xlow , xhigh ) ;
      }
    //
    if ( std::abs ( n + 1 ) <= 0.10 ) // almost log
      {
        const double lyl    = std::log ( std::abs ( xlow  ) ) ;
        const double lyh    = std::log ( std::abs ( xhigh ) ) ;
        //
        const double d   = n + 1 ;
        //
        double logl    = lyl ;
        double logh    = lyh ;
        double dd      = 1 ;
        double invfact = 1 ;
        //
        double result = logl - logh ;
        for ( unsigned short i = 2 ; i < 100 ; ++i )
          {
            invfact  /= i   ;
            logl     *= lyl ;
            logh     *= lyh ;
            dd       *= d   ;
            //
            const double term = ( logl - logh ) * dd * invfact ;
            if ( 5 <= i && ( !term || s_zero ( term ) || s_equal ( result , result + term ) ) ) 
              {
                result   += term ;
                break ;
              } 
            result   += term ;
          }
        //
        return result ;
      }
    ///
    if ( ( ( xlow <= 0 ) || ( xhigh <= 0 ) ) && n + 1 < 0 )
      { return std::numeric_limits<double>::quiet_NaN() ; }   // return NaN 
    /// regular  case
    return ( std::pow ( xhigh  , n + 1 ) - std::pow ( xlow  , n + 1 ) ) / ( n + 1 ) ;
  }
  // ==========================================================================
}
// ============================================================================
/** trivial Cavalieri's integral 
 *  \f[ I = \int\limits_{x_{low}}^{x_{high}} \left( ax + b \right)^n dx   \f]
 *  @see https://en.wikipedia.org/wiki/Cavalieri%27s_quadrature_formula
 *  - for n<0 , it is assumed that function is continuous between x_low and xhigh 
 */
// ============================================================================
double Ostap::Math::cavalieri
( const double n     ,
  const double xlow  ,
  const double xhigh ,
  const double a     ,
  const double b     )
{
  if      ( !n                          ) { return                        xhigh - xlow   ; }
  else if ( !a  &&  ( 0 < b || 0 <= n ) ) { return std::pow ( b , n ) * ( xhigh - xlow ) ; }  
  else if ( s_equal ( xlow , xhigh )    ) { return 0 ; }
  //
  if ( !a ) { return std::numeric_limits<double>::quiet_NaN() ; }   // return NaN
  //
  // y -> a * x + b
  return _cavalieri_  ( n , a * xlow + b , a * xhigh + b ) / a ;
  //
}
// ======================================================================
/** trivial Cavalieri's integral 
 *  \f[ I = \int\limits_{x_{low}}^{x_{high}} \left( ax + b \right)^n dx   \f]
 *  @see https://en.wikipedia.org/wiki/Cavalieri%27s_quadrature_formula
 */
// ======================================================================
double Ostap::Math::cavalieri
( const int    n     ,
  const double xlow  ,
  const double xhigh , 
  const double a     ,
  const double b     ) 
{
  if      ( !n                      ) { return                         xhigh - xlow   ; }
  else if ( !a                      ) { return  std::pow ( b , n ) * ( xhigh - xlow ) ; }  
  else if ( s_equal ( xlow ,xhigh ) ) { return 0 ; }
  //
  // y -> a * x + b
  return _cavalieri_  ( n , a * xlow + b , a * xhigh + b ) / a ;
}
// ============================================================================
/*  Student's t-CDF 
 *  \f[ f(t;\nu) = \left\{
 *  \begin{array}{ll}
 *   1-\frac{1}{2}I_{x(t}}\left(\frac{\nu}{2}, \frac{1}{2}\right)   
 *   & \mathrm{for}~t\ge0 \\
 *  \frac{1}{2}I_{x(t}}\left(\frac{\nu}{2}, \frac{1}{2}\right)   
 *   & \mathrm{for}~t\<0
 *  \end{array} \right. f]
 *  where \f$ x(t) = \frac{\nu}{t^2+\nu}\f$ and 
 *  \f$I_{x}(a,b)\f$ is incomplete beta function; 
 *  @param  t t-value 
 *  @param  nu parameter nu , $\nu>0$
 */
// ============================================================================
double Ostap::Math::student_cdf 
( const double t  , 
  const double nu ) 
{
  const double anu = std::abs ( nu ) ; // NB!!
  //
  const double xt    = anu / ( t * t + anu ) ;
  const double value = 0.5 * gsl_sf_beta_inc ( 0.5 * anu , 0.5 , xt ) ;
  return t >= 0 ? 1 - value : value ;
}
// ============================================================================
/* Normalized incomplete Beta function  
 *  \f$ f(\alpha_1,\alpha_2, z ) = 
 *      I_z( \alpha_1, \alpha_2 ) = 
 *      \frac{\Beta_z(\alpha_1,\alpha_2}}
 *           {\Beta  (\alpha_1,\alpha_2}
 *  - \f$ 0<\alpha_1\f$
 *  - \f$ 0<\alpha_2\f$ 
 *  - \f$ 0<z<1\f$
 */
// ============================================================================
double Ostap::Math::beta_inc 
( const double alpha1 , 
  const double alpha2 , 
  const double z      ) 
{
  //
  if ( alpha1 <= 0       ) { return std::numeric_limits<double>::quiet_NaN () ; }
  if ( alpha2 <= 0       ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  if ( s_zero  ( z     ) ) { return 0 ; }
  if ( s_equal ( z , 1 ) ) { return 1 ; }
  //
  if ( z < 0 || 1 < z    ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  if ( 0.99 < alpha1       && 0.99 < alpha2 &&
       alpha1 < 50         && alpha2 < 50   &&
       isushort ( alpha1 ) && alpha2 < 50     )
    { return beta_inc ( (unsigned short) round ( alpha1 ) ,
                        (unsigned short) round ( alpha2 ) , z ) ; }
  // =============================================================================  
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_beta_inc_e ( alpha1 , alpha2 , z , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_beta_inc_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  return result.val ;         
}
// ============================================================================
/*  Normalized incomplete Beta function  
 *  \f$ f ( \alpha_1,\alpha_2, z ) = 
 *      I_z( \alpha_1, \alpha_2 ) = 
 *      \frac{\Beta_z(\alpha_1,\alpha_2}}
 *           {\Beta  (\alpha_1,\alpha_2}
 *  - \f$ 0<z<1\f$
 *  - \f$ 0<\alpha_1\f$
 *  - \f$ 0<\alpha_2\f$ 
 */
// ============================================================================
double Ostap::Math::beta_inc 
( const unsigned short alpha1 , 
  const unsigned short alpha2 , 
  const double         z      )
{
  if ( !alpha1 || !alpha2 ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  if ( s_zero  ( z     ) ) { return 0 ; }
  if ( s_equal ( z , 1 ) ) { return 1 ; }
  //
  if ( z < 0 || 1 < z    ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  if ( 51 <= alpha1 || 51 <= alpha2 ) { return beta_inc ( 1.0 * alpha1 , 1.0 * alpha2 , z ) ; }
  //
  const unsigned short k = alpha1 - 1 ;
  const unsigned short m = alpha2 - 1 ;
  const unsigned short N { static_cast<unsigned short> ( k + m ) } ;
  //
  const Ostap::Math::Bernstein::Basic bb { k , N } ;
  const Ostap::Math::Bernstein        B  { bb    } ;
  //
  return B.integral ( 0 , z ) * ( k + m + 1 ) ;
}
// ============================================================================
/** Derivatibe of the normalized incomplete Beta function  
 *  \f$ f ( \alpha_1,\alpha_2, z ) = 
 *      I_z( \alpha_1, \alpha_2 ) = 
 *      \frac{\Beta_z(\alpha_1,\alpha_2}}
 *           {\Beta  (\alpha_1,\alpha_2}
 *  - \f$ 0<z<1\f$
 *  - \f$ 0<\alpha_1\f$
 *  - \f$ 0<\alpha_2\f$ 
 */
// ============================================================================
double Ostap::Math::dbeta_inc 
( const double alpha1 , 
  const double alpha2 , 
  const double z      )
{
  //
  if ( alpha1 <= 0       ) { return std::numeric_limits<double>::quiet_NaN () ; }
  if ( alpha2 <= 0       ) { return std::numeric_limits<double>::quiet_NaN () ; }
  if ( z < 0 || 1 < z    ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  if ( 0.99 < alpha1       && 0.99 < alpha2 &&
       alpha1 < 50         && alpha2 < 50   &&
       isushort ( alpha1 ) && alpha2 < 50     )
    { return dbeta_inc ( (unsigned short) round ( alpha1 ) ,
                         (unsigned short) round ( alpha2 ) , z ) ; }
  // ==========================================================================
  if ( s_zero  ( z     ) && 1 < alpha1 ) { return  0 ; }
  if ( s_equal ( z , 1 ) && 1 < alpha2 ) { return  0 ; }
  // ==========================================================================
  const double result = std::pow ( z , alpha1 - 1 ) * std::pow ( 1 - z , alpha2 - 1 ) ;
  return result / beta ( alpha1 , alpha2 ) ;
  // ==========================================================================
}
// ============================================================================
/*  Derivative of the normalized incomplete Beta function  
 *  \f$ f ( \alpha_1,\alpha_2, z ) = 
 *      I_z( \alpha_1, \alpha_2 ) = 
 *      \frac{\Beta_z(\alpha_1,\alpha_2}}
 *           {\Beta  (\alpha_1,\alpha_2}
 *  - \f$ 0<z<1\f$
 *  - \f$ 0<\alpha_1\f$
 *  - \f$ 0<\alpha_2\f$ 
 */
// ============================================================================
double Ostap::Math::dbeta_inc 
( const unsigned short alpha1 , 
  const unsigned short alpha2 , 
  const double         z      )
{
  if ( !alpha1 || !alpha2 ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  if ( z < 0 || 1 < z     ) { return std::numeric_limits<double>::quiet_NaN () ; }
  //
  if ( 51 <= alpha1 || 51 <= alpha2 ) { return dbeta_inc ( 1.0 * alpha1 , 1.0 * alpha2 , z ) ; }
  //
  const unsigned short k = alpha1 - 1 ;
  const unsigned short m = alpha2 - 1 ;
  const unsigned short N { static_cast<unsigned short> ( k + m ) } ;
  //
  const Ostap::Math::Bernstein::Basic bb { k , N } ;
  const Ostap::Math::Bernstein        B  { bb    } ;
  //
  return B ( z ) * ( k + m + 1 ) ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  typedef unsigned long long ULL ;
  // ==============================================================================
  // calculate Pochhammer symbol
  // ==============================================================================
  inline long double _pochhammer_ ( const long double x , const unsigned short N ) 
  {
    if      (  0 == N ) { return 1 ; }
    else if (  1 == N ) { return x ; }
    else if (  2 == N ) { return Ostap::Math::Pochhammer_< 2>::evaluate ( x ) ; }
    else if (  3 == N ) { return Ostap::Math::Pochhammer_< 3>::evaluate ( x ) ; }
    else if (  4 == N ) { return Ostap::Math::Pochhammer_< 4>::evaluate ( x ) ; }
    else if (  5 == N ) { return Ostap::Math::Pochhammer_< 5>::evaluate ( x ) ; }
    else if (  6 == N ) { return Ostap::Math::Pochhammer_< 6>::evaluate ( x ) ; }
    else if (  7 == N ) { return Ostap::Math::Pochhammer_< 7>::evaluate ( x ) ; }
    else if (  8 == N ) { return Ostap::Math::Pochhammer_< 8>::evaluate ( x ) ; }
    else if (  9 == N ) { return Ostap::Math::Pochhammer_< 9>::evaluate ( x ) ; }
    else if ( 10 == N ) { return Ostap::Math::Pochhammer_<10>::evaluate ( x ) ; }
    else if ( 11 == N ) { return Ostap::Math::Pochhammer_<11>::evaluate ( x ) ; }
    else if ( 12 == N ) { return Ostap::Math::Pochhammer_<12>::evaluate ( x ) ; }
    else if ( 13 == N ) { return Ostap::Math::Pochhammer_<13>::evaluate ( x ) ; }
    else if ( 14 == N ) { return Ostap::Math::Pochhammer_<14>::evaluate ( x ) ; }
    else if ( 15 == N ) { return Ostap::Math::Pochhammer_<15>::evaluate ( x ) ; }
    else if ( 16 == N ) { return Ostap::Math::Pochhammer_<16>::evaluate ( x ) ; }
    //
    // more   specific treatment 
    //
    if ( s_zero ( x )  ) { return 0 ; }  // RETURN 
    // avoid too negative values 
    if ( x < 0.5L - N ) 
    { return _pochhammer_ ( std::abs ( x ) - N + 1 , N ) * ( N % 2 ? -1 : 1 ) ; }
    //
    const double s_delta        = 1.e-8 ;
    const bool   use_dimidation = 
      ( 1 - N - s_delta < x && x < s_delta ) && 
      ( std::abs ( x - Ostap::Math::round ( x ) ) < s_delta ) ;
    //
    // use the dimidation formula
    if ( 96 >= N || use_dimidation ) 
    {
      //
      const unsigned short K2 =          N / 2       ;
      const unsigned short K1 =  N % 2 ? K2 + 1 : K2 ;
      //
      return 
        std::ldexp ( _pochhammer_ ( std::ldexp ( x     , -1 ) , K1 ) , K1 ) * 
        std::ldexp ( _pochhammer_ ( std::ldexp ( x + 1 , -1 ) , K2 ) , K2 ) ;
      //
    }
    /// use the generic formula 
    return std::exp ( std::lgamma ( x + N ) - std::lgamma ( x ) ) ;
  }
  // ==========================================================================
  typedef std::pair<long double, long double> RESULT ;
  // ==========================================================================
  inline RESULT _pochhammer2_ ( const long double x , const unsigned short N ) 
  {
    //
    if      (  0 == N ) { return std::make_pair ( 1.L , 0.L ) ; }
    else if (  1 == N ) { return std::make_pair ( x   , 1.L ) ; }
    else if (  2 == N ) { return Ostap::Math::Pochhammer_< 2>::value_with_derivative ( x ) ; }
    else if (  3 == N ) { return Ostap::Math::Pochhammer_< 3>::value_with_derivative ( x ) ; }
    else if (  4 == N ) { return Ostap::Math::Pochhammer_< 4>::value_with_derivative ( x ) ; }
    else if (  5 == N ) { return Ostap::Math::Pochhammer_< 5>::value_with_derivative ( x ) ; }
    else if (  6 == N ) { return Ostap::Math::Pochhammer_< 6>::value_with_derivative ( x ) ; }
    else if (  7 == N ) { return Ostap::Math::Pochhammer_< 7>::value_with_derivative ( x ) ; }
    else if (  8 == N ) { return Ostap::Math::Pochhammer_< 8>::value_with_derivative ( x ) ; }
    else if (  9 == N ) { return Ostap::Math::Pochhammer_< 9>::value_with_derivative ( x ) ; }
    else if ( 10 == N ) { return Ostap::Math::Pochhammer_<10>::value_with_derivative ( x ) ; }
    else if ( 11 == N ) { return Ostap::Math::Pochhammer_<11>::value_with_derivative ( x ) ; }
    else if ( 12 == N ) { return Ostap::Math::Pochhammer_<12>::value_with_derivative ( x ) ; }
    else if ( 13 == N ) { return Ostap::Math::Pochhammer_<13>::value_with_derivative ( x ) ; }
    else if ( 14 == N ) { return Ostap::Math::Pochhammer_<14>::value_with_derivative ( x ) ; }
    else if ( 15 == N ) { return Ostap::Math::Pochhammer_<15>::value_with_derivative ( x ) ; }
    else if ( 16 == N ) { return Ostap::Math::Pochhammer_<16>::value_with_derivative ( x ) ; }
    //    
    // avoid too negative values 
    if ( x < 0.5L - N ) 
    {  
      const RESULT r = _pochhammer2_ ( std::abs ( x ) - N + 1 , N ) ;
      const int    s =  ( N % 2 ? -1 : 1 ) ;
      return std::make_pair ( s * r.first , -1 * s * r.second ) ;
    }
    //
    const double s_delta        = 1.e-8 ;
    const bool   use_dimidation = 
      ( 1 - N - s_delta < x && x < s_delta ) && 
      ( std::abs ( x - Ostap::Math::round ( x ) ) < s_delta ) ;
    //
    // use the dimidation formula
    if ( 96 >= N || use_dimidation ) 
    {
      //
      const unsigned short K2 =          N / 2       ;
      const unsigned short K1 =  N % 2 ? K2 + 1 : K2 ;
      //
      const RESULT r1 = _pochhammer2_ ( std::ldexp ( x     , -1 ) , K1 ) ;
      const RESULT r2 = _pochhammer2_ ( std::ldexp ( x + 1 , -1 ) , K2 ) ;
      //
      return std::make_pair ( std::ldexp ( r1.first  * r2.first  , N     ) ,                         
                              std::ldexp ( r1.first  * r2.second , N - 1 ) + 
                              std::ldexp ( r1.second * r2.first  , N - 1 ) ) ;
    }
    ///
    /// use the generic formula 
    const  double p = std::exp ( std::lgamma ( x + N ) - std::lgamma ( x ) ) ;
    //
    return std::make_pair ( p , 
                            p * ( Ostap::Math::psi ( x + N ) - Ostap::Math::psi ( x ) ) );
  }
  // ==========================================================================
  inline long double __pochhammer__ ( const long double x , const unsigned short N ) 
  {
    return  
      0 == N       ? 1 :
      1 == N       ? x :
      s_zero ( x ) ? 0 :
      ( 0.5L - N ) < x ? 
      _pochhammer_ ( x                      , N ) :
      _pochhammer_ ( std::abs ( x ) - N + 1 , N ) * ( N % 2 ? -1 : 1 ) ;
  }
  // ==========================================================================
}
// ============================================================================
/*  Pochhammer symbol, aka rising factorial 
 *  \f[ P(x,n) = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
 *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
 *  @see Ostap::Math::rising_factorial
 */
// ============================================================================
double Ostap::Math::pochhammer ( const double x , const unsigned short n ) 
{ return __pochhammer__ ( x , n ) ; }
// ============================================================================
/*  Rising  factorial, aka Pochhammer's symbol   
 *  \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
 *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
 *  @see Ostap::Math::pochhammer 
 *  @see Ostap::Math::falling_factorial
 */
// ============================================================================
double Ostap::Math::rising_factorial  ( const double x , const unsigned short n ) 
{ return __pochhammer__ ( x , n ) ; }
// ============================================================================
/*  Falling factorial, aka Pochhammer's symbol   
 *  \f[ (x)_n = \Pi^{k-1}_{k=0}  (x - k) \f] 
 *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
 *  @see Ostap::Math::rising_factorial
 */
// ============================================================================
double Ostap::Math::falling_factorial ( const double x , const unsigned short n ) 
{ return __pochhammer__ ( -1 * x , n ) * ( n % 2 ? -1 : 1 ) ; }
// ============================================================================
/*  Pochhammer symbol, aka "rising factorial" and its derivative 
 *  \f[ P(x,n) = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
 *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
 *  @see Ostap::Math::rising_factorial
 *  @see Ostap::Math::pochhammer
 */
// ============================================================================
std::pair<double,double> Ostap::Math::pochhammer_with_derivative 
( const double         x , 
  const unsigned short n ) { return _pochhammer2_ ( x , n ) ; }
// ============================================================================



// ============================================================================
// Elliptic integrals 
// ============================================================================



// ============================================================================
/*  Complete elliptic integral \f$ K(k) \f$  
 *  \[ K(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_K_gsl ( const double k ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_Kcomp_e ( k , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_Kcomp_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/**  Complete elliptic integral \f$ E(k) \f$  
 *  \[ E(k) \equiv E ( \frac{\pi}{2}, k ) \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_E_gsl ( const double k   ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_Ecomp_e ( k , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_Ecomp_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/* Complete elliptic integral \f$ K(k) \f$  
 *  \[ K(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_K ( const double k ) { return elliptic_Km ( k * k ) ; } 
// ============================================================================
/*  Complete elliptic integral \f$ E(k) \f$  
 *  \[ E(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_E ( const double k ) { return elliptic_Em ( k * k ) ; } 
// ============================================================================
/*  Complete elliptic integral \f$ K[m] \f$  
 *  \f[ K[m] = K(k) = 
 *     \frac{\pi}{2 \mathrm{agm} (1, \sqrt{1-k^2}} = 
 *     \frac{\pi}{2 \mathrm{agm} (1, \sqrt{1-m}  } = 
 *     \frac{\pi}{2 \mathrm{agm} (1, \sqrt{ m^{\prime}}} \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_Km 
( const double m   )
{
  if ( s_zero ( m ) )  { return 0.5 * M_PI ; }
  if ( ( 1 < m ) || ( m < 0 ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  const long double sq_mprime = std::sqrt ( 1.0L - m ) ;
  return 0.5 * M_PI / agm ( 1.0 , sq_mprime) ;
}
// ============================================================================
/*  Complete elliptic integral \f$ E[m] \f$ as function of parameter m  
 *  \[ E(,) \equiv F ( \frac{\pi}{2}, m ) \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  \f[ E(m) \equiv 2 R_{G}(0, 1 - m , 1 ) \f]
 *  @see Eq. (55) in arXiv:math/9409227
 */
// ============================================================================
double Ostap::Math::elliptic_Em 
( const double m )
{
  if      ( s_zero  ( m     ) ) { return 0.5 * M_PI ; }
  else if ( s_equal ( m , 1 ) ) { return 1 ; }
  if ( ( 1 < m ) || ( m < 0 ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  return 2 * carlson_RG ( 0 , 1 - m , 1 ) ;
}
// ============================================================================
/* Trigonometric form of incomplete elliptic integral \f$ F(\phi,k) \f$
 *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \dfrac{ d \psi }{\sqrt{1-k^2 \sin^2 \psi }}\f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_F
( const double phi , 
  const double k   )
{ return elliptic_Fm ( phi , k * k ) ; }
// ============================================================================
/*   Trigonometric form of incomplete elliptic integral \f$ E(\phi,k) \f$
 *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \sqrt{1-k^2 \sin^2 \psi } d \psi \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_E
( const double phi , 
  const double k   ) 
{ return elliptic_Em ( phi , k * k ) ; }
// ============================================================================
/* Trigonometric form of incomplete elliptic integral \f$ F(\phi,m) \f$
 *  \f[ F(\phi,m) \equiv \int_{0}^{\phi} \dfrac{ d \psi }{\sqrt{1-m \sin^2 \psi }}\f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Eq. (59) in arXiv:math/9409227
 */
// ============================================================================
double Ostap::Math::elliptic_Fm
( const double phi , 
  const double m   )
{
  if ( s_zero ( phi ) ) { return 0 ; }
  if ( ( 1 < m ) || ( m < 0 ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  if ( std::abs ( phi ) > 0.5 * M_PI )
    {
      const double nc = std::floor ( phi / M_PI + 0.5 ) ;
      const double kk = elliptic_Km ( m ) ;      
      return elliptic_Fm ( phi - nc * M_PI , m ) + 2 * nc * kk ;	
    }
  //
  const long double sinphi = std::sin ( phi ) ;
  const long double cosphi = std::cos ( phi ) ;
  //
  return sinphi * carlson_RF ( cosphi * cosphi , 1.0 - m * sinphi * sinphi , 1 ) ;
}
// ========================================================================
/* Trigonometric form of incomplete elliptic integral \f$ E(\phi,m) \f$
 *  \f[ F(\phi,m) \equiv \int_{0}^{\phi} \sqrt{1-m \sin^2 \psi } d \psi \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Eq. (60) in arXiv:math/9409227
 */
// ============================================================================
double Ostap::Math::elliptic_Em
( const double phi , 
  const double m   )
{
  if ( s_zero ( phi ) ) { return 0 ; }
  if ( ( 1 < m ) || ( m < 0 ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  if ( std::abs ( phi ) > 0.5 * M_PI )
    {
      const double nc = std::floor ( phi / M_PI + 0.5 ) ;
      const double ee = elliptic_Em ( m ) ;      
      return elliptic_Em ( phi - nc * M_PI , m ) + 2 * nc * ee ;	
    }
  //  
  const long double sinphi = std::sin ( phi ) ;
  const long double cosphi = std::cos ( phi ) ;
  const long double sp2    = sinphi * sinphi  ;
  const long double cp2    = cosphi * cosphi  ;
  //
  const double rf = carlson_RF ( cp2 , 1 - m * sp2 , 1 ) ;
  const double rd = carlson_RD ( cp2 , 1 - m * sp2 , 1 ) ;
  //
  return sinphi * ( rf - m * sp2 * rd / 3 ) ;
}
// ============================================================================
/*  Trigonometric form of incomplete elliptic integral \f$ F(\phi,k) \f$
 *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \dfrac{ d \psi }{\sqrt{1-k^2 \sin^2 \phi }}\f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_F_gsl ( const double phi , const double k   ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_F_e ( phi , k , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_F_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/* Trigonometric form of incomplete elliptic integral \f$ E(\phi,k) \f$
 *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \sqrt{1-k^2 \sin^2 \phi } d \psi \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 */
// ============================================================================
double Ostap::Math::elliptic_E_gsl ( const double phi , const double k   ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_E_e ( phi , k , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_E_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/*  Jacobi zeta function \f$ Z(\beta k) \f$
 *  \f[ K(k) Z( \beta , k ) = K(k) E(\beta, k ) - E(k) F(\beta,k) \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  http://functions.wolfram.com/EllipticIntegrals/JacobiZeta/introductions/IncompleteEllipticIntegrals/ShowAll.html
 */
// ============================================================================
double Ostap::Math::elliptic_Z ( const double beta  , const double k   ) 
{
  const double K_k  = elliptic_K ( k ) ;
  const double E_k  = elliptic_E ( k ) ;
  const double E_bk = elliptic_E ( beta , k ) ;
  const double F_bk = elliptic_F ( beta , k ) ;
  //
  return E_bk - E_k * F_bk / K_k ;
}
// ============================================================================
/*  Product of Jacobi zeta function \f$ Z(\beta,k) \f$
 *  and complete elliptic integral \f$ K(k) \f$
 *  \f[ K(k) Z( \beta , k ) = K(k) E(\beta, k ) - E(k) F(\beta,k) \f] 
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227 Eq. (63) 
 */
// ============================================================================
double Ostap::Math::elliptic_KZ ( const double beta  , const double k   ) 
{
  const double sinbeta = std::sin ( beta ) ;
  const double cosbeta = std::cos ( beta ) ;
  const double alpha   = 1.0L - k * k * sinbeta * sinbeta ;
  //
  const double r = carlson_RJ  ( 0 , 1 - k * k , 1 , alpha ) ;
  return k * k * sinbeta * cosbeta * std::sqrt ( alpha ) * r / 3 ;
}
// ============================================================================
/*  difference in complete elliptic integrals  \f$ K(k) \f$ and \f$ E(k) \f$
 *  \f[ K(k) - E(k) = \frac{k^2}{3}R_D\left(0,1-k^2,1\right)\f],
 *  where \f$ R_D(x,y,z)\f$ is a symmetric Carlson form 
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227  Eq.(57)
 */
// ============================================================================
double Ostap::Math::elliptic_KmE ( const double k   )
{
  // ==========================================================================
  // see https://arxiv.org/abs/math/9409227  Eq. (57)
  return k * k * carlson_RD ( 0 , 1 - k * k , 1 ) / 3 ; 
  // ==========================================================================
}
// ============================================================================
/* elliptic \f$ \Pi(\alpha^2,k)\f$ function 
 *  - \f$ alpha^2 < 1 \f$ 
 *  - \f$ k      < 1 \f$ 
 *  \f[ \Pi(\alpha^2, k) - K(k) = 
 *   \frac{1}{3}\alpha^2 R_J( 0, 1-k^2, 1 , 1 - \alpha^2) \f] 
 */ 
// ============================================================================
double Ostap::Math::elliptic_PI
( const double alpha2 , 
  const double k      ) 
{
  return elliptic_K ( k ) + 
    alpha2 * carlson_RJ ( 0 , 1 - k * k , 1 , 1 - alpha2 ) / 3 ;
}
// ============================================================================
/* elliptic \f$ \Pi(\alpha^2,k) - K(k) \f$ function 
 *  \f[ \Pi(\alpha^2, k) - K(k) \equiv  
 *   \frac{1}{3}\alpha^2 R_J( 0, 1-k^2, 1 , 1 - \alpha^2) \f] 
 *  - \f$ alpha^2 < 1 \f$ 
 *  - \f$ k      < 1 \f$ 
 */ 
// ============================================================================
double Ostap::Math::elliptic_PImK  
( const double alpha2 , 
  const double k      )
{ return alpha2 * carlson_RJ ( 0 , 1 - k * k , 1 , 1 - alpha2 ) / 3 ; }
// ========================================================================



// ============================================================================
// Symmetric Carlson forms 
// ============================================================================


// ============================================================================
/* Symmetric Carlson form 
 *  \f[ R_F(x,y,z) = \int_0^{+\infty} 
 *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2}dt \f]
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RF 
( const double x , 
  const double y , 
  const double z ) 
{
  //
  const double A0 = ( x + y + z ) / 3 ;
  // 
  double xm = x  ;
  double ym = y  ;
  double zm = z  ;
  double Am = A0 ;
  //
  const double dx = A0 - x ;
  const double dy = A0 - y ;
  const double dz = A0 - z ;
  //
  static const double s_R = std::pow ( 3 * 2 * GSL_DBL_EPSILON , -1.0/8 ) ;
  const double Q  = s_R * std::max ( std::abs ( dx ) , std::max ( std::abs ( dy ) , std::abs ( dz ) ) ) ;
  //
  double       Qn = Q  ;
  double       p4 = 1  ;
  for ( unsigned int m = 0 ; m < 100  ; ++m ) 
  {
    if ( 0 < m && Qn < std::abs ( Am ) ) { break ; }
    //
    const double xm_sq = std::sqrt ( xm ) ;
    const double ym_sq = std::sqrt ( ym ) ;
    const double zm_sq = std::sqrt ( zm ) ;
    //
    const double lm = xm_sq * ym_sq + xm_sq * zm_sq + ym_sq * zm_sq ;
    //
    const double Am1 = 0.25 * ( Am + lm ) ;
    const double xm1 = 0.25 * ( xm + lm ) ;
    const double ym1 = 0.25 * ( ym + lm ) ;
    const double zm1 = 0.25 * ( zm + lm ) ;
    //
    xm = xm1 ;
    ym = ym1 ;
    zm = zm1 ;
    Am = Am1 ;
    Qn /= 4  ;    
    p4 /= 4  ;
  }
  //
  const double iA4 = p4 / Am  ;
  const double X   = dx * iA4 ;
  const double Y   = dy * iA4 ;
  const double Z   = - ( X + Y );
  //
  const double XY = X * Y ;
  const double Z2 = Z * Z ;
  //
  const double E2 = XY - Z2 ;
  const double E3 = XY * Z  ;
  //
  static const double s_c1 = -1. /  10 ;
  static const double s_c2 = +1. /  14 ;
  static const double s_c3 = +1. /  24 ;
  static const double s_c4 = -3. /  44 ;
  static const double s_c5 = -5. / 208 ;
  static const double s_c6 = +3. / 104 ;
  static const double s_c7 = +1. /  16 ;
  //
  return ( 1 
           + s_c1 * E2 
           + s_c2 * E3 
           + s_c3 * E2 * E2 
           + s_c4 * E2 * E3 
           + s_c5 * E2 * E2 * E2  
           + s_c6 * E3 * E3 
           + s_c7 * E2 * E2 * E3 ) / std::sqrt ( Am ) ;
}
// ============================================================================
/* Symmetric Carlson form 
 *  \f[ R_F(x,y,z) = \int_0^{+\infty} 
 *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2}dt \f]
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RF_gsl  
( const double x , 
  const double y , 
  const double z ) 
{
  //
  // if ( s_zero ( x ) ) { return carlson_RF ( y , z ) ; }
  // if ( s_zero ( y ) ) { return carlson_RF ( x , z ) ; }
  // if ( s_zero ( z ) ) { return carlson_RF ( x , y ) ; }
  //
  // use GSL 
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_RF_e 
    ( x , y , z , GSL_PREC_DOUBLE , &result ) ;
  //
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_RF_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  return result.val ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_J(x,y,z,p) = \int_0^{+\infty} 
 *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2} (t+p) dt \f]
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RJ
( const double x , 
  const double y , 
  const double z , 
  const double p ) 
{
  //
  const double A0 = 0.20 * ( x + y + z + 2 * p ) ;
  // 
  double xm = x  ;
  double ym = y  ;
  double zm = z  ;
  double pm = p  ;
  double Am = A0 ;
  //
  const double dx = A0 - x ;
  const double dy = A0 - y ;
  const double dz = A0 - z ;
  const double dp = A0 - p ;
  //
  const double delta = ( p - x ) * ( p - y ) * ( p - z ) ;
  //
  static const double s_R =std::pow ( 0.25 * 2 * GSL_DBL_EPSILON , -1.0/8 ) ;
  const double Q  = s_R * std::max ( std::abs ( dx ) , 
                                     std::max ( std::abs ( dy ) , 
                                                std::max ( std::abs ( dz ) , 
                                                           std::abs ( dp ) ) ) ) ;
  //
  unsigned int n       = 0 ;
  double       Qn      = Q ;
  double       p4      = 1 ;
  // 
  double       result  = 0 ;
  double       sm      = 1 ;
  //
  for ( unsigned int m = 0 ; m < 100  ; ++m ) 
  {
    if ( 0 < m && Qn < std::abs ( Am ) ) { break ; }
    //
    const double xm_sq = std::sqrt ( xm ) ;
    const double ym_sq = std::sqrt ( ym ) ;
    const double zm_sq = std::sqrt ( zm ) ;
    const double pm_sq = std::sqrt ( pm ) ;
    //
    const double lm = xm_sq * ym_sq + xm_sq * zm_sq + ym_sq * zm_sq ;
    const double dm =  ( pm_sq + xm_sq ) * ( pm_sq + ym_sq ) * ( pm_sq + zm_sq ) ;
    //
    if ( 0 == m ) { sm = 0.5 * dm ; }
    //
    const double em = p4 * p4 * p4 * delta / ( dm * dm ) ;
    if ( -1.5 < em && em < -0.5 &&  m < 3 )
    {
      const double b  = 2 * pm_sq * ( pm + xm_sq * ( ym_sq  +  zm_sq) + ym_sq * zm_sq ) / dm ;
      result += 6 * p4 * carlson_RC ( 1 , b  ) / dm ;
    }
    else 
    { result += 6 * p4 * carlson_RC ( 1 , 1 + em ) / dm  ; }
    
    //
    const double Am1 = 0.25 * ( Am + lm ) ;
    const double xm1 = 0.25 * ( xm + lm ) ;
    const double ym1 = 0.25 * ( ym + lm ) ;
    const double zm1 = 0.25 * ( zm + lm ) ;
    const double pm1 = 0.25 * ( pm + lm ) ;
    //
    xm = xm1 ;
    ym = ym1 ;
    zm = zm1 ;
    pm = pm1 ;
    Am = Am1 ;
    ++n      ;
    Qn /= 4  ;
    p4 /= 4  ;
    //
    const double rm = sm * ( 1 + std::sqrt  ( 1 + p4 * delta / ( sm * sm ) ) ) ;
    sm = 0.5 * ( dm * rm - p4 * p4 * delta ) / ( dm + p4 * rm ) ;
    //
  }
  //
  const double iA4 =  p4 / Am  ;
  const double X   =  dx * iA4 ;
  const double Y   =  dy * iA4 ;
  const double Z   =  dz * iA4 ;
  const double P   = -0.5 * ( X + Y + Z );
  //
  const double XY  = X * Y ;
  const double XZ  = X * Z ;
  const double YZ  = Y * Z ;
  const double XYZ = XY * Z ;
  //
  const double Z2  = Z  * Z ;
  const double P2  = P  * P ;
  const double P3  = P2 * P ;
  //
  const double E2  = XY + XZ + YZ - 3 * P2 ;
  const double E3  = XYZ + 2 * E2 * P + 4 * P3 ;
  const double E4  = ( 2 * XYZ + E2 * P + 3 * P3 ) * P ;
  const double E5  = XYZ * P2 ;
  //
  static const double s_c1  = -3.  /  14 ;
  static const double s_c2  = +1.  /   6 ;
  static const double s_c3  = +9.  /  88 ;
  static const double s_c4  = -3.  /  22 ;
  static const double s_c5  = -9.  /  52 ;
  static const double s_c6  = +3.  /  26 ;
  //
  static const double s_c7  = -1.  /  16 ;
  static const double s_c8  = +3.  /  40 ;
  static const double s_c9  = +3.  /  20 ;
  static const double s_c10 = +45. / 272 ;
  static const double s_c11 = -9.  /  68 ;
  //
  const double res = iA4 / std::sqrt ( Am ) * ( 1 
                                                + s_c1  * E2 
                                                + s_c2  * E3 
                                                + s_c3  * E2 * E2 
                                                + s_c4  * E4 
                                                + s_c5  * E2 * E3 
                                                + s_c6  * E5   
                                                //
                                                + s_c7  * E2 * E2 * E2 
                                                + s_c8  * E3 * E3  
                                                + s_c9  * E2 * E4  
                                                + s_c10 * E2 * E2 * E3   
                                                + s_c11 * ( E3 * E4 + E2 * E5 ) ) ;
  //
  // const double add = 3 * carlson_RC ( 1 , 1 + p4 * delta / ( sm * sm ) ) / sm  ;
  //
  return res + result ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_J(x,y,z,p) = \int_0^{+\infty} 
 *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2} (t+p) dt \f]
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RJ_gsl 
( const double x , 
  const double y , 
  const double z , 
  const double p ) 
{
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_RJ_e 
    ( x , y , z , p , GSL_PREC_DOUBLE , &result ) ;
  //
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_RJ_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_C(x,y) = R_F(x,y,y) = \int_0^{+\infty} 
 *   (t+x)^{-1/2}(t+y)^{-1}dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 *  For negative y, Cauchy principal value it returned 
 *  \f[ R_C(x,-y) = \left(\frac{x}{x+y}\right)R_C(x+y,y), 0 < y \f] 
 */
// ============================================================================
double Ostap::Math::carlson_RC
( const double x , 
  const double y ) 
{
  //
  if ( y < 0 ) 
  { return s_zero ( x ) ? 0.0 : std::sqrt ( x / ( x - y ) ) * carlson_RC ( x - y , -y ) ; }
  //
  const double A0 = ( x + 2 * y ) / 3 ;
  // 
  double xm = x  ;
  double ym = y  ;
  double Am = A0 ;
  //
  const double dx = A0 - x ;
  //
  static const double s_R =std::pow ( 3 * 2 * GSL_DBL_EPSILON , -1.0/8 ) ;
  const double Q  = s_R * std::abs ( dx ) ;
  //
  double Qn = Q ;
  double p4 = 1 ;
  //
  for ( unsigned int m = 0 ; m < 100 ; ++m ) 
  {
    if ( 0 < m && Qn < std::abs ( Am ) ) { break ; }
    //
    const double xm_sq = std::sqrt ( xm ) ;
    const double ym_sq = std::sqrt ( ym ) ;
    //
    const double lm = 2 * xm_sq * ym_sq + ym ;
    //
    const double xm1 =  0.25 * ( xm + lm ) ;
    const double ym1 =  0.25 * ( ym + lm ) ;
    const double Am1 =  0.25 * ( Am + lm ) ;
    //
    xm = xm1 ;
    ym = ym1 ;
    Am = Am1 ;
    //
    Qn /= 4  ;
    p4 /= 4  ;
  }
  //
  const double s = p4 * ( y - A0 ) / Am ;
  //
  static  const std::array<long double,8> s_poly 
  { { 9.0L/8 , 159.0L/208 , 9.0L/22 , 3.0L/8 , 1.0L/7 , 3.0L/10 , 0.0 , 1.0L } } ;
  //
  const double result = 
    Ostap::Math::Clenshaw::monomial_sum ( s_poly.begin () , s_poly.end   () , s ).first  ;
  //
  return result / std::sqrt ( Am )  ;  
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_C(x,y) = R_F(x,y,y) = \int_0^{+\infty} 
 *   (t+x)^{-1/2}(t+y)^{-1}dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 *  For negative y, Cauchy principal value it returned 
 *  \f[ R_C(x,-y) = \left(\frac{x}{x+y}\right)R_C(x+y,y), 0 < y \f] 
 */
// ============================================================================
double Ostap::Math::carlson_RC_gsl
( const double x , 
  const double y ) 
{
  //
  if ( y < 0 ) 
  { return s_zero ( x ) ? 0.0 : std::sqrt ( x / ( x - y ) ) * carlson_RC_gsl ( x - y , -y ) ; }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_RC_e 
    ( x , y , GSL_PREC_DOUBLE , &result ) ;
  //
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_RC_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_D (x,y,z) = R_J(x,y,z,z) = \int_0^{+\infty} 
 *  \left[  (t+x)(t+y)\right]^{-1/2} (t+z)^{-3/2} dt \f]
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RD_gsl  
( const double x , 
  const double y , 
  const double z ) 
{
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_ellint_RD_e 
    ( x , y , z , GSL_PREC_DOUBLE , &result ) ;
  //
  if ( ierror ) 
  {
    //
    gsl_error ( "Error from gsl_sf_ellint_RD_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
  }
  //
  return result.val ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_D (x,y,z) = R_J(x,y,z,z) = \int_0^{+\infty} 
 *  \left[  (t+x)(t+y)\right]^{-1/2} (t+z)^{-3/2} dt \f]
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RD
( const double x , 
  const double y , 
  const double z ) 
{
  //
  const double A0 = 0.20 * ( x + y + 3 * z ) ;
  // 
  double xm = x  ;
  double ym = y  ;
  double zm = z  ;
  double Am = A0 ;
  //
  const double dx = A0 - x ;
  const double dy = A0 - y ;
  const double dz = A0 - z ;
  //
  static const double s_R =std::pow ( 0.25 * 2 * GSL_DBL_EPSILON , -1.0/8 ) ;
  const double Q  = s_R * std::max ( std::abs ( dx ) , std::max ( std::abs ( dy ) , std::abs ( dz ) ) ) ;
  //
  double result   = 0  ;
  unsigned int n  = 0  ;
  double       Qn = Q  ;
  double       p4 = 1  ;
  for ( unsigned int m = 0 ; m < 100  ; ++m ) 
  {
    if ( 0 < m && Qn < std::abs ( Am ) ) { break ; }
    //
    const double xm_sq = std::sqrt ( xm ) ;
    const double ym_sq = std::sqrt ( ym ) ;
    const double zm_sq = std::sqrt ( zm ) ;
    //
    const double lm = xm_sq * ym_sq + xm_sq * zm_sq + ym_sq * zm_sq ;
    //
    result += 3.0 * p4 / ( zm_sq * ( zm + lm ) ) ;
    //
    const double Am1 = 0.25 * ( Am + lm ) ;
    const double xm1 = 0.25 * ( xm + lm ) ;
    const double ym1 = 0.25 * ( ym + lm ) ;
    const double zm1 = 0.25 * ( zm + lm ) ;
    //
    xm = xm1 ;
    ym = ym1 ;
    zm = zm1 ;
    Am = Am1 ;
    ++n      ;
    Qn /= 4  ;
    p4 /= 4  ;
    
  }
  //
  const double iA4 = p4 / Am  ;
  const double X   = dx * iA4 ;
  const double Y   = dy * iA4 ;
  const double Z   = - ( X + Y ) / 3   ;
  //
  const double XY = X * Y ;
  const double Z2 = Z * Z ;
  //
  const double E2 =       XY - 6 * Z2 ;
  const double E3 = ( 3 * XY - 8 * Z2 ) * Z  ;
  const double E4 = 3 * ( XY -     Z2 ) * Z2 ;
  const double E5 = XY * Z2 * Z ;
  //
  static const double s_c1  = -3.  /  14 ;
  static const double s_c2  = +1.  /   6 ;
  static const double s_c3  = +9.  /  88 ;
  static const double s_c4  = -3.  /  22 ;
  static const double s_c5  = -9.  /  52 ;
  static const double s_c6  = +3.  /  26 ;
  static const double s_c7  = -1.  /  16 ;
  static const double s_c8  = +3.  /  40 ;
  static const double s_c9  = +3.  /  20 ;
  static const double s_c10 = +45. / 272 ;
  static const double s_c11 = -9.  /  68 ;
  //
  result += iA4 / std::sqrt ( Am ) * ( 1 
                                       + s_c1  * E2 
                                       + s_c2  * E3 
                                       + s_c3  * E2 * E2 
                                       + s_c4  * E4 
                                       + s_c5  * E2 * E3 
                                       + s_c6  * E5   
                                       + s_c7  * E2 * E2 * E2 
                                       + s_c8  * E3 * E3  
                                       + s_c9  * E2 * E4  
                                       + s_c10 * E2 * E2 * E3   
                                       + s_c11 * ( E3 * E4 + E2 * E5 ) ) ;
  //
  return result ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_G (x,y,z) = \rac{1}{4} \int_0^{+\infty} 
 *  \left[  (t+x)(t+y)\right]^{-1/2} 
 *   \left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right)  dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RG 
( const double x , 
  const double y , 
  const double z ) 
{
  //
  if ( s_zero  ( z ) ) { return carlson_RG ( z , x , y ) ; }
  //
  return ( z * carlson_RF ( x , y , z ) 
           - ( x - z ) * ( y - z ) * carlson_RD ( x , y , z ) / 3   
           + std::sqrt ( x * y / z ) ) / 2 ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_G (x,y,z) = \rac{1}{4} \int_0^{+\infty} 
 *  \left[  (t+x)(t+y)\right]^{-1/2} 
 *   \left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right)  dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RG_gsl 
( const double x , 
  const double y , 
  const double z ) 
{
  //
  if ( s_zero  ( z ) ) { return carlson_RG_gsl ( z , x , y ) ; }
  //
  return ( z * carlson_RF_gsl ( x , y , z ) 
           - ( x - z ) * ( y - z ) * carlson_RD_gsl ( x , y , z ) / 3   
           + std::sqrt ( x * y / z ) ) / 2 ;
}
// ===========================================================================
/*  Symmetric Carlson form 
 *  \f[ R_F(x,y) = R_F(x,y,0)\f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RF 
( const double x , 
  const double y ) 
{
  if ( x <= 0 || y <= 0 ) 
  {
    //
    gsl_error ( "Invalid arguments for Ostap::Math::carlson_RF(2)" , __FILE__ , __LINE__ , GSL_EDOM ) ;
    return std::numeric_limits<double>::quiet_NaN();
    //
  }
  else if (  s_equal ( x , y ) ) { return 0.5 * M_PI / std::sqrt ( x ) ; }
  //
  double xm = std::sqrt ( x ) ;
  double ym = std::sqrt ( y ) ;
  /// precision 
  static const double s_R = 2 * std::sqrt ( 2 * GSL_DBL_EPSILON ) ;
  //
  for ( unsigned int n = 0 ; n < 10000 ; ++n ) 
  {
    if ( 0 < n && std::abs ( xm - ym ) < s_R * std::abs ( xm ) ) 
    { return M_PI / ( xm + ym ) ; }
    const double xm1 = 0.5 *     ( xm + ym )     ;
    const double ym1 = std::sqrt ( xm * ym ) ;
    xm  = xm1 ;
    ym  = ym1 ;  
  }
  //
  gsl_error ( "Too many iterations for Ostap::Math::carlson_RF(2)" , __FILE__ , __LINE__ , GSL_EMAXITER ) ;
  return M_PI / ( xm + ym ) ;
  //
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_G(x,y) = R_G(x,y,0)\f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RG 
( const double x , 
  const double y ) 
{
  if ( x <= 0 || y <= 0 ) 
  {
    //
    gsl_error ( "Invalid arguments for Ostap::Math::carlson_RG(2)" , __FILE__ , __LINE__ , GSL_EDOM ) ;
    return std::numeric_limits<double>::quiet_NaN();
    //
  }
  //
  double xm = std::sqrt ( x ) ;
  double ym = std::sqrt ( y ) ;
  /// precision 
  static const double s_R = 2 * std::sqrt ( 2 * GSL_DBL_EPSILON ) ;
  //
  double result = ( xm + ym ) * ( xm + ym ) ;
  //
  for ( unsigned int n = 0 ; n < 10000 ; ++n ) 
  {
    if ( 0 < n && std::abs ( xm - ym ) < s_R * std::abs ( xm ) ) 
    { return 0.125 * result * M_PI / ( xm + ym ) ; }
    //
    const double xm1 = 0.5 *     ( xm + ym )     ;
    const double ym1 = std::sqrt ( xm * ym ) ;
    xm      = xm1 ;
    ym      = ym1 ;
    result -= std::pow ( xm - ym , 2 ) * std::pow ( 2.0 , n + 1 ) ; 
  }
  //
  gsl_error ( "Too many iterations for Ostap::Math::carlson_RG(2)" , __FILE__ , __LINE__ , GSL_EMAXITER ) ;
  return 0.125 * result * M_PI / ( xm + ym ) ;
  //
}




// ============================================================================
// Symmetric Carlson forms via plain integrtaion
// (just for cross-checks, validation and debugging
// ============================================================================
namespace 
{
  // =========================================================================
  /// helper structure to calcualte RF via the plain integration
  struct RF 
  {
    // ======================================================================
    RF ( const double x , const double y , const double z ) 
      : m_x ( x ) , m_y ( y ) , m_z ( z ) {}
    // ======================================================================
    // RF 
    double operator () ( const double t ) const 
    { return 0.5 / std::sqrt ( ( t + m_x ) * ( t + m_y ) * ( t + m_z ) ) ; }
    // ======================================================================
    double m_x ;
    double m_y ;
    double m_z ;    
    // ======================================================================
  };  
  // =========================================================================
  /// helper structure to calcualte RJ via the plain integration
  struct RJ 
  {
    // ======================================================================
    RJ ( const double x , const double y , const double z , const double p ) 
      : m_x ( x ) , m_y ( y ) , m_z ( z ) , m_p ( p ) {}
    // ======================================================================
    // RJ 
    double operator () ( const double t ) const 
    { return 1.5 / std::sqrt ( ( t + m_x ) * ( t + m_y ) * ( t + m_z ) ) / ( t + m_p )  ; }
    // ======================================================================
    double m_x ;
    double m_y ;
    double m_z ;    
    double m_p ;    
    // ======================================================================
  };  
  // =========================================================================
  /// helper structure to calcualte RC via the plain integration
  struct RC
  {
    // ======================================================================
    RC ( const double x , const double y ) 
      : m_x ( x ) , m_y ( y ) {}
    // ======================================================================
    // RC 
    double operator () ( const double t ) const 
    { return 0.5 / ( std::sqrt ( t + m_x ) * ( t + m_y ) ) ; }
    // ======================================================================
    double m_x ;
    double m_y ;
    // ======================================================================
  };  
  // =========================================================================
  /// helper structure to calcualte RD via the plain integration
  struct RD 
  {
    RD ( const double x , const double y , const double z ) 
      : m_x ( x ) , m_y ( y ) , m_z ( z ) {}
    // ======================================================================
    // RD 
    double operator () ( const double t ) const 
    { return 1.5 / std::sqrt ( ( t + m_x ) * ( t + m_y ) * ( t + m_z ) ) / ( t + m_z )  ; }
    // ======================================================================
    double m_x ;
    double m_y ;
    double m_z ;    
    // ======================================================================
  };  
  // =========================================================================
  /// helper structure to calcualte RG via the plain integration
  struct RG 
  {
    // ======================================================================
    RG ( const double x , const double y , const double z ) 
      : m_x ( x ) , m_y ( y ) , m_z ( z ) {}
    // RG 
    double operator () ( const double t ) const 
    {
      double sum = m_x / ( t + m_x ) + m_y / ( t + m_y ) + m_z / ( t + m_z ) ;
      return 0.25 * t / std::sqrt ( ( t + m_x ) * ( t + m_y ) * ( t + m_z ) ) * sum ; 
    }
    // ======================================================================
    double m_x ;
    double m_y ;
    double m_z ;    
    // ======================================================================
  };  
  // =========================================================================
  /// statis integrtaion workspace 
  const Ostap::Math::WorkSpace s_workspace { 100000 } ;
  // =========================================================================
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_F(x,y,z) = \int_0^{+\infty} 
 *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2}dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RF_int 
( const double x , 
  const double y , 
  const double z ) 
{
  //
  static const Ostap::Math::GSL::Integrator1D<RF> s_integrator {} ;
  static char s_message[] = "Integral(carlson_RF)" ;
  //
  const double  a1 = 0.1 ;
  const double  a2 = 20  ;
  //
  const RF rf { x , y , z } ;
  const auto F = s_integrator.make_function ( &rf ) ;
  //
  int    ierror   = 0   ;
  double result1  = 0.0 ;
  double error1   = 1.0 ;
  std::tie ( ierror , result1 , error1 ) = s_integrator.qagiu_integrate
    ( &F     ,
      a2     ,                       // low limit 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION        ,          // absolute precision
      s_RPRECISION        ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;  
  //
  double result2  = 0.0 ;
  double error2   = 1.0 ;
  std::tie ( ierror , result2 , error2 ) = s_integrator.qagp_integrate
    ( &F     ,
      0      ,  a1              ,    // integrtaion limits 
      std::vector<double>()     ,    // no internal singular points 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION        ,          // absolute precision
      s_RPRECISION        ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  double result3  = 0.0 ;
  double error3   = 1.0 ;
  std::tie ( ierror , result3 , error3 ) = s_integrator.qag_integrate
    ( &F     ,
      a1     , a2         ,          // integration limits 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION        ,          // absolute precision
      s_RPRECISION        ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //

  return result1 + result2 + result3 ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_J(x,y,z,p) = \int_0^{+\infty} 
 *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2} (t+p) dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RJ_int
( const double x , 
  const double y , 
  const double z , 
  const double p ) 
{
  //
  static const Ostap::Math::GSL::Integrator1D<RJ> s_integrator {} ;
  static char s_message[] = "Integral(carlson_RJ)" ;
  //
  const double  a1 = 0.1 ;
  const double  a2 = 20  ;
  //
  const RJ rj { x , y , z , p } ;
  const auto F = s_integrator.make_function ( &rj ) ;
  //
  int    ierror   = 0   ;
  double result1  = 0.0 ;
  double error1   = 1.0 ;
  std::tie ( ierror , result1 , error1 ) = s_integrator.qagiu_integrate
    ( &F     ,
      a2     ,                       // low limit 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;  
  //
  double result2  = 0.0 ;
  double error2   = 1.0 ;
  std::tie ( ierror , result2 , error2 ) = s_integrator.qagp_integrate
    ( &F     ,
      0      ,  a1              ,    // integrtaion limits 
      std::vector<double>()     ,    // no internal singular points 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  double result3  = 0.0 ;
  double error3   = 1.0 ;
  std::tie ( ierror , result3 , error3 ) = s_integrator.qag_integrate
    ( &F     ,
      a1     , a2         ,          // integration limits 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result1 + result2 + result3 ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_C(x,y) = R_F(x,y,y) = \int_0^{+\infty} 
 *   (t+x)^{-1/2}(t+y)^{-1}dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 *  For negative y, Cauchy principal value it returned 
 *  \f[ R_C(x,-y) = \left(\frac{x}{x+y}\right)R_C(x+y,y), 0 < y \f] 
 */
// ============================================================================
double Ostap::Math::carlson_RC_int
( const double x , 
  const double y ) 
{
  if ( y < 0 ) 
  { return s_zero ( x ) ? 0.0 : std::sqrt ( x / ( x - y ) ) * carlson_RC_int ( x - y , -y ) ; }
  //
  static const Ostap::Math::GSL::Integrator1D<RC> s_integrator {} ;
  static char s_message[] = "Integral(carlson_RC)" ;
  //
  const double  a1 = 0.1 ;
  const double  a2 = 20  ;
  //
  const RC rc { x , y } ;
  const auto F = s_integrator.make_function ( &rc ) ;
  //
  int    ierror   = 0   ;
  double result1  = 0.0 ;
  double error1   = 1.0 ;
  std::tie ( ierror , result1 , error1 ) = s_integrator.qagiu_integrate
    ( &F     ,
      a2     ,                       // low limit 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;  
  //
  double result2  = 0.0 ;
  double error2   = 1.0 ;
  std::tie ( ierror , result2 , error2 ) = s_integrator.qagp_integrate
    ( &F     ,
      0      ,  a1              ,    // integrtaion limits 
      std::vector<double>()     ,    // no internal singular points 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  double result3  = 0.0 ;
  double error3   = 1.0 ;
  std::tie ( ierror , result3 , error3 ) = s_integrator.qag_integrate
    ( &F     ,
      a1     , a2         ,          // integration limits 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result1 + result2 + result3 ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_D (x,y,z) = R_J(x,y,z,z) = \int_0^{+\infty} 
 *  \left[  (t+x)(t+y)\right]^{-1/2} (t+z)^{-3/2} dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RD_int  
( const double x , 
  const double y , 
  const double z ) 
{
  //
  static const Ostap::Math::GSL::Integrator1D<RD> s_integrator {} ;
  static char s_message[] = "Integral(carlson_RD)" ;
  //
  const double  a1 = 0.1 ;
  const double  a2 = 20  ;
  //
  const RD rd { x , y , z } ;
  const auto F = s_integrator.make_function ( &rd ) ;
  //
  int    ierror   = 0   ;
  double result1  = 0.0 ;
  double error1   = 1.0 ;
  std::tie ( ierror , result1 , error1 ) = s_integrator.qagiu_integrate
    ( &F     ,
      a2     ,                       // low limit 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;  
  //
  double result2  = 0.0 ;
  double error2   = 1.0 ;
  std::tie ( ierror , result2 , error2 ) = s_integrator.qagp_integrate
    ( &F     ,
      0      ,  a1              ,    // integrtaion limits 
      std::vector<double>()     ,    // no internal singular points 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  double result3  = 0.0 ;
  double error3   = 1.0 ;
  std::tie ( ierror , result3 , error3 ) = s_integrator.qag_integrate
    ( &F     ,
      a1     , a2         ,          // integration limits 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result1 + result2 + result3 ;
}
// ============================================================================
/*  Symmetric Carlson form 
 *  \f[ R_G (x,y,z) = \rac{1}{4} \int_0^{+\infty} 
 *  \left[  (t+x)(t+y)\right]^{-1/2} 
 *   \left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right)  dt \f]
 *  @see https://en.wikipedia.org/wiki/Elliptic_integral
 *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
 *                Numerical Algorithms, 10, 1995,  13
 *  @see https://doi.org/10.1007/BF02198293
 *  @see https://arxiv.org/abs/math/9409227
 */
// ============================================================================
double Ostap::Math::carlson_RG_int  
( const double x , 
  const double y , 
  const double z ) 
{
  //
  static const Ostap::Math::GSL::Integrator1D<RG> s_integrator {} ;
  static char s_message[] = "Integral(carlson_RG)" ;
  //
  const double  a1 = 0.1 ;
  const double  a2 = 20  ;
  //
  const RG rg { x , y , z } ;
  const auto F = s_integrator.make_function ( &rg ) ;
  //
  int    ierror   = 0   ;
  double result1  = 0.0 ;
  double error1   = 1.0 ;
  std::tie ( ierror , result1 , error1 ) = s_integrator.qagiu_integrate
    ( &F     ,
      a2     ,                       // low limit 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;  
  //
  double result2  = 0.0 ;
  double error2   = 1.0 ;
  std::tie ( ierror , result2 , error2 ) = s_integrator.qagp_integrate
    ( &F     ,
      0      ,  a1              ,    // integrtaion limits 
      std::vector<double>()     ,    // no internal singular points 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  double result3  = 0.0 ;
  double error3   = 1.0 ;
  std::tie ( ierror , result3 , error3 ) = s_integrator.qag_integrate
    ( &F     ,
      a1     , a2         ,          // integration limits 
      workspace ( s_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      s_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //

  return result1 + result2 + result3 ;
}
// ============================================================================



// ============================================================================
// Jacobi elliptic functions
// ============================================================================

// ============================================================================
/* Elliptic amplitude \f$ \mathrm{am}(u,m)=\phi\f$, where 
 *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
 *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
 */
// ============================================================================
double Ostap::Math::am
( const double u ,
  const double m )
{
  //
  if      ( s_zero  ( u    )  ) { return      0   ; } 
  else if ( s_zero  ( m     ) ) { return      u   ; }
  else if ( s_zero  ( 1 - m ) ) { return gd ( u ) ; }
  else if ( s_equal ( m , 1 ) ) { return gd ( u ) ; }
  else if (  m < 0 || 1 < m   ) { return std::numeric_limits<double>::quiet_NaN(); }
  else if ( u < 0             ) { return - am ( -u , m ) ; }
  //
  // some reduction
  if ( 0.5 * M_PI < std::abs ( u ) )
    {
      const double K = elliptic_Km ( m ) ;      
      if ( 2 * K < std::abs ( u ) )
	{
	  double ur, nc ;
	  std::tie ( ur , nc ) = reduce ( u , -2 * K , +2 * K ) ;
	  if ( !s_zero ( nc ) ) { return am ( ur , m ) + nc * 2 * M_PI ; }
	}
    }
  //
  typedef std::array<long double,3>     ITEM  ;
  typedef std::vector<ITEM>             ITEMS ;
  typedef ITEMS::const_reverse_iterator CRI   ;
  //
  constexpr unsigned short NR = 50 ;
  ITEMS items{} ; items.reserve ( NR ) ;
  //
  long double a = 1.0L ;
  long double b = std::sqrt ( 1.0L - m ) ;
  long double c = a - b ;
  //
  for  ( unsigned short n = 0 ; n < NR ; ++n )
    {
      long double a1 = 0.5L *    ( a + b ) ;
      long double b1 = std::sqrt ( a * b ) ;
      long double c1 = 0.5L *    ( a - b ) ;
      a = a1 ;
      b = b1 ;
      c = c1 ;
      //
      items.push_back ( ITEM { a , b , c } ) ;
      //
      if      ( 3 <= n && s_zero  ( c1 * 1000 ) ) { break ; }
      else if ( 5 <= n && s_equal ( a , b     ) ) { break ; }
    }
  const unsigned short NN = items.size() ;
  //
  long double phi = a * u * std::pow ( 2 , NN ) ;
  //
  for ( CRI it = items.rbegin() ; it != items.rend() ; ++it )
    {
      a = std::get<0>(*it) ;
      c = std::get<2>(*it) ;
      phi = 0.5L * ( phi + std::asin ( c * std::sin ( phi ) / a ) ) ;
    }
  //
  return phi ;
}
// ============================================================================
/*  Sine elliptic amplitude \f$ \mathrm{sn} (u,m)=\sin \mathrm{am} ( u, m ) \f$, where 
 *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
 *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
 */
// ============================================================================
double Ostap::Math::sn
( const double u ,
  const double m )
{
  //
  if      ( s_zero  ( u     ) ) { return             u   ; }
  else if ( s_zero  ( m     ) ) { return std::sin  ( u ) ; }
  else if ( s_zero  ( 1 - m ) ) { return std::tanh ( u ) ; }
  else if ( s_equal ( m , 1 ) ) { return std::tanh ( u ) ; }
  else if ( u < 0             ) { return - sn ( -u , m ) ; }
  else if (  m < 0 || 1 < m   ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  // reduce the argument
  if ( 0.5 * M_PI < u ) 
    {
      const long double KK = elliptic_Km ( m ) ;
      if      ( s_equal ( u , 2 * KK ) ) { return 2 * KK - u ; }
      else if ( s_equal ( u , 4 * KK ) ) { return u - 4 * KK ; }
      else if ( 4 * KK < u )
	{      
	  const double uk = u - 4 * KK * ( std::floor ( u * 0.25L / KK ) + 1 ) ;
	  return sn ( uk , m ) ;
	}
      else if ( 2 * KK < u )
	{      
	  const double uk = 4 * KK - u ;
	  return - sn ( uk , m ) ;
	}
      else if ( KK < u )
	{
	  const double uk = 2 * KK - u ;
	  return sn ( uk , m ) ;      
	}
    }
  //
  typedef std::array<long double,3>     ITEM  ;
  typedef std::vector<ITEM>             ITEMS ;
  typedef ITEMS::const_reverse_iterator CRI   ;
  //
  constexpr unsigned short NR = 50 ;
  ITEMS items{} ; items.reserve ( NR ) ;
  //
  long double a = 1.0L ;
  long double b = std::sqrt ( 1.0L - m ) ;
  long double c = a - b ;
  //
  for  ( unsigned short n = 0 ; n < NR ; ++n )
    {
      long double a1 = 0.5L *    ( a + b ) ;
      long double b1 = std::sqrt ( a * b ) ;
      long double c1 = 0.5  *    ( a - b ) ;
      //
      a = a1 ;
      b = b1 ;
      c = c1 ;
      //
      items.push_back ( ITEM { a , b , c } ) ;
      //
      if      ( 3 <= n && s_zero  ( c1 * 1000 ) ) { break ; }
      else if ( 5 <= n && s_equal ( a , b     ) ) { break ; }
      if ( 5 <= n && s_zero  ( b * 1000 ) ) { break ; }
    }
  const unsigned short NN = items.size() ;
  //
  long double y = a / std::sin ( a * u ) ;
  //
  for ( CRI it = items.rbegin() ; it != items.rend() ; ++it )
    {
      a = std::get<0>(*it) ;
      c = std::get<2>(*it) ;
      y = y + a * c / y    ;
      
    }
  //
  return 1.0 / y ;
}
// ============================================================================
/*  Sine elliptic amplitude \f$ \mathrm{sn} (u,m)=\sin \mathrm{am} ( u, m ) \f$, where 
 *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
 *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
 */
// ============================================================================
double Ostap::Math::sn_
( const double u ,
  const double m )
{
  //
  if      ( s_zero  ( u     ) ) { return             u   ; }
  else if ( s_zero  ( m     ) ) { return std::sin  ( u ) ; }
  else if ( s_zero  ( 1 - m ) ) { return std::tanh ( u ) ; }
  else if ( s_equal ( m , 1 ) ) { return std::tanh ( u ) ; }
  else if ( u < 0             ) { return - sn ( -u , m ) ; }
  else if (  m < 0 || 1 < m   ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  // reduce the argument
  if ( 0.5 * M_PI < u ) 
    {
      const long double KK = elliptic_Km ( m ) ;
      if      ( s_equal ( u , 2 * KK ) ) { return 2 * KK - u ; }
      else if ( s_equal ( u , 4 * KK ) ) { return u - 4 * KK ; }
      else if ( 4 * KK < u )
	{      
	  const double uk = u - 4 * KK * ( std::floor ( u * 0.25L / KK ) + 1 ) ;
	  return sn ( uk , m ) ;
	}
      else if ( 2 * KK < u )
	{      
	  const double uk = 4 * KK - u ;
	  return - sn ( uk , m ) ;
	}
      else if ( KK < u )
	{
	  const double uk = 2 * KK - u ;
	  return sn ( uk , m ) ;      
	}
    }
  //
  typedef std::array<long double,2>     ITEM  ;
  typedef std::vector<ITEM>             ITEMS ;
  typedef ITEMS::const_reverse_iterator CRI   ;
  //
  constexpr unsigned short NR = 50 ;
  ITEMS items{} ; items.reserve ( NR ) ;
  //
  const long double kp = std::sqrt ( 1.0L - m ) ;
  
  long double a = u ;
  long double b = ( 1.0L - kp ) / ( 1.0L + kp ) ;
  //
  items.push_back ( ITEM { a , b } ) ;
  for ( unsigned short n = 0 ; n < NR ; ++n )
    {
      a = a / ( 1.0L + b ) ;
      const long double bb = std::sqrt ( 1.L - b * b ) ;
      b = ( 1.0L - bb ) / ( 1.0L + bb ) ;
      items.push_back ( ITEM { a , b } ) ;
      if ( 5 <= n && s_zero  ( b * 1000 ) ) { break ; }
    }
  
  long double y = std::sin ( a ) ;  
  for ( CRI it = items.rbegin() ; it != items.rend() ; ++it )
    {
      b = std::get<1>(*it) ;
      y = y * ( 1.0L + b ) / ( 1.0L + b * y *y ) ;
    }
  //
  return y ;
}
// ============================================================================
/*  Cosine elliptic amplitude \f$ \mathrm{sn} (u,m)=\cos \mathrm{am} ( u, m ) \f$, where 
 *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
 *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
 */
// ============================================================================
double Ostap::Math::cn
( const double u ,
  const double m )
{
  //
  if      ( s_zero  ( u     ) ) { return             1   ; }
  else if ( s_zero  ( m     ) ) { return std::cos  ( u ) ; }
  else if ( s_zero  ( 1 - m ) ) { return      sech ( u ) ; }
  else if ( s_equal ( m , 1 ) ) { return      sech ( u ) ; }
  else if ( u < 0             ) { return    cn ( u , m ) ; }
  else if (  m < 0 || 1 < m   ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  // reduce the argument
  if ( 0.5 * M_PI < u ) 
    {
      const long double KK = elliptic_Km ( m ) ;
      if      ( s_equal ( u ,     KK ) ) { return KK - u     ; }
      else if ( s_equal ( u , 2 * KK ) ) { return -1         ; }
      else if ( s_equal ( u , 3 * KK ) ) { return u - 3 * KK ; }
      else if ( s_equal ( u , 4 * KK ) ) { return  1         ; }
      else if ( 4 * KK < u )
       	{      
       	  const double uk = u - 4 * KK * ( std::floor ( u * 0.25L / KK ) ) ;
       	  return cn ( uk , m ) ;
       	}
      else if ( 2 * KK < u )
       	{      
       	  const double uk = 4 * KK - u ;
       	  return cn ( uk , m ) ;
       	}
      else if ( KK < u )
       	{
      	  const double uk = 2 * KK - u ;
       	  return - cn ( uk , m ) ;      
       	}
    }
  //
  return std::cos ( am ( u , m ) ) ;
}
// ============================================================================
/*  Elliptic delta amplitude \f$ \mathrm{sn} (u,m)=\frac{d}{du} \mathrm{am} ( u, m ) \f$, where 
 *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
 *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
 */
// ============================================================================
double Ostap::Math::dn
( const double u ,
  const double m )
{
  //
  if      ( s_zero  ( u     ) ) { return        1   ; }
  else if ( s_zero  ( m     ) ) { return        1   ; }
  else if ( s_zero  ( 1 - m ) ) { return sech ( u ) ; }
  else if ( s_equal ( m , 1 ) ) { return sech ( u ) ; }
  else if ( u < 0             ) { return cn ( u , m ) ; }
  else if (  m < 0 || 1 < m   ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  const double snu = sn ( u , m ) ;
  return std::sqrt ( 1.0 - m * snu * snu ) ;
}
// ============================================================================
/*  elliptic function 
 *  \f[ \matmrm{sc}\,(u,m) = \frac{ \mathrm{sn} ( u, m) } { \mathrm{cn} ( u , m ) } \f] 
 *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
 */ 
// ============================================================================
double Ostap::Math::sc
( const double u ,
  const double m )
{
  if      ( s_zero  ( u     ) ) { return             u   ; }
  else if ( s_zero  ( m     ) ) { return std::tan  ( u ) ; }
  else if ( s_zero  ( 1 - m ) ) { return std::sinh ( u ) ; }
  else if ( s_equal ( m , 1 ) ) { return std::sinh ( u ) ; }
  else if ( u < 0             ) { return - sc ( -u , m ) ; }
  else if (  m < 0 || 1 < m   ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  // reduce the argument
  if ( 0.5 * M_PI < u ) 
    {
      const long double KK = elliptic_Km ( m ) ;
      if ( 4 * KK < u  ) 
	{      
	  const double uk = u - 4 * KK * ( std::floor ( u * 0.25L / KK ) + 1 ) ;
	  return sn ( uk , m ) ;
	}
    } ;
  //
  const double a = am ( u , m ) ;
  //
  return std::tan ( a ) ;
}
// ============================================================================

// ============================================================================
/* Dilogarithm function (real case) 
 *  \f$ Li_2(x) = - Re \int\limits_0^{x}\draf{\log ( 1-s) } {s} ds  \f$ 
 */
// ============================================================================
double Ostap::Math::dilog ( const double x )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_dilog_e ( x , &result) ;
  if ( ierror ) 
    {
      //
      gsl_error ( "Error from gsl_sf_dilog_e function" , __FILE__ , __LINE__ , ierror ) ;
      if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
      //
    }
  return result.val ;
}
// ============================================================================
/* Dilogarithm function (complex case) 
 *  \f$ Li_2(x) = - \int\limits_0^{x}\draf{\log ( 1-s) } {s} ds  \f$ 
 */
// ============================================================================
std::complex<double>
Ostap::Math::dilog ( const std::complex<double>& z ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result_re ;
  gsl_sf_result result_im ;
  const int ierror = gsl_sf_complex_dilog_e ( std::abs ( z ) ,
                                              std::arg ( z ) ,
                                              &result_re     ,
                                              &result_im     ) ;
  if ( ierror ) 
    {
      //
      gsl_error ( "Error from gsl_sf_dilog_e function" , __FILE__ , __LINE__ , ierror ) ;
      if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
      //
    }
  return std::complex<double>( result_re.val , result_im.val ) ;
}
// ============================================================================


// ============================================================================
/*  Riemann's Zeta function \f$ n\ne 1\f$:
 *  \f$ \zeta ( n ) = \sum_k k^{-n}\f$ 
 */
// ============================================================================
double Ostap::Math::zeta ( const int n  )
{
  if ( 1 == n ) { return std::numeric_limits<double>::quiet_NaN() ; }
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_zeta_int_e ( n , &result) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_zeta_int_e function" , __FILE__ , __LINE__ , ierror ) ;
      if ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
    }
  return result.val ;
}
// ============================================================================
/*  Riemann's Zeta function \f$ s\ne 1\f$:
 *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
 */
// ============================================================================
double Ostap::Math::zeta ( const double s   )
{
  if ( 1 == s || s_equal ( s , 1.0 ) )
    { return std::numeric_limits<double>::quiet_NaN() ; }
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_zeta_e ( s , &result) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_zeta_e function" , __FILE__ , __LINE__ , ierror ) ;
      if ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
    }
  return result.val ;
}
// ============================================================================
/* Riemann's Zeta function minus 1 \f$ n\ne 1\f$:
 *  \f$ f ( n ) = \zeta ( n ) - 1 = -1 + \sum_k k^{-n}\f$ 
 */
// ============================================================================
double Ostap::Math::zetam1 ( const int n  )
{
  if ( 1 == n ) { return std::numeric_limits<double>::quiet_NaN() ; }
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_zetam1_int_e ( n , &result) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_zetam1_int_e function" , __FILE__ , __LINE__ , ierror ) ;
      if ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
    }
  return result.val ;
}
// ============================================================================
/*  Riemann's Zeta function minus 1 \f$ s \ne 1\f$:
 *  \f$ f ( s ) = \zeta ( s ) - 1 = -1 + \sum_k k^{-s}\f$ 
 */
// ============================================================================
double Ostap::Math::zetam1 ( const double s   )
{
  if ( 1 == s || s_equal ( s , 1.0 ) )
    { return std::numeric_limits<double>::quiet_NaN() ; }
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_zetam1_e ( s , &result) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_zetam1_e function" , __FILE__ , __LINE__ , ierror ) ;
      if ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
    }
  return result.val ;
}
// ============================================================================

// ============================================================================
/* Hurwitz Zeta function 
 *  \f$ zeta ( s , q ) = \sum_k  ( k + q )^{-s}\f$
 *  - \f$ 1 < s \f$
 *  - \f$ 0 < q \f$
 */
// ============================================================================
double Ostap::Math::hurwitz 
( const double s ,
  const double q )
{
  if ( s <= 1 || q <= 0 )
    { return std::numeric_limits<double>::quiet_NaN() ; }
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_hzeta_e ( s , q , &result) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_hzeta_e function" , __FILE__ , __LINE__ , ierror ) ;
      if ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
    }
  return result.val ;
}
// ============================================================================


// ============================================================================
/* Dirichlet's Eta function 
 *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
 */
// ============================================================================
double Ostap::Math::eta ( const int n  )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_eta_int_e ( n , &result) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_eta_int_e function" , __FILE__ , __LINE__ , ierror ) ;
      if ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
    }
  return result.val ;
}
// ============================================================================
/*  Dirichlet's Eta function 
 *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
 */
// ============================================================================
double Ostap::Math::eta ( const double s   )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_eta_e ( s , &result) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_eta_e function" , __FILE__ , __LINE__ , ierror ) ;
      if ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN() ; }
    }
  return result.val ;
}
// ============================================================================

// ============================================================================
/* Harmonic number 
 *  \f$ H_n = \sum_{k=1}^{n}  \frac{1}{k} \f$ 
 */
// ============================================================================
double Ostap::Math::harmonic ( const unsigned int n )
{
  long double result = 0 ;
  if ( n <= 100 )
    {
      for ( unsigned int k = 1 ; k <= n ; ++k ) { result += 1.0L / k ; }
      return result ; 
    }
  return M_EULER + psi ( 1.0 + n )  ;
}
// ============================================================================
/* Harmonic number (generalized) 
 *  \f$ H(x) = \gamma + \psi ( 1 + x ) \f$ 
 */
// ============================================================================
double Ostap::Math::harmonic ( const double x  )
{ return M_EULER + psi ( 1.0 + x ) ; }

// ============================================================================
/* Mill's ratio for normal distribution
 *  - \f$ m (x) = \frac{1 - \Phi(x)}{\phi(x)}\f$  
 *  - \f$ m(x) =  \sqrt{ \frac{\pi}{2} } \mathrm{e}^{ \frac{x^2}{2}} \mathrm{erfc}{ \frac{x}{\sqrt{2}}}\f$ 
 *  - \f$ m(x) =  \sqrt{ \frac{\pi}{2} } \mathrm{erfcx}{ \frac{x}{\sqrt{2}}}\f$ 
 *  @see https://en.wikipedia.org/wiki/Mills_ratio
 *  @see Ostap::Math::erfcx 
 */
// ============================================================================
double Ostap::Math::mills_normal ( const double x ) 
{ return s_SQRTPIHALF * Ostap::Math::erfcx ( s_SQRT2i * x ) ; }
// ============================================================================
/* Product of the Gaussian PDF and Millt's ratio 
 *  \f$ f(a,b) = \phi(a) R(b) \f$
 */
// ============================================================================
double Ostap::Math::gauss_mills
( const double a , 
  const double b ) 
{
  return 
    0 <= b ? 
    Ostap::Math::gauss_pdf ( a )             * Ostap::Math::mills_normal ( b ) : 
    std::exp ( 0.5 * ( b - a ) * ( b + a ) ) * Ostap::Math::erfc ( b * s_SQRT2i ) * 0.5 ; 
}
// ============================================================================
/*  get the intermediate polynom \f$ g_l (x)\f$ used for the calculation of 
 *  the angular-momentum Blatt-Weisskopf centrifugal-barrier factor 
 *  @see S.U.Chung "Formulas for Angular-Momentum Barrier Factors", BNL-QGS-06-01
 *  @see https://physique.cuso.ch/fileadmin/physique/document/2015_chung_brfactor1.pdf
 *
 *  The complex-valued polynomials \f$ g_l(x) \f$  with integer 
 *  coefficients can be written as 
 *  \f[ g_l(x) = \sum_{k=0}^{l} a_{lk}(-ix)^{l-k}\f] with
 *  \f$ a_{lk} = \frac{(l+k)!}{2^k k! (l-k)!}\f$ and \f$a_{l0}=1\f$.
 *
 *  It satisfies the recurrent relation
 *  \f[ g_{l+1}(x) = (2l+1)g_l(x) -  x^2 g_{l-1}(x)\f] 
 *  with the initial values of \f$ g_0(x) \equiv 1 \f$ 
 *  and \f$ g_1(x) \equiv -ix + 1\f$.
 *  This recurrense relation is used for the actual calculation.
 *
 *  @param  x  the value of scaled relative momentum 
 *  @param  l  the orbital momentum 
 *  @return the value of \f$ g_l(x) \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-10-13
 *  @see Ostap::Math::barrier_factor 
 *  @see Ostap::Math::barrier_absg 
 */
// ============================================================================
std::complex<double>
Ostap::Math::barrier_g 
( const double       x , 
  const unsigned int l ) 
{
  //
  if      ( 0 == l ) {  return                        1        ; }
  else if ( 1 == l ) {  return std::complex<double> ( 1 , -x ) ; }
  //
  // real part 
  const long double re_g = Ostap::Math::Clenshaw::term 
    ( x  , 
      l  , 
      [] ( const unsigned int    k    , 
           const long double  /* t */ ) -> long double { return  2 * k+1 ; } , 
      [] ( const unsigned int /* k */ , 
           const long double     t    ) -> long double { return -t * t   ; } , 
      [] ( const long double  /* t */ ) -> long double { return  1       ; } , 
      [] ( const long double  /* t */ ) -> long double { return  1       ; } ) ;
  //
  // imaginary part 
  const long double im_g = Ostap::Math::Clenshaw::term 
    ( x  , 
      l  , 
      [] ( const unsigned int    k    , 
           const long double  /* t */ ) -> long double { return  2 * k+1 ; } , 
      [] ( const unsigned int /* k */ , 
           const long double     t    ) -> long double { return -t * t   ; } , 
      [] ( const long double  /* t */ ) -> long double { return  0       ; } , 
      [] ( const long double     t    ) -> long double { return -t       ; } ) ;
  //
  return std::complex<double> ( re_g , im_g ) ;
}
// ============================================================================


// ============================================================================
/* regular Bessel function of the first kind
 *  \f$ J_n(x)\f$ 
 * @see https://en.wikipedia.org/wiki/Bessel_function
 * @see gsl_sf_bessel_J0
 * @see gsl_sf_bessel_J1
 * @see gsl_sf_bessel_Jn
 */
// ============================================================================
double Ostap::Math::bessel_Jn
( const int    n , 
  const double x ) 
{
  gsl_sf_result result ;
  const int ierror = 
    ( 0 == n ) ? gsl_sf_bessel_J0_e ( x , &result ) :
    ( 1 == n ) ? gsl_sf_bessel_J1_e ( x , &result ) :
    gsl_sf_bessel_Jn_e ( n , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( 0 == n ? "Error from gsl_sf_bessel_J0_e" : 
                1 == n ? "Error from gsl_sf_bessel_J1_e" : 
                "Error from gsl_sf_bessel_Jn_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  //
  return result.val ;  
}
// ============================================================================
/* Irregular Bessel function of the first kind
 *  \f$ Y_n(x)\f$ 
 * @see https://en.wikipedia.org/wiki/Bessel_function
 * @see gsl_sf_bessel_Y0
 * @see gsl_sf_bessel_Y1
 * @see gsl_sf_bessel_Yn
 */
// ============================================================================
double Ostap::Math::bessel_Yn
( const int    n , 
  const double x ) 
{
  gsl_sf_result result ;
  const int ierror = 
    ( 0 == n ) ? gsl_sf_bessel_Y0_e ( x , &result ) :
    ( 1 == n ) ? gsl_sf_bessel_Y1_e ( x , &result ) :
    gsl_sf_bessel_Yn_e ( n , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( 0 == n ? "Error from gsl_sf_bessel_Y0_e" : 
                1 == n ? "Error from gsl_sf_bessel_Y1_e" : 
                "Error from gsl_sf_bessel_Yn_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  //
  return result.val ;  
}


// ============================================================================
/** Regular Bessel function 
 *  \f$ J_{\nu}(x) \f$ 
 *  @see https://en.wikipedia.org/wiki/Bessel_function
 *  @see gsl_sf_bessel_Jnu_e 
 */
// ============================================================================
double Ostap::Math::bessel_Jnu 
( const double nu , 
  const double x  ) 
{ 
  //
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round ( nu ) ;
    return Ostap::Math::bessel_Jn ( n , x ) ;
  }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_bessel_Jnu_e ( nu , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_bessel_Jnu_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}

// ============================================================================
/** Itregular Bessel function
 *  \f$ Y_{\nu}(x) \f$ 
 *  @see https://en.wikipedia.org/wiki/Bessel_function
 *  @see gsl_sf_bessel_Ynu_e 
 */
// ============================================================================
double Ostap::Math::bessel_Ynu 
( const double nu , 
  const double x  ) 
{ 
  //
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round ( nu ) ;
    return Ostap::Math::bessel_Yn ( n , x ) ;
  }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_bessel_Ynu_e ( nu , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_bessel_Ynu_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}

// ============================================================================
/*  modified Bessel function of the fist kind  
 *  \f$ I_n(x) \f$ for \f$ x>0 \f$
 * @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1 
 * @see gsl_sf_bessel_I0
 * @see gsl_sf_bessel_I1
 * @see gsl_sf_bessel_I2
 */
// ============================================================================
double Ostap::Math::bessel_In ( const int n , const double x ) 
{
  gsl_sf_result result ;
  const int ierror = 
    ( 0 == n ) ? gsl_sf_bessel_I0_e ( x , &result ) :
    ( 1 == n ) ? gsl_sf_bessel_I1_e ( x , &result ) :
    gsl_sf_bessel_In_e ( n , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( 0 == n ? "Error from gsl_sf_bessel_I0_e" : 
                1 == n ? "Error from gsl_sf_bessel_I1_e" : 
                "Error from gsl_sf_bessel_In_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  //
  return result.val ;  
}

// ============================================================================
/*  modified Bessel function of the second kind 
 *  \f$ K_n(x) \f$ for \f$ x>0 \f$
 *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
 *  @see gsl_sf_bessel_K0_e 
 *  @see gsl_sf_bessel_K1_e 
 *  @see gsl_sf_bessel_Kn_e 
 */
// ============================================================================
double Ostap::Math::bessel_Kn 
( const int    n , 
  const double x ) 
{ 
  gsl_sf_result result ;
  const int ierror = 
    ( 0 == n ) ? gsl_sf_bessel_K0_e ( x , &result ) :
    ( 1 == n ) ? gsl_sf_bessel_K1_e ( x , &result ) :
    gsl_sf_bessel_Kn_e ( n , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( 0 == n ? "Error from gsl_sf_bessel_K0_e" : 
                1 == n ? "Error from gsl_sf_bessel_K1_e" : 
                "Error from gsl_sf_bessel_Kn_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================
/** modified Bessel function of the second kind  
 *  \f$ I_{\nu}(x) \f$ for \f$ x>0, \nu>0 \f$
 *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
 *  @see gsl_sf_bessel_Inu_e 
 */
// ============================================================================
double Ostap::Math::bessel_Inu 
( const double nu , 
  const double x  ) 
{ 
  //
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round ( nu ) ;
    return Ostap::Math::bessel_In ( n , x ) ;
  }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_bessel_Inu_e ( nu , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_bessel_Inu_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================
/** modified Bessel function of the second kind  
 *  \f$ K_{\nu}(x) \f$ for \f$ x>0, \nu>0 \f$
 *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
 *  @see gsl_sf_bessel_Knu_e 
 */
// ============================================================================
double Ostap::Math::bessel_Knu 
( const double nu , 
  const double x  ) 
{ 
  //
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round ( nu ) ;
    return Ostap::Math::bessel_Kn ( n , x ) ;
  }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_bessel_Knu_e ( std::abs ( nu ) , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_bessel_Knu_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}



// ============================================================================
/** scaled modified Bessel function of the second kind 
 *  \f$ \mathrm{e}^x K_n(x) \f$ for \f$ x>0 \f$
 *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
 *  @see gsl_sf_bessel_K0_scaled_e 
 *  @see gsl_sf_bessel_K1_scaled_e 
 *  @see gsl_sf_bessel_Kn_scaled_e 
 */
// ============================================================================
double Ostap::Math::bessel_Kn_scaled 
( const int    n , 
  const double x ) 
{ 
  gsl_sf_result result ;
  const int ierror = 
    ( 0 == n ) ? gsl_sf_bessel_K0_scaled_e ( x , &result ) :
    ( 1 == n ) ? gsl_sf_bessel_K1_scaled_e ( x , &result ) :
    gsl_sf_bessel_Kn_scaled_e ( n , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( ( 0 == n ) ? "Error from gsl_sf_bessel_K0_scaled_e" : 
                ( 1 == n ) ? "Error from gsl_sf_bessel_K1_scaled_e" : 
                "Error from gsl_sf_bessel_Kn_scaled_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================
/** scaled modified Bessel function of the second kind 
 *  \f$ \mathrm{e}^x K_{\nu}(x) \f$ for \f$ x>0, \nu>0 \f$
 *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
 *  @see gsl_sf_bessel_Knu_scaled_e 
 */
// ============================================================================
double Ostap::Math::bessel_Knu_scaled 
( const double nu , 
  const double x  ) 
{ 
  //
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round ( nu ) ;
    return Ostap::Math::bessel_Kn_scaled ( n , x ) ;
  }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_bessel_Knu_scaled_e ( std::abs ( nu ) , x , &result )  ;
  //
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_bessel_Knu_scaled_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================



// ============================================================================
/*  derivative for the  regular Bessel function of the first kind
 *  - \f$ J_0^{\prime}(x) =  0 J_1(s)\f$
 *  - \f$ J_{n}^{\prime}(x) = \frac{1}{2}\left(J_{n-1}(x) - J_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_Jn
 */
// ============================================================================
double Ostap::Math::der_bessel_Jn
( const int    n , 
  const double x )
{
  return
    0 == n ? -Ostap::Math::bessel_Jn ( 1 , x ) :
    0.5 * ( Ostap::Math::bessel_Jn ( n - 1 , x ) -
            Ostap::Math::bessel_Jn ( n + 1 , x ) ) ; 
}
// ============================================================================
/*  derivative for the  regular Bessel function of the first kind
 *  - \f$ J_0^{\prime}(x) =  0 J_1(s)\f$
 *  - \f$ J_{n}^{\prime}(x) = \frac{1}{2}\left(J_{n-1}(x) - J_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_Jnu
 */
// ============================================================================
double Ostap::Math::der_bessel_Jnu
( const double nu , 
  const double x  )
{
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round  ( nu ) ;
    return Ostap::Math::der_bessel_Jn ( n , x ) ;
  }
  return
    s_zero ( nu ) ? -Ostap::Math::bessel_Jn ( 1 , x ) :
    0.5 * ( Ostap::Math::bessel_Jnu ( nu - 1 , x ) -
            Ostap::Math::bessel_Jnu ( nu + 1 , x ) ) ; 
}
// ============================================================================
/*  derivative for the irregular Bessel function of the first kind 
 *  - \f$ Y_0^{\prime}(x) =  - Y_1(s)\f$
 *  - \f$ Y_{n}^{\prime}(x) = \frac{1}{2}\left(Y_{n-1}(x) - Y_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_Yn
 */
// ============================================================================
double Ostap::Math::der_bessel_Yn
( const int    n , 
  const double x )
{
  return
    0 == n ? -Ostap::Math::bessel_Yn ( 1 , x ) :
    0.5 * ( Ostap::Math::bessel_Yn ( n - 1 , x ) -
            Ostap::Math::bessel_Yn ( n + 1 , x ) ) ; 
}
// ============================================================================
/*  derivative for the irregular Bessel function of the first kind 
 *  - \f$ Y_0^{\prime}(x) =  - Y_1(s)\f$
 *  - \f$ Y_{n}^{\prime}(x) = \frac{1}{2}\left(Y_{n-1}(x) - Y_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_Yn
 */
// ============================================================================
double Ostap::Math::der_bessel_Ynu
( const double nu , 
  const double x  )
{
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round  ( nu ) ;
    return Ostap::Math::der_bessel_Yn ( n , x ) ;
  }
  return
    s_zero ( nu ) ? -Ostap::Math::bessel_Yn ( 1 , x ) :
    0.5 * ( Ostap::Math::bessel_Ynu ( nu - 1 , x ) -
            Ostap::Math::bessel_Ynu ( nu + 1 , x ) ) ; 
}
// ============================================================================
/*  derivative for the  modified Bessel function
 *  - \f$ I_0^{\prime}(x) =  I_1(s)\f$
 *  - \f$ I_{n}^{\prime}(x) = \frac{1}{2}\left(I_{n-1}(x) + I_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_In
 */
// ============================================================================
double Ostap::Math::der_bessel_In
( const int    n , 
  const double x )
{
  return
    0 == n ? Ostap::Math::bessel_In ( 1 , x ) :
    0.5 * ( Ostap::Math::bessel_In ( n - 1 , x ) +
            Ostap::Math::bessel_In ( n + 1 , x ) ) ; 
}
// ============================================================================
/*  derivative for the  modified Bessel function
 *  - \f$ I_0^{\prime}(x) =  I_1(s)\f$
 *  - \f$ I_{n}^{\prime}(x) = \frac{1}{2}\left(I_{n-1}(x) + I_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_Inu
 */
// ============================================================================
double Ostap::Math::der_bessel_Inu
( const double nu , 
  const double x  )
{
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round  ( nu ) ;
    return Ostap::Math::der_bessel_In ( n , x ) ;
  }
  return
    s_zero ( nu ) ?  Ostap::Math::bessel_In ( 1 , x ) :
    0.5 * ( Ostap::Math::bessel_Inu ( nu - 1 , x ) +
            Ostap::Math::bessel_Inu ( nu + 1 , x ) ) ; 
}
// ============================================================================
/*  derivative for the  modified Bessel function
 *  - \f$ K_0^{\prime}(x) = -K_1(s)\f$
 *  - \f$ K_{n}^{\prime}(x) = -\frac{1}{2}\left(K_{n-1}(x) + K_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_Kn
 */
// ============================================================================
double Ostap::Math::der_bessel_Kn
( const int    n , 
  const double x )
{
  return
    0 == n ? - Ostap::Math::bessel_Kn ( 1 , x ) :
    -0.5 * ( Ostap::Math::bessel_Kn ( n - 1 , x ) +
             Ostap::Math::bessel_Kn ( n + 1 , x ) ) ; 
}
// ============================================================================
/*  derivative for the  modified Bessel function
 *  - \f$ K_0^{\prime}(x) = -K_1(s)\f$
 *  - \f$ K_{n}^{\prime}(x) = -\frac{1}{2}\left(K_{n-1}(x) + K_{n+1}(x)\right)\f$ 
 *  @see Ostap::Math::bessel_Knu
 */
// ============================================================================
double Ostap::Math::der_bessel_Knu
( const double nu , 
  const double x  )
{
  if ( isint ( nu ) ) 
  {
    const int n = Ostap::Math::round  ( nu ) ;
    return Ostap::Math::der_bessel_Kn ( n , x ) ;
  }
  return
    s_zero ( nu ) ? - Ostap::Math::bessel_Kn ( 1 , x ) :
    - 0.5 * ( Ostap::Math::bessel_Knu ( nu - 1 , x ) +
              Ostap::Math::bessel_Knu ( nu + 1 , x ) ) ; 
}






// ============================================================================
/*  Laguerre polynomila of non-integher order 
 *  \f$ L_{q}(x) = {}_1F_1(-q; 1; x ) \f$, where 
 *  \f$ {}_1F_1(-1; 1; x ) \f$ is a confluent hypergeometrical function 
 */
// ============================================================================
double Ostap::Math::laguerre_q ( const double q , const double x ) 
{
  gsl_sf_result result ;
  const int ierror = gsl_sf_hyperg_1F1_e ( -q , 1.0 , x , &result ) ;
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_hyperg_1F1_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}

// ============================================================================
/* Helpful function \f$ H_a(a,u_1,u_2)\f$ for the relativistic Voigt profile
 * 
 * The relativistic Voigt profile \f$ V_2(m;\mu,\Gamma,\sigma) \f$  is
 *  \f$ V_2(m; \mu,\Gamma,\sigma) \equiv  S_2(m;\mu,\Gamma)\ast G(\delta m;\sigma)\f$ 
 *  where 
 *  - \f$ S_2 = \frac{1}{\pi}\frac{\mu\Gamma}{ (m^2-\mu^2)^2 + \mu^2\Gamma^2 } \f$    
 *  - \f$ G(\delta m ; \sigma) = \frac{1}{\sqrt{2\pi\sigma^2}} 
 *     \mathrm{e}^{-\frac{1}{2} \left( \frac{\delta m }{\sigma} \right)^2} \f$$
 *  
 *  \f$ V_2(m; \mu,\Gamma,\sigma = \frac{H_2(a,u_1,u_2)}{2\sqrt{\pi}\sigma^2}\f$, where 
 *  - \f$ u_1 = \frac{m-\mu }{\sqrt{2}\sigma} \f$
 *  - \f$ u_2 = \frac{m+\mu }{\sqrt{2}\sigma} \f$
 *  - \f$ a   = \frac{\mu\Gamma}{2\sigma^2} \f$
 *
 *  \f[ H_2(a,u_1,u_2) = 
 *   \frac{a}{\pi} \int_{-\infty}{+\infty}  
 *    \frac{dt}{  (u_1-t)^2(u_2-t0^2 + a^2 } \f] 
 * @see Kycia, Radoslaw A.; Jadach, Stanislaw. 
 *      "Relativistic Voigt profile for unstable particles in high energy physics". 
 *      Journal of Mathematical Analysis and Applications. 463 (2): 1040–1051 
 *      arXiv:1711.09304. doi:10.1016/j.jmaa.2018.03.065. ISSN 0022-247X.
 * @see https://arxiv.org/abs/1711.09304
 */
// ============================================================================
double Ostap::Math::H2 ( const double a  , 
                         const double u1 , 
                         const double u2 ) 
{
  if ( a  < 0 ) { return H2 ( std::abs ( a ) , u1 , u2 ) ; }
  return 0 ;
}
// ============================================================================
/* \f$ \left| \frac{\Gamma(x+iy)}{\Gamma(x)} \right|^2 \f$ for 
 *  \f$  x> 0\f$.
 *  
 *  \f$ \left| \right|^2 = \left| \frac{1}{F(-iy,iy,x,1)}\right| \f$, where 
 *  \f$ F(a,b,c;z)\f$ is a hypergeometrical function.
 *
 *  This expression appears in normalization for  the Pearson Type IV pdf 
 *  @see J. Heinrich, "A guide to the Pearson Type IV distribution", 
 *       CDF/MEMO/STATISTICS/PUBLIC/6820, 2004 
 *  @see http://www-cdf.fnal.gov/physics/statistics/notes/cdf6820_pearson4.pdf
 */
// ============================================================================
double Ostap::Math::pearsonIV_g2 ( const double x , const double y ) 
{
  if ( s_zero ( y ) ) { return 1 ; }
  //
  if ( x <= 0 || s_zero ( x ) ) 
  {
    // use GSL error reporting/handling system 
    Ostap::Math::GSL::GSL_Error_Handler sentry ;
    gsl_error ( "Invalid x for Ostap::Math::pearsonIV" , __FILE__ , __LINE__ , GSL_EDOM ) ;
    return std::numeric_limits<double>::quiet_NaN() ;
  }
  //
  const double y2   = y * y ;
  const double xmin = std::max ( 2 * y2 , 10.0 ) ;
  //
  double xx = x ;
  double r  = 1 ;
  while ( xx < xmin ) 
  {
    const double t = y / xx++ ;
    r *= 1 + t * t ;
  }
  double p = 1 ;
  double s = 1 ;
  double f = 0 ;
  while ( p > s * DBL_EPSILON ) 
  {
    p *= y2   + f * f ;
    p /= xx++ *   ++f ;
    s += p ;
  }
  return 1 / ( r * s ) ;
}
// ============================================================================

// ============================================================================
/*  simple infinitely smooth and finite function 
 *  \f$ f(x) = \mathrm{e}^{ - \frac{1}{1-x^2}}\f$ 
 *   for \f$ \left| x \right| \< 1\f$. else 0.
 */
// ============================================================================
double Ostap::Math::hat ( const double x ) 
{ return 1 <= std::abs ( x )  ?  0.0 : std::exp ( - 1.0 / ( 1.0 - x * x ) ) ; }
// ============================================================================
/*  Sinc function 
 *  \f$ f(x) = \frac{ \sin x }{x}  \f$ 
 *  @see https://en.wikipedia.org/wiki/Sinc_function
 */
// ============================================================================
double Ostap::Math::sinc ( const double x ) 
{
  const double absx = std::abs ( x ) ;
  // 
  if ( 1.e-3 < absx ) { return std::sin ( absx ) / absx ; } 
  //
  const long double x2     = absx * absx ;
  //
  long double       term   = - x2 / 6 ;
  long double       result = 1 ;
  unsigned long     k      = 3 ;
  //
  while ( std::abs ( term ) > 0.5 * s_EPSILON  ) 
  {
    result += term ;
    term   *= - x2 / ( ( k + 1 ) * ( k + 2 ) ) ;
    k      += 2 ;
  }
  //
  return result ;
}

// ========================================================================
/*  \f$ f(x) = \frac{ \sin x }{x}  \f$ 
 *  @see https://en.wikipedia.org/wiki/Sinc_function
 */
// ========================================================================
double Ostap::Math::sin_x  ( const double x ) 
{ return Ostap::Math::sinc ( x ) ; }
// ========================================================================
/*  \f$ f(x) = \frac{ \sinh x }{x}  \f$ 
 */
// ========================================================================
double Ostap::Math::sinh_x ( const double x ) 
{
  const double absx = std::abs ( x ) ;
  if ( 1.e-3 < absx ) { return std::sinh ( absx ) / absx ; } 
  //
  const long double x2     = absx * absx ;
  //
  long double       term   = x2 / 6 ;
  long double       result = 1 ;
  unsigned long     k      = 3 ;
  //
  while ( std::abs ( term ) > 0.5 * s_EPSILON  ) 
  {
    result += term ;
    term   *= x2 / ( ( k + 1 ) * ( k + 2 ) ) ;
    k      += 2 ;
  }
  //
  return result ;
}

// ============================================================================
/*  \f$ f(x) = \frac{ \asinh x }{x}  \f$ 
 */
// ========================================================================
double Ostap::Math::asinh_x ( const double x ) 
{
  const double absx = std::abs ( x ) ;
  if ( 1.e-3 < absx ) { return std::asinh ( absx ) / absx ; } 
  //
  const long double x2     = absx * absx ;
  //
  long double       term   = - x2 / 2 ;
  long double       result =    1 ;
  unsigned long     k      =    3 ;
  //
  while ( std::abs ( term ) > 0.5 * s_EPSILON  ) 
  {
    result += term / k             ;
    term   *= - x2 * k / ( k + 1 ) ;
    k      += 2 ;
  }
  //
  return result ;
}
// ============================================================================
/*  \f$ f(x) = \frac{ \log ( 1 + x ) }{x}  \f$ 
 *  precise for small x 
 */
// ============================================================================
double Ostap::Math::log1p_x ( const double x ) 
{
  //
  const double absx = std::abs ( x ) ;
  if ( 1.e-3 < absx ) { return std::log1p ( x ) / x ; } 
  //
  long double       term   = -x ;
  long double       result =  1 ;
  unsigned long     k      =  1 ;
  //
  while ( k * std::abs ( term ) > 0.5 * s_EPSILON  ) 
  {
    result += term / ( k + 1 ) ;
    term   *= -x ;
    k      +=  1 ;
  }
  //
  return result ;
}
// ============================================================================
/*  \f$ f(x) = \frac{ \mathrm{e}^{x} -1 }{x}  \f$ 
 *  precise for small x 
 */
// ============================================================================
double Ostap::Math::expm1_x ( const double x ) 
{
  const double absx = std::abs ( x ) ;
  if ( 1.e-3 < absx ) { return std::expm1 ( x ) / x  ; } 
  //
  long double       term   = x / 2 ;
  long double       result = 1 ;
  unsigned long     k      = 2 ;
  //
  while ( std::abs ( term ) > 0.5 * s_EPSILON  ) 
  {
    result += term ;
    term   *= x  / ( k + 1 ) ;
    k      += 1 ;
  }
  //
  return result ;
}
// ============================================================================

// ============================================================================
/** Fourrier-image of the finite atomic function <code>up</code> 
 *  \f$ \hat{up}(p) = \Pi_{k=1}^{\infty} \frac{ \sin p 2^{-k}}{ p 2^{-k}}\f$ 
 */
// ============================================================================
double Ostap::Math::up_F ( const double p ) 
{
  long double   result = 1     ;
  long double   arg    = p / 2 ;
  for  (; ;)
  {
    const double term = Ostap::Math::sinc ( arg ) ;  
    if      ( s_zero ( term   ) ) { return 0.0 ; }
    else if ( s_zero ( result ) ) { return 0.0 ; }
    else if ( std::abs ( term - 1.0L ) <= 1.2 * s_EPSILON ) 
    { return result * term ;}
    //
    result *= term ;
    arg    /= 2 ;
  }
  return result ;
}
// ========================================================================
/** Fourrier-image of the finite atomic function <code>fup_N</code> 
 *  \f$ \hat{fup_N}(p) = 
 *   \left(  \frac{\sin x /2 }{ x/2 }  \right )^N 
 * \Pi_{k=1}^{\infty} \frac{ \sin p 2^{-k}}{ p 2^{-k}}\f$ 
 */
// ========================================================================
double Ostap::Math::fupN_F 
( const unsigned short N , 
  const double         p ) 
{
  if ( 0 == N ) { return up_F ( p ) ; }
  long double   result = std::pow ( Ostap::Math::sinc ( p / 2 ) , N ) ;
  long double   arg    = p / 2 ;
  for  (; ;)
  {
    const double term = Ostap::Math::sinc ( arg ) ;  
    if      ( s_zero ( term   ) ) { return 0.0 ; }
    else if ( s_zero ( result ) ) { return 0.0 ; }
    else if ( std::abs ( term - 1.0L ) <= 1.2 * s_EPSILON ) 
    { return result * term ;}
    //
    result *= term ;
    arg    /= 2 ;
  }
  return result ;
}
// ============================================================================


// ============================================================================
/*  get Arithmetic-geometric mean 
 *  @see https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean
 *  @param x  x-value, \f$ x \ge 0 \f$   
 *  @param y  y-value, \f$ y \ge 0 \f$    
 *  @return Arithmetic-geometric mean 
 */
// ============================================================================
double Ostap::Math::agm
( const double x ,
  const double y )
{
  //
  if ( !std::isfinite ( x ) || x < 0 ) 
  { return std::numeric_limits<double>::quiet_NaN() ; }
  //
  if ( !std::isfinite ( y ) || y < 0 ) 
  { return std::numeric_limits<double>::quiet_NaN() ; }
  //
  if ( 0 == x || 0 == y ) { return 0 ; }
  //
  long double a = x ;
  long double b = y ;
  //
  for ( unsigned short i = 0 ; i < 50 ; ++i )
  {
    long double aa = 0.5L *    ( a + b ) ;
    long double bb = std::sqrt ( a * b ) ;
    //
    if ( 3 < i && s_equal ( aa , bb ) ) { return 0.5L * ( aa + bb ) ; }
    a = aa ;
    b = bb ;
  }
  return 0.5L * ( a + b ) ;  
}
// ============================================================================
/*  get Geometric-harmonic  mean         
 *  @see https://en.wikipedia.org/wiki/Geometric%E2%80%93harmonic_mean
 *  @param x  x-value, \f$ x > 0 \f$   
 *  @param y  y-value, \f$ y > 0 \f$    
 *  @return Geometric-harmonic  mean         
 */
// ============================================================================
double Ostap::Math::ghm 
( const double x ,
  const double y )
{
  //
  if ( !std::isfinite ( x ) || x <= 0 ) 
  { return std::numeric_limits<double>::quiet_NaN() ; }
  //
  if ( !std::isfinite ( y ) || y <= 0 ) 
  { return std::numeric_limits<double>::quiet_NaN() ; }
  //
  return x * y / agm ( x , y ) ;
}
// ============================================================================
/* simple geometric mean for two real numbers (just for completeness)
 *  @param a the first number 
 *  @param b the second number 
 *  @return geometric mean 
 */
// ============================================================================
double
Ostap::Math::geometric_mean 
( const double a , 
  const double b ) 
{
  return 
    s_zero ( a ) ? 0.0 :
    s_zero ( b ) ? 0.0 : std::sqrt ( a * b ) ;
}
// ============================================================================
/*  simple geometric mean for two complex numbers 
 *  branch is chosen according to this:
 *  @see https://www.johndcook.com/blog/2022/07/30/complex-agm/#:~:text=Complex%20AGM,%E2%88%9A(an%20bn)
 *  @param a the first number 
 *  @param b the second number 
 *  @return geometric mean 
 */
// ============================================================================
std::complex<double>
Ostap::Math::geometric_mean 
( const std::complex<double>& a , 
  const std::complex<double>& b ) 
{
  //
  static const Ostap::Math::Zero     < std::complex<double> >  s_czero  {} ;
  static const Ostap::Math::Equal_To < std::complex<double> >  s_cequal {} ;
  //
  if      ( s_czero  ( a     ) ) { return 0.0 ; }
  else if ( s_czero  ( b     ) ) { return 0.0 ; }
  else if ( s_cequal ( a , b ) ) { return a   ; }
  //
  const double d1 = std::norm ( a + b ) ;
  const double d2 = std::norm ( a - b ) ;
  //
  const std::complex<double> r { std::sqrt ( a * b ) } ;
  //
  if      (           d2 < d1                              ) { return r ; }
  else if ( s_equal ( d2 , d2 ) && 0 < std::imag ( b / a ) ) { return r ; }
  //
  return -r ;
}
// ============================================================================
/** get Arithmetic-geometric mean for complex numbers 
 *  @see https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean
 *  @param x  x-value, \f$ x \ge 0 \f$   
 *  @param y  y-value, \f$ y \ge 0 \f$    
 *  @return Arithmetic-geometric mean 
 *  @see https://www.johndcook.com/blog/2022/07/30/complex-agm/#:~:text=Complex%20AGM,%E2%88%9A(an%20bn)
 */
// ============================================================================
std::complex<double> 
Ostap::Math::agm
( const std::complex<double>& x , 
  const std::complex<double>& y ) 
{
  //
  static const Ostap::Math::IsReal   < std::complex<double> > s_isreal {} ;
  static const Ostap::Math::Zero     < std::complex<double> > s_czero  {} ;
  static const Ostap::Math::Equal_To < std::complex<double> > s_cequal {} ;
  //
  if ( s_isreal ( x ) && s_isreal ( y ) ) { return agm ( x.real () , y.real () ) ; }
  //
  if      ( s_czero  ( x     ) ) { return 0.0 ; }
  else if ( s_czero  ( y     ) ) { return 0.0 ; }
  else if ( s_cequal ( x , y ) ) { return x   ; }
  //
  std::complex<double> a = x ;
  std::complex<double> b = y ;
  //
  for ( unsigned short i = 0 ; i < 100 ; ++i )
  {
    std::complex<double> aa = 0.5 *     ( a + b ) ;
    std::complex<double> bb = std::sqrt ( a * b ) ;
    //
    if ( s_czero  ( aa      ) ) { return 0.0 ; }
    if ( s_cequal ( aa , bb ) ) { return 0.5 * ( aa + bb ) ; }
    //
    const double d1 = std::norm ( a + b ) ;
    const double d2 = std::norm ( a - b ) ;
    //
    if      (           d2 < d1 )                 
    {
      a =  aa ;
      b =  bb ;
    }
    else if ( s_equal ( d1 , d2 ) && 0 < std::imag ( bb / aa ) ) 
    {
      a =  aa ;
      b =  bb ;
    }
    else 
    {
      a =  aa ;
      b = -bb ;
    }
    //
  }
  return 0.5 * ( a + b ) ; 
}
// ============================================================================

// ============================================================================
/*  simple harmonic mean for two real numbers (just for completeness)
 *  @param a the first number 
 *  @param b the second number 
 *  @return harmonic mean 
 */
// ============================================================================
double
Ostap::Math::harmonic_mean 
( const double a , 
  const double b ) 
{ return
    s_zero ( a ) ? 0.0 :
    s_zero ( b ) ? 0.0 :
    2.0 / ( 1. / a + 1. / b ) ; }
// ============================================================================
/*  simple harmonic mean for two real numbers (just for completeness)
 *  @param a the first number 
 *  @param b the second number 
 *  @return harmonic mean 
 */
// ============================================================================
std::complex<double>
Ostap::Math::harmonic_mean 
( const std::complex<double>& a , 
  const std::complex<double>& b ) 
{
  return
    s_zero ( a.real() ) && s_zero ( a.imag() ) ? 0.0 :
    s_zero ( b.real() ) && s_zero ( b.imag() ) ? 0.0 :
    2.0 / ( 1. / a + 1. / b ) ;
}
// ============================================================================
/*  Heronian mean for two real numbers (just for completeness)
 *  @see https://en.wikipedia.org/wiki/Heronian_mean 
 *  @param a the first number 
 *  @param b the second number 
 *  @return heronian mean 
 */
// ============================================================================
double
Ostap::Math::heronian_mean 
( const double a , 
  const double b )
{
  return ( a + std::sqrt ( a * b ) + b  ) / 3 ;
}
// ============================================================================
/** Power mean for two real numbers 
 *  @see https://en.wikipedia.org/wiki/Power_mean 
 *  @param x-parameter: 0<=x<=0.5
 *  @param a the first number 
 *  @param b the second number 
 *  @return power mean 
 */
// ============================================================================
double
Ostap::Math::power_mean 
( const double p , 
  const double a ,
  const double b )
{
  return
    s_zero ( p ) ? geometric_mean ( a , b ) :
    std::pow ( 0.5 * ( std::pow ( a , p ) + std::pow ( b , p ) ) , 1.0 / p ) ;
}
// ============================================================================
/*  Heinz mean for two real numbers (just for completeness)
 *  @see https://en.wikipedia.org/wiki/Heinz_mean 
 *  \f$ 0 \le x \l1 1/2 \f$ 
 *  @param x-parameter: 0<=x<=0.5
 *  @param a the first number 
 *  @param b the second number 
 *  @return Heinz mean 
 */
// ============================================================================
double
Ostap::Math::heinz_mean 
( const double x , 
  const double a ,
  const double b )
{
  return
    ( x   <= 0 ) || s_zero  ( x       ) ? 0.5 *     ( a + b ) : 
    ( 0.5 <= x ) || s_equal ( x , 0.5 ) ? std::sqrt ( a * b ) : 
    0.5 * ( std::pow ( a ,     x ) * std::pow ( b ,  1 - x ) +
	    std::pow ( a , 1 - x ) * std::pow ( b ,      x ) ) ;
}
// ============================================================================
/*  Lehmer mean for two real numbers (just for completeness)
 *  @see https://en.wikipedia.org/wiki/Lehmer_mean 
 *  - \f$ p \rigtharrow - \infty\f$ : minimal value
 *  - \f$ p = 0 \f$   : harmonic mean 
 *  - \f$ p = 1/2 \f$ : geometric mean 
 *  - \f$ p = 1 \f$   : ariphmetic mean 
 *  - \f$ p = 2 \f$   : contraharmonic mean 
 *  - \f$ p \rigtharrow + \infty\f$ : maximal value
 *  @param p-parameter
 *  @param a the first number 
 *  @param b the second number 
 *  @return Lehmer mean 
 */
// ============================================================================
double
Ostap::Math::lehmer_mean 
( const double p , 
  const double a ,
  const double b )
{
  return
    s_zero  ( p       ) ? harmonic_mean   ( a , b ) : 
    s_equal ( p , 0.5 ) ? geometric_mean  ( a , b ) :
    s_equal ( p , 1   ) ? arithmetic_mean ( a , b ) : 
    ( std::pow ( a , p     ) + std::pow ( b , p     ) ) /
    ( std::pow ( a , p - 1 ) + std::pow ( b , p - 1 ) ) ; 
}
// ============================================================================



// ============================================================================
/*  Gudermannian function 
 *  @see https://en.wikipedia.org/wiki/Gudermannian_function
 *  @param x argument 
 *  @return value of Gudermannian function 
 */
// ============================================================================
double Ostap::Math::gd     ( const double x ) 
{ return 2 * std::atan  ( std::tanh ( 0.5L * x ) ) ; }
// ============================================================================
/*  inverse Gudermannian function 
 *  @see https://en.wikipedia.org/wiki/Gudermannian_function
 *  @param x argument 
 *  @return value of inverse Gudermannian function 
 */
// ============================================================================
double Ostap::Math::gd_inv ( const double x ) 
{ return
    std::abs ( x ) < 0.5 * M_PI ? 
                     2 * std::atanh ( std::tan ( 0.5L * x ) ) : 
    std::numeric_limits<double>::quiet_NaN() ;
}                     
// ============================================================================


// ============================================================================
/*  Airy function Ai 
 *  @see https://en.wikipedia.org/wiki/Airy_function
 *  @parameter x argument 
 *  @return value of Airy fuinction Ai 
 */
// ============================================================================
double Ostap::Math::Ai ( const double x ) 
{ 
  gsl_sf_result result ;
  const int ierror = gsl_sf_airy_Ai_e ( x , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_airy_Ai_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================
/*  Airy function Bi 
 *  @see https://en.wikipedia.org/wiki/Airy_function
 *  @parameter x argument 
 *  @return value of Airy fuinction Bi 
 */
// ============================================================================
double Ostap::Math::Bi ( const double x ) 
{
  gsl_sf_result result ;
  const int ierror = gsl_sf_airy_Bi_e ( x , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_airy_Bi_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================
/*  derivative Airy function Ai 
 *  @see https://en.wikipedia.org/wiki/Airy_function
 *  @parameter x argument 
 *  @return value of derivative of Airy fuinction Ai 
 */
// ============================================================================
double Ostap::Math::der_Ai ( const double x ) 
{ 
  gsl_sf_result result ;
  const int ierror = gsl_sf_airy_Ai_deriv_e ( x , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_airy_Ai_deriv_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================
/*  derivative Airy function Bi 
 *  @see https://en.wikipedia.org/wiki/Airy_function
 *  @parameter x argument 
 *  @return value of derivative Airy fuinction Bi 
 */
// ============================================================================
double Ostap::Math::der_Bi ( const double x ) 
{
  gsl_sf_result result ;
  const int ierror = gsl_sf_airy_Bi_deriv_e ( x , GSL_PREC_DOUBLE , &result ) ;
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_airy_Bi_deriv_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ============================================================================


// ============================================================================
/* Ramanudjan' sum
 *  \f$ R(x) = \sum_{n=1}^{+\infty} \frac{ (-1)^{n-1}x^n}{n!2^{n-1}
 *      \sum_{k=0}{ floor (\frac{n-1}{2})} \frac{1}{2k+1}\f$ 
 *  Helpe sum function to calcaltuon integral logarithm or integral exponent 
 */
// ============================================================================
double Ostap::Math::ramanudjan_sum ( const double x )
{
  double result =  0 ;
  double term   = -2 ;
  double sum    =  0 ;
  //
  static const unsigned short s_USMX = std::numeric_limits<unsigned short>::max () ;
  //
  const double        ax = std::abs ( x ) ;
  const unsigned int  N  =
    ( ax < s_USMX ) ? std::max ( 100u , (unsigned int) ax ) : s_USMX ;
  //
  for ( unsigned short n = 1 ; n < N ;  ++n  )
    {
      if ( 1 == n % 2 ) { sum  += 1. / n ; } 
      term *= - x / ( n * 2 ) ;
      const double delta = term  * sum  ;
      if ( s_zero ( delta ) || s_equal ( result , result + delta ) )
        { return result + delta ; }                              // RETURN 
      result += delta ;
    }
  //
  return result ;
}
// ============================================================================
/* get the Logarithmic integral function 
 *   \f$ li(x) \equiv \int\limits_{0}^{x} \frac{dt}{\log t} \f$ 
 *   for \f$ 0 < x \f$ 
 *  @see https://en.wikipedia.org/wiki/Logarithmic_integral_function
 */
// ============================================================================
double Ostap::Math::li ( const double x )
{
  //
  if      ( 0 == x || s_zero ( x ) ) { return 0 ; }
  else if ( x <  0 )                 { return std::numeric_limits<double>::quiet_NaN() ; }
  //
  const double lnx = std::log  ( x ) ;
  const double sqx = std::sqrt ( x ) ;
  //
  double result = s_EULER + std::log ( std::abs ( lnx ) ) ;
  //
  return result + sqx * ramanudjan_sum ( lnx ) ; 
}
// ============================================================================
/*  get the Logarithmic integral function 
 *   \f$ Li(x) \equiv \int\limits_{2}^{x} \frac{dt}{\log t} \f$ 
 *   for \f$ 0 < x \f$ 
 *  - \f$ Li(x) = li(x) - li(2) \f$ 
 *  @see https://en.wikipedia.org/wiki/Logarithmic_integral_function
 */
// ============================================================================
double Ostap::Math::Li  ( const double x )
{
  static const double s_li2 = li(2.0) ;
  return li ( x ) - s_li2 ;
}
// ============================================================================
/* Get the Exponential  integral function 
 * \f$ Ei(x) \equiv \int\limits_{-\infty}^{x} \frac{e^t}{t}dt \f$ 
 *  @see https://en.wikipedia.org/wiki/Exponential_integral
 */
// ============================================================================
double Ostap::Math::Ei  ( const double x )
{
  return s_EULER + std::log ( x ) + std::exp ( 0.5 * x ) * ramanudjan_sum ( x ) ;
  // return s_EULER + std::log ( x ) - Ein ( -x ) ;
}
// ============================================================================
/** get the Exponential  integral function  Ein(x) 
 *   \f$ Ein(x) \equiv \int\limits_{0}^{x} \frac{1-e^t}{t}dt \f$ 
 *  @see https://en.wikipedia.org/wiki/Exponential_integral
 */
// ============================================================================
double Ostap::Math::Ein  ( const double x )
{ return - std::exp ( -0.5 * x ) * ramanudjan_sum ( -x ) ; }
// ============================================================================
/*  get the Exponential  integral function \f$ E_1(x)\f$ 
 *  \f$ E_1(x) = -\gamma - \log x + Ein(x) \f$ 
 *  @see https://en.wikipedia.org/wiki/Exponential_integral
 */
// =============================================================================
double Ostap::Math::E1  ( const double x )
{ return -s_EULER - std::log ( x ) + Ein( x ) ; }
// =============================================================================


namespace
{
  // ===========================================================================
  double pade_f ( const double x )
  {
    static const std::array<long double,11> A = {
      +1.0L                       ,
      +7.44437068161936700618e+2  , 
      +1.96396372895146869801e+5  ,
      +2.37750310125431834034e+7  ,
      +1.43073403821274636888e+9  ,
      +4.33736238870432522765e+10 ,
      +6.40533830574022022911e+11 ,
      +4.20968180571076940208e+12 ,
      +1.00795182980368574617e+13 ,
      +4.94816688199951963482e+12 ,
      -4.94701168645415959931e+11 } ;
    static const std::array<long double,10> B = {
      +1.0L                       ,
      +7.46437068161927678031e+2  ,
      +1.97865247031583951450e+5  ,
      +2.41535670165126845144e+7  ,
      +1.47478952192985464958e+9  , 
      +4.58595115847765779830e+10 ,
      +7.08501308149515401563e+11 , 
      +5.06084464593475076774e+12 ,
      +1.43468549171581016479e+13 ,
      +1.11535493509914254097e+13 } ;
    //
    const long double x2 = 1.0 / ( x * x ) ;
    return
      Ostap::Math::Clenshaw::monomial_sum ( A.rbegin() , A.rend() , x2 ).first /
      Ostap::Math::Clenshaw::monomial_sum ( B.rbegin() , B.rend() , x2 ).first / x ;
  }
  // ================================================================================
  double pade_g ( const double x )
  {
    //
    static const std::array<long double,11> A = {
      +1.0L                    ,
      +8.1359520115168615e+2   ,
      +2.35239181626478200e+5  , 
      +3.12557570795778731e+7  ,
      +2.06297595146763354e+9  , 
      +6.83052205423625007e+10 ,
      +1.09049528450362786e+12 , 
      +7.57664583257834349e+12 ,
      +1.81004487464664575e+13 , 
      +6.43291613143049485e+12 ,
      -1.36517137670871689e+12 } ;
    //
    static const std::array<long double,10> B = { {
        +1.0L ,
        +8.19595201151451564e+2  ,
        +2.40036752835578777e+5  , 
        +3.26026661647090822e+7  ,
        +2.23355543278099360e+9  ,
        +7.87465017341829930e+10 ,
        +1.39866710696414565e+12 , 
        +1.17164723371736605e+13 ,
        +4.01839087307656620e+13 , 
        +3.99653257887490811e+13 } } ;
    //
    const long double x2 = 1.0 / ( x * x ) ;   
    return
      Ostap::Math::Clenshaw::monomial_sum ( A.rbegin() , A.rend() , x2 ).first /
      Ostap::Math::Clenshaw::monomial_sum ( B.rbegin() , B.rend() , x2 ).first * x2 ;
  }
  // ===========================================================================
}
// =============================================================================
/*  get the integral sine 
 *  \f$ Si(x) = \int\limits_{0}^{x} \frac{\sin t}{t}dt \f$ 
 *  @see https://en.wikipedia.org/wiki/Trigonometric_integral
 */
// =============================================================================
double Ostap::Math::Si   ( const double x )
{
  //
  static const std::array<long double,8> A = {
    +1.0L                    , 
    -4.54393409816329991e-2  ,
    +1.15457225751016682e-3  ,
    -1.41018536821330254e-5  ,
    +9.43280809438713025e-8  ,
    -3.53201978997168357e-10 ,
    +7.08240282274875911e-13 ,
    -6.05338212010422477e-16 } ;
  //
  static const std::array<long double,7> B = {
    +1.0L ,
    +1.01162145739225565e-2  ,
    +4.99175116169755106e-5  ,
    +1.55654986308745614e-7  ,
    +3.28067571055789734e-10 , 
    +4.5049097575386581e-13  ,
    +3.21107051193712168e-16 } ;
  //
  if ( std::abs ( x ) <= 4 )
    {
      const long double x2 = x * x ;
      return
        Ostap::Math::Clenshaw::monomial_sum ( A.rbegin() , A.rend() , x2 ).first /
        Ostap::Math::Clenshaw::monomial_sum ( B.rbegin() , B.rend() , x2 ).first * x ;
    }
  //
  const double ax = std::abs ( x ) ;
  const double result = M_PI * 0.5
    - pade_f ( ax ) * std::cos ( ax )
    - pade_g ( ax ) * std::sin ( ax ) ;
  ///
  return 0 < x ? result : -result ;
}
// =============================================================================
/* get the integral cosine  
 *  \f$ Cin(x) = \int\limits_{0}^{x} \frac{1 - \cos t}{t}dt \f$ 
 *  @see https://en.wikipedia.org/wiki/Trigonometric_integral
 */
// =============================================================================
double Ostap::Math::Cin  ( const double x )
{
  //
  static const std::array<long double,7> A = {
    -0.25L ,
    +7.51851524438898291e-3  , 
    -1.27528342240267686e-4  ,
    +1.05297363846239184e-6  ,
    -4.68889508144848019e-9  ,
    +1.06480802891189243e-11 , 
    -9.93728488857585407e-15 } ;
  static const std::array<long double,8> B = {
    +1.0L , 
    +1.1592605689110735e-2   ,
    +6.72126800814254432e-5  , 
    +2.55533277086129636e-7  ,
    +6.97071295760958946e-10 , 
    +1.38536352772778619e-12 ,
    +1.89106054713059759e-15 , 
    +1.39759616731376855e-18 } ;
  //
  if ( std::abs ( x ) <= 4 )
    {
      const double x2 = x * x ;
      return -1 * 
        Ostap::Math::Clenshaw::monomial_sum ( A.rbegin() , A.rend() , x2 ).first /
        Ostap::Math::Clenshaw::monomial_sum ( B.rbegin() , B.rend() , x2 ).first * x2 ;
    }
  //
  const double ax     = std::abs ( x ) ;
  const double result = pade_f ( ax ) * std::sin ( ax ) - pade_g ( ax ) *  std::cos ( ax ) ; 
  //
  return s_EULER + std::log ( ax ) - result ; 
}
// =============================================================================
/** get Clausen function \f$ Cl_2 \f$ using GSL 
 *  @param x argument 
 *  @return value of clausen fnction \f$ Cl_2 \f$
 */
// =============================================================================
double Ostap::Math::clausen ( const double x )
{
  // ===========================================================================
  gsl_sf_result result ;
  const int ierror = gsl_sf_clausen_e ( x , &result ) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_clausen_e" , __FILE__ , __LINE__ , ierror ) ;
      if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN(); }
    }
  
  return result.val ;      
}
// ======================================================================

// =============================================================================
/* get the integral cosine  
 *  \f$ Ci(x) = - \int\limits_{x}^{+\infrty} \frac{\cos t}{t}dt \f$ 
 *  for \f$ 0 < x \f$ 
 *  @see https://en.wikipedia.org/wiki/Trigonometric_integral
 */
// =============================================================================
double Ostap::Math::Ci   ( const double x ) 
{
  if ( x <= 0 || s_zero ( x ) ) { return std::numeric_limits<double>::quiet_NaN(); }
  return s_EULER + std::log ( x ) - Cin ( x ) ;
}
// ============================================================================
/* get the Lambert W_0 function for \f$- \frac{1}{e} < x \f$ 
 *  @see https://en.wikipedia.org/wiki/Lambert_W_function
 */
// ============================================================================
double Ostap::Math::lambert_W0 ( const double x )
{
  //
  static const double s_mnmn = -1 / M_E ;
  if (x <= s_mnmn ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_lambert_W0_e ( x , &result ) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_lambert_W0_e" , __FILE__ , __LINE__ , ierror ) ;
      if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN(); }
    }
  return result.val ;
}
// ============================================================================
/*  get the Lambert W_0 function for \f$- \frac{1}{e} < x < 0  \f$ 
 *  @see https://en.wikipedia.org/wiki/Lambert_W_function
 */
// ============================================================================
double Ostap::Math::lambert_Wm1 ( const double x )
{
  //
  static const double s_mnmn = -1 / M_E ;
  if ( x <= s_mnmn | 0 <= x ) { return std::numeric_limits<double>::quiet_NaN(); }
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_lambert_Wm1_e ( x , &result ) ;
  if ( ierror ) 
    {
      gsl_error ( "Error from gsl_sf_lambert_Wm1_e" , __FILE__ , __LINE__ , ierror ) ;
      if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
        { return std::numeric_limits<double>::quiet_NaN(); }
    }
  return result.val ;
}


// ========================================================================
namespace
{
  // ======================================================================
  const Ostap::Math::Interpolation::Table s_BRING_table {
    {{ 0.302430000000, 
       0.317052447762, 
       0.331796375994, 
       0.346678757435, 
       0.361718076876, 
       0.376934395700, 
       0.392349416423, 
       0.407986547232, 
       0.423870966521, 
       0.440029687433, 
       0.456491622400, 
       0.473287647677, 
       0.490450667886, 
       0.508015680550, 
       0.526019840639, 
       0.544502525100, 
       0.563505397403, 
       0.583072472075, 
       0.603250179244, 
       0.624087429173, 
       0.645635676800, 
       0.667948986280, 
       0.691084095521, 
       0.715100480722, 
       0.740060420914, 
       0.766029062500, 
       0.793074483790, 
       0.821267759542, 
       0.850683025503, 
       0.881397542944, 
       0.913491763200, 
       0.947049392211, 
       0.982157455060, 
       1.018906360509, 
       1.057389965541, 
       1.097705639900, 
       1.139954330625, 
       1.184240626594, 
       1.230672823058, 
       1.279362986187, 
       1.330427017600, 
       1.383984718911, 
       1.440159856263, 
       1.499080224872, 
       1.560877713561, 
       1.625688369300, 
       1.693652461748, 
       1.764914547789, 
       1.839623536070, 
       1.917932751542, 
       2.000000000000 }} ,
    {{ -0.30000000, 
       -0.31400000, 
       -0.32800000, 
       -0.34200000, 
       -0.35600000, 
       -0.37000000, 
       -0.38400000, 
       -0.39800000, 
       -0.41200000, 
       -0.42600000, 
       -0.44000000, 
       -0.45400000, 
       -0.46800000, 
       -0.48200000, 
       -0.49600000, 
       -0.51000000, 
       -0.52400000, 
       -0.53800000, 
       -0.55200000, 
       -0.56600000, 
       -0.58000000, 
       -0.59400000, 
       -0.60800000, 
       -0.62200000, 
       -0.63600000, 
       -0.65000000, 
       -0.66400000, 
       -0.67800000, 
       -0.69200000, 
       -0.70600000, 
       -0.72000000, 
       -0.73400000, 
       -0.74800000, 
       -0.76200000, 
       -0.77600000, 
       -0.79000000, 
       -0.80400000, 
       -0.81800000, 
       -0.83200000, 
       -0.84600000, 
       -0.86000000, 
       -0.87400000, 
       -0.88800000, 
       -0.90200000, 
       -0.91600000, 
       -0.93000000, 
       -0.94400000, 
       -0.95800000, 
       -0.97200000, 
       -0.98600000, 
       -1.00000000 }} , true } ;
  // 
  const Ostap::Math::FloaterHormann s_BRING { s_BRING_table , 4 } ;
  // ======================================================================
  static const std::array<long double,5> s_BR1 { -1.0L , +1.0L , -5.0L , +35.0L , -285.0L } ; 
  // ======================================================================
}
// ========================================================================
/*  Bring radical or ultra-radical, the real solution of
 *  the equation \f$ x^5+x+a=0 \f$ 
 *  @see https://en.wikipedia.org/wiki/Bring_radical
 */
// ========================================================================
double Ostap::Math::bring   ( const double x )
{
  if      ( s_zero ( x ) ) { return 0                         ; }
  else if ( x < 0        ) { return -bring ( std::abs ( x ) ) ; }
  //
  const long double xx = x       ;
  const long double x2 = xx * xx ;
  const long double x4 = x2 * x2 ;
  //
  const long double x1 = - std::pow ( xx , 0.2 ) ;
  const long double ax = std::abs ( xx ) ;
  long double x0 =
    ax < 0.5 ?
    xx * Ostap::Math::Clenshaw::monomial_sum ( s_BR1.rbegin() , s_BR1.rend() , x4 ).first :
    ( s_BRING.xmin() <= ax && ax <= s_BRING.xmax() ? s_BRING ( x ) : 
      x1
      - 0.20       / std::pow ( x1 ,  3 )
      - 0.08       / std::pow ( x1 ,  7 ) 
      + 4.0/(25*3) / std::pow ( x1 , 11 ) ) ;
  
  if ( 5 < std::abs ( x0 ) ) 
    {
      // few fixed point iterations 
      for ( unsigned short i = 0 ; i < 4 ; ++ i )
	{ x0 = - std::pow ( xx + x0 , 0.2 ) ; }
    }
  /// Halley's method 
  auto fun  = [ xx ] ( const long double z ) -> long double
  { return      std::pow ( z , 5 ) + z + xx ; } ;
  auto der1 = [ xx ] ( const long double z ) -> long double
  { return  5 * std::pow ( z , 4 ) + 1    ; } ;
  auto der2 = [ xx ] ( const long double z ) -> long double
  { return 20  *std::pow ( z , 3 )        ; } ;
  
  long double r = x0 ;
  for ( unsigned short n ; n < 100 ; ++n )
    {
      const long double fn = fun  ( r ) ;
      if ( s_zero ( fn ) ) { return r ; }      // RETURN
      
      const long double d1 = der1 ( r ) ;
      const long double d2 = der2 ( r ) ;
      
      const long double fd = fn/d1 ;
      const long double dr = - fd * ( 1 - 0.5 * fd * d2 / d1 ) ;
      //
      if ( s_zero ( dr ) || s_equal ( r , r + dr ) ) { return r + dr ; }
      r += dr ;
    }
  //
  return r ;
}
// ============================================================================
/*  Bring radical or ultra-radical, the real solution of
 *  the equation \f$ x^5+x+a=0 \f$ 
 *  @see https://en.wikipedia.org/wiki/Bring_radical
 */
// ============================================================================
std::complex<double> Ostap::Math::bring   ( const std::complex<double>& x )
{
  const double im = std::imag ( x ) ;
  if ( s_zero ( im ) ) { return bring ( std::real ( x ) ) ; }
  const double re = std::real ( x ) ;
  static const std::complex<double> s_I { 0 , 1 } ;
  if ( s_zero ( re ) ) { return s_I * bring ( im ) ;  }
  //
  const double ax = std::abs ( x ) ;
  //
  const std::complex<double> xx = x       ;
  const std::complex<double> x2 = xx * xx ;
  const std::complex<double> x4 = x2 * x2 ;
  //
  const std::complex<double> x1 = - std::pow ( xx , 0.2 ) ;
  //
  std::complex<double> x0 =
    ax   < 0.5 ?
    xx * Ostap::Math::Clenshaw::monomial_sum ( s_BR1.rbegin() , s_BR1.rend() , x4 ).first :
    ( ax < 1.5 ? x1 :
      x1
      - 0.20       / std::pow ( x1 ,  3 )
      - 0.08       / std::pow ( x1 ,  7 ) 
      + 4.0/(25*3) / std::pow ( x1 , 11 ) ) ;  
  //
  if ( 5 < std::abs ( x0 ) ) 
    {
      // few fixed point iterations 
      for ( unsigned short i = 0 ; i < 4 ; ++ i )
	{ x0 = - std::pow ( xx + x0 , 0.2 ) ; }
    }
  //  
  /// Halley's method 
  auto fun  = [ xx ] ( const std::complex<double>& z ) -> std::complex<double> 
  { return       std::pow ( z , 5 ) + z + xx ; } ;
  auto der1 = [ xx ] ( const std::complex<double>& z ) -> std::complex<double> 
  { return  5. * std::pow ( z , 4 ) + 1.     ; } ;
  auto der2 = [ xx ] ( const std::complex<double>& z ) -> std::complex<double>
  { return 20.  *std::pow ( z , 3 )          ; } ;
  
  std::complex<double> r = x0 ;
  for ( unsigned short n ; n < 100 ; ++n )
    {
      const std::complex<double> fn = fun  ( r ) ;
      
      if ( s_czero ( fn ) ) { return r ; }      // RETURN
      
      const std::complex<double> d1 = der1 ( r ) ;
      const std::complex<double> d2 = der2 ( r ) ;
      
      const std::complex<double> fd = fn/d1 ;
      const std::complex<double> dr = - fd * ( 1.0 - 0.5 * fd * d2 / d1 ) ;

      if ( s_czero ( dr ) || s_cequal ( r , r + dr ) ) { return r + dr ; }
      r += dr ;
    }
  //
  return r ;
 
}
// ============================================================================




// ============================================================================
/*  complete Fermi-Dirac integral 
 *  \f$ F_j(x) = \frac{1}{\Gamma(j+1)}\int^{+\infty}_{0} 
 *   \frac{t^j}{ \exp ( t-x) + 1 } dt \f$ 
 *  @param j parameter
 *  @param x argument 
 *  @return return value of complete Fermi-Dirac integral 
 */
// ============================================================================
double Ostap::Math::fermi_dirac 
( const unsigned short j , 
  const double         x ) 
{
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_fermi_dirac_int_e ( j , x , &result ) ;
  if ( ierror ) 
  {
    gsl_error ( "Error from gsl_sf_fermi_diract_int_e" , __FILE__ , __LINE__ , ierror ) ;
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
  }
  return result.val ;
}
// ========================================================================

// ============================================================================
/* smoothstep (polynomial) function
 *  @see https://en.wikipedia.org/wiki/Smoothstep
 *  Transition function for \f$ 0 \le x \le 1\f$ 
 *  Transition function for \f$ 0 \le x \le 1\f$ 
 *  - \f$ f(x)=0\f$ for  \f$ x \le 0 \f$
 *  - \f$ f(x)=1\f$ for  \f$ x \ge 1 \f$
 *  - \f$ f(x)\f$ if a \f$ 2n+1 \f$ polynomial fuction inbetween      
 *  @param x variable
 *  @param n index, polynomial of order \f$ 2n+1 \f$     
 */
// ============================================================================
double Ostap::Math::smooth_step
( const double         x , 
  const unsigned short n ) 
{
  //
  static const std::array<int,4> s_3 { -20  ,   +70 , -84    , +35                             } ;
  static const std::array<int,5> s_4 { +70  ,  -315 , +540   , -420   , +126                   } ;
  static const std::array<int,6> s_5 { -252 , +1386 , -3080  , +3465  , -1980  , +462          } ;
  static const std::array<int,7> s_6 { +924 , -6006 , +16380 , -24024 , +20020 , -9009 , +1716 } ;
  //
  return 
    x <= 0 ? 0 : 
    x >= 1 ? 1 : 
    0 == n ? x : // simple clump function 
    1 == n ? x * x *     (            -  2 * x +  3 ) : 
    2 == n ? x * x * x * (  6 * x * x - 15 * x + 10 ) :
    3 == n ? std::pow ( x , 4 ) * Ostap::Math::Clenshaw::monomial_sum ( s_3.begin () , s_3.end() , x ).first : 
    4 == n ? std::pow ( x , 5 ) * Ostap::Math::Clenshaw::monomial_sum ( s_4.begin () , s_4.end() , x ).first : 
    5 == n ? std::pow ( x , 6 ) * Ostap::Math::Clenshaw::monomial_sum ( s_5.begin () , s_5.end() , x ).first : 
    6 == n ? std::pow ( x , 7 ) * Ostap::Math::Clenshaw::monomial_sum ( s_6.begin () , s_6.end() , x ).first : 
    // use n = 7 for other cases too.. 
    std::pow ( x , 7 ) * Ostap::Math::Clenshaw::monomial_sum ( s_6.begin () , s_6.end() , x ).first ; 
}

// ============================================================================
namespace 
{
  // ==========================================================================
  inline double  _smooth_psi_
  ( const double x )
  { return x <= 0 || s_zero ( x ) ? 0.0 : std::exp ( -1 / x ) ; }
  // ==========================================================================
  inline double  _smooth_phi_
  ( const double x ) 
  {
    const double psi1 = _smooth_psi_ (     x ) ;
    const double psi2 = _smooth_psi_ ( 1 - x ) ;
    return psi1 / ( psi1 + psi2 ) ;
  }
  // ==========================================================================
}  
// ========================================================================
/*  smooth transition function 
 *   \f[ \phix() = left\{ 
 *   \begin{array}{ll}
 *     0 &  x\le b \                            \
 *     1 &  x\ge b \\ 
 *    smooth & 
 *    \end{array} \rigth. \f] 
 */
// ========================================================================
double Ostap::Math::smooth_transtion 
( const double x , 
  const double a ,
  const double b ) 
{
  const double xmin = std::min ( a , b ) ;
  const double xmax = std::max ( a , b ) ;
  //
  return
    x <= xmin ? 0.0 :
    x >= xmax ? 1.0 :
    ::_smooth_phi_ ( ( x - xmin ) / ( xmax - xmin ) ) ;
}
// ============================================================================
/* get quantile function for standard normal distribution
 *  @see http://en.wikipedia.org./wiki/Probitq
 *  @param alpha argument    \f$  0<\alpha<1 \f$  
 *  @return quantile value 
 */
// ============================================================================
double Ostap::Math::probit
( const double alpha  )
{
  return
    s_zero  ( alpha     ) ? -std::numeric_limits<double>::max () :
    s_equal ( alpha , 1 ) ?  std::numeric_limits<double>::max () : 
    alpha < 0             ?  std::numeric_limits<double>::quiet_NaN () :
    alpha > 1             ?  std::numeric_limits<double>::quiet_NaN () :
    gsl_cdf_ugaussian_Pinv ( alpha ) ; 
}
// ============================================================================
//                                                                      The END 
// ============================================================================
