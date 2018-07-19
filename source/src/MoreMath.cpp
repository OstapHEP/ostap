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
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_log.h"
#include "gsl/gsl_sf_psi.h"
// ============================================================================
// LHCbMath
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Choose.h"
// ============================================================================
// Local
// ============================================================================
#include "GSL_sentry.h"
#include "Faddeeva.hh"
#include "gauss.h"
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
  const Ostap::Math::Equal_To<double> s_equal{} ;       // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>     s_zero {} ;       // zero for doubles
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
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
    if ( ierror == GSL_ERANGE   ||    // output range error, e.g. exp(1e100)
         ierror == GSL_EINVAL   ||    // invalid argument supplied by user
         ierror == GSL_EUNDRFLW ||    // underflow
         ierror == GSL_EOVRFLW  ||    // overflow
         ierror == GSL_ELOSS    ||    // loss of accuracy
         ierror == GSL_EROUND    )    // failed because of roundoff error
    {}
    else
    {
      gsl_error ( "Error from kummer function" ,
                  __FILE__ , __LINE__ , ierror ) ;
    } 
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
 *  @return the value of the coplmex error function 
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
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
    if ( ierror == GSL_ERANGE   ||    // output range error, e.g. exp(1e100)
         ierror == GSL_EINVAL   ||    // invalid argument supplied by user
         ierror == GSL_EUNDRFLW ||    // underflow
         ierror == GSL_EOVRFLW  ||    // overflow
         ierror == GSL_ELOSS    ||    // loss of accuracy
         ierror == GSL_EROUND    )    // failed because of roundoff error
    {}
    else
    {
      gsl_error ( "Error from gsl_sf_gammainv_e" ,
                  __FILE__ , __LINE__ , ierror ) ;
    } 
  }
  //
  return result.val ;
}
// ============================================================================
/*  compute psi function 
 *  \$f f(x) = \frac{d}{dx}\ln \Gamma(x)\f$
 *  @return the value of psi function 
 */
// ============================================================================G
double Ostap::Math::psi ( const double x ) 
{
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  gsl_sf_result result ;
  const int ierror = gsl_sf_psi_e ( x , &result ) ;
  if ( ierror ) 
  {
    //
    if      ( ierror == GSL_EDOM     ) // input domain error, e.g sqrt(-1)
    { return std::numeric_limits<double>::quiet_NaN(); }
    //
    if ( ierror == GSL_ERANGE   ||    // output range error, e.g. exp(1e100)
         ierror == GSL_EINVAL   ||    // invalid argument supplied by user
         ierror == GSL_EUNDRFLW ||    // underflow
         ierror == GSL_EOVRFLW  ||    // overflow
         ierror == GSL_ELOSS    ||    // loss of accuracy
         ierror == GSL_EROUND    )    // failed because of roundoff error
    {}
    else
    {
      gsl_error ( "Error from gsl_sf_psi_e" ,
                  __FILE__ , __LINE__ , ierror ) ;
    } 
  }
  //
  return result.val ;
}
// ============================================================================



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


// ============================================================================
/*  get the standard gaussian pdf 
 *  @see https://en.wikipedia.org/wiki/Normal_distribution
 *  @param x x-value  
 *  @param sigma sigma (width)  
 *  @param mu    mu (location)
 *  @return the value of gaussian pdf 
 */
// ============================================================================
double Ostap::Math::gauss_pdf ( const double x     ,
                                const double mu    ,
                                const double sigma )
{
  static const double s_norm = 1.0/std::sqrt( 2.0 * M_PI ) ;
  const double dx = ( x  - mu ) / std::abs ( sigma ) ;
  return s_norm * std::exp ( -0.5 * dx * dx ) / std::abs ( sigma ) ;
}
// ============================================================================
/*  get the standard gaussian cdf 
 *  @see https://en.wikipedia.org/wiki/Normal_distribution
 *  \f$ f(x) = \frac{1}{2} \left( 1 + erf ( \frac{x} { \sqrt{2} } ) \right) \f$ 
 *  @param x x-value  
 *  @return the value of gaussian cdf 
 */
// ============================================================================
double Ostap::Math::gauss_cdf ( const double x     ,
                                const double mu    ,
                                const double sigma )
{
  //
  static const double s_sqrt2 = std::sqrt( 2.0 ) ;
  const double y = ( x - mu ) / ( s_sqrt2 * std::abs ( sigma ) ) ;
  return 0.5 * ( 1 + std::erf ( y ) ) ;
}
// ============================================================================
#include <iostream>
namespace 
{
  // ==========================================================================
  typedef unsigned long long ULL ;
  // ==========================================================================
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
// The END 
// ============================================================================
