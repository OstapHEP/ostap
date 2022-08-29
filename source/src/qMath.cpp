// ============================================================================
// Include files
// ============================================================================ \
// STD&STL 
// ============================================================================
#include <cmath> 
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/qMath.h"
// ============================================================================
// Local
// ============================================================================
#include "local_math.h"
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for fiction from Ostap/qMath.h
 *  Collection of functions related qto Tsallis statistics 
 *  @see https://en.wikipedia.org/wiki/Tsallis_statistics
 *  @see Umarov, Sabir; Tsallis, Constantino; Steinberg, Stanly (2008). 
 *      "On a q-Central Limit Theorem Consistent with Nonextensive 
 *      Statistical Mechanics" Milan J. Math. Birkhauser Verlag. 76: 307–328. 
 *  @see doi:10.1007/s00032-008-0087-y. S2CID 55967725..
 *  @see https://doi.org/10.1007%2Fs00032-008-0087-y
 *  @date 2022-08-28 
 *  @author anya Belyaev Ivan/Belyaev@itep.ru
 */
// ============================================================================
/*  q-sum of two variables in Tsallis statistics 
 *  \f$ x + \oplus_q y = x + y + (1-q)xy\f$ 
 */
// ============================================================================
double Ostap::Math::tsallis_qsum        
( const double x ,
  const double y , 
  const double q )
{ return 1 == q || s_equal ( 1 , q ) ? x + y : x + y + ( 1.0L - q ) * x * y ; }
// ============================================================================
/*  q-subtraction of two variables in Tsallis statistics 
 *  \f$ x + \ominus_q y = \frac{x-y}{1+(1-q)y}\f$ 
 */
// ============================================================================
double Ostap::Math::tsallis_qsubtraction
( const double x ,
  const double y , 
  const double q )
{ return 1 == q || s_equal ( 1 , q ) ? x - y : ( x - y ) / ( 1 + ( 1.0L - q ) * y ) ; }
// ============================================================================
/* q-product of two variables in Tsallis statistics 
 *  \f$ x + \otimes_q y = \left{ x^{1-1} + y^{1-q} -1 \right]_+^{\frac{1}{1-q}}\f$ 
 */
// ============================================================================
double Ostap::Math::tsallis_qproduct 
( const double x ,
  const double y , 
  const double q )
{
  if ( 1 == q || s_equal ( q , 1 ) ) { return x * y ; } 
  const long double xl  = x ;
  const long double yl  = y ;
  const long double arg =  std::pow ( xl , 1.0L - q ) + std::pow ( yl , 1.0L - q ) - 1.0L ;
  return 0 <= arg ? std::pow ( arg , 1 / ( 1.0L - q ) ) : 0 ;
}
// =======================================================================
/*  q-division of two variables in Tsallis statistics 
 *  \f$ x + \oslash_q y = \left{ x^{1-1} - y^{1-q} +1 \right]_+^{\frac{1}{1-q}}\f$ 
 */
// =======================================================================
double Ostap::Math::tsallis_qdivision     
( const double x ,
  const double y , 
  const double q )
{
  if ( 1 == q || s_equal ( q , 1 ) ) { return x / y ; }
  const long double xl = x ;
  const long double yl = y ;
  const long double arg =  std::pow ( xl , 1.0L - q ) - std::pow ( yl , 1.0L - q ) + 1 ;
  return 0 <= arg ? std::pow ( arg , 1 / ( 1.0L - q ) ) : 0 ;
}
// ============================================================================
/* q-exponent in Tsallis statistics 
 *  \f$ e_q(x) = \left[1+(1-q)x\right]_+^{\frac{1}{1-q}\f$ 
 */
// ============================================================================
double Ostap::Math::tsallis_qexp         
( const double x ,
  const double q ) 
{
  if ( 1 == q || s_equal ( 1 , q ) ) { return std::exp ( x ) ; }
  const long double arg = 1 + ( 1.0L - q ) * x ;
  return 0 <= arg ? std::pow ( arg , 1 / ( 1.0L - q ) ) : 0 ; 
}
// ============================================================================
/* q-logarithm in Tsallis statistics 
 *  \f$ \log_q(x) = \frac{x^{1-q}-1}{1-q}\f$ 
 */
// ============================================================================
double Ostap::Math::tsallis_qlog  
( const double x ,
  const double q ) 
{
  if ( 1 == q || s_equal ( 1 , q ) ) { return std::log ( x ) ; }
  const long double xl = x ;
  return ( std::pow ( xl ,  1.0L - q ) - 1 ) / ( 1.0L - q ) ;
}
// ============================================================================
/** (unnormalized) q-gaussian in Tsallis statistics
 *  \f$ \hat{G}_q(x, \beta , q) = e_q( - \left| \beta\right| x^2 ) \f$ 
 *  - for \f$ q=1\f$, it correspoind to Gaussian 
 *  - for \f$ q<1\f$, it it a function with the finite support  
 *  - for \f$ 1<q\f$, it it a variant of generalized Student't distribution
 *  - for \f$ q=2\f$, it it a Cauchy distribution 
 *  - for \f$ 3\le q\f$, it cannot be normalzed 
 *  @see Ostap::Math::tsallis_qexp 
 */
// ============================================================================
double Ostap::Math::tsallis_qgaussian_u   
( const double x     ,
  const double beta  , 
  const double q     ) 
{
  const double absbeta = std::abs ( beta ) ;
  const double arg     = - absbeta * x * x ;
  return Ostap::Math::tsallis_qexp ( arg ) ;
}
// ============================================================================
/* (normalized) q-gaussian in Tsallis statistics for ( q < 3 ) 
 *  \f$ G_q(x, \beta,q) = 
 *   \frac{\sqrt{beta}}{C_q} e_q( - \left| \beta\right| x^2 ) \f$ 
 *  - for \f$ q<1\f$, it it a function with a finite support  
 *  - for \f$ q=1\f$, it is a Gaussian 
 *  - for \f$ q=2\f$, it it a Cauchy distribution 
 *  - for \f$ 1<q<3\f$, it it a variant of generalized Student't distribution
 */
// ============================================================================
double Ostap::Math::tsallis_qgaussian   
( const double x     ,
  const double beta  , 
  const double q     ) 
{
  Ostap::Assert ( q < 3 , 
                  "Invalid value of q (It must be <3)" , 
                  "Ostap::Math::tsallis_qgaussian"  ) ;
  //
  const double absbeta = std::abs ( beta ) ;
  const double arg     = - absbeta * x * x ;
  const double c1      = std::sqrt ( absbeta / M_PI ) ;
  //
  const double result  = c1 * tsallis_qexp ( arg , q ) ;
  //
  if ( 1 == q || s_equal ( 1 , q ) ) { return result ; }
  //
  if ( 1 >  q ) 
  { 
    const long double g1 = std::lgamma  ( 1.0             / ( 1.0L - q ) ) ;
    const long double g2 = std::lgamma  ( 0.5 * ( 3 - q ) / ( 1.0L - q ) ) ;
    const long double cq = 0.5 * std::exp ( g2  + std::log ( 3.0L - q ) +
                                            0.5 * std::log ( 1.0L - q ) - g1 ) ;
    return result * cq ; 
    
  }
  // 1 < q < 3
  const long double g1 = std::lgamma ( 1.0                / ( q - 1.0L ) ) ;
  const long double g2 = std::lgamma ( 0.5 * ( 3.0L - q ) / ( q - 1.0L ) ) ;
  const long double cq = std::exp    ( g1 + 0.5 * std::log  ( q - 1.0L ) - g2 ) ;
  //
  return result * cq ;
  //
}
// ============================================================================
/* (normalized) q-gaussian in Tsallis statistics for \f$ q < 3 \f$, \f$0<\sigma\f$  
 *  \f$ G_q(x,\mu, \sigma, q) = \frac{1}{\sigma} G_q (\frac{ x-\mu}{\sigma} , \frac{1}{2}, q \f$ 
 *  - for \f$ q<1\f$, it it a function with a finite support  
 *  - for \f$ q=1\f$, it is a Gaussian 
 *  - for \f$ q=2\f$, it it a Cauchy distribution 
 *  - for \f$ 1<q<3\f$, it it a variant of generalized Student't distribution
 */
// ============================================================================
double Ostap::Math::tsallis_qgaussian
( const double x     ,
  const double mu    , 
  const double sigma , 
  const double q     ) 
{ return tsallis_qgaussian ( ( x - mu ) / sigma , 0.5 , q ) / std::abs ( sigma ) ; }
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
