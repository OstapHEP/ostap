// ===========================================================================
// Inclode files
// ===========================================================================
/// STD & STL
// ===========================================================================
#include <cmath>
#include <complex>
#include <map>
#include <numeric>
#include <tuple>
// ===========================================================================
// GSL 
// ===========================================================================
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_zeta.h"
// ===========================================================================
// Ostap
// ===========================================================================
#include "Ostap/Zeta.h"
#include "Ostap/Choose.h"
#include "Ostap/Gamma.h"
#include "Ostap/Bernulli.h"
#include "Ostap/ToStream.h"
// ===========================================================================
// Local
// ===========================================================================
#include "syncedcache.h"
#include "GSL_sentry.h"
#include "local_math.h"
#include "local_gsl.h"
#include "status_codes.h"
// ===========================================================================
namespace
{
  // =========================================================================
  /// Use GSL to calculate zeta(n)
  double GSL_zeta ( const int n ) 
  {
    if      (  1 == n ) { return std::numeric_limits<double>::quiet_NaN() ; }
    else if (  0 == n ) { return -0.5         ; }
    else if ( -1 == n ) { return -1.0/12      ; }  
    // trivial zeroes ?
    else if (  n < -1 && 0 == n % 2 ) { return 0 ; }
    //
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
    // 
    return result.val ;
  }
  // =========================================================================
  double GSL_zeta ( const double s   )
  {
    /// integer argument ? 
    if ( Ostap::Math::isint ( s ) )
    {
       const int nn = Ostap::Math::round ( s ) ;
       return GSL_zeta ( nn ) ;
    }
    //
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
  // ===========================================================================
  double GSL_zetam1 ( const int n  )
  {
    if      (  1 == n ) { return std::numeric_limits<double>::quiet_NaN() ; }
    else if (  0 == n ) { return -0.5         -1 ; }
    else if ( -1 == n ) { return -1.0/12      -1 ; }  
    // trivial zeroes ?
    else if (  n < -1 && 0 == n % 2 ) { return -1 ; }
    //
    // (2) use GSL: 
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
   //
   return result.val ;
  }
  // =========================================================================
  double GSL_zetam1 ( const double s   )
  {
    /// integer argument ? 
    if ( Ostap::Math::isint ( s ) )
    {
      const int nn = Ostap::Math::round ( s ) ;
      return ::GSL_zetam1 ( nn ) ;
    }
    //
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
  // =========================================================================
  /* Dirichlet's Eta function 
   *  \f$ \eta ( s ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
   */
  double GSL_eta ( const int n  )
  {
    /// trvial zeroes of zeta-function:
    if ( n < -1 && 0 == n % 2 ) { return 0 ; }
    //  
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
    //
    return result.val ;
  }
  // ===========================================================================
  /* Dirichlet's Eta function 
   *  \f$ \eta ( s ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
   */
  double GSL_eta ( const double s )
  {
    /// integer argument ? 
    if ( Ostap::Math::isint ( s ) )
    {
       const int nn = Ostap::Math::round ( s ) ;
       return GSL_eta ( nn ) ;
    }
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
  // =========================================================================
  /* Hurwitz Zeta function 
   *  \f$ zeta ( s , q ) = \sum_k  ( k + q )^{-s}\f$
   *  - \f$ 1 < s \f$
   *  - \f$ 0 < q \f$
   */
  double GSL_hurwitz 
  ( const double s ,
    const double q )
  {
    if ( s <= 1 || q <= 0 ) { return std::numeric_limits<double>::quiet_NaN() ; }
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
  // =========================================================================
} // =========================================================================
// ===========================================================================
/*  Riemann's Zeta function \f$ n\ne 1\f$:
 *  \f$ \zeta ( n ) = \sum_k k^{-n}\f$ 
 */
// ============================================================================
double Ostap::Math::zeta ( const int n  )
{
  if      (  1 == n ) { return std::numeric_limits<double>::quiet_NaN() ; }
  else if (  0 == n ) { return -0.5         ; }
  else if ( -1 == n ) { return -1.0/12      ; }  
  // trivial zeroes ?
  else if (  n < -1 && 0 == n % 2 ) { return 0 ; }
  //
  typedef int                   KEY    ;
  typedef long double           RESULT ; 
  typedef std::map<KEY,RESULT>  MAP    ;
  typedef SyncedCache<MAP>      CACHE  ;
  /// the cache
  static CACHE                 s_CACHE {} ; // the cache
  static const std::size_t s_MAX_CACHE { 5000 } ; 
  //
  const KEY key { n } ;
  //
  { // (1) check a value already calculated 
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //
  /// (2) get the result 
  const RESULT result = ::GSL_zeta ( n ) ; 
  //
  { // (3) add calculated value into the cache 
    CACHE::Lock lock { s_CACHE.mutex () } ;
    if ( s_MAX_CACHE < s_CACHE->size() ) { s_CACHE->clear() ; }
    s_CACHE->insert ( std::make_pair ( key , result ) ) ;
    return result ; 
  }  
  //
}
// ============================================================================
/*  Riemann's Zeta function \f$ s\ne 1\f$:
 *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
 */
// ============================================================================
double Ostap::Math::zeta ( const double s   )
{
  /// integer argument ? 
  if ( Ostap::Math::isint ( s ) )
  {
    const int nn = Ostap::Math::round ( s ) ;
    return Ostap::Math::zeta ( nn ) ;
  }
  //
  return ::GSL_zeta ( s ) ; 
}
// ===========================================================================
/*  Riemann's Zeta function \f$ s\ne 1\f$:
 *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
 */
// ============================================================================
long double Ostap::Math::zeta ( const long double s  )
{ return zeta ( static_cast<double> ( s ) ) ; }

// ============================================================================
/* Riemann's Zeta function minus 1 \f$ n\ne 1\f$:
 *  \f$ f ( n ) = \zeta ( n ) - 1 = -1 + \sum_k k^{-n}\f$ 
 */
// ============================================================================
double Ostap::Math::zetam1 ( const int n  )
{
  if      (  1 == n ) { return std::numeric_limits<double>::quiet_NaN() ; }
  else if (  0 == n ) { return -0.5         -1 ; }
  else if ( -1 == n ) { return -1.0/12      -1 ; }  
  // trivial zeroes ?
  else if (  n < -1 && 0 == n % 2 ) { return -1 ; }
  //
  typedef int                   KEY    ;
  typedef double                RESULT ; 
  typedef std::map<KEY,RESULT>  MAP    ;
  typedef SyncedCache<MAP>      CACHE  ;
  /// the cache 
  static CACHE                 s_CACHE {} ; // the cache
  static const std::size_t s_MAX_CACHE { 10000 } ; 
  //
  const KEY key { n } ;
  //
  { // (1) check a value already calculated 
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }  
  //
  // (2) calculate the value 
  const RESULT result = ::GSL_zetam1 ( key ) ;
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
/*  Riemann's Zeta function minus 1 \f$ s \ne 1\f$:
 *  \f$ f ( s ) = \zeta ( s ) - 1 = -1 + \sum_k k^{-s}\f$ 
 */
// ============================================================================
double Ostap::Math::zetam1 ( const double s   )
{
  /// integer argument ? 
  if ( Ostap::Math::isint ( s ) )
  {
    const int nn = Ostap::Math::round ( s ) ;
    return Ostap::Math::zetam1 ( nn ) ;
  }
  //
  return ::GSL_zetam1 ( s ) ; 
}
// ============================================================================
/*  Riemann's Zeta function minus 1 \f$ s \ne 1\f$:
 *  \f$ f ( s ) = \zeta ( s ) - 1 = -1 + \sum_k k^{-s}\f$ 
 */
// ============================================================================
long double Ostap::Math::zetam1 ( const long  double s   )
{ return zetam1 ( static_cast<double> ( s ) ) ;  }
// ============================================================================

// ============================================================================
/* Dirichlet's Eta function 
 *  \f$ \eta ( s ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
 */
// ============================================================================
double Ostap::Math::eta ( const int n  )
{
  /// trvial zeroes of zeta-function:
  if ( n < -1 && 0 == n % 2 ) { return 0 ; }
  //
  typedef int                   KEY    ;
  typedef double                RESULT ;
  typedef std::map<KEY,RESULT>  MAP    ;
  typedef SyncedCache<MAP>      CACHE  ;
  /// the cache
  static CACHE                 s_CACHE {} ; // the cache
  static const std::size_t s_MAX_CACHE { 10000 } ;
  //
  const KEY key { n } ;
  //
  { // (1) check a value already calculated 
    CACHE::Lock lock { s_CACHE.mutex () } ;
    auto it = s_CACHE->find ( key ) ;
    if ( s_CACHE->end () != it ) { return it->second ; }   // RETURN 
  }
  //  
  //(2)  use GSL: 
  const RESULT result = ::GSL_eta ( key ) ;
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
/*  Dirichlet's Eta function 
 *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
 */
// ============================================================================
double Ostap::Math::eta ( const double s   )
{
  /// integer argument ? 
  if ( Ostap::Math::isint ( s ) )
  {
    const int nn = Ostap::Math::round ( s ) ;
    return Ostap::Math::eta ( nn ) ;
  }
  //
  return ::GSL_eta ( s ) ; 
}
// ===========================================================================
/*  Dirichlet's Eta function 
 *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
 */
// ============================================================================
long double Ostap::Math::eta ( const long double s   )
{ return eta ( static_cast<double> ( s ) ) ;  }

// ========================================================================
/** Dirichet's beta function
 *  \f[ \beta ( s ) 
 *   \equiv \sum_{n=0}^{+\innfty}\frac{(-1)^n}{(2n+1)^s} \]
 *    \f[ \beta ( s ) = \frac{1}{\Gamma(s)}
 * \int_0^{+\infty}\frac{x^{s-1}\mathrm{e}^{-x}}{1+\mathrm{e}^{-2x}}dx \f]
 *   \f[ \beta ( s ) = 4^{-1}\left( \zeta ( z , \frac{1}{4} - \zeta ( s , frac{3}{4} ) \right) \f] 
 *  where \f$\zeta(a,b)\f$ is Hurwitz's zeta function
 */
// ===========================================================================
double Ostap::Math::dirichlet_beta
( const int n ) 
{ 
  static const unsigned short NN = 10 ;
  static const long double s_beta_d [NN] =  { 
  0.5L ,  
  0.7853981633974483096156608L ,
  0.9159655941772190150546035L ,
  0.9689461462593693804836348L ,
  0.9889445517411053361084226L , 
  0.9961578280770880640063194L , 
  0.9986852222184381354416008L , 
  0.9995545078905399094963465L ,
  0.9998499902468296563380671L ,
  0.9999496841872200898213589L } ;
 //
 if      ( 0 <= n && n < NN     ) { return s_beta_d [ n ] ;}
 else if ( n < 0  && 1 == n % 2 ) { return 0 ; }
 // 
 return dirichlet_beta ( 1.0 * n ) ;
}
// ===========================================================================
/* Dirichet's beta function
 *  \f[ \beta ( s ) 
 *   \equiv \sum_{n=0}^{+\innfty}\frac{(-1)^n}{(2n+1)^s} \]
 *    \f[ \beta ( s ) = \frac{1}{\Gamma(s)}
 * \int_0^{+\infty}\frac{x^{s-1}\mathrm{e}^{-x}}{1+\mathrm{e}^{-2x}}dx \f]
 *   \f[ \beta ( s ) = 4^{-1}\left( \zeta ( z , \frac{1}{4} - \zeta ( s , frac{3}{4} ) \right) \f] 
 *  where \f$\zeta(a,b)\f$ is Hurwitz's zeta function
 */
// ==========================================================================   
double Ostap::Math::dirichlet_beta
( const double x ) 
{ 
  if ( x < 9.001 && Ostap::Math::isint ( x ) )
  {
    const int n = Ostap::Math::round ( x ) ; 
    if      ( 0 <= n && n <= 9     ) { return dirichlet_beta ( n ) ; }
    else if ( n <  0 && 1 == n % 2 ) { return 0 ; }
  }
  return std::pow ( 4.0 , -x ) * ( hurwitz_zeta ( x , 0.25 )  - hurwitz_zeta ( x , 0.75 ) ) ; 
}
// ===========================================================================
/* Dirichet's beta function
 *  \f[ \beta ( s ) 
 *   \equiv \sum_{n=0}^{+\innfty}\frac{(-1)^n}{(2n+1)^s} \]
 *    \f[ \beta ( s ) = \frac{1}{\Gamma(s)}
 * \int_0^{+\infty}\frac{x^{s-1}\mathrm{e}^{-x}}{1+\mathrm{e}^{-2x}}dx \f]
 *   \f[ \beta ( s ) = 4^{-1}\left( \zeta ( z , \frac{1}{4} - \zeta ( s , frac{3}{4} ) \right) \f] 
 *  where \f$\zeta(a,b)\f$ is Hurwitz's zeta function
 */
// ===========================================================================  
long double Ostap::Math::dirichlet_beta
( const long  double x ) 
{ return dirichlet_beta ( static_cast<double> ( x ) ) ; }
// ===========================================================================
  
// ===========================================================================
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
  if ( s <= 1 || q <= 0 ) { return std::numeric_limits<double>::quiet_NaN() ; }
  return ::GSL_hurwitz ( s , q ) ; 
}
// ============================================================================
/* Hurwitz Zeta function 
 *  \f$ zeta ( s , q ) = \sum_k  ( k + q )^{-s}\f$
 *  - \f$ 1 < s \f$
 *  - \f$ 0 < q \f$
 */
// ============================================================================
long double Ostap::Math::hurwitz
( const long double s ,
  const long double q ) 
{ return hurwitz ( static_cast<double> ( s ) ,  static_cast<double> ( q ) ) ;  }
// ===========================================================================

// ===========================================================================
// Complex arguments 
// ===========================================================================
namespace 
{
  // =========================================================================
  typedef std::complex<double>       DC ;
  typedef std::complex<long double> LDC ;
  // =========================================================================
}
// ===========================================================================
/*  Dirichlet's Eta function 
 *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
 */
// ============================================================================
std::complex<double>
Ostap::Math::eta ( const  std::complex<double>& s )
{
  const double x = s.real () ;
  const double y = s.imag () ;
  // real?
  if ( !y || s_zero ( y ) || s_equal ( x + y , x ) ) { return eta ( x ) ; }
  //
  const std::complex<long double> ss { s } ;
  return std::complex<double> ( eta ( ss ) ) ;
}
// ===========================================================================
#include <iostream> 
namespace
{
  // =========================================================================
  /// Laurent expansing for zeta in vicinty of 1
  std::complex<long double>
  laurent_zeta
  ( const  std::complex<long double>& s )
  {
    /// Types :
    typedef std::complex<long double>        TYPE      ;
    static const Ostap::Math::Equal_To<TYPE> xequal {} ;
    static const Ostap::Math::Zero<TYPE>     xzero  {} ; 
    //
    const TYPE ds { s - 1.0L } ;
    //
    /// pole term + 0th order term
    TYPE           result = 1.0L / ds + 1.0L * Ostap::Math::stieltjes ( 0 ) ;
    TYPE           term   = 1.0L ;
    /// number of consequitive "small" terms 
    unsigned short nSmall = 0 ;    
    for  ( unsigned short k = 1 ; k <= Ostap::Math::N_STIELTJES ; ++k )
    {
      // collect factorials 
      term             *= ds   / ( 1.0L * k )  ;
      const TYPE delta  = term * ( 1.0L * Ostap::Math::sign ( k ) * Ostap::Math::stieltjes ( k ) )  ;
      //
      /// the term is zero, or very small or does not change the sum 
      if (  is_not ( delta ) || xzero ( delta ) || xequal ( result + delta , result ) ) { ++nSmall     ; }
      else                                                                              {   nSmall = 0 ; }
      //
      result += delta ;
      /// requre at least teo consequitive small terms 
      if ( 2 <= nSmall )  { break ; }    
    }
    //
    return result ;      
  }
  // =========================================================================
  // Borwain's algorithm for the complex eta 
  // =========================================================================
  template <unsigned short N>
  std::pair<std::complex<long double>, double> 
  borwain_eta
  ( const  std::complex<long double>& s )
  {
    /// Types: 
    typedef std::complex<long double>    TYPE   ;
    ///
    typedef std::array<long double,N+1>  EARRAY ;
    /// the precomputed array 
    static const EARRAY ek { [] () -> EARRAY 
    {
      /// array of scaled binomial coefficients 
      EARRAY ea ;
      for ( unsigned short k = 0 ; k <= N ; ++ k )
      { ea [ k ] = Ostap::Math::choose2 ( N , k ) ; }
      ///
      for ( unsigned short k = 0 ; k <= N ; ++ k )
      { ea [ k ] = std::accumulate ( ea.begin() +  k , ea.end() , 0.0L ) ;}
      ///
      return ea ;
    } () } ;
    //
    std::complex<long double> result = 0 ;
    for ( unsigned short k = 1 ; k <= N ; ++ k )
      {
      const long double c      = 1.0L * Ostap::Math::sign ( k + 1 ) ; 
      const auto        delta1 = c    * std::pow ( 1.0L /   k       , s  ) ;
      const auto        delta2 = c    * std::pow ( 1.0L / ( k + N ) , s  ) * ek [ k ] ;
      //
      result += delta1 ; // the firstt sum 
      result += delta2 ; // the second sum 
    }
    //
    const        float sigma   = s.real () ;
    const        float t       = s.imag () ;
    static const float log8    = std::log ( 8.0f ) ;
    static const float log4    = std::log ( 4.0f ) ;
    static const float pi_half = s_pi_2  ; 
    //
    if ( 0 < s.real () )
    {
      const float e1 = std::log ( 1 + std::abs ( t / sigma  ) ) + std::abs ( t ) * pi_half - N * log8 ; 
      return std::make_pair ( result , std::exp ( e1 ) ) ;
    }    
    // if ( 1.0 -N  < s.real () )
    const float e2 = std::abs ( sigma ) * log4 - N * log8 
      - std::real ( Ostap::Math::lgamma ( s ) ) ;
    //
    return std::make_pair ( result , std::exp ( e2 )  ) ;
    //
  }
  // =========================================================================
  // Borwain's algorithm for the complex eta 
  // =========================================================================
  template <unsigned short N>
  std::pair<std::complex<long double>, double> 
  borwain_zeta
  ( const  std::complex<long double>& s )
  {
    /// (1) calcualate eta-function  
    std::complex<long double> eta_value {} ;
    double                    eta_error {} ; 
    std::tie ( eta_value , eta_error ) =  borwain_eta<N> ( s ) ;
    ///
    const std::complex<long double> scale  { 1.0L - std::pow ( 2.0L , 1.0L - s ) };    
    return std::make_pair ( eta_value / scale , eta_error / std::abs ( scale ) ) ;
  }
  // =========================================================================
  // Euler-Maclaurin algorithm for the the complex zeta 
  // =========================================================================
  template <unsigned short N>
  std::pair<std::complex<long double>,double> 
  emaclaurin_zeta
  ( const std::complex<long double>& s ) 
  {
    /// 2M is an index for the maximal involved Bernulli's number 
    constexpr unsigned short M { Ostap::Math::N_BERNULLI_MAX2 } ;

    /// Types:
    
    typedef std::complex<long double>        TYPE      ;
    static const Ostap::Math::Equal_To<TYPE> xequal {} ;
    static const Ostap::Math::Zero<TYPE>     xzero  {} ; 

    /// useful constant N**(1-s)
    const TYPE NS1 = std::pow ( 1.0L * N , 1.0L - s ) ;
    const TYPE ds  = s - 1.0L ;
    TYPE  result   = NS1 * ( 1.0L / ( 2 * N ) + 1.0L / ds ) ;

    /// (1) direct summation of N-terms [1,N)
    for ( unsigned short k = 1 ; k < N ; ++k )
    { result     += std::pow ( 1.0L / k , s ) ;  } 
    
    /// (2) Bernulli's term
    TYPE term                 = s * NS1 / ( 1.0L * N * N ) ;
    /// absolte value of "previous" step 
    long double     previous = std::numeric_limits<long double>::max() ;
    /// number of consequitive "small" terms
    unsigned short nSmall    = 0 ;
    /// maximal index  before break   (used for error estimation)
    unsigned short rr        = 0 ;
    for ( unsigned short   r = 1 ; r < M ; ++r )
    {
      // 
      rr = r + 1 ; 
      /// Bernulli's index 
      const unsigned short r2 = 2 * r ;
      /// Bernulli's factor: B_2k/(2k)! 
      const long double    B_fct = Ostap::Math::bernulli_k ( r2 ) ;
      /// current trm in te sum 
      const TYPE           delta = B_fct * term ; 
      /// Term is "zero" or such small that does not change the sum ?
      if ( is_not ( delta ) || xzero ( delta ) || xequal ( result  + delta , result ) )  { ++nSmall ; }
      /// update the result 
      result += delta ;      
      /// check the magnitude of the current term 
      const long double abs_delta = std::abs ( delta ) ;
      /// 
      if ( 1 <= nSmall || previous < abs_delta ) { break ; }       // BREAK IS HERE!!!
      //
      /// advance "previous" value 
      previous = abs_delta ; 
      //
      /// update the term 
      term *= ( s + ( 1.0L * r2 - 1.0L )  ) * ( s + 1.0L * r2 ) / ( 1.0L * N * N ) ; 
    }
    ///
    /// Estimate the error term :
    const double sigma = s.real() ;
    const double error = std::abs ( term ) * Ostap::Math::bernulli_k ( 2 * rr ) / ( sigma + rr - 1 ) ;
    return std::make_pair ( result , error ) ;
  };
  // ========================================================================
}
// ===========================================================================
/*  Dirichlet's Eta function 
 *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
 */
// ===========================================================================
std::complex<long double>
Ostap::Math::eta ( const  std::complex<long double>& s )
{
  const long double sigma = s.real () ;
  const long double t     = s.imag () ;
  /// real case?
  if ( !t || s_zero ( t ) || s_equal ( sigma + t , sigma ) ) { return eta ( sigma ) ; }
  //
  
  // (1) estimate the number of terms to get the desired precision with Borwain's method 
  constexpr    float eps     = 1.e-16 ;
  static const float log_eps = std::log ( eps  ) ;
  static const float pi_half = s_pi_2  ;
  static const float log4    = std::log ( 4.0f ) ;
  
  //
  const float tf = t ;
  const float sf = sigma ;
  //
  float n = 0 ;
  if ( 0 < sigma )
  {
    n =
      + std::log ( 1.0f + std::abs ( tf / sf ) )    
      + std::log ( std::abs ( tf ) ) * pi_half 
      - log_eps ; 
  }
  else
  {
    n =
      + std::abs ( sf ) * log4 
      - std::log ( float ( std::abs ( gamma ( s ) ) ) ) 
      - log_eps ;
    //      
  }
  //
  
  /// try to use the functional equation
  if ( sigma <= 0  || s_zero ( sigma ) )
  {    
    return ( std::pow ( 2 , s ) - 1.0L ) * std::pow ( s_1_pi , s ) * gamma ( s ) *
      std::cos ( s_pi_2 * s ) * eta ( 1.0L - s ) ;
  }
  
  if      ( n <  10  && ( 1 - n ) < sigma ) { return borwain_eta<10>   ( s ).first ; }
  else if ( n <  20  && ( 1 - n ) < sigma ) { return borwain_eta<20>   ( s ).first ; }
  else if ( n <  30  && ( 1 - n ) < sigma ) { return borwain_eta<30>   ( s ).first ; }
  else if ( n <  40  && ( 1 - n ) < sigma ) { return borwain_eta<40>   ( s ).first ; }
  else if ( n <  50  && ( 1 - n ) < sigma ) { return borwain_eta<50>   ( s ).first ; }
  else if ( n <  70  && ( 1 - n ) < sigma ) { return borwain_eta<70>   ( s ).first ; }
  else if ( n <  80  && ( 1 - n ) < sigma ) { return borwain_eta<80>   ( s ).first ; }
  else if ( n < 100  && ( 1 - n ) < sigma ) { return borwain_eta<100>  ( s ).first ; }
  else if ( n < 120  && ( 1 - n ) < sigma ) { return borwain_eta<120>  ( s ).first ; }
  else if ( n < 150  && ( 1 - n ) < sigma ) { return borwain_eta<150>  ( s ).first ; }

  //
  return borwain_eta<200>  ( s ).first ; 
}
// ===========================================================================
/*  Riemann's Zeta function \f$ s \ne 1\f$:
 *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
 */
// ===========================================================================
std::complex<double> Ostap::Math::zeta ( const std::complex<double>& s  ) 
{
   const double sigma = s.real () ;
   const double tau   = s.imag () ;
   if ( !tau || s_zero ( tau ) ) { return zeta ( sigma ) ; }
   //
   const std::complex<long double> ss { s } ;
   return std::complex<double> ( zeta ( ss ) ) ;   
}
// ===========================================================================
/*  Riemann's Zeta function \f$ s \ne 1\f$:
 *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
 */
// ==========================================================================
#include <iostream> 
std::complex<long double> 
Ostap::Math::zeta 
( const std::complex<long double>& s  ) 
{
  // 
  const long double sigma = s.real () ;
  const long double tau   = s.imag () ;
  /// real case ? 
  if ( !tau || s_zero ( tau ) || s_equal ( tau + sigma , sigma ) ) { return zeta ( sigma ) ; }
  
  /// Types :
  typedef std::complex<long double>        TYPE      ;
  static const Ostap::Math::Equal_To<TYPE> xequal {} ;
  static const Ostap::Math::Zero<TYPE>     xzero  {} ; 
  
  // ======================================================================
  // (1) try Laurent expansion at s= 1 
  //      Usually it it not the best approach, but here it should be cheap,
  //      since all involved Stielties' constants are cached.
  // ======================================================================
  // The actual convergency range is:  |1-s|<2pi
  const TYPE ds { s - 1.0L } ;  
  if ( std::abs ( ds ) <= 1.5 * s_pi ) { return ::laurent_zeta ( s ) ; }

  /// use the reflection formula if real part is non-positive 
  if ( sigma <= 0 || s_zero ( sigma ) ) { return chi ( s ) * zeta ( 1.0L - s ) ; }

  /// check the number of terms in Borwain's algorithm for eta 
  const unsigned int n =  15 + std::ceil ( 1.2 * std::abs ( tau ) ) ;

  if      ( n <   5  && ( 1 - n ) < sigma ) { return borwain_zeta<5>   ( s ).first ; }
  else if ( n <  10  && ( 1 - n ) < sigma ) { return borwain_zeta<10>  ( s ).first ; }
  else if ( n <  20  && ( 1 - n ) < sigma ) { return borwain_zeta<20>  ( s ).first ; }
  else if ( n <  30  && ( 1 - n ) < sigma ) { return borwain_zeta<30>  ( s ).first ; }
  else if ( n <  40  && ( 1 - n ) < sigma ) { return borwain_zeta<40>  ( s ).first ; }
  else if ( n <  50  && ( 1 - n ) < sigma ) { return borwain_zeta<50>  ( s ).first ; }
  else if ( n <  70  && ( 1 - n ) < sigma ) { return borwain_zeta<70>  ( s ).first ; }
  else if ( n <  80  && ( 1 - n ) < sigma ) { return borwain_zeta<80>  ( s ).first ; }
  else if ( n < 100  && ( 1 - n ) < sigma ) { return borwain_zeta<100> ( s ).first ; }
  else if ( n < 120  && ( 1 - n ) < sigma ) { return borwain_zeta<120> ( s ).first ; }
  else if ( n < 150  && ( 1 - n ) < sigma ) { return borwain_zeta<150> ( s ).first ; }
  else if ( n < 200  && ( 1 - n ) < sigma ) { return borwain_zeta<200> ( s ).first ; }

  auto r10  = emaclaurin_zeta<10>  ( s ) ; 
  auto r15  = emaclaurin_zeta<15>  ( s ) ; 
  auto r20  = emaclaurin_zeta<20>  ( s ) ; 
  auto r25  = emaclaurin_zeta<25>  ( s ) ; 
  auto r50  = emaclaurin_zeta<50>  ( s ) ; 
  auto r75  = emaclaurin_zeta<75>  ( s ) ; 
  auto r100 = emaclaurin_zeta<100> ( s ) ; 

  /** 
  std::cerr << " EM[10]  : " ;
  Ostap::Utils::toStream ( r10 , std::cerr ) << std::endl ;
  std::cerr << " EM[15]  : " ;
  Ostap::Utils::toStream ( r15 , std::cerr ) << std::endl ;
  std::cerr << " EM[20]  : " ;
  Ostap::Utils::toStream ( r20 , std::cerr ) << std::endl ;
  std::cerr << " EM[50]  : " ;
  Ostap::Utils::toStream ( r25 , std::cerr ) << std::endl ; 
  std::cerr << " EM[75]  : " ;
  Ostap::Utils::toStream ( r75 , std::cerr ) << std::endl ;
  std::cerr << " EM[100] : " ;
  Ostap::Utils::toStream ( r100 , std::cerr ) << std::endl ;
  */

  
  return r100.first ;
  
  /// Use Borwain' method for precise calculation of eta:
  Ostap::Assert ( false , 
		  "zeta is not yet implelmented for complex argument!" , 
		  "Ostap::Math::zeta" , 
                   NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
  
   //
  return 0 ; 
}
// ===========================================================================

// ===========================================================================
/* helper chi-function:
 * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{2} \Gamma (1-s), \f]
 * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s ) f$  
 */
// ============================================================================
double Ostap::Math::chi ( const double s ) 
{ return chi ( s * 1.0L ) ; }
// ========================================================================
/* helper chi-function:
 * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{s} \Gamma (1-s), \f]
 * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s \f$  
 */
// ======================================================================== 
long double Ostap::Math::chi ( const long double s ) 
{
  if ( isint ( s ) )
  {
    const int n = round ( s ) ;
    // poles 
    if      ( n >= 1 && 1 == n % 2 ) { return std::numeric_limits<double>::quiet_NaN () ; }
    // finite value 
    else if ( n >= 2 && 0 == n % 2 ) { return std::pow ( s_2pi , n ) * 0.5L * igamma ( n )  * sign ( n / 2 ) ; }
    // zeroes 
    else if ( n <= 0 && 0 == n % 2 ) { return 0 ; }
  }
  //
  return std::pow ( s_2pi , s ) * s_1_pi * std::sin ( s_pi * s / 2 ) *  gamma ( 1.0L - s ) ;  
}
// ========================================================================
/* helper chi-function:
 * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{s} \Gamma (1-s), \f]
 * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s \f$  
 */
// ===========================================================================
std::complex<double> 
Ostap::Math::chi 
( const std::complex<double>&  s ) 
{ 
    //
    const double x = s.real () ;
    const double y = s.imag () ;
    if ( !y ||  s_zero ( y ) ) { return chi ( x ) ; }
    //
    const LDC ss { s } ; 
    return DC ( chi ( ss ) ) ;  
}
// ==========================================================================
/* helper chi-function:
 * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{s} \Gamma (1-s), \f]
 * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s \f$  
 */
// ===========================================================================
std::complex<long double> 
Ostap::Math::chi 
( const std::complex<long double>&  s ) 
{ 
    //
    const long double x = s.real () ;
    const long double y = s.imag () ;
    if ( !y ||  s_zero ( y ) ) { return chi ( x ) ; }
    //
    return std::pow ( s_2pi , s ) * s_1_pi * std::sin ( s_pi * s * 0.5L ) * gamma ( 1.0L - s ) ;  
} 
  

// ===========================================================================
//                                                                     The END
// ===========================================================================

