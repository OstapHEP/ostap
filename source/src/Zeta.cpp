// ===========================================================================
// Inclode files
// ===========================================================================
/// STD & STL
// ===========================================================================
#include <cmath>
#include <complex>
#include <map>
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
    typedef std::complex<double>       DC ;
    typedef std::complex<long double> LDC ;
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
   const LDC ss { s } ;
   return DC ( zeta ( ss ) ) ;   
}
// ===========================================================================
/*  Riemann's Zeta function \f$ s \ne 1\f$:
 *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
 */
// ==========================================================================
std::complex<long double> 
Ostap::Math::zeta 
( const std::complex<long double>& s  ) 
{
   const long double sigma = s.real () ;
   const long double tau   = s.imag () ;
   if ( !tau || s_zero ( tau ) ) { return zeta ( sigma ) ; }
    
   // 
   static const Ostap::Math::Equal_To<LDC> xequal {} ;
   static const Ostap::Math::Zero<LDC>     xzero  {} ; 
   //
    
   // ======================================================================== 
   /// (1) try the Laurent expansion at s= 1 
   //      Usually it it not the best approach, but here it should be cheap,
   //      since all involved Stielties' constants are cached.
   // ========================================================================
   const  LDC ds { s - 1.0L } ;
   if ( std::abs ( ds ) <= 1 * s_pi ) // the actual convergency range is |1-s|<2pi
   {
     LDC            result = 1.0L / ds + 1.0L * stieltjes ( 0 ) ; // pole term + 0 order term  
     LDC            term   = 1.0L ;
     unsigned short nSmall = 0 ; 
     for  ( unsigned short k = 1 ; k <= N_STIELTJES ; ++k )
     {
        // collect factorials 
        term            *= ds   / ( 1.0L * k )  ;
        const LDC delta  = term * ( 1.0L * sign ( k ) * stieltjes ( k ) )  ;
        //
        if (   xzero ( delta ) || xequal ( result + delta , result ) ) { ++nSmall     ; }
        else                                                           {   nSmall = 0 ; }
        //
        result += delta ;
        //
        if ( 2 <= nSmall )  { break ; }    
     }
     //
     return result ;
   }
   //
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
double Ostap::Math::chi ( const double s ) { return chi ( s * 1.0L ) ; }
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
        else if ( n >= 1 && 0 == n % 2 ) { return std::pow ( s_2pi , n ) * 0.5L * igamma ( n )  * sign ( n / 2 ) ; }
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

