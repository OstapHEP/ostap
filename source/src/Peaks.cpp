// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <array>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Peaks.h"
#include "Ostap/Tails.h"
#include "Ostap/StatusCode.h"
#include "Ostap/MoreMath.h"
#include "Ostap/qMath.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Polynomials.h"
// ============================================================================
//  Local
// ============================================================================
#include "local_math.h"
#include "local_gsl.h"
#include "local_hash.h"
#include "gauss.h"      
#include "Integrator1D.h"      
#include "syncedcache.h"  // the cache 
#include "status_codes.h" // the cache 
// ============================================================================
/** @file 
 *  implmentation file for classes from the file Ostap/Peaks.h
 *  @author 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-04-19
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /** evaluate the helper function \f[ f = \frac{\sinh (x) }{x} \f]
  *  it allows to calculate Novosibirsk's function in efficient and regular way
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-05-23
  */
  double x_sinh ( const double x , double precision = s_APRECISION )
  {
    //
    if      ( s_equal   ( x , 0 )    ) { return 1 ; }
    else if ( std::abs ( x ) < 0.1  )
    {
      double result = 1.0  ;
      double delta  = x    ;
      //
      precision = std::abs ( precision ) ;
      precision = std::min ( precision , std::abs ( s_APRECISION_TAIL ) ) ;
      unsigned int n = 1 ;
      //
      do
      {
        delta  *= x * x / ( ( n + 1 )  * ( n + 2 ) ) ;
        result += delta ;
        n      += 2 ;
      }
      while ( std::abs ( delta ) > 0.1 * precision && n < 10000 ) ;
      //
      return result ;
    }
    //
    if ( 100 < std::abs ( x )  ) { return  s_INFINITY ; }
    //
    // the generic evaluation
    //
    return std::sinh ( x ) / x ;
  }
  // ==========================================================================
  // Crystal Ball & Co 
  // ==========================================================================
  /** @var s_TRUNC
   *  trunkating parameter for CrystalBall-functions
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double      s_TRUNC = 15.0 ;
  // ==========================================================================
  /** evaluate very simple power-law integral
   *  
   *  \f[ I = \int_{x_{low}}^{x_{high}} \left( \frac{A}{B+Cx}\right)^{N} dx \f]
   * 
   *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
   */
  double tail_integral 
  ( const double A    , 
    const double B    , 
    const double C    , 
    const double N    , 
    const double low  , 
    const double high )
  {
    //
    // few really very simple cases:
    if      ( s_equal ( N , 0 ) ) { return high - low ; }
    else if ( s_equal ( A , 0 ) ) { return 0 ; }
    else if ( s_equal ( C , 0 ) ) { return std::pow ( A / B , N ) * ( high - low ) ; }
    //
    // again the trivial cases 
    if ( s_equal ( low , high ) ) { return 0 ; }
    else if      ( low > high   ) { return -1 * tail_integral ( A , B , C , N , high , low ) ; }
    //
    //  y = (B+C*x)/A
    //
    const double y_low  = ( B + C * low  ) / A ;
    const double y_high = ( B + C * high ) / A ;
    //
    // the special case 
    if ( s_equal ( N , 1 ) ) { return A / C * my_log ( y_high / y_low ) ; } // RETURN 
    //
    // the regular case 
    return A / C * ( std::pow ( y_high , 1 - N ) - 
                     std::pow ( y_low  , 1 - N ) ) / ( 1 - N ) ;
  }
  // ========================================================================== 
  // Atlas/Zeus& Co 
  // ==========================================================================
  /** @var S_ATLAL 
   *  magic constant - integral for Atlas function 
   *  @see Ostap::Math::Atlas 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-08-21
   */   
  const double s_ATLAS = 3.052369876253939 ;
  // ==========================================================================
  // Sinh-asinh 
  // ==========================================================================
  inline double shash
  ( const double x   , 
    const double eps , 
    const double dlt ) 
  {
    const double y = eps + dlt * std::asinh ( x ) ;
    return 
      (     GSL_LOG_DBL_MAX < y ) ?    s_INFINITY :
      ( -1* GSL_LOG_DBL_MAX > y ) ? -1*s_INFINITY : std::sinh ( y ) ;
  }
  // // ==========================================================================
  // /en-t'
  // // ==========================================================================
  // inline double student_cdf (  const double t , const double nu ) 
  // {
  //   const double xt    = nu / ( t * t + nu ) ;
  //   const double value = 0.5 * gsl_sf_beta_inc ( 0.5 * nu , 0.5 , xt ) ;
  //   return t >= 0 ? 1 - value : value ;
  // }
  // ==========================================================================
  /** product of Gaussian PDF and Mill's ratio 
   *  \f$ f(x; a ) = \phi ( x ) R( a + x ) \f$ 
   */
  // inline double gauss_mills 
  // ( const double x , 
  //  const double a ) 
  // { 
  //   return Ostap::Math::gauss_pdf ( x ) * Ostap::Math::mills_normal ( a + x ) ; 
  // } ;
  // ==========================================================================
}
// ============================================================================

// ============================================================================
// BifurcatedGauss
// ============================================================================
/*  constructor from all parameters
 *  @param peak    the peak posiion
 *  @param sigma_L (left sigma)
 *  @param sigma_R (right sigma)
 */
// ============================================================================
  Ostap::Math::BifurcatedGauss::BifurcatedGauss
  ( const double peak   ,
    const double sigmaL ,
    const double sigmaR )
  : m_peak    ( peak )
  , m_sigmaL  ( -1 )
  , m_sigmaR  ( -1 )
  , m_kappa   (  0 )
  , m_psi     (  0 )  
{
  setSigma ( sigmaL , sigmaR ) ;
}
// ============================================================================
/*  constructor from all parameters
 *  @param peak    the peak posiion
 *  @param sigma   avegrate sigma 
 */
// ============================================================================
Ostap::Math::BifurcatedGauss::BifurcatedGauss
( const double peak   ,
  const double sigma  )
: BifurcatedGauss ( peak , sigma , sigma ) 
{}
// ============================================================================
/*  constructor from all parameters
 *  @param peak    gaussian 
 */
// ============================================================================
Ostap::Math::BifurcatedGauss::BifurcatedGauss
( const Ostap::Math::Gauss& gauss ) 
: BifurcatedGauss ( gauss.peak() , gauss.sigma() , gauss.sigma() ) 
{}
// =====================================================================
// evaluate Bifurcated Gaussian
// ============================================================================
double Ostap::Math::BifurcatedGauss::evaluate ( const double x ) const
{
  const double dx = 
    ( x < m_peak ) ?
    ( x - m_peak ) / m_sigmaL :
    ( x - m_peak ) / m_sigmaR ;
  //
  const double norm = s_SQRTPIHALF * ( m_sigmaL + m_sigmaR ) ;
  //
  return std::exp ( -0.5 * dx * dx ) / norm ;
}
// ============================================================================
/*  log-derivative \f$ \frac{ f^\prime}{f}  \f$
 *  Useful to attach the radiative tail to ensure 
 *  the continuity of the function and the 1st derivatibve 
 */
// ============================================================================
double Ostap::Math::BifurcatedGauss::dFoF ( const double x ) const
{
  const double dx = 
    ( x < m_peak ) ?
    ( x - m_peak ) / ( m_sigmaL * m_sigmaL ) :
    ( x - m_peak ) / ( m_sigmaR * m_sigmaR ) ;
  //
  return -dx ; 
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::BifurcatedGauss::integral () const { return 1 ; }
// ============================================================================
//  get CDF 
// ============================================================================
double Ostap::Math::BifurcatedGauss::cdf ( const double x ) const 
{
  // left half-gaussian
  if ( x <= m_peak )
    {
      const double sigma = sigmaL () ;
      const double sf    = s_SQRT2i / sigma  ;
      const double nf    = sigma    / ( sigmaL() + sigmaR() ) ;
      const double b     = ( x - m_peak ) * sf ;
      return std::erfc ( -b ) * nf ; // RETURN
    }
  //
  const double bias = sigmaL() / ( sigmaL() + sigmaR() ) ;
  return bias + integral ( m_peak , x ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::BifurcatedGauss::integral
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }                         // RETURN
  if (           low > high   ) { return - integral ( high , low ) ; } // RETURN
  //
  // left half-gaussian
  if       ( high <= m_peak )
    {
      const double sigma = sigmaL () ;
      const double sf    = s_SQRT2i / sigma  ;
      const double nf    = sigma    / ( sigmaL() + sigmaR() ) ;
      const double a     = ( low  - m_peak ) * sf ;
      const double b     = ( high - m_peak ) * sf ;
      return  ( std::erf ( b ) -  std::erf ( a ) ) * nf ; // RETURN
    }
  // right half-gaussian
  else if ( low >= m_peak )
    {
      const double sigma = sigmaR () ;
      const double sf    = s_SQRT2i / sigma  ;
      const double nf    = sigma    / ( sigmaL() + sigmaR() ) ;
      const double a     = ( low  - m_peak ) * sf ;
      const double b     = ( high - m_peak ) * sf ;
      return  ( std::erf ( b ) -  std::erf ( a ) ) * nf ; // RETURN
    }
  // split into two intervals 
  return
    integral ( low    , m_peak ) +
    integral ( m_peak , high   ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::BifurcatedGauss::tag () const 
{
  static const std::string s_name = "BiFurcatedGauss" ;
  return Ostap::Utils::hash_combiner ( s_name   , 
                                       m_peak   , 
                                       m_sigmaL , 
                                       m_sigmaR ) ; 
}
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setSigma
 ( const double valueL ,
   const double valueR )
{
  const double valueL_ = std::abs ( valueL ) ;
  const double valueR_ = std::abs ( valueR ) ;
  //
  if ( s_equal ( m_sigmaL , valueL_ ) && 
       s_equal ( m_sigmaR , valueR_ ) ) { return false ; }
  //
  Ostap::Assert ( valueL_ && valueR_ ,
		  "Parameters 'sigmaL/R' must be non-zero"    ,
		  "Ostap::Math::BifurcatedGauss::setSigma" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //
  m_sigmaL = valueL_ ;
  m_sigmaR = valueR_ ;
  //
  m_kappa  = ( m_sigmaL - m_sigmaR ) / ( m_sigmaL + m_sigmaR ) ;
  m_psi    = m_kappa ? std::atanh ( m_kappa ) : 0.0 ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setSigmaL ( const double value )
{ 
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( m_sigmaL , value_ ) ) { return false ; }
  return setSigma ( value , m_sigmaR ) ;
} 
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setSigmaR ( const double value )
{ 
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( m_sigmaL , value_ ) ) { return false ; }
  return setSigma ( m_sigmaL , value ) ; 
}
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setPeak( const double value )
{
  if ( s_equal ( m_peak , value ) ) { return false ; }
  m_peak   = value  ;
  //
  return true ;
}
// ============================================================================
/*  set asymmetry keeping average sigma untouched
 *  \f$ \left| \kappa \right| < 1 \f$ 
 */
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setKappa
( const double value )
{
  // =========================================================================
  Ostap::Assert ( std::abs ( value ) < 1                   ,
		  "Parameter 'kappa' must be |kappa|<1"    ,
		  "Ostap::Math::BifurcatedGauss::setKappa" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //
  if ( s_equal ( value , m_kappa ) ) { return false ; }
  //  
  const double s = sigma () ; 
  //
  m_kappa  = value ;
  m_psi    = m_kappa ? std::atanh  ( m_kappa ) : 0.0 ; 
  // 
  m_sigmaL = s * ( 1 + m_kappa ) ;
  m_sigmaR = s * ( 1 - m_kappa ) ;
  //
  return true ;
}
// ============================================================================
/*  set asymmetry keeping average sigma untouched
 *  \f$ \left| \kappa \right| < 1 \f$ 
 */
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setPsi
( const double value )
{ 
  if ( s_equal ( m_psi , value ) ) { return false ; } 
  //
  const double s = sigma () ; 
  //
  m_psi   = value                ;
  m_kappa = m_psi ? std::tanh ( m_psi ) : 0.0  ;
  //
  m_sigmaL = s * ( 1 + m_kappa ) ;
  m_sigmaR = s * ( 1 - m_kappa ) ; 
  //
  return true ;
}
// ============================================================================
/* constructor from all parameters
 *  @param peak     the peak posiion
 *  @param sigmaL   the sigma for first component
 *  @param fraction the fraction of first component 
 *  @param scale    the ratio of sigmas for second and first components
 */
// ============================================================================
Ostap::Math::DoubleGauss::DoubleGauss
( const double peak     ,
  const double sigma    , 
  const double fraction , 
  const double scale    ) 
  : m_peak     ( peak               ) 
  , m_sigma    ( std::abs ( sigma ) )
  , m_fraction ( std::min ( std::max ( fraction , 0.0 ), 1.0 ) ) 
  , m_scale    ( std::abs ( scale ) )
{}
// ============================================================================
bool Ostap::Math::DoubleGauss::setPeak   ( const double value ) 
{
  if ( s_equal ( value , m_peak ) ) { return false ; }
  m_peak = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::DoubleGauss::setSigma ( const double value ) 
{
  const double value_ =  std::abs (value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::DoubleGauss::setSigma"     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //    
  m_sigma = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::DoubleGauss::setFraction ( const double value ) 
{
  const double value_ = std::min ( std::max ( value , 0.0 ), 1.0 ) ;
  if ( s_equal ( value_ , m_fraction ) ) { return false ; }
  m_fraction = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::DoubleGauss::setScale ( const double value ) 
{
  const double value_ =  std::abs (value ) ;
  if ( s_equal ( value_ , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'scale' must be non-zero"     ,
		  "Ostap::Math::DoubleGauss::setScale"     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //    
  m_scale = value_ ;
  return true ;
}
// ============================================================================
// evaluate  double Gaussian
// ============================================================================
double Ostap::Math::DoubleGauss::pdf ( const double x ) const 
{
  const double mu       = m_peak     ;
  const double sigma    = m_sigma    ;
  const double scale    = m_scale    ;
  const double fraction = m_fraction ;
  //
  const double sigma2   =  scale *  sigma ;
  //
  const double dx1      = ( x - mu ) / sigma  ;
  const double dx2      = ( x - mu ) / sigma2 ;
  //
  const double  f1 = fraction ;
  const double  f2 = 1 - f1   ;
  //
  static const double s_norm = 1.0 / std::sqrt ( 2.0 * M_PI  ) ;
  //
  return 
    s_norm * ( f1 * std::exp ( -0.5 * dx1 * dx1 ) / sigma  +
               f2 * std::exp ( -0.5 * dx2 * dx2 ) / sigma2 ) ;  
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::DoubleGauss::integral 
( const double xmin , 
  const double xmax ) const 
{
  const double mu       = m_peak     ;
  const double sigma    = m_sigma    ;
  const double scale    = m_scale    ;
  const double fraction = m_fraction ;
  //
  const double sigma2   =  scale *  sigma ;
  //
  const double  f1 = fraction ;
  const double  f2 = 1 - f1 ;
  //
  static const double s_isqrt2 = 1.0 / std::sqrt ( 2.0 ) ;
  //
  const double ixscale1 = s_isqrt2 / sigma  ;
  const double ixscale2 = s_isqrt2 / sigma2 ;
  //
  const double r1 = 
    std::erf ( ( xmax - mu ) * ixscale1 ) - 
    std::erf ( ( xmin - mu ) * ixscale1 ) ;
  //
  const double r2 = 
    std::erf ( ( xmax - mu ) * ixscale2 ) - 
    std::erf ( ( xmin - mu ) * ixscale2 ) ;
  //
  return 0.5 * ( f1 * r1 + f2 * r2  ) ;
} 
// ============================================================================
// get cdf 
// ============================================================================
double Ostap::Math::DoubleGauss::cdf ( const double x )  const 
{
  const double mu       = m_peak     ;
  const double sigma    = m_sigma    ;
  const double scale    = m_scale    ;
  const double fraction = m_fraction ;
  //
  const double sigma2   =  scale *  sigma ;
  //
  const double  f1 = fraction ;
  const double  f2 =  1 - f1  ;
  //
  static const double s_isqrt2 = 1.0 / std::sqrt ( 2.0 ) ;
  //
  const double ixscale1 = s_isqrt2 / sigma  ;
  const double ixscale2 = s_isqrt2 / sigma2 ;
  //
  const double r1 = std::erf ( ( x - mu ) * ixscale1 ) ;
  const double r2 = std::erf ( ( x - mu ) * ixscale2 ) ;
  //
  return  0.5 * ( f1 * ( r1 + 1  ) + f2 * ( r2  + 1 ) ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::DoubleGauss::tag () const 
{
  static const std::string s_name = "DoubleGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_peak , m_sigma , m_fraction , m_scale ) ; 
}
// ============================================================================

// ============================================================================
// Gauss
// ============================================================================
/*  constructor from all parameters
 *  @param peak    the peak posiion
 *  @param sigma_L (left sigma)
 *  @param sigma_R (right sigma)
 */
// ============================================================================
Ostap::Math::Gauss::Gauss
( const double peak   ,
  const double sigma  )
  : m_peak  ( peak )
  , m_sigma ( std::abs ( sigma ) )
  //
{
 Ostap::Assert ( m_sigma ,
                  "Invalid parameter `sigma` : must be non-zero!" ,
                  "Ostap::Math::Gauss" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
// evaluate Bifurcated Gaussian
// ============================================================================
double Ostap::Math::Gauss::evaluate ( const double x ) const
{ return Ostap::Math::gauss_pdf ( x , m_peak , m_sigma ) ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Gauss::integral () const { return 1 ; }
// ============================================================================
//  get CDF 
// ============================================================================
double Ostap::Math::Gauss::cdf ( const double x ) const 
{ return Ostap::Math::gauss_cdf ( x , m_peak , m_sigma ) ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Gauss::integral
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }                         // RETURN
  return Ostap::Math::gauss_int ( low , high , m_peak , m_sigma ) ;
}
// ============================================================================
bool Ostap::Math::Gauss::setSigma ( const double value )
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( m_sigma , value_ ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Gauss::setSigma"           ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //      
  m_sigma = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Gauss::setPeak ( const double value )
{
  if ( s_equal ( m_peak , value ) ) { return false ; }
  m_peak   = value  ;
  return true ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Gauss::tag () const 
{ 
  static const std::string s_name = "Gauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_peak , m_sigma ) ; 
}
// ============================================================================
/* get the logarithmic derivative 
 * \f$ \frac{ f^\prime}{f}  \f$
 *  usefut to attach the radiative tail to ensure 
 *  the continuity of the function and the 1st derivatibve 
 */  
// ===========================================================================
double Ostap::Math::Gauss::dFoF
( const double x ) const 
{ return - ( x - m_peak ) / ( m_sigma * m_sigma ) ; } 

// ============================================================================
/*  constructor from all agruments 
 *  @param mu     location/peak posiiton 
 *  @param alpha  "scale" parameter 
 *  @param beta   "shape" parameter 
 */
// ============================================================================
Ostap::Math::GenGaussV1::GenGaussV1 
( const double mu    ,
  const double alpha , 
  const double beta  ) 
  : m_mu     ( mu                 ) 
  , m_alpha  ( std::abs ( alpha ) ) 
  , m_beta   ( std::abs ( beta  ) ) 
  , m_gbeta1 ( 0 )
  , m_gbeta2 ( 0 )
{
  //
  setBeta ( beta ) ;
  //
}
// ============================================================================
bool Ostap::Math::GenGaussV1::setMu        ( const double value ) 
{
  if ( s_equal ( value , m_mu) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV1::setAlpha     ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV1::setBeta     ( const double value ) 
{
  //
  const double value_ = std::max ( std::abs ( value ) , 1.5/GSL_SF_GAMMA_XMAX ) ;
  //
  if ( s_equal ( value_ , m_beta ) ) { return false ; }
  //
  m_beta = value_ ;
  //
  if   ( beta() * GSL_SF_GAMMA_XMAX < 6 ) 
  { 
    m_gbeta1  = 0 ;
    m_gbeta2  = std::lgamma ( 3 / beta() ) ;    
    m_gbeta2 -= std::lgamma ( 1 / beta() ) ;
    m_gbeta2  = my_exp      ( m_gbeta2   ) ;
  }
  else 
  { 
    m_gbeta1 = 1. / std::tgamma ( 1 / beta() )            ;
    m_gbeta2 =      std::tgamma ( 3 / beta() ) * m_gbeta1 ;
  }
  //
  return true ;
}
// ============================================================================
double Ostap::Math::GenGaussV1::pdf ( const double x ) const 
{
  //
  const double delta  = std::abs ( x - m_mu )         ;
  const double delta1 =            delta  / m_alpha   ;
  const double delta2 = std::pow ( delta1 , m_beta  ) ;
  //
  if ( delta2 > 60 || 0 == m_gbeta1 || beta() * GSL_SF_GAMMA_XMAX < 4 ) 
  {
    double result  = std::log    ( 0.5 * beta() / alpha() ) ;
    result        -= delta2 ;
    result        -= std::lgamma ( 1            / beta()  ) ;
    return my_exp ( result ) ;
  }
  //
  double result   = 0.5 * beta() / alpha() ;
  result         *= my_exp   ( -delta2 ) ;
  result         *= m_gbeta1  ;
  //
  return result ;
}
// ============================================================================
double Ostap::Math::GenGaussV1::cdf ( const double x ) const 
{
  //
  const double delta  = std::abs ( x - m_mu )         ;
  const double delta1 =            delta  / m_alpha   ;
  const double delta2 = std::pow ( delta1 , m_beta  ) ;
  //
  const double c = 0.5 * gsl_sf_gamma_inc_P ( 1/beta() , delta2 ) ;
  //
  return x < m_mu ?  0.5 - c : 0.5 + c ;
}
// ============================================================================ 
double Ostap::Math::GenGaussV1::integral ( const double low  , 
                                           const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================ 
double  Ostap::Math::GenGaussV1::variance () const 
{ return alpha() * alpha() *             m_gbeta2   ; }
// ============================================================================
double  Ostap::Math::GenGaussV1::sigma    () const 
{ return alpha()           * std::sqrt ( m_gbeta2 ) ; }
// ============================================================================
double  Ostap::Math::GenGaussV1::kurtosis () const 
{
  double result  =   std::lgamma ( 5 / beta() ) ;
  result        +=   std::lgamma ( 1 / beta() ) ;
  result        -= 2*std::lgamma ( 3 / beta() ) ;
  //
  return my_exp ( result ) - 3 ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::GenGaussV1::tag () const 
{ 
  static const std::string s_name = "GenGaussV1" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_alpha , m_beta ) ; 
}
// ============================================================================

// ============================================================================
/* constructor from all agruments 
 *  @param xi     location/peak posiiton 
 *  @param alpha  "scale" parameter 
 *  @param kappa  "shape" parameter 
 */
// ============================================================================
Ostap::Math::GenGaussV2::GenGaussV2 
( const double xi    ,
  const double alpha , 
  const double kappa )  // kappa=0 correponds to gaussian  
  : m_xi      ( xi                 ) 
  , m_alpha   ( std::abs ( alpha ) ) 
  , m_kappa   (            kappa   ) 
{
  //
  setKappa ( kappa ) ;
  //
}
// ============================================================================
bool Ostap::Math::GenGaussV2::setXi        ( const double value ) 
{
  if ( s_equal ( value , m_xi) ) { return false ; }
  m_xi = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV2::setAlpha     ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV2::setKappa  ( const double value ) 
{
  //
  double value_ = value ;
  //
  if ( s_equal ( value_ , 0       ) ) { value_ = 0 ; }
  if ( s_equal ( value_ , m_kappa ) ) { return false ; }
  //
  m_kappa = value_;
  //
  return true ;
}
// ============================================================================
double  Ostap::Math::GenGaussV2::y ( const double x ) const
{
  if  ( s_equal( m_kappa , 0 ) ) { return ( x - xi() ) / alpha() ; }
  //
  const double delta =  - ( x - xi () ) * kappa() / alpha () ;
  //
  return 
    delta > 1 ?
    -std::log   ( 1 + delta ) / kappa() :
    -std::log1p (     delta ) / kappa() ;  
}
// ============================================================================
double Ostap::Math::GenGaussV2::pdf ( const double x ) const 
{
  //
  if      ( s_equal ( m_kappa , 0  ) ) {}
  else if ( m_kappa * x >= m_kappa * m_xi + m_alpha ) { return 0 ; } // cover both cases(?)
  //
  const double y_   = y ( x ) ;
  //
  const double gau  = my_exp ( -0.5 * y_ * y_ ) / s_SQRT2PI ;
  //
  return gau / ( alpha() - kappa() * ( x - xi () ) ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::cdf ( const double x ) const 
{
  //
  if      ( s_equal ( m_kappa , 0 ) ) {}
  else if ( kappa () > 0 && ( m_kappa * x >= m_kappa * m_xi + m_alpha ) ) { return 1 ; }
  else if ( kappa () < 0 && ( m_kappa * x >= m_kappa * m_xi + m_alpha ) ) { return 0 ; }
  //
  const double y_   = y ( x ) ;
  //
  const double e_   = std::erf ( y_ * s_SQRT2i ) ;
  //
  return 0.5 * ( 1 + e_ ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::integral ( const double low  , 
                                           const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::mean () const 
{
  if ( s_equal ( kappa () , 0 ) ) { return xi () ; }
  //
  const double k2 = 0.5 * kappa() * kappa() ;
  //
  return xi() - 0.5 * alpha() * kappa() * exprel ( k2 ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::variance () const 
{
  if ( s_equal ( kappa() , 0 ) ) { return alpha () * alpha() ; }
  //
  const double k2 = kappa() * kappa() ;
  //
  return alpha () * alpha() * std::exp ( k2 ) * exprel ( k2 ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::sigma () const 
{ return std::sqrt ( variance() ) ; }
// ============================================================================
double Ostap::Math::GenGaussV2::skewness () const 
{
  const double k2     = kappa () * kappa() ;
  //
  const double a1     = exprel (     k2 ) ;
  const double a3     = exprel ( 3 * k2 ) ;
  //
  const double a      = std::pow ( a1 , 1.5 ) ;
  //
  const double result = 3 *  ( a1 - a3 ) / a ;
  //
  return kappa() * result ;
  //  
}
// ============================================================================
double Ostap::Math::GenGaussV2::kurtosis () const 
{
  //
  const double ek2 = my_exp ( kappa() * kappa() ) ;
  //
  return  Ostap::Math::POW ( ek2 , 4 )  
    + 2 * Ostap::Math::POW ( ek2 , 3 ) 
    + 3 * Ostap::Math::POW ( ek2 , 2 ) - 6 ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::GenGaussV2::tag () const 
{ 
  static const std::string s_name = "GenGaussV2" ;
  return Ostap::Utils::hash_combiner ( s_name , m_xi , m_alpha , m_kappa ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor from all agruments 
 *  @param xi     location/peak posiiton 
 *  @param omega  "scale" parameter 
 *  @param alpha  "shape" parameter 
 */
// ============================================================================
Ostap::Math::SkewGauss::SkewGauss
( const double xi    ,
  const double omega , 
  const double alpha )  // alpha=0 correponds to gaussian  
  : m_xi      ( xi                 ) 
  , m_omega   ( std::abs ( omega ) ) 
  , m_alpha   (            alpha   ) 
{}
// ============================================================================
bool Ostap::Math::SkewGauss::setXi        ( const double value ) 
{
  if ( s_equal ( value , m_xi) ) { return false ; }
  m_xi = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SkewGauss::setOmega     ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_omega ) ) { return false ; }
  m_omega = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SkewGauss::setAlpha  ( const double value ) 
{
  //
  if ( s_equal ( value , m_alpha ) ) { return false ; }
  m_alpha = value ;
  if ( s_equal ( 0     , m_alpha ) ) { m_alpha = 0  ; }
  //
  return true ;
}
// ============================================================================
double Ostap::Math::SkewGauss::pdf ( const double x ) const 
{
  const double y = ( x - m_xi ) / m_omega ;
  return 2* gauss_pdf ( y ) * gauss_cdf ( m_alpha * y ) / m_omega ;
}
// ============================================================================
double Ostap::Math::SkewGauss::cdf ( const double x ) const 
{
  const double y = ( x - m_xi ) / m_omega ;
  return gauss_cdf ( y ) - 2 * owen ( y , m_alpha ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::integral ( const double low  , 
                                          const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::mean () const 
{
  static const double s_c = std::sqrt ( 2.0L / M_PI ) ;
  const double delta = m_alpha / std::sqrt ( 1 + m_alpha * m_alpha ) ;
  return m_xi + m_omega * delta * s_c ;
}
// ============================================================================
double Ostap::Math::SkewGauss::variance () const 
{
  const double delta = m_alpha / std::sqrt( 1 + m_alpha * m_alpha ) ;
  return m_omega * m_omega * ( 1 - 2 * delta * delta / M_PI ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::skewness () const 
{
  static const double s_c1 =  ( 4.0- M_PI ) / 2.0 ;
  static const double s_c2 =  std::sqrt ( 2.0 / M_PI ) ;
  //
  const double delta = m_alpha / std::sqrt( 1 + m_alpha * m_alpha ) ;
  return s_c1 * std::pow ( delta * s_c2 , 3 ) /
    std::pow ( 1 - 2 * delta * delta / M_PI , 1.5  ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::kurtosis () const 
{
  static const double s_c1 =  2.0 * ( M_PI - 3 ) ;
  static const double s_c2 =  std::sqrt ( 2.0 / M_PI ) ;
  //
  const double delta = m_alpha / std::sqrt ( 1 + m_alpha * m_alpha ) ;
  return s_c1 * std::pow ( delta * s_c2 , 4 ) / std::pow ( 1 - 2 * delta * delta / M_PI , 2 ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::sigma  () const 
{ return std::sqrt ( variance () ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::SkewGauss::tag () const 
{ 
  static const std::string s_name = "SkewGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_xi , m_omega , m_alpha ) ; 
}
// ============================================================================

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Math::ExGauss::ExGauss
( const double mu       , 
  const double varsigma , 
  const double k        ) 
  : m_mu ( mu ) 
  , m_varsigma ( std::abs ( varsigma ) )
  , m_k  ( 0  )
  , m_mk ( 0  )
{
  setK ( k ) ;
}
// ============================================================================
double Ostap::Math::ExGauss::evaluate          ( const double x ) const 
{
  //
  const double z     = ( x - m_mu ) / m_varsigma ;
  const bool k_zero  = s_zero   ( m_k ) ;
  //
  const double kk    = std::abs ( m_k ) ;
  //
  return 
    k_zero  ? Ostap::Math::gauss_pdf   ( z ) / m_varsigma : 
    m_k > 0 ? Ostap::Math::gauss_mills ( z , 1.0/kk - z ) / ( kk * m_varsigma ) :
    m_k < 0 ? Ostap::Math::gauss_mills ( z , 1.0/kk + z ) / ( kk * m_varsigma ) :
    Ostap::Math::gauss_pdf   ( z ) / m_varsigma  ;
}
// ============================================================================
bool Ostap::Math::ExGauss::setMu        ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::ExGauss::setVarsigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_varsigma , avalue ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::ExGauss::setK ( const double value ) 
{
  if ( s_equal ( m_k , value ) ) { return false ; }
  //
  if ( s_zero ( value ) ) 
  {
    m_k  = 0 ;
    m_mk = 0 ;
  }
  else 
  {
    m_k = value ;
    if    ( std::abs ( value ) < 1.e-4 ) { m_mk = m_k ; }
    else 
    {
      //
      const double kk = 1.0 / m_k ;
      static const double s_C2 = std::sqrt ( 2.0 / M_PI ) ;
      const double aa = s_SQRT2 * Ostap::Math::erfcxinv ( s_C2 / std::abs ( kk ) ) ;
      //
      m_mk = 0 < m_k ? kk - aa : kk + aa ;
    }
  } 
  return true ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::ExGauss::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::ExGauss::integral   
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low )        { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::ExGauss::cdf ( const double x ) const 
{
  const double z     = ( x - m_mu ) / m_varsigma ;
  const bool k_zero  = s_zero ( m_k ) ;
  //
  const double gauss = Ostap::Math::gauss_cdf ( z ) ;
  const double kk    = std::abs ( m_k ) ;
  //
  return 
    k_zero  ? gauss  :
    m_k > 0 ? gauss - Ostap::Math::gauss_mills ( z , 1.0 / kk - z ) :
    m_k < 0 ? gauss + Ostap::Math::gauss_mills ( z , 1.0 / kk + z ) :
    gauss ;
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::ExGauss::mean        () const 
{ return m_mu + m_k * m_varsigma ; }
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::ExGauss::variance    () const 
{ return m_varsigma * m_varsigma * ( 1 + m_k * m_k ) ; }
// ============================================================================
// RMS value 
// ============================================================================
double Ostap::Math::ExGauss::rms        () const
{ return std::sqrt ( variance () )  ; }
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::ExGauss::skewness    () const 
{ return cumulant ( 3 ) / std::pow ( cumulant ( 2 ) , 1.5 ) ; }
// ============================================================================
// (excess) kurtosis  
// ============================================================================
double Ostap::Math::ExGauss::kurtosis    () const 
{
  const double k4 = cumulant ( 4 ) ;
  const double k2 = cumulant ( 2 ) ;
  const double s2 = variance ()  ;
  //
  return ( k4 + 3 * k2 * k2 ) / ( s2 * s2 ) - 3 ;
}
// ============================================================================
// get cumulant 
// ============================================================================
double Ostap::Math::ExGauss::cumulant    ( const unsigned short r ) const 
{
  return 
    0 == r ? 0.0         :
    1 == r ? mean     () :
    2 == r ? variance () :
    s_zero ( m_k ) ? 0.0 :
    std::tgamma ( r ) * std::pow ( m_k * m_varsigma , r ) ;  
}
// ============================================================================
// get the mode
// ============================================================================
double Ostap::Math::ExGauss::mode () const 
{ return m_mu + m_varsigma * m_mk ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::ExGauss::tag () const 
{ 
  static const std::string s_name = "ExGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma , m_k  ) ; 
}
// ============================================================================

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Math::ExGauss2::ExGauss2
( const double mu       , // the mode 
  const double varsigma , 
  const double k        )
  : m_emg ( mu , varsigma , k ) 
{
  setMu ( mu  ) ;  
}
// ============================================================================
// set new mode 
// ============================================================================
bool Ostap::Math::ExGauss2::setMu ( const double value ) 
{ 
  return m_emg.setMu ( value - m_emg.delta () ) ;
}
// ============================================================================
// set new k
// ============================================================================
bool Ostap::Math::ExGauss2::setK ( const double value ) 
{
  const double m1      = m_emg.mode () ;
  const bool   changed = m_emg.setK   ( value ) ;
  if ( !changed ) { return false ; }
  //
  const double m2    = m_emg.mode () ;
  if ( !s_equal ( m1 , m2 ) ) { setMu ( m1    ) ; }
  //
  return true ;
}
// ============================================================================
// set new sigma 
// ============================================================================
bool Ostap::Math::ExGauss2::setVarsigma ( const double value ) 
{
  const double m1      = m_emg.mode () ;
  const bool   changed = m_emg.setVarsigma ( value ) ;
  if ( !changed ) { return false ; }
  //
  const double m2      = m_emg.mode () ;
  if ( !s_equal ( m1 , m2 ) ) { setMu ( m1 ) ; }
  //
  return true ;
}
// ============================================================================
// evaluate it 
// ============================================================================
double Ostap::Math::ExGauss2::evaluate ( const double x ) const 
{ return m_emg ( x ) ; }
// ============================================================================
// get the CDF
// ============================================================================
double Ostap::Math::ExGauss2::cdf 
( const double x ) const 
{ return m_emg.cdf ( x ) ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::ExGauss2::integral   () const { return 1 ; }
// ============================================================================
/// get the integral between low and high limits
// ============================================================================
double Ostap::Math::ExGauss2::integral
( const double low  ,
  const double high ) const 
{ return m_emg.integral ( low , high ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::ExGauss2::tag () const 
{ 
  static const std::string s_name = "ExGauss2" ;
  return Ostap::Utils::hash_combiner ( s_name , mu () , varsigma () , k () ) ; 
}
// ============================================================================

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Bukin2::Bukin2
( const double mu        ,
  const double varsigmaA , 
  const double varsigmaB , 
  const double kA        , 
  const double kB        , 
  const double phi       )
  : m_A   ( mu , varsigmaA , kA )
  , m_B   ( mu , varsigmaB , kB )
  , m_phi ( 0   )
  , m_fA  ( 0.5 ) 
  , m_fB  ( 0.5 ) 
{
  setPhi ( phi ) ;
}
// ============================================================================
// set new value for parameter mu 
// ============================================================================
bool Ostap::Math::Bukin2::setMu ( const double value )
{
  const bool changedA = m_A.setMu ( value ) ;
  const bool changedB = m_B.setMu ( value ) ;
  return changedA || changedB ;
}
// ============================================================================
// set new value for parameter phi 
// ============================================================================
bool Ostap::Math::Bukin2::setPhi ( const double value )
{
  if ( s_equal ( value  , m_phi  ) ) { return false ; }
  //
  m_phi = value ;
  //
  static const double pi4 = M_PI * 0.25 ;
  const        double s   = std::sin ( value + pi4 ) ;
  m_fA = s * s  ;
  m_fB = 1.0 - m_fA ;
  return true ;
}
// ============================================================================
// evaluate the function
// ============================================================================
double  Ostap::Math::Bukin2::evaluate ( const double x  ) const
{ return m_fA * m_A ( x ) + m_fB * m_B ( x  ) ; }
// ============================================================================
// get the mean value 
// ============================================================================
double  Ostap::Math::Bukin2::mean () const
{ return m_fA * m_A.mean() + m_fB * m_B.mean() ; }
// ============================================================================
// get the integral 
// ============================================================================
double  Ostap::Math::Bukin2::integral () const { return 1 ; }
// ============================================================================
// get the cdf  
// ============================================================================
double  Ostap::Math::Bukin2::cdf ( const double x  ) const
{ return m_fA * m_A.cdf ( x ) + m_fB * m_B.cdf ( x ) ; }
// ============================================================================
// get the  integral 
// ============================================================================
double  Ostap::Math::Bukin2::integral
( const double low  ,
  const double high ) const
{
  return
    m_fA * m_A.integral ( low , high ) +
    m_fB * m_B.integral ( low , high ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Bukin2::tag () const 
{ 
  static const std::string s_name = "Bukin2" ;
  return Ostap::Utils::hash_combiner ( s_name    , mu () , m_phi ,  
                                       varsigmaA() , kA () ,  
                                       varsigmaB() , kB () ) ;
}
// ============================================================================

// ============================================================================
// constructor 
// ============================================================================
Ostap::Math::NormalLaplace::NormalLaplace 
( const double mu       ,
  const double varsigma ,
  const double kL       , 
  const double kR       ) 
  : m_mu ( mu ) 
  , m_varsigma ( std::abs ( varsigma ) ) 
  , m_kL       ( std::abs ( kL       ) ) 
  , m_kR       ( std::abs ( kR       ) ) 
{}
// ============================================================================
bool Ostap::Math::NormalLaplace::setMu        ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::NormalLaplace::setVarsigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_varsigma , avalue ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::NormalLaplace::setKL ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_kL , avalue ) ) { return false ; }
  m_kL = s_zero ( avalue ) ? 0.0 : avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::NormalLaplace::setKR ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_kR , avalue ) ) { return false ; }
  m_kR = s_zero ( avalue ) ? 0.0 : avalue ;
  return true ;
}
// ============================================================================
double Ostap::Math::NormalLaplace::evaluate ( const double x ) const 
{
  //
  const double z = ( x - m_mu ) / m_varsigma ;
  const bool l_zero = s_zero ( m_kL ) ;
  const bool r_zero = s_zero ( m_kR ) ;
  //
  return 
    l_zero && r_zero ? Ostap::Math::gauss_pdf   ( z ) / m_varsigma :
    l_zero           ? Ostap::Math::gauss_mills ( z , 1/m_kR - z ) / ( m_kR * m_varsigma ) :
    r_zero           ? Ostap::Math::gauss_mills ( z , 1/m_kL + z ) / ( m_kL * m_varsigma ) :
    ( Ostap::Math::gauss_mills ( z , 1/m_kR - z ) +
      Ostap::Math::gauss_mills ( z , 1/m_kL + z ) ) / ( ( m_kL + m_kR ) * m_varsigma ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::NormalLaplace::integral   () const { return 1.0 ; } 
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::NormalLaplace::integral  
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low )        { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::NormalLaplace::cdf ( const double x ) const 
{
  const double z = ( x - m_mu ) / m_varsigma ;
  const bool l_zero = s_zero ( m_kL ) ;
  const bool r_zero = s_zero ( m_kR ) ;
  //
  const double gauss = Ostap::Math::gauss_cdf ( z ) ;
  //
  return 
    l_zero && r_zero ? gauss : 
    l_zero           ? gauss - Ostap::Math::gauss_mills ( z , 1 / m_kR - z ) :
    r_zero           ? gauss + Ostap::Math::gauss_mills ( z , 1 / m_kL + z ) :
    gauss  - ( Ostap::Math::gauss_mills ( z , 1 / m_kR - z ) * m_kR - 
               Ostap::Math::gauss_mills ( z , 1 / m_kL + z ) * m_kL ) / ( m_kL + m_kR ) ;
}
// ============================================================================
// get cumulant 
// ============================================================================
double Ostap::Math::NormalLaplace::cumulant    ( const unsigned short r ) const 
{
  return 
    0 == r ? 0.0         :
    1 == r ? mean     () :
    2 == r ? variance () :
    std::tgamma ( r ) * ( std::pow ( m_kR * m_varsigma , r )  +
                          std::pow ( m_kL * m_varsigma , r ) ) ; 
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::NormalLaplace::mean        () const 
{  return m_mu + m_varsigma * ( m_kR - m_kL ) ; }
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::NormalLaplace::variance    () const 
{  return m_varsigma * m_varsigma * ( 1 + m_kR * m_kR  + m_kL * m_kL ) ; }
// ============================================================================
// RMS
// ============================================================================
double Ostap::Math::NormalLaplace::rms        () const 
{ return std::sqrt ( variance () )  ; }
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::NormalLaplace::skewness    () const 
{ return cumulant ( 3 ) / std::pow ( cumulant ( 2 ) , 1.5 ) ; }
// ============================================================================
// (excess) kurtosis  
// ============================================================================
double Ostap::Math::NormalLaplace::kurtosis    () const 
{
  const double k4 = cumulant ( 4 ) ;
  const double k2 = cumulant ( 2 ) ;
  const double s2 = variance ()  ;
  //
  return ( k4 + 3 * k2 * k2 ) / ( s2 * s2 ) - 3 ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::NormalLaplace::tag () const 
{ 
  static const std::string s_name = "NormalLaplace" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma , m_kL , m_kR ) ; 
}
// ============================================================================

// ============================================================================
// Bukin
// ============================================================================
/*  constructor from all parameters
 *  @param peak  the peak position
 *  @param sigma the effective sigma, defined as FWHM/2.35
 *  @param xi    the asymmetry parameter
 *  @param rhoL  the left  tail paramter
 *  @param rhoR  the right tail paramter
 */
// ============================================================================
Ostap::Math::Bukin::Bukin
( const double peak   ,
  const double sigma  ,
  const double xi     ,
  const double rho_L  ,
  const double rho_R  )
//
  : m_peak      ( M_PI + peak  )
  , m_sigma     ( M_PI + sigma )
  , m_xi        ( M_PI + xi    )
  , m_rho_L     ( M_PI + rho_L )
  , m_rho_R     ( M_PI + rho_R )
//
  , m_x1        ( M_PI )
  , m_x2        ( M_PI )
//
  , m_workspace ()
{
  //
  setXi    ( xi    ) ; // must be the first
  //
  setPeak  ( peak  ) ;
  //
  setSigma ( sigma ) ;
  //
  setRhoL  ( rho_L ) ;
  //
  setRhoR  ( rho_R ) ;
  //
}
// ============================================================================
bool Ostap::Math::Bukin::setPeak  ( const double value )
{
  if ( s_equal ( value , m_peak ) ) { return false ; }
  //
  m_peak   = value ;
  //
  const double xi_ = m_xi / std::sqrt ( 1 + m_xi * m_xi ) ;
  m_x1     = m_peak + m_sigma * s_Bukin * ( xi_ - 1 ) ;
  m_x2     = m_peak + m_sigma * s_Bukin * ( xi_ + 1 ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setSigma ( const double value )
{
  //
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Bukin::setSigma"           ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //        
  m_sigma  = value_ ;
  //
  const double xi_ = m_xi/std::sqrt ( 1 + m_xi * m_xi ) ;
  m_x1 = m_peak + m_sigma * s_Bukin * ( xi_ - 1 ) ;
  m_x2 = m_peak + m_sigma * s_Bukin * ( xi_ + 1 ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setXi    ( const double value )
{
  // no need for update
  if ( s_equal ( value , m_xi ) ) { return false ; } // no need for update
  //
  m_xi     = value ;
  //
  const double xi      = m_xi    ;
  const double xi2     =   xi*xi ;
  const double xi2sqrt = std::sqrt ( 1 + xi2 ) ;
  //
  const double alpha = 2 * xi  * xi2sqrt / s_Bukin   ;
  const double beta  = 2 * xi  * ( xi - xi2sqrt )    ;
  // well, it is actually alpha/beta:
  const double ab    = xi2sqrt / ( xi - xi2sqrt ) / s_Bukin ;
  //
  m_A     = alpha             ;
  //
  m_B2    = 1/Ostap::Math::log1p_x ( beta )  ;
  m_B2   *= m_B2              ;
  m_B2   *= ab*ab             ;
  //
  //
  const double delta = xi + xi2sqrt - 1 ;
  const double tail  =
    0.5 * s_Bukin * xi2sqrt * ( 1 + xi + xi2sqrt ) / ( xi + xi2sqrt ) / Ostap::Math::log1p_x  ( delta ) ;
  //
  // left  tail parameter
  //
  m_L  = tail ;
  m_L /= ( xi2sqrt - xi ) ;
  m_L /= ( xi2sqrt - xi ) ;
  //
  // right tail parameter
  //
  m_R  = tail ;
  m_R /= ( xi2sqrt + xi ) ;
  m_R /= ( xi2sqrt + xi ) ;
  //
  // central region
  //
  const double xi_ = m_xi / xi2sqrt ;
  m_x1 = m_peak + m_sigma * s_Bukin * ( xi_ - 1 ) ;
  m_x2 = m_peak + m_sigma * s_Bukin * ( xi_ + 1 ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setRhoL  ( const double value )
{
  //
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_rho_L  ) ) { return false ; }
  //
  m_rho_L    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setRhoR  ( const double value )
{
  //
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( value_ , m_rho_R  ) ) { return false ; }
  //
  m_rho_R    = value_ ;
  //
  return true ;
}
// ============================================================================
// evaluate Bukin's function
// ============================================================================
double Ostap::Math::Bukin::pdf ( const double x ) const
{
  //
  //  left tail :
  //
  if       ( m_x1 >= x )
  {
    const double dx  = x - m_x1               ;
    const double dx2 = dx / ( m_peak - m_x1 ) ;
    return  0.5 * my_exp (   m_L * dx / m_sigma  - m_rho_L * m_rho_L * dx2 * dx2 ) ;
  }
  //
  // right tail :
  //
  else if ( m_x2 <=  x )
  {
    const double dx  = x - m_x2               ;
    const double dx2 = dx / ( m_peak - m_x2 ) ;
    return 0.5 * my_exp ( - m_R * dx / m_sigma  - m_rho_R * m_rho_R * dx2 * dx2 ) ;
  }
  //
  // central region
  //
  const double dx   = ( x - m_peak ) / m_sigma ;
  const double A    = Ostap::Math::log1p_x ( m_A * dx ) ;
  //
  return my_exp ( - s_ln2 * dx * dx * A * A * m_B2 ) ;
}
// =========================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::Bukin::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0        ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low  ) ; } // RETURN
  //
  // split into reasonable sub-intervals
  //
  if ( low < m_x1    && m_x1   < high )
  { return integral (  low , m_x1   ) + integral ( m_x1   , high ) ; }
  if ( low < m_x2    && m_x2   < high )
  { return integral (  low , m_x2   ) + integral ( m_x2   , high ) ; }
  if ( low < m_peak  && m_peak < high )
  { return integral (  low , m_peak ) + integral ( m_peak , high ) ; }
  //
  const bool in_tail = 
    ( high < m_x1 - 5 * std::abs ( m_x2 - m_x1 ) )  || 
    ( low  > m_x2 + 5 * std::abs ( m_x2 - m_x1 ) ) ; 
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Bukin> s_integrator {} ;
  static char s_message[] = "Integral(Bukin)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     ,  
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size () ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::Bukin::integral () const
{
  //
  // Tails
  //
  static const Ostap::Math::GSL::Integrator1D<Bukin> s_integrator {} ;
  static char s_message1[] = "Integral(Bukin/left)"  ;
  static char s_message2[] = "Integral(Bukin/right)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  //
  int    ierror1  =  0 ;
  double result1  =  1 ;
  double error1   = -1 ;
  std::tie ( ierror1 , result1 , error1 ) = s_integrator.qagil_integrate
    ( tag () , 
      &F     , 
      m_x1   ,                       // low edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION_TAIL    ,          // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message1          , 
      __FILE__ , __LINE__ ) ;
  //
  int    ierror2  =  0 ;
  double result2  =  1 ;
  double error2   = -1 ;
  std::tie ( ierror2 , result2 , error2 ) = s_integrator.qagiu_integrate
    ( tag () , 
      &F     , 
      m_x2   ,                       // high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION_TAIL    ,          // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message2          , 
      __FILE__ , __LINE__ ) ;
  //
  return result1 + result2 + integral ( m_x1 , m_x2 ) ;
  //
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Bukin::tag () const 
{ 
  static const std::string s_name = "Bukin" ;
  return Ostap::Utils::hash_combiner ( s_name , m_peak , m_sigma , m_xi , m_rho_L , m_rho_R ) ; 
}
// ============================================================================

// ============================================================================
// Novosibirsk function
// ============================================================================
/*  constructor from all parameters
 *  @param m0    the peak posiion
 *  @param sigma the effective sigma
 *  @param tau   the tail paramter
 */
// ============================================================================
Ostap::Math::Novosibirsk::Novosibirsk
( const double m0    ,
  const double sigma ,
  const double tau   )
  : m_m0        ( m0                   )
  , m_sigma     ( std::fabs ( sigma )  )
  , m_tau       (             tau      )
    //
  , m_workspace ()
{
  //
  m_lambda = x_sinh ( m_tau * s_Novosibirsk ) ;
  //
}
// ============================================================================
// set parameter m0
// ============================================================================
bool Ostap::Math::Novosibirsk::setM0    ( const double value )
{
  if ( s_equal ( m_m0 ,  value ) ) { return false ; }
  m_m0 = value ;
  return true ;
}
// ============================================================================
// set parameter sigma
// ============================================================================
bool Ostap::Math::Novosibirsk::setSigma ( const double value )
{
  const double avalue = std::fabs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( avalue                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Novosibirsk::setSigma"     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //        
  m_sigma    = value ;
  return true ;
}
// ============================================================================
// set parameter tau
// ============================================================================
bool Ostap::Math::Novosibirsk::setTau ( const double value )
{
  if ( s_equal ( value , m_tau ) ) { return false ; }
  m_tau      = value ;
  m_lambda   = x_sinh ( m_tau * s_Novosibirsk ) ;
  return true ;
}
// ============================================================================
// evaluate Novosibirsk's function
// ============================================================================
double Ostap::Math::Novosibirsk::pdf  ( const double x ) const
{
  const double dx     = ( x - m_m0 ) / m_sigma ;
  const double arg    = m_lambda * dx * m_tau ;
  if ( arg <= -1 || s_equal ( arg , -1 ) ) { return 0 ; } // RETURN
  const double l      = Ostap::Math::log1p_x  ( arg ) * m_lambda * dx ;
  const double result = l * l ; // + m_tau * m_tau ;
  //
  return  my_exp ( -0.5 * result ) * s_SQRT2PIi / m_sigma ;
}
// =========================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::Novosibirsk::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                  0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low   ) ; } // RETURN
  //
  // split into reasonable sub intervals
  //
  if ( low <  m_m0  && m_m0 < high ) { return integral ( low , m_m0 ) + integral ( m_m0 , high ) ; }
  //
  {
    const double x1 = m_m0 +  3 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 -  3 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = m_m0 +  5 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 -  5 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }  
  //
  {
    const double x1 = m_m0 + 10 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 - 10 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = m_m0 + 15 * m_sigma ;
    if ( 0 < m_tau && low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 - 15 * m_sigma ;
    if ( 0 > m_tau && low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  const double x1     = m_m0 - 15 * m_sigma  ;
  const double x2     = m_m0 + 15 * m_sigma  ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Novosibirsk> s_integrator {} ;
  static char s_message[] = "Integral(Novosibirsk)" ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high             ,    // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,           // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::Novosibirsk::integral () const
{
  //
  if ( s_zero ( m_tau ) ) { return 1 ; }
  //
  const double tau1 = std::max ( 1.0 , std::abs ( m_tau ) ) ;
  const double tau2 = 1 ;
  //
  const double x_low  = m_m0 - ( 0 <= m_tau ?  5 * tau2 : 15 * tau1 ) * m_sigma ;
  const double x_high = m_m0 + ( 0 <= m_tau ? 15 * tau1 :  5 * tau2 ) * m_sigma ;
  //
  // use GSL to evaluate the tails:
  //
  static const Ostap::Math::GSL::Integrator1D<Novosibirsk> s_integrator {} ;
  static char s_message1[] = "Integral(Novosibirsk/left)"  ;
  static char s_message2[] = "Integral(Novosibirs/right)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  //
  int    ierror1  =  0 ;
  double result1  =  1 ;
  double error1   = -1 ;
  std::tie ( ierror1 , result1 , error1 ) = s_integrator.qagil_integrate
    ( tag () , 
      &F     , 
      x_low  ,                       // low edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION_TAIL   ,          // absolute precision
      s_APRECISION_TAIL   ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message1          , 
      __FILE__ , __LINE__ ) ;
  //
  int    ierror2  =  0 ;
  double result2  =  1 ;
  double error2   = -1 ;
  std::tie ( ierror2 , result2 , error2 ) = s_integrator.qagiu_integrate
    ( tag () , 
      &F     , 
      x_high ,                       // high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION_TAIL    ,          // absolute precision
      s_RPRECISION_TAIL    ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message2          , 
      __FILE__ , __LINE__ ) ;
  //
  return result1 + result2 + integral ( x_low ,  x_high ) ;
  //  
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Novosibirsk::tag () const 
{ 
  static const std::string s_name = "Novosibirsk" ;
  return Ostap::Utils::hash_combiner ( s_name , m_m0 , m_sigma , m_tau ) ; 
}
// ============================================================================

// ============================================================================
// Crystal Ball & Co
// ============================================================================

// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBall::CrystalBall
( const double m0     ,
  const double sigma  ,
  const double alpha  ,
  const double n      )
  : CrystalBall ( Ostap::Math::Gauss ( m0    , sigma ) ,
		  Ostap::Math::Tail  ( alpha , n     ) )
{}
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBall::CrystalBall
( const Ostap::Math::Gauss& core  ,  
  const double              alpha ,
  const double              n     )
  : CrystalBall ( core , Ostap::Math::Tail  ( alpha , n  ) )
{}
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBall::CrystalBall
( const Ostap::Math::Gauss& core  ,  
  const Ostap::Math::Tail&  tail  )
  : m_core   ( core )
  , m_tail   ( tail )
{}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBall::pdf ( const double x ) const
{
  //
  const double xl = xL () ; 
  
  // Gaussian ?
  if ( xl <= x ) { return m_core ( x ) ; }
  
  // Power-law tail
  
  // f( xL )
  const double F    = m_core      ( xl ) ; // f (xl)   
  
  // f'/f(xL) 
  const double dFoF = m_core.dFoF ( xl ) ; // f'/f    (xl) 
  //
  return m_tail ( x , xl , F , dFoF ) ;
}
// ============================================================================
/*  quantify the effect of tail, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over the function 
 */
// ============================================================================
double Ostap::Math::CrystalBall::non_gaussian
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB =        integral ( xlow , xhigh ) ;
  const double I_G  = m_core.integral ( xlow , xhigh ) ; 
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBall::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return   0; } // RETURN
  else if (           low > high   ) { return - integral ( high , low ) ; } // RETURN
  // 
  const double xl = xL() ; 
  //
  // split into the proper sub-intervals: peak-region & tail-region 
  if ( low < xl && xl < high ) { return integral ( low , xl ) + integral ( xl , high ) ; }
  //
  // Gaussian region ?
  if ( xl <= low ) { return m_core.integral ( low , high ) ; }
  //
  // function value at normalization point  x = xL 
  const double F     = m_core      ( xl ) ;
  // log-derivative at normalization point  x = xL
  const double dFoF  = m_core.dFoF ( xl ) ;
  //
  return m_tail.integral ( low  ,
			   high ,
			   xl   ,
			   F    ,
			   dFoF ) ;  
}
// ============================================================================
/*  get the integral from the the negative to positive infinity 
 *  @attention +infinity is returned for <code>n=0(N=1)</code>
 */
// ============================================================================
double Ostap::Math::CrystalBall::integral() const
{
  const double nn = N() ;
  if ( nn <= 1 || s_equal ( nn , 1.0 ) ) { return s_POSINF ; }
  //
  const double xl = xL () ;
  ///
  return m_tail.integral ( xl , xl , m_core ( xl ) , m_core.dFoF ( xl ) )
    + ( 1 - m_core.cdf ( xl ) ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBall::tag () const 
{
  static const std::string s_name = "CrystalBall" ;
  return Ostap::Utils::hash_combiner ( s_name , m_core.tag () , m_tail.tag () ) ; 
}
// ============================================================================

// ============================================================================
// Needham function
// ============================================================================
/* constructor from all parameters
 *  @param m0     m0       parameter
 *  @param sigma  sigma    parameter
 *  @param a0     c0       parameter
 *  @param a1     c1       parameter
 *  @param a2     c2       parameter
 */
// ============================================================================
Ostap::Math::Needham::Needham
( const double m0    ,
  const double sigma ,
  const double c0    ,
  const double c1    ,
  const double c2    ,
  const double n     , 
  const double amin  ) 
  : m_cb   ( m0 , sigma , 1 , 0 ) // Ostap::Math:CrystalBall
  , m_c0   ( -1 )
  , m_c1   ( -1 )
  , m_c2   ( -1 )
  , m_amin ( std::abs ( amin )  )
{
  setC      ( c0 , c1 , c2 ) ;
  m_cb.setN ( n            ) ;
  Ostap::Assert ( 0 < m_amin  && m_amin < 1 ,
		  "Parameter 'amin' must be 0<amin<1"    ,
		  "Ostap::Math::Needham" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
}
// ============================================================================
bool Ostap::Math::Needham::setSigma
( const double value )
{
  if ( !m_cb.setSigma ( value ) ) { return false ; }
  //
  return m_cb.setAlpha ( alpha ( sigma () ) ) ;
}
// ============================================================================
// set all three values togather
// ============================================================================
bool Ostap::Math::Needham::setC
( const double c0 ,
  const double c1 ,
  const double c2 )
{
  const double c0_ = std::abs ( c0  ) ;
  const double c1_ = std::abs ( c1  ) ;
  const double c2_ =            c2    ;
  //
  if ( s_equal ( c0_ , m_c0 ) &&
       s_equal ( c1_ , m_c1 ) &&
       s_equal ( c2_ , m_c2 ) ) { return false ; }
  //
  m_c0 = c0_ ;
  m_c1 = c1_ ;
  m_c2 = c2_ ;
  //
  return m_cb.setAlpha ( alpha ( sigma () ) ) ;
}
// ===========================================================================
// evaluate Needham's function
// ===========================================================================
double Ostap::Math::Needham::pdf ( const double x ) const
{ return m_cb ( x ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Needham::tag () const 
{
  static const std::string s_name = "Needham" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       m_cb.tag () ,
				       m_c0        ,
				       m_c1        ,
				       m_c2        ,
				       m_amin      ) ; 
}
// ============================================================================
// show alpha as function of sigma 
// ============================================================================
double Ostap::Math::Needham::alpha
( const double sigma ) const 
{ return Ostap::Math::needham_alpha ( sigma  , 
                                      m_c0   , 
                                      m_c1   , 
                                      m_c2   , 
                                      m_amin ) ; }
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBallRightSide::CrystalBallRightSide
( const double m0    ,
  const double sigma ,
  const double alpha ,
  const double n     )
  : m_core ( m0    , sigma )
  , m_tail ( alpha , n     ) 
{}
// ============================================================================
// constructor from gaussian and tail parameeter 
// ============================================================================
Ostap::Math::CrystalBallRightSide::CrystalBallRightSide
( const Ostap::Math::Gauss& core  , 
  const double              alpha , 
  const double              n     )
  : m_core ( core      )
  , m_tail ( alpha , n ) 
{}
// ============================================================================
// constructor from gaussian and tail parameeter 
// ============================================================================
Ostap::Math::CrystalBallRightSide::CrystalBallRightSide
( const Ostap::Math::Gauss& core  ,
  const Ostap::Math::Tail&  tail  )
  : m_core ( core )
  , m_tail ( tail )
{}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallRightSide::pdf ( const double x ) const
{
  //
  const double xr = xR () ; 
  
  // Gaussian ?
  if ( x <= xr ) { return m_core ( x ) ; }
  
  // Power-law tail
  
  // f( xr )
  const double F    = m_core ( xr ) ; // f    (xr)   
  
  // f'/f(xr) 
  const double dFoF = m_core.dFoF ( xr )     ; // f'/f (xr) 
  //
  return m_tail ( x , xr , F , dFoF ) ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallRightSide::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return   0; } // RETURN
  else if (           low > high   ) { return - integral ( high , low ) ; } // RETURN
  // 
  const double xr = xR() ; 
  //
  // split into the proper sub-intervals: peak-region & tail-region 
  if ( low < xr && xr < high ) { return integral ( low , xr ) + integral ( xr , high ) ; }
  //
  // Gaussian region ?
  if ( high <= xr ) { return m_core.integral ( low , high ) ; }
  //
  // function value at normalization point  x = xr 
  const double F     = m_core ( xr ) ;
  // log-derivative at normalization point  x = xr
  const double dFoF  = m_core.dFoF ( xr ) ;
  //
  return m_tail.integral ( low  ,
			   high ,
			   xr   ,
			   F    ,
			   dFoF ) ;
}
// ============================================================================
/*  get the integral from the the negative to positive infinity 
 *  @attention +infinity is returned for <code>n=0(N=1)</code>
 */
// ============================================================================
double Ostap::Math::CrystalBallRightSide::integral() const
{
  const double nn = N() ;
  if ( nn <= 1 || s_equal ( nn , 1.0 ) ) { return s_POSINF ; }
  //
  const double xr = xR () ;
  ///
  return m_tail.integral ( xr , xr , m_core ( xr ) , m_core.dFoF ( xr ) ) + m_core.cdf ( xr ) ;
}
// ============================================================================
  
// =========================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallRightSide::tag () const 
{  
  static const std::string s_name = "CrystalBallRightSide" ;
  return Ostap::Utils::hash_combiner ( s_name , m_core.tag() , m_tail.tag ()  ) ;
}
// ============================================================================
/** quantify the effect of tail, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 */
// ============================================================================
double Ostap::Math::CrystalBallRightSide::non_gaussian
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB =        integral ( xlow , xhigh ) ;
  const double I_G  = m_core.integral ( xlow , xhigh ) ; 
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================


// ============================================================================
/* constructor from all parameters
 *  @param m0     m0          parameter
 *  @param sigma  sigma       parameter
 *  @param alphaL alpha_L     parameter
 *  @param nL     n_L         parameter  (N-1 for "standard" definition)
 *  @param alphaR alpha_R parameter
 *  @param nR     n_R         parameter  (N-1 for "standard" definition)
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSided::CrystalBallDoubleSided
( const double m0      ,
  const double sigma   ,
  const double alphaL  ,
  const double nL      ,
  const double alphaR  ,
  const double nR      )
  : m_core  ( m0     , sigma )
  , m_left  ( alphaL , nL    )
  , m_right ( alphaR , nR    )
{}
// ============================================================================
/* constructor from core and tail  parameters
 *  @param core core Gaussian 
 *  @param alphaL  alpha_L     parameter
 *  @param nL      n_L         parameter  (N-1 for "standard" definition)
 *  @param alphaR  alpha_R parameter
 *  @param nR      n_R         parameter  (N-1 for "standard" definition)
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSided::CrystalBallDoubleSided
( const Ostap::Math::Gauss& core    , 
  const double              alphaL  ,
  const double              nL      ,
  const double              alphaR  ,
  const double              nR      )
  : m_core  ( core )
  , m_left  ( alphaL , nL    )
  , m_right ( alphaR , nR    )
{}
// ========================================================================
/*  constructor from all components  
 *  @param core core Gaussian 
 *  @param left  left  tail
 *  @param right right tail
 */
// ========================================================================
Ostap::Math::CrystalBallDoubleSided::CrystalBallDoubleSided
( const Ostap::Math::Gauss&     core  ,
  const Ostap::Math::LeftTail&  left  ,
  const Ostap::Math::RightTail& right )
  : m_core  ( core  )
  , m_left  ( left  )
  , m_right ( right )
{}
// ========================================================================
/*  constructor from core & tail components  
 *  @param cb CrystalBall function (left tail) 
 *  @param right right tail
 */
// ========================================================================
Ostap::Math::CrystalBallDoubleSided::CrystalBallDoubleSided
( const Ostap::Math::CrystalBall& cb    , 
  const Ostap::Math::RightTail  & right )
  : CrystalBallDoubleSided ( cb.core()  , cb.tail_left() , right )
{}
// ========================================================================
/* constructor from core & tail components  
 *  @param cb CrystalBallRightSide function (right tail) 
 *  @param left  left  tail
 */
// ========================================================================
Ostap::Math::CrystalBallDoubleSided::CrystalBallDoubleSided
( const Ostap::Math::CrystalBallRightSide& cb   , 
  const Ostap::Math::LeftTail            & left )
  : CrystalBallDoubleSided ( cb.core()  , left , cb.tail_right () )
{}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::pdf ( const double x ) const
{
  //
  // Left tail ? 
  const double xl =  xL () ;
  if ( x < xl )
    {
      // f( xL )
      const double F    = m_core ( xl ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xl )     ; // f'/f    (xl) 
      return m_left ( x , xl , F , dFoF ) ;
    } 
  //
  // right tail ?
  const double xr =  xR () ;
  if ( xr < x )
    {      
      // f( xL )
      const double F    = m_core ( xr ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xr )     ; // f'/f    (xl) 
      return m_right ( x , xr , F , dFoF ) ;      
    }
  //
  // core gaussian 
  return m_core ( x ) ; 
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low ) ; }
  //
  // split at xL 
  const double xl = xL () ;
  if ( low < xl && xl < high ) { return integral ( low , xl ) + integral ( xl , high ) ; }
  
  // split at x 
  const double xr = xR () ;
  if ( low < xr && xr < high ) { return integral ( low , xr ) + integral ( xr , high ) ; }
    
  // left tail
  if      ( high <= xl )
    {
      // f(xl)
      const double F    = m_core ( xl ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xl )     ; // f'/f    (xl)
      //
      return m_left.integral ( low , high , xl , F , dFoF ) ;
    }
  // right tail 
  else if ( xr <= low  )
    {
      // f(xr)
      const double F    = m_core ( xr ) ; // f    (xr)         
      // f'/f(xr) 
      const double dFoF = m_core.dFoF ( xr )     ; // f'/f (xr)
      //
      return m_right.integral ( low , high , xr , F , dFoF ) ;
    }
  // core Gaussian 
  return m_core.integral ( low , high ) ;
}
// ============================================================================
/*  Get the integral from the the negative to positive infinity 
 *  @attention +infinity is returned for <code>n=0(N=1)</code>
 */
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::integral() const
{
  const double nl = NL () ;
  if ( nl <= 1 || s_equal ( nl , 1.0 ) ) { return s_POSINF ; }
  const double nr = NR () ;
  if ( nr <= 1 || s_equal ( nr , 1.0 ) ) { return s_POSINF ; }
  //
  const double xl = xL () ;
  const double xr = xR () ;
  ///
  return m_core.integral ( xl , xr ) 
    + m_left .integral ( xl , xl , m_core ( xl ) , m_core.dFoF ( xl ) ) 
    + m_right.integral ( xr , xr , m_core ( xr ) , m_core.dFoF ( xr ) ) ;
}
// ============================================================================
/* quantify the effect of tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 */
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::non_gaussian
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB =        integral ( xlow , xhigh ) ;
  const double I_G  = m_core.integral ( xlow , xhigh ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallDoubleSided::tag () const 
{ 
  static const std::string s_name = "CrystalBallDoubleSided" ;
  return Ostap::Utils::hash_combiner ( s_name         ,
				       m_core .tag () ,
				       m_left .tag () ,
				       m_right.tag () ) ;
}
// ============================================================================


// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBallA::CrystalBallA
( const double m0     ,
  const double sigmaL ,
  const double sigmaR ,
  const double alpha  ,
  const double n      )
  : m_core ( m0    , sigmaL , sigmaR )
  , m_tail ( alpha , n ) 
{}
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBallA::CrystalBallA
( const Ostap::Math::BifurcatedGauss& core  ,  
  const double                        alpha ,
  const double                        n     )
  : m_core ( core )
  , m_tail ( alpha , n ) 
{}
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBallA::CrystalBallA
( const Ostap::Math::BifurcatedGauss& core  ,  
  const Ostap::Math::Tail&            tail  )
  : m_core   ( core )
  , m_tail   ( tail )
{}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallA::pdf ( const double x ) const
{
  //
  const double xl = xL () ;
  
  // Gaussian ?
  if ( xl <= x ) { return m_core ( x ) ; }
  
  // Power-law tail
  
  // f( xL )
  const double F    = m_core      ( xl ) ; // f (xl)   
  
  // f'/f(xL) 
  const double dFoF = m_core.dFoF ( xl ) ; // f'/f    (xl) 
  //
  return m_tail ( x , xl , F , dFoF ) ;
}
// ============================================================================
/*  quantify the effect of tail, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over the function 
 */
// ============================================================================
double Ostap::Math::CrystalBallA::non_gaussian
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB =        integral ( xlow , xhigh ) ;
  const double I_G  = m_core.integral ( xlow , xhigh ) ; 
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallA::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return   0; } // RETURN
  else if (           low > high   ) { return - integral ( high , low ) ; } // RETURN
  // 
  const double xl = xL() ; 
  //
  // split into the proper sub-intervals: peak-region & tail-region 
  if ( low < xl && xl < high ) { return integral ( low , xl ) + integral ( xl , high ) ; }
  //
  // Gaussian region ?
  if ( xl <= low ) { return m_core.integral ( low , high ) ; }
  //
  // function value at normalization point  x = xL 
  const double F     = m_core      ( xl ) ;
  // log-derivative at normalization point  x = xL
  const double dFoF  = m_core.dFoF ( xl ) ;
  //
  return m_tail.integral ( low  ,
			   high ,
			   xl   ,
			   F    ,
			   dFoF ) ;  
}
// ============================================================================
/*  get the integral from the the negative to positive infinity 
 *  @attention +infinity is returned for <code>n=0(N=1)</code>
 */
// ============================================================================
double Ostap::Math::CrystalBallA::integral() const
{
  const double nn = N() ;
  if ( nn <= 1 || s_equal ( nn , 1.0 ) ) { return s_POSINF ; }
  //
  const double xl = xL () ;
  ///
  return m_tail.integral ( xl , xl , m_core ( xl ) , m_core.dFoF ( xl ) )
    + ( 1 - m_core.cdf ( xl ) ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallA::tag () const 
{
  static const std::string s_name = "CrystalBallA" ;
  return Ostap::Utils::hash_combiner ( s_name , m_core.tag () , m_tail.tag () ) ; 
}
// ============================================================================


// ============================================================================
/* constructor from all parameters
 *  @param m0     m0          parameter
 *  @param sigmaL left sigma       parameter
 *  @param sigmaR right sigma       parameter
 *  @param alphaL alpha_L     parameter
 *  @param nL     n_L         parameter  (N-1 for "standard" definition)
 *  @param alphaR alpha_R parameter
 *  @param nR     n_R         parameter  (N-1 for "standard" definition)
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSidedA::CrystalBallDoubleSidedA
( const double m0      ,
  const double sigmaL  ,
  const double sigmaR  ,
  const double alphaL  ,
  const double nL      ,
  const double alphaR  ,
  const double nR      )
  : m_core  ( m0     , sigmaL , sigmaR  )
  , m_left  ( alphaL , nL     )
  , m_right ( alphaR , nR     )
{}
// ============================================================================
/* constructor from core and tail  parameters
 *  @param core core Gaussian 
 *  @param alphaL  alpha_L     parameter
 *  @param nL      n_L         parameter  (N-1 for "standard" definition)
 *  @param alphaR  alpha_R parameter
 *  @param nR      n_R         parameter  (N-1 for "standard" definition)
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSidedA::CrystalBallDoubleSidedA
( const Ostap::Math::BifurcatedGauss& core    , 
  const double                        alphaL  ,
  const double                        nL      ,
  const double                        alphaR  ,
  const double                        nR      )
  : m_core  ( core )
  , m_left  ( alphaL , nL    )
  , m_right ( alphaR , nR    )
{}
// ========================================================================
/*  constructor from all components  
 *  @param core core Gaussian 
 *  @param left  left  tail
 *  @param right right tail
 */
// ========================================================================
Ostap::Math::CrystalBallDoubleSidedA::CrystalBallDoubleSidedA
( const Ostap::Math::BifurcatedGauss& core  ,
  const Ostap::Math::LeftTail&        left  ,
  const Ostap::Math::RightTail&       right )
  : m_core  ( core  )
  , m_left  ( left  )
  , m_right ( right )
{}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedA::pdf ( const double x ) const
{
  //
  // Left tail ? 
  const double xl =  xL () ;
  if ( x < xl )
    {
      // f( xL )
      const double F    = m_core ( xl ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xl )     ; // f'/f    (xl) 
      return m_left ( x , xl , F , dFoF ) ;
    } 
  //
  // right tail ?
  const double xr =  xR () ;
  if ( xr < x )
    {      
      // f( xL )
      const double F    = m_core ( xr ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xr )     ; // f'/f    (xl) 
      return m_right ( x , xr , F , dFoF ) ;      
    }
  //
  // core gaussian 
  return m_core ( x ) ; 
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedA::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low ) ; }
  //
  // split at xL 
  const double xl = xL () ;
  if ( low < xl && xl < high ) { return integral ( low , xl ) + integral ( xl , high ) ; }
  
  // split at x 
  const double xr = xR () ;
  if ( low < xr && xr < high ) { return integral ( low , xr ) + integral ( xr , high ) ; }
    
  // left tail
  if      ( high <= xl )
    {
      // f(xl)
      const double F    = m_core ( xl ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xl )     ; // f'/f    (xl)
      //
      return m_left.integral ( low , high , xl , F , dFoF ) ;
    }
  // right tail 
  else if ( xr <= low  )
    {
      // f(xr)
      const double F    = m_core ( xr ) ; // f    (xr)         
      // f'/f(xr) 
      const double dFoF = m_core.dFoF ( xr )     ; // f'/f (xr)
      //
      return m_right.integral ( low , high , xr , F , dFoF ) ;
    }
  // core Gaussian 
  return m_core.integral ( low , high ) ;
}
// ============================================================================
/*  Get the integral from the the negative to positive infinity 
 *  @attention +infinity is returned for <code>n=0(N=1)</code>
 */
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedA::integral() const
{
  const double nl = NL () ;
  if ( nl <= 1 || s_equal ( nl , 1.0 ) ) { return s_POSINF ; }
  const double nr = NR () ;
  if ( nr <= 1 || s_equal ( nr , 1.0 ) ) { return s_POSINF ; }
  //
  const double xl = xL () ;
  const double xr = xR () ;
  ///
  return m_core.integral ( xl , xr ) 
    + m_left .integral ( xl , xl , m_core ( xl ) , m_core.dFoF ( xl ) ) 
    + m_right.integral ( xr , xr , m_core ( xr ) , m_core.dFoF ( xr ) ) ;
}
// ============================================================================
/* quantify the effect of tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 */
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedA::non_gaussian
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB =        integral ( xlow , xhigh ) ;
  const double I_G  = m_core.integral ( xlow , xhigh ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallDoubleSidedA::tag () const 
{ 
  static const std::string s_name = "CrystalBallDoubleSidedA" ;
  return Ostap::Utils::hash_combiner ( s_name         ,
				       m_core .tag () ,
				       m_left .tag () ,
				       m_right.tag () ) ;
}
// ============================================================================


// ============================================================================
/* constructor from all parameters
 *  @param m0     m0          parameter
 *  @param sigmaL left sigma       parameter
 *  @param sigmaR right sigma       parameter
 *  @param alphaL alpha_L     parameter
 *  @param nL     n_L         parameter  (N-1 for "standard" definition)
 *  @param alphaR alpha_R parameter
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSidedE::CrystalBallDoubleSidedE
( const double m0      ,
  const double sigmaL  ,
  const double sigmaR  ,
  const double alphaL  ,
  const double nL      ,
  const double alphaR  )
  : m_core  ( m0     , sigmaL , sigmaR  )
  , m_left  ( alphaL , nL     )
  , m_right ( alphaR )
{}
// ============================================================================
/* constructor from core and tail  parameters
 *  @param core core Gaussian 
 *  @param alphaL  alpha_L     parameter
 *  @param nL      n_L         parameter  (N-1 for "standard" definition)
 *  @param alphaR  alpha_R parameter
 *  @param nR      n_R         parameter  (N-1 for "standard" definition)
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSidedE::CrystalBallDoubleSidedE
( const Ostap::Math::BifurcatedGauss& core    , 
  const double                        alphaL  ,
  const double                        nL      ,
  const double                        alphaR  ) 
  : m_core  ( core )
  , m_left  ( alphaL , nL    )
  , m_right ( alphaR )
{}
// ========================================================================
/*  constructor from all components  
 *  @param core core Gaussian 
 *  @param left  left  tail
 *  @param right right tail
 */
// ========================================================================
Ostap::Math::CrystalBallDoubleSidedE::CrystalBallDoubleSidedE
( const Ostap::Math::BifurcatedGauss& core  ,
  const Ostap::Math::LeftTail&        left  ,
  const Ostap::Math::RightExpTail&    right )
  : m_core  ( core  )
  , m_left  ( left  )
  , m_right ( right )
{}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedE::pdf ( const double x ) const
{
  //
  // Left tail ? 
  const double xl =  xL () ;
  if ( x < xl )
    {
      // f( xL )
      const double F    = m_core ( xl ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xl )     ; // f'/f    (xl) 
      return m_left ( x , xl , F , dFoF ) ;
    } 
  //
  // right tail ?
  const double xr =  xR () ;
  if ( xr < x )
    {      
      // f( xL )
      const double F    = m_core ( xr ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xr )     ; // f'/f    (xl) 
      return m_right ( x , xr , F , dFoF ) ;      
    }
  //
  // core gaussian 
  return m_core ( x ) ; 
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedE::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low ) ; }
  //
  // split at xL 
  const double xl = xL () ;
  if ( low < xl && xl < high ) { return integral ( low , xl ) + integral ( xl , high ) ; }
  
  // split at x 
  const double xr = xR () ;
  if ( low < xr && xr < high ) { return integral ( low , xr ) + integral ( xr , high ) ; }
    
  // left tail
  if      ( high <= xl )
    {
      // f(xl)
      const double F    = m_core ( xl ) ; // f (xl)         
      // f'/f(xL) 
      const double dFoF = m_core.dFoF ( xl )     ; // f'/f    (xl)
      //
      return m_left.integral ( low , high , xl , F , dFoF ) ;
    }
  // right tail 
  else if ( xr <= low  )
    {
      // f(xr)
      const double F    = m_core ( xr ) ; // f    (xr)         
      // f'/f(xr) 
      const double dFoF = m_core.dFoF ( xr )     ; // f'/f (xr)
      //
      return m_right.integral ( low , high , xr , F , dFoF ) ;
    }
  // core Gaussian 
  return m_core.integral ( low , high ) ;
}
// ============================================================================
/*  Get the integral from the the negative to positive infinity 
 *  @attention +infinity is returned for <code>n=0(N=1)</code>
 */
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedE::integral() const
{
  const double nl = NL () ;
  if ( nl <= 1 || s_equal ( nl , 1.0 ) ) { return s_POSINF ; }
  //
  const double xl = xL () ;
  const double xr = xR () ;
  ///
  return m_core.integral ( xl , xr ) 
    + m_left .integral ( xl , xl , m_core ( xl ) , m_core.dFoF ( xl ) ) 
    + m_right.integral ( xr , xr , m_core ( xr ) , m_core.dFoF ( xr ) ) ;
}
// ============================================================================
/* quantify the effect of tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 */
// ============================================================================
double Ostap::Math::CrystalBallDoubleSidedE::non_gaussian
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB =        integral ( xlow , xhigh ) ;
  const double I_G  = m_core.integral ( xlow , xhigh ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallDoubleSidedE::tag () const 
{ 
  static const std::string s_name = "CrystalBallDoubleSidedE" ;
  return Ostap::Utils::hash_combiner ( s_name         ,
				       m_core .tag () ,
				       m_left .tag () ,
				       m_right.tag () ) ;
}
// ============================================================================


// ============================================================================
// Apollonios
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 *  @param b     b-parameter 
 */
// ============================================================================
Ostap::Math::Apollonios::Apollonios
( const double m0     ,
  const double sigmaL ,
  const double sigmaR ,
  const double beta   )
  : m_m0         (  0 )
  , m_sigmaL     (  1 )
  , m_sigmaR     (  1 )
  , m_beta       (  1 )
  , m_workspace  () 
{
  //
  setM0     ( m0     ) ;
  setSigmaL ( sigmaL ) ;
  setSigmaR ( sigmaR ) ;
  setBeta   ( beta   ) ;
  //
}
// ============================================================================
bool  Ostap::Math::Apollonios::setM0 ( const double value )
{
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  m_m0       = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setSigmaL ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigmaL ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigmaL' must be non-zero"    ,
		  "Ostap::Math::Aplollonious::setSigmaL"   ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //        
  m_sigmaL    = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setSigmaR ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigmaR ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigmaR' must be non-zero"    ,
		  "Ostap::Math::Aplollonious::setSigmaR"   ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //          
  m_sigmaR    = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setSigma
( const double valueL , 
  const double valueR ) 
{
  const bool m1 = setSigmaL ( valueL ) ;
  const bool m2 = setSigmaR ( valueR ) ;
  return m1 || m2 ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setBeta ( const double value )
{
  //
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_beta ) ) { return false ; }
  //
  m_beta    = value_  ;
  //
  if ( s_equal ( m_beta , 0 ) ) { m_beta = 0 ; }
  if ( s_equal ( m_beta , 1 ) ) { m_beta = 1 ; }
  //
  return true ;
}
// ============================================================================
/*  log-derivative \f$ \frac{ f^\prime}{f}  \f$
 *  Useful to attach the radiative tail to ensure 
 *  the continuity of the function and the 1st derivatibve 
 */
// ============================================================================
double Ostap::Math::Apollonios::dFoF ( const double x ) const
{
  const double dx = 
    ( x < m_m0 ) ?
    ( x - m_m0 ) / m_sigmaL :
    ( x - m_m0 ) / m_sigmaR ;
  //
  const double h2 = std::hypot ( s_SQRT2 , m_beta ) ;
  const double hx = std::hypot ( dx      , m_beta ) ;
  //
  return -h2 * dx / ( hx * ( dx < 0 ? m_sigmaL : m_sigmaR ) ) ;
}
// ============================================================================
//  evaluate Apollonios' function
// ============================================================================
double Ostap::Math::Apollonios::pdf ( const double x ) const
{
  //
  const double dx = 
    ( x < m_m0 ) ?
    ( x - m_m0 ) / m_sigmaL :
    ( x - m_m0 ) / m_sigmaR ;
  //
  // the peak
  //
  const double h2 = std::hypot ( s_SQRT2 , m_beta ) ;
  const double hx = std::hypot ( dx      , m_beta ) ;
  // 
  return std::exp ( h2 * ( m_beta - hx ) ) * s_SQRT2PIi / sigma()  ;  
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Apollonios::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  /// split at the maximum
  if ( low < m_m0 && m_m0 < high ) 
  { return integral ( low , m_m0 ) + integral ( m_m0 , high ) ; }
  //
  /// split into reasonable intervals of 2-sigma width
  const unsigned int N = 6 ;
  for ( unsigned int n = 2 ; n <= N ; n += 2 )
  {
    const double xr = m_m0 + n * m_sigmaR ;
    if ( low < xr && xr < high ) 
      { return integral ( low , xr ) + integral ( xr , high ) ; }
    //
    const double xl = m_m0 - n * m_sigmaL ;
    if ( low < xl && xl < high ) 
      { return integral ( low , xl ) + integral ( xl , high ) ; }      
    //
  }
  //
  const double xR = m_m0 + N * m_sigmaR ;
  const double xL = m_m0 - N * m_sigmaL ;
  //
  const double in_tail = ( low >= xR || high <= xL ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Apollonios> s_integrator {} ;
  static char s_message[] = "Integral(Apollonios)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    ,                                     // low integration edge
      high   ,                                     // high integration edge
      workspace ( m_workspace )                  , // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()                         , // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Apollonios::tag () const 
{ 
  static const std::string s_name = "Apollonios" ;
  return Ostap::Utils::hash_combiner ( s_name , m_m0 , m_sigmaL , m_sigmaR , m_beta ) ; 
}
// ============================================================================

// ============================================================================
// ApolloniosL 
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 *  @param b     b-parameter 
 */
// ============================================================================
Ostap::Math::ApolloniosL::ApolloniosL
( const double m0     ,
  const double sigmaL ,
  const double sigmaR ,
  const double beta   , 
  const double alpha  ,
  const double n      )
  : ApolloniosL ( Ostap::Math::Apollonios ( m0    , sigmaL , sigmaR , beta ) ,
		  Ostap::Math::LeftTail   ( alpha , n ) )
{}
// ============================================================================
/*  constructor from two parameters
 *  @param core core component 
 *  @param tail left-tail component 
 */
// ============================================================================
Ostap::Math::ApolloniosL::ApolloniosL
( const Ostap::Math::Apollonios& core ,
  const Ostap::Math::Tail&       tail )
  : m_core ( core )
  , m_tail ( tail )
{}
// ============================================================================
//  evaluate Apollonios' function
// ============================================================================
double Ostap::Math::ApolloniosL::pdf ( const double x ) const
{
  const double xl = xL () ;
  //
  /// Apollonious core?
  if ( xl <= x ) { return m_core ( x ) ; }
  //
  /// Power-law tail
  const double F    = m_core      ( xl ) ;
  const double dFoF = m_core.dFoF ( xl ) ;
  //
  return m_tail ( x , xl , F , dFoF ) ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::ApolloniosL::integral
( const double low ,
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  // split at the meeting point 
  const double xl = xL () ;
  if ( low < xl && xl < high )
  { return integral ( low , xl ) + integral ( xl  , high ) ; }
  //
  /// core Apollonious region ?
  if ( xl <= low ) { return m_core.integral ( low , high ) ; }
  //
  /// Power-law tail
  const double F    = m_core      ( xl ) ;
  const double dFoF = m_core.dFoF ( xl ) ;
  //
  return m_tail.integral ( low  ,
			   high ,
			   xl   ,
			   F    ,
			   dFoF ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::ApolloniosL::tag () const 
{ 
  static const std::string s_name = "ApolloniosL" ;
  return Ostap::Utils::hash_combiner ( s_name , m_core.tag() , m_tail.tag() ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param sigma \f$\sigma\f$-parameter 
 */
// ============================================================================
Ostap::Math::Atlas::Atlas
( const double mean  ,
  const double sigma ) 
  : m_mean      ( mean               ) 
  , m_sigma     ( std::abs ( sigma ) )    
  , m_workspace ()
{}
// ============================================================================
// get variance:  very good numerical approximation d
// ============================================================================
double Ostap::Math::Atlas::variance () const { return 3 * m_sigma * m_sigma ; }
// ============================================================================
// get rms :  very good numerical approximation 
// ============================================================================
double Ostap::Math::Atlas::rms      () const { return s_SQRT3     * m_sigma ; }
// ============================================================================
bool Ostap::Math::Atlas::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean ) ) { return false ; }
  m_mean = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Atlas::setSigma ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Atlas::setSigma"           ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //          
  m_sigma = value_ ;
  return true ;
}
// ============================================================================
// evaluate atlas function 
// ============================================================================
double Ostap::Math::Atlas::pdf        ( const double x ) const 
{
  const double dx = std::abs  ( x - m_mean ) / m_sigma ;
  if ( s_zero ( dx ) ) { return 1 ; }                        // return 1 
  const double x2 = std::pow ( dx , 1.0 + 1 / ( 1 + 0.5 * dx ) ) ;
  return std::exp ( -0.5 * x2 ) / ( s_ATLAS * m_sigma )  ;
}
// ============================================================================
double Ostap::Math::Atlas::integral
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low ,high ) ) { return 0 ; }
  else if ( low > high            ) { return -integral ( high , low ) ; }
  //
  // split 
  if ( low < m_mean && m_mean < high ) 
  { return integral ( low , m_mean ) + integral ( m_mean , high ) ; }
  //
  const double left  = m_mean - 5 * m_sigma ;  
  if ( low < left   &&  left < high ) 
  { return integral ( low , left   ) + integral ( left   , high ) ; }
  //
  const double right = m_mean + 5 * m_sigma ;  
  if ( low < right  && right  < high ) 
  { return integral ( low , right  ) + integral ( right  , high ) ; }
  //
  const bool in_tail = ( high <= left || low >= right ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Atlas> s_integrator {} ;
  static char s_message[] = "Integral(Atlas)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// overall integral, not exact but precise enough...
// ============================================================================
double Ostap::Math::Atlas::integral () const { return 1 ; }
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::Atlas::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Atlas::tag () const 
{ 
  static const std::string s_name = "Atlas" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_sigma ) ; 
}
// ============================================================================


// ============================================================================
// Sech 
// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param sigma \f$\sigma\f$-parameter
 */
// ============================================================================
Ostap::Math::Sech::Sech  
( const double mean  ,
  const double sigma ) 
  : m_mean  (             mean    ) 
  , m_sigma (  std::abs ( sigma ) )
{}
// ============================================================================
// evaluate sech function 
// ============================================================================
double Ostap::Math::Sech::pdf ( const double x ) const 
{
  const long double y = ( x - m_mean ) * M_PI_2 / m_sigma ;
  return 
    GSL_LOG_DBL_MAX < std::abs ( y )  ? 0 : 
    0.5 / ( m_sigma * std::cosh ( y ) ) ;
}
// ============================================================================
bool Ostap::Math::Sech::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean  ) ) { return false ; }
  m_mean  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Sech::setSigma ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma  ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Sech::setSigma"            ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //            
  m_sigma  = value_ ;
  return true ;
}
// ============================================================================
// get integral from low to high 
// ============================================================================
double Ostap::Math::Sech::integral 
( const double low  ,
  const double high ) const 
{ return s_equal ( low , high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// get integral from -infinity to + infinity 
// ============================================================================
double Ostap::Math::Sech::integral () const { return 1 ; } 
// ============================================================================
// evaluate CDF function 
// ============================================================================
double Ostap::Math::Sech::cdf ( const double x ) const 
{
  const long double y = ( x - m_mean ) * M_PI_2 / m_sigma ;
  return
    ( GSL_LOG_DBL_MAX < y ) ? 1 :
    ( GSL_LOG_DBL_MIN > y ) ? 0 : 
    std::atan (  std::exp ( y ) ) / M_PI_2 ;
}
// ============================================================================
// get quantile (0<p<1)
// ============================================================================
double Ostap::Math::Sech::quantile ( const double p ) const 
{ return 
    0 >= p || s_zero  ( p     ) ? -s_INFINITY :
    1 <= p || s_equal ( p , 1 ) ? +s_INFINITY : 
    m_mean + m_sigma * 2 / M_PI * std::log( std::tan ( M_PI * p /  2 ) ) ; }
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::Sech::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================ 
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Sech::tag () const 
{ 
  static const std::string s_name = "Sech" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_sigma ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param alpha \f$\alpha\f$-parameter 
 *  @param beta  \f$\beta\f$-parameter 
 */
// ============================================================================
Ostap::Math::Losev::Losev
( const double mu    , 
  const double alpha , 
  const double beta  ) 
  : m_mu        ( mu                 ) 
  , m_alpha     ( std::abs ( alpha ) ) 
  , m_beta      ( std::abs ( beta  ) ) 
  , m_norm      ( -1 )
  , m_workspace () 
{}
// ============================================================================
bool Ostap::Math::Losev::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu  = value ;
  return true ;
}
// =============================================================================
bool Ostap::Math::Losev::setAlpha ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_alpha  ) ) { return false ; }
  m_alpha  = v  ;
  m_norm   = -1 ;
  return true ;
}
// =============================================================================
bool Ostap::Math::Losev::setBeta ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_beta  ) ) { return false ; }
  m_beta   = v  ;
  m_norm   = -1 ;
  return true ;
}
// =============================================================================
// the mode of the distribution 
// =============================================================================
double Ostap::Math::Losev::mode () const 
{ return m_mu + std::log ( m_alpha / m_beta ) / ( m_alpha + m_beta ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Losev::tag () const 
{ 
  static const std::string s_name = "Losev" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_alpha , m_beta ) ; 
}
// =============================================================================
// evaluate the function 
// =============================================================================
double Ostap::Math::Losev::pdf ( const double x ) const 
{
  if ( m_norm <= 0 ) 
  {
    const double sumab = m_alpha + m_beta ;
    m_norm = sumab * std::sin ( M_PI * m_beta / sumab ) / M_PI ;
  }
  //
  const double dx = x - m_mu ;
  return 0 <= dx  ? 
    m_norm * std::exp ( -m_beta  * dx ) / ( 1 + std::exp  ( - ( m_alpha + m_beta ) * dx ) ) :
    m_norm * std::exp (  m_alpha * dx ) / ( 1 + std::exp  (   ( m_alpha + m_beta ) * dx ) ) ;  
}
// ============================================================================
/*  get the integral between low and high values 
 *  \f$ \int_{low}^{high}f(x) dx\f$
 */
// ============================================================================
double Ostap::Math::Losev::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low ,high ) ) { return 0 ; }
  else if ( low > high            ) { return -integral ( high , low ) ; }
  //  
  // split 
  const double left  = m_mu - 6 * m_alpha;  
  if ( low < left && left < high ) 
  { return integral ( low , left ) + integral ( left   , high ) ; }
  //
  const double right = m_mu + 6 * m_beta ;  
  if ( low < right && right < high ) 
  { return integral ( low , right ) + integral ( right , high ) ; }
  //
  const bool in_tail = ( high <= left || low >= right ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Losev> s_integrator {} ;
  static char s_message[] = "Integral(Losev)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// Logistic
// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param sigma \f$\sigma\f$-parameter
 */
// ============================================================================
Ostap::Math::Logistic::Logistic
( const double mean  ,
  const double sigma ) 
  : m_mean  (             mean    ) 
  , m_sigma (  std::abs ( sigma ) )
{}
// ============================================================================
// evaluate sech function 
// ============================================================================
double Ostap::Math::Logistic::pdf ( const double x ) const 
{
  const double s = m_sigma * s_SQRT3overPI ;
  const long double y = ( x - m_mean ) / ( 2 * s ) ;
  if  ( GSL_LOG_DBL_MAX < std::abs( y ) ) { return 0 ; }
  const long double c = std::cosh ( y ) ;
  return 0.25 / c / c / s  ;
}
// ============================================================================
bool Ostap::Math::Logistic::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean  ) ) { return false ; }
  m_mean  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Logistic::setSigma ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma  ) ) { return false ; }
  //  
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Logistic::setSigma"        ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //          
  m_sigma  = value_ ;
  return true ;
}
// ============================================================================
// get integral from low to high 
// ============================================================================
double Ostap::Math::Logistic::integral 
( const double low  ,
  const double high ) const 
{ return s_equal ( low , high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// get integral from -infinity to + infinity 
// ============================================================================
double Ostap::Math::Logistic::integral () const { return 1 ; } 
// ============================================================================
// evaluate CDF function 
// ============================================================================
double Ostap::Math::Logistic::cdf ( const double x ) const 
{
  const double s = m_sigma * s_SQRT3overPI ;
  const long double y = ( x - m_mean ) / ( 2 * s ) ;
  return 0.5 * ( 1 + std::tanh ( y ) ) ;
}
// ============================================================================
// get parameter s 
// ============================================================================
double Ostap::Math::Logistic::s() const 
{ return m_sigma * s_SQRT3overPI ; }
// ============================================================================
// quantile function  (0<p<1)
// ============================================================================
double Ostap::Math::Logistic::quantile ( const double p ) const 
{ return
    0 >= p || s_zero  ( p     ) ? -s_INFINITY :
    1 <= p || s_equal ( p , 1 ) ? +s_INFINITY : 
    m_mean + m_sigma * s_SQRT3overPI * std::log ( p / ( 1 - p ) ) ; }
// ===========================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::Logistic::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Logistic::tag () const 
{ 
  static const std::string s_name = "Logistic" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_sigma ) ;
}
// ============================================================================

// ============================================================================
/*  constructor with all parameters
 *  @param mu    \f$\mu  \f$-parameter
 *  @param sigma \f$sigma\f$-parameter
 *  @param alpha \f$alpha\f$-parameter
 *  @param beta  \f$beta\f$-parameter
 */
// ============================================================================
Ostap::Math::GenLogisticIV::GenLogisticIV
( const double mu    ,
  const double sigma , 
  const double alpha , 
  const double beta  ) 
  : m_mu    (            mu      ) 
  , m_sigma ( std::abs ( sigma ) ) 
  , m_alpha ( std::abs ( alpha ) ) 
  , m_beta  ( std::abs ( beta  ) ) 
  , m_tilda_mu (  0 ) 
  , m_tilda_s  (  1 ) 
  , m_norm     ( -1 )
{
  setAlpha ( alpha ) ;
  setBeta  ( beta  ) ;
}
// ============================================================================
// set parameter mu 
// ============================================================================
bool Ostap::Math::GenLogisticIV::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;  
}
// ============================================================================
// set parameter sigma
// ============================================================================
bool Ostap::Math::GenLogisticIV::setSigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  //  
  Ostap::Assert ( avalue                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::GenLogisticIV::setSigma"   ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //            
  m_sigma = avalue ;
  return true ;  
}
// ============================================================================
// set parameter alpha
// ============================================================================
bool Ostap::Math::GenLogisticIV::setAlpha ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_alpha ) && 0 < m_norm ) { return false ; }
  //
  m_alpha    = avalue ;
  //
  m_tilda_mu =             Ostap::Math::psi ( m_alpha     ) - Ostap::Math::psi ( m_beta     )   ;
  m_tilda_s  = std::sqrt ( Ostap::Math::psi ( m_alpha , 1 ) + Ostap::Math::psi ( m_beta , 1 ) ) ;
  m_norm     = 1.0/Ostap::Math::beta ( m_alpha , m_beta ) ;
  //
  return true ;  
}
// ============================================================================
// set parameter beta 
// ============================================================================
bool Ostap::Math::GenLogisticIV::setBeta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_beta ) && 0 < m_norm ) { return false ; }
  //
  m_beta     = avalue ;
  //
  m_tilda_mu =             Ostap::Math::psi ( m_alpha     ) - Ostap::Math::psi ( m_beta     )   ;
  m_tilda_s  = std::sqrt ( Ostap::Math::psi ( m_alpha , 1 ) + Ostap::Math::psi ( m_beta , 1 ) ) ;
  m_norm     = 1.0/Ostap::Math::beta ( m_alpha , m_beta ) ;
  //
  return true ;  
}
// ============================================================================
// get the "standard" Generalized Type IV Logistic distribution 
// ============================================================================
double  Ostap::Math::GenLogisticIV::std_type4 ( const double t ) const 
{
  //
  // const double tt = std::tanh ( 0.5 * t ) ;
  // const double s1 = 0.5 * ( 1 + tt ) ;
  // const double s2 = 0.5 * ( 1 - tt ) ;
  // return m_norm * std::pow ( s1 , m_alpha ) * std::pow ( s2 , m_beta ) ;
  //
  // better numeraicla properties :
  return 0 <= t ? 
    m_norm * std::exp ( - m_beta  * t ) / std::pow ( 1 + std::exp  ( -t ) , m_alpha + m_beta ) :
    m_norm * std::exp (   m_alpha * t ) / std::pow ( 1 + std::exp  (  t ) , m_alpha + m_beta ) ;  
}
// ============================================================================
// get the helper variable y 
// ============================================================================
double Ostap::Math::GenLogisticIV::y ( const double z ) const 
{ return m_mu      + m_sigma    * ( z - m_tilda_mu ) / m_tilda_s ; }
// ============================================================================
// get the helper variable z 
// ============================================================================
double Ostap::Math::GenLogisticIV::z ( const double y ) const 
{ return m_tilda_mu + m_tilda_s * ( y - m_mu       ) / m_sigma ; }
// ============================================================================
//  evaluate the function
// ============================================================================
double  Ostap::Math::GenLogisticIV::evaluate ( const double x ) const 
{
  const double r = m_tilda_s / m_sigma ;
  const double t = m_tilda_mu + r * ( x - m_mu ) ;
  return std_type4 ( t ) * r ;
}
// ============================================================================
// integral from -infinity to +infinity
// ============================================================================
double Ostap::Math::GenLogisticIV::integral () const { return 1 ; }
// ============================================================================
// get integral from low to high
// ============================================================================
double Ostap::Math::GenLogisticIV::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  // split into reasonable subintervals 
  //
  if ( low < m_mu && m_mu < high ) 
  { return integral ( low , m_mu ) + integral ( m_mu , high ) ; }
  //
  {
    const double x1 = m_mu +  3 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mu -  3 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }  
  //
  {
    const double x1 = m_mu +  6 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mu -  6 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }  
  //
  {
    const double x1 = m_mu + 10 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mu - 10 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }  
  //
  const double x1     = m_mu - 10 * m_sigma  ;
  const double x2     = m_mu + 10 * m_sigma  ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenLogisticIV> s_integrator {} ;
  static char s_message[] = "Integral(GenLogisticIV)" ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high             ,    // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,           // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ===========================================================================
// get the mode
// ============================================================================
double Ostap::Math::GenLogisticIV::mode() const 
{ return y ( std::log ( m_alpha / m_beta ) ) ; }
// ===========================================================================
// get the skewness
// ============================================================================
double Ostap::Math::GenLogisticIV::skewness () const 
{ return cumulant ( 3 ) / std::pow ( m_sigma , 3 ) ; }
// ===========================================================================
// get the (excess) kurtosis
// ============================================================================
double Ostap::Math::GenLogisticIV::kurtosis () const 
{
  const double mu4 = cumulant ( 4 ) + 3 * std::pow ( variance () , 2 ) ;
  return mu4 / std::pow ( m_sigma , 4 ) - 3 ; 
}
// ===========================================================================
// get the cumulant 
// ============================================================================
double Ostap::Math::GenLogisticIV::cumulant
( const unsigned short k ) const 
{ 
  return 
    0 == k ? 0.0         :
    1 == k ? mean     () :
    2 == k ? variance () :
    ( Ostap::Math::psi ( m_alpha , k - 1 )
      + ( 0 == k % 2 ? 1 : -1 ) * Ostap::Math::psi ( m_beta , k - 1 ) ) 
    * std::pow ( m_sigma / m_tilda_s , k ) ;
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::GenLogisticIV::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ===========================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::GenLogisticIV::tag () const 
{ 
  static const std::string s_name = "GenLogisticIV" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_sigma , m_alpha , m_beta ) ;
}
// ============================================================================

// ============================================================================
// Student-T 
// ============================================================================
/*  constructor from mass, resolution and "n"-parameter 
 *  @param M     mass 
 *  @param sigma width parameter
 *  @param n     n-parameter  ( actually  N=N(n) ) 
 */
// ============================================================================
Ostap::Math::StudentT::StudentT 
( const double mass  , 
  const double scale ,
  const double n     ) 
  //
  : m_M     (      std::abs ( mass  ) )
  , m_scale (      std::abs ( scale ) )
  , m_n     ( -1 )
  , m_nu    ( -1 )    
  , m_norm  ( -1 ) 
{
  setN ( n ) ;  
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setM ( const double x )
{
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_M ) ) { return false ; }
  m_M = v ;
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setScale ( const double value )
{
  const double avalue = std::abs ( value ) ;
  //
  Ostap::Assert ( avalue  ,
		  "Parameter 'scale/sigma' must be non-zero" ,
		  "Ostap::Math::StudentT::setScale"          ,
		  INVALID_PARAMETER , __FILE__ , __LINE__    ) ;
  //
  if ( s_equal ( avalue , m_scale ) ) { return false ; }
  m_scale = avalue ;
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setN ( const double value )
{
  //
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_n ) && 0 < m_nu && 0 < m_norm ) { return false ; }
  //
  m_n    = avalue ;
  //
  // ATTENTION HERE!
  m_nu   = nu ( m_n ) ; // ATTENTION!! nu=nu(n)! 
  //
  m_norm = 1 / ( Ostap::Math::beta ( 0.5 , 0.5 * m_nu ) * std::sqrt ( m_nu ) )  ;
  //
  return true ;
}
// =========================================================================
// get the expression nu=nu(n) 
// =========================================================================
double Ostap::Math::StudentT::nu ( const double n )
{
  /// ATTENTION:  nu=nu(n)!
  return std::hypot ( 2 , n ) ;
}
// ==========================================================================
double Ostap::Math::StudentT::pdf ( const double x ) const
{
  const double y = ( x - m_M  ) / m_scale ;
  const double f = std::pow ( 1 + y * y / m_nu , -0.5 * ( m_nu + 1 ) ) ;
  return m_norm * f / m_scale  ; 
}
// ============================================================================
double Ostap::Math::StudentT::cdf ( const double y ) const
{
  const double  t    = ( y - m_M ) / m_scale ;
  return Ostap::Math::student_cdf ( t , m_M ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::StudentT::integral() const
{
  if ( m_nu <= 1 || s_equal ( m_nu , 1 ) ) { return std::numeric_limits<double>::quiet_NaN () ; }
  return 1 ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::StudentT::integral
( const double low  , 
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::StudentT::tag () const 
{ 
  static const std::string s_name = "StudentT" ;
  return Ostap::Utils::hash_combiner ( s_name , m_M , m_scale , m_n ) ; 
}
// ============================================================================
// variance
// ============================================================================
double Ostap::Math::StudentT::variance   () const
{
  if      ( m_nu <= 1 ) { return std::numeric_limits<double>::quiet_NaN () ; }
  else if ( m_nu <= 2 || s_equal ( m_nu , 2 ) )
    { return std::numeric_limits<double>::infinity () ;}
  //
  return m_scale * m_scale * m_nu / ( m_nu - 2 ) ;
}
// ============================================================================
// RMS 
// ============================================================================
double Ostap::Math::StudentT::rms         () const
{
  if      ( m_nu <= 1 ) { return std::numeric_limits<double>::quiet_NaN () ; }
  else if ( m_nu <= 2 || s_equal ( m_nu , 2 ) )
    { return std::numeric_limits<double>::infinity () ;}
  //
  return m_scale * std::sqrt ( m_nu / ( m_nu - 2 ) ) ;
}
// ============================================================================
//  get (excess) kurtosis 
// ============================================================================
double Ostap::Math::StudentT::kurtosis   () const 
{
  if      ( m_nu <= 2 ) { return std::numeric_limits<double>::quiet_NaN () ; }
  else if ( m_nu <= 4 || s_equal ( m_nu , 4 ) )
    { return std::numeric_limits<double>::infinity () ;}
  //
  return 6 / ( m_nu - 4 ) ;
}
// ============================================================================
// Bifurcated Student-T 
// ============================================================================
/*  constructor from mass, resolution and "n"-parameter 
 *  @param M     mass 
 *  @param sigma width parameter
 *  @param N     n-parameter  ( actually  n=1+|N| ) 
 */
// ============================================================================
Ostap::Math::BifurcatedStudentT::BifurcatedStudentT 
( const double mass   , 
  const double sigmaL ,
  const double sigmaR ,
  const double nL     , 
  const double nR     ) 
//
  : m_M     (      std::abs ( mass   ) )
  , m_sL    (      std::abs ( sigmaL ) )
  , m_sR    (      std::abs ( sigmaR ) )
  , m_nL    ( -1 )
  , m_nR    ( -1 )
  , m_nuL   ( -1 )
  , m_nuR   ( -1 )
  , m_normL ( -1 ) 
  , m_normR ( -1 ) 
{
  setNL ( nL ) ;  
  setNR ( nR ) ;  
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setM ( const double x )
{
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_M ) ) { return false ; }
  m_M = v ;
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setSigmaL ( const double x )
{
  const double v = std::abs ( x ) ;
  //
  Ostap::Assert ( v  ,
		  "Parameter 'sigmaL' must be non-zero"        ,
		  "Ostap::Math::BifurcatedStudentT::setSigmaL" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__      ) ;
  //
  if ( s_equal ( v , m_sL ) ) { return false ; }
  m_sL = v ;
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setSigmaR ( const double x )
{
  const double v = std::abs ( x ) ;
  //
  Ostap::Assert ( v  ,
		  "Parameter 'sigmaR' must be non-zero"        ,
		  "Ostap::Math::BifurcatedStudentT::setSigmaR" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__      ) ;
  //
  if ( s_equal ( v , m_sR ) ) { return false ; }
  m_sR = v ;
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setNL ( const double value )
{
  //
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_nL ) && 0 < m_nuL && 0 < m_normL ) { return false ; }
  //
  m_nL    = avalue ;
  //
  // ATTENTION HERE!
  m_nuL   = Ostap::Math::StudentT::nu ( m_nL ) ; // ATTENTION!! nu=nu(n)! 
  //
  m_normL  = 1 / ( Ostap::Math::beta ( 0.5 , 0.5 * m_nuL ) * std::sqrt ( m_nuL ) )  ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setNR ( const double value )
{
  //
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_nR ) && 0 < m_nuR && 0 < m_normR ) { return false ; }
  //
  m_nR    = avalue ;
  //
  // ATTENTION HERE!
  m_nuR   = Ostap::Math::StudentT::nu ( m_nR ) ; // ATTENTION!! nu=nu(n)! 
  //
  m_normR = 1 / ( Ostap::Math::beta ( 0.5 , 0.5 * m_nuR ) * std::sqrt ( m_nuR ) )  ;
  //
  return true ;
}
// ==========================================================================
double Ostap::Math::BifurcatedStudentT::pdf ( const double x ) const
{
  //
  const double y = ( x <= m_M ) ? ( x - m_M ) / m_sL : ( x - m_M ) / m_sR  ;
  //
  const double f = ( x <= m_M ) ? 
    std::pow ( 1 + y * y / m_nuL ,  -0.5 * ( m_nuL + 1 ) ) :
    std::pow ( 1 + y * y / m_nuR ,  -0.5 * ( m_nuR + 1 ) ) ;
  //
  const double n_1 = m_normL       / m_sL ;
  const double n_2 = m_normR       / m_sR ;
  const double n_t = 2 * n_1 * n_2 / ( n_1 + n_2 ) ;
  //
  return n_t * f ; 
}
// ============================================================================
double Ostap::Math::BifurcatedStudentT::cdf ( const double y ) const
{
  //
  const double n_1 = m_normL / m_sL ;
  const double n_2 = m_normR / m_sR ;
  //
  if ( y <= m_M ) 
  {
    const double  t    = ( y - m_M ) / m_sL ;
    return     2 * n_2 / ( n_1 + n_2 ) * Ostap::Math::student_cdf (  t , m_nuL ) ;  
  }
  //
  const   double  t    = ( y - m_M ) / m_sR ;
  return   1 - 2 * n_1 / ( n_1 + n_2 ) * Ostap::Math::student_cdf ( -t , m_nuR ) ;  
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::BifurcatedStudentT::integral() const 
{
  if ( m_nuL <= 1 || s_equal ( m_nuL , 1 ) ) { return std::numeric_limits<double>::quiet_NaN () ; }
  if ( m_nuR <= 1 || s_equal ( m_nuR , 1 ) ) { return std::numeric_limits<double>::quiet_NaN () ; }
  return 1 ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::BifurcatedStudentT::integral
( const double low  , 
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::BifurcatedStudentT::tag () const 
{  
  static const std::string s_name = "BiFurcatedStudentT" ;
  return Ostap::Utils::hash_combiner ( s_name , m_M , m_sL , m_sR  , m_nL , m_nR ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor from all parameters 
 *  @param mu    location parameter 
 *  @param sigma width/scale parameter 
 *  @param n     n-parameter 
 *  @param kappa asymmetry parameter 
 */
// ============================================================================
Ostap::Math::PearsonIV::PearsonIV
( const double mu       , 
  const double varsigma , 
  const double n        , 
  const double kappa     )
  : m_mu       ( mu ) 
  , m_varsigma ( std::abs ( varsigma ) )
  , m_n        ( std::abs ( n        ) ) 
  , m_kappa    ( kappa )
  , m_C        ( -1 ) 
{
  setN ( n ) ;
}
// ===========================================================================
// get value of the function 
// ===========================================================================
double Ostap::Math::PearsonIV::evaluate ( const double x ) const 
{
  const double y = ( x - m_mu ) / m_varsigma ;
  const double s = m_C * std::pow ( 1 + y * y , - m() ) / m_varsigma ;
  return s_zero ( m_kappa ) ? s : s * std::exp ( -m_kappa * std::atan ( y ) ) ;
}
// ===========================================================================
// set location parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setMu ( const double value )
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ===========================================================================
// set width/scale parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setVarsigma ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_varsigma , avalue ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ===========================================================================
// set n-parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setN ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_n , avalue ) && 0 < m_C ) { return false ; }
  m_n = avalue ;
  //
  // m_C = std::norm ( Ostap::Math::gamma ( std::complex<double> ( m() , 0.5 * nu () ) ) /
  //                   Ostap::Math::gamma ( m() ) ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  m_C = Ostap::Math::pearsonIV_g2 ( m() , 0.5 * nu() ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  return true ;
}   
// ===========================================================================
// set asymmetry parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setKappa ( const double value )
{
  if ( s_equal ( m_kappa , value ) && 0 < m_C ) { return false ; }
  m_kappa = value ;
  //
  // m_C     = std::norm ( Ostap::Math::gamma ( std::complex<double> ( m() , 0.5 * nu () ) ) /
  //                       Ostap::Math::gamma ( m() ) ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  m_C = Ostap::Math::pearsonIV_g2 ( m() , 0.5 * nu() ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  return true ;
}  
// ===========================================================================
// get the integral
// ===========================================================================
double Ostap::Math::PearsonIV::integral () const { return 1 ; }
// ===========================================================================
// get the integral between low and high limits
// ===========================================================================
double Ostap::Math::PearsonIV::integral
( const double low  ,
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low ) { return - integral ( high , low ) ; }
  //
  const bool symmetric = s_zero ( nu () ) ;
  //
  const double m0      = ( 1 < m() ) ? 0.5 * ( mode () + mean () ) : mode () ;
  const double width   = ( 2 * m() <= 3 ) ? 0.5 * 
    std::max ( m_varsigma , 0.5 * infection_width () ) : rms () ;
  //
  { // split at mode 
    const double x1 = mode () ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = x1 + 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    const double x3 = x1 - 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x3 ) + integral ( x3 , high ) ; }
  }
  //
  if ( !symmetric && ( 1 < m () ) ) 
  {
    // split at mean 
    const double x1 = mean ()  ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = x1 + 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    const double x3 = x1 - 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x3 ) + integral ( x3 , high ) ; }
  }
  //
  // more splits 
  { 
    const double x2 = m0 - 6 * width  ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    const double x3 = m0 + 6 * width  ;
    if ( low < x3 && x3 < high ) { return integral ( low , x3 ) + integral ( x3 , high ) ; }
  }
  //
  if ( !symmetric && ( 0 < nu() ) )
  {
    const double xx = m0 - 12 * width  ;
    if ( low < xx && xx < high ) { return integral ( low , xx ) + integral ( xx , high ) ; }
  }
  //
  if ( !symmetric && ( 0 > nu() ) )
  {
    const double xx = m0 + 12 * width  ;
    if ( low < xx && xx < high ) { return integral ( low , xx ) + integral ( xx , high ) ; }
  }
  //
  //
  const bool in_tail = ( high <= m0 - 10 * width ) || ( low >=  m0 + 10 * width ) ;             
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PearsonIV> s_integrator {} ;
  static char s_message[] = "Integral(PEarsonIV)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
             ( tag () , 
               &F     , 
               low    , high  ,               // low & high edges
               workspace ( m_workspace ) ,    // workspace
               in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
               in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
               m_workspace.size()   ,          // size of workspace
               s_message           , 
               __FILE__ , __LINE__ ) ;
  //
  return result ; 
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::PearsonIV::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// mode 
// ============================================================================
double Ostap::Math::PearsonIV::mode () const // mode of the distribution 
{ return m_mu - 0.5 * nu() * a() / m() ; }
// ============================================================================
// mean
// ============================================================================
double Ostap::Math::PearsonIV::mean () const // mean of the distribution 
{
  return 
    s_zero ( nu() ) ? m_mu                                  :
    1 < m ()        ? m_mu - 0.5 * a() * nu() / ( m() - 1 ) :
    std::copysign ( std::numeric_limits<double>::infinity() , -nu() ) ;
}
// ============================================================================
// (central) moment 
// ============================================================================
double Ostap::Math::PearsonIV::moment ( const unsigned short k ) const 
{
  if      ( 0 == k ) { return 1 ; }
  else if ( 1 == k ) { return 0 ; }
  //
  const bool odd = ( 1 == k % 2 )  ;
  if ( odd && s_zero ( nu () ) ) { return 0 ; }
  //
  if ( r () + 1 <= k ) 
  {
    return odd ?
      std::copysign ( std::numeric_limits<double>::infinity () , -nu() ) :
      std::numeric_limits<double>::infinity () ;
  }
  //
  const double r2  = std::pow ( r  () , 2 ) ;
  const double nu2 = std::pow ( nu () , 2 ) ;
  //
  double m2 = 1  ;
  double m1 = 0  ;
  double m  = m1 ;
  for  ( unsigned short kk = 2 ; kk <= k ; ++kk ) 
  {
    const double c = a () * ( kk - 1 ) / ( r2 * ( r () - ( kk - 1 ) ) ) ;
    m  = -2 * nu() * r () * m1 + a () * ( r2 + nu2 ) * m2 ;
    m *= c  ;
    m2 = m1 ;
    m1 = m  ;
  }
  return m ;
}
// ============================================================================
// variance                        (for m>3.5)
// ============================================================================
double Ostap::Math::PearsonIV::variance () const  // variance of the distribution 
{ return 2 * m () <= 3 ? std::numeric_limits<double>::infinity () : moment ( 2 ) ; }
// ============================================================================
// rms                              (for m>3.5)
// ============================================================================
double Ostap::Math::PearsonIV::rms     () const  // rms of the distribution 
{ return 2 * m () <= 3 ? std::numeric_limits<double>::infinity () : std::sqrt ( moment ( 2 ) ) ; }
// ============================================================================
// skewness ( for m>2) 
// ============================================================================
double Ostap::Math::PearsonIV::skewness () const  // skewness 
{
  return 
    s_zero ( nu () ) ? 0.0 :
    m () <= 2 ? std::copysign ( std::numeric_limits<double>::infinity () , -nu() ) :
    moment ( 3 ) / std::pow ( moment ( 2 ) , 1.5 ) ;
}
// ============================================================================
// (excessive) kurtosis ( for m>5/2) 
// ============================================================================
double Ostap::Math::PearsonIV::kurtosis () const  // (excessive) kurtosis 
{
  return 
    2 * m () <= 5 ? std::numeric_limits<double>::infinity () :
    moment ( 4 ) / std::pow ( moment ( 2 ) , 2) - 3 ;
}
// ============================================================================
// beta1 parameter of Pearson family (m>2) 
// ============================================================================
double Ostap::Math::PearsonIV::beta1 () const  // beta1 parameter of Pearson family 
{ 
  return 
    s_zero ( nu () ) ? 0.0 :
    m () <= 2 ? std::numeric_limits<double>::infinity () :
    std::pow ( moment ( 3 ) , 2 ) / std::pow ( moment ( 2 ) , 3 ) ;
}
// ============================================================================
// beta2 parameter of Pearson family (m>2) 
// ============================================================================
double Ostap::Math::PearsonIV::beta2 () const  // beta2 parameter of Pearson family 
{ 
  return 
    2 * m () <= 5 ? std::numeric_limits<double>::infinity () :
    moment ( 4 ) / std::pow ( moment ( 2 ) , 2) ;
}
// ============================================================================
/* distance between two infection points:
 *  distance between two points with \f$ f^{\prime\prime}=0\f$.
 *  the twp points are equidstance fro mthe mode 
 */
// ============================================================================
double Ostap::Math::PearsonIV::infection_width () const 
{
  return a () / m() * std::sqrt ( ( 4 * std::pow ( m() , 2 ) + 
                                    std::pow ( nu () , 2 ) ) / ( 2 * m() + 1 ) ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PearsonIV::tag () const 
{  
  static const std::string s_name = "PearsonIV" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma , m_n , m_kappa ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor with all parameters
 *  @param location \f$\mu\f$-parameter       \f$-\inf<\mu<+\inf\f$
 *  @param scale    \f$\sigma\f$-parameter    \f$0<\sigma\f$
 *  @param epsilon  \f$\epsilon\f$-parameter  \f$-\inf<\epsilon<+\inf\f$
 *  @param delta    \f$\delta\f$-parameter    \f$0<\epsilon<+\inf\f$
 */
// ============================================================================
Ostap::Math::SinhAsinh::SinhAsinh 
( const double location  ,
  const double scale     , 
  const double epsilon   , 
  const double delta     )
  : m_mu       (            location   ) 
  , m_sigma    ( std::abs ( scale    ) ) 
  , m_epsilon  (            epsilon    ) 
  , m_delta    ( std::abs ( delta    ) ) 
{}
// ============================================================================
bool Ostap::Math::SinhAsinh::setMu      ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SinhAsinh::setSigma      ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma  ) ) { return false ; }
  //
  Ostap::Assert ( value_                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::SinhAsinh::setSigma"       ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //              
  m_sigma  = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SinhAsinh::setEpsilon ( const double value ) 
{
  if ( s_equal ( value , m_epsilon  ) ) { return false ; }
  m_epsilon  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SinhAsinh::setDelta ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_delta  ) ) { return false ; }
  m_delta  = value_ ;
  return true ;
}
// ============================================================================
// evaluate sinhasinh-distributions
// ============================================================================
double Ostap::Math::SinhAsinh::pdf ( const double x ) const 
{
  //
  const double y = ( x - mu () ) / sigma()  ;
  const double z = shash ( y , epsilon() , delta() )  ;
  //
  const double r = s_SQRT2PIi * delta() 
    // * std::sqrt  ( ( 1.0 + z * z ) / ( 1.0 + y * y )  ) 
    *    std::hypot ( 1 , z )          / std::hypot ( 1 , y )  
    *    my_exp ( -0.5 * z * z ) ;
  //
  return  r / sigma() ;
}
// ============================================================================
// evaluate sinhasinh cimulative distribution
// ============================================================================
double Ostap::Math::SinhAsinh::cdf ( const double x ) const 
{
  //
  const double y = ( x - mu() ) / sigma()  ;
  const double z = shash ( y , epsilon () , delta () )  ;
  //
  return gsl_cdf_ugaussian_P ( z ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::SinhAsinh::integral
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::SinhAsinh::integral () const { return 1 ; } 
// ======================================================================
// get median
// ======================================================================
double Ostap::Math::SinhAsinh::median () const 
{ return m_mu - m_sigma * std::sinh ( m_epsilon / m_delta ) ; }
// ============================================================================
// get the mean for the distgribution
// ============================================================================
double Ostap::Math::SinhAsinh::mean   () const 
{
  const double d1 = 0.5 * ( 1 - m_delta ) / m_delta ;
  const double d2 = 0.5 * ( 1 - m_delta ) / m_delta ;
  //
  static const double s_const1 = std::pow ( std::exp ( 1 ) , 0.25 ) / std::sqrt ( 8 *M_PI ) ;
  //
  const double a = std::sinh ( m_epsilon / m_delta ) * s_const1 *
    ( Ostap::Math::bessel_Knu ( d1 , 0.25 ) + 
      Ostap::Math::bessel_Knu ( d2 , 0.25 ) ) ;
  //
  return m_mu - m_sigma * a ;
}
// ============================================================================
// get the variance for the distribution
// ============================================================================
double Ostap::Math::SinhAsinh::variance () const
{
  //
  const double d1 = 0.5 * ( 1 + m_delta ) / m_delta ;
  const double d2 = 0.5 * ( 1 - m_delta ) / m_delta ;
  //
  static const double s_const1 = std::pow ( std::exp ( 1 ) , 0.25 ) / std::sqrt ( 8 *M_PI ) ;
  //
  const double a = std::sinh (     m_epsilon / m_delta ) * s_const1 *
    ( Ostap::Math::bessel_Knu ( d1 , 0.25 ) + 
      Ostap::Math::bessel_Knu ( d2 , 0.25 ) ) ;
  //
  const double p1 = 0.5 * ( 2 + m_delta ) / m_delta ;
  const double p2 = 0.5 * ( 2 - m_delta ) / m_delta ;
  //
  static const double s_const2 = s_const1 / 2 ;
  //
  const double b = std::cosh ( 2 * m_epsilon / m_delta ) * s_const2 *
    ( Ostap::Math::bessel_Knu ( p1 , 0.25 ) + 
      Ostap::Math::bessel_Knu ( p2 , 0.25 ) ) ;
  //
  return m_sigma * m_sigma * ( b - a * a - 0.5 ) ;
}
// ============================================================================
// get the RMS for the distribution
// ============================================================================
double Ostap::Math::SinhAsinh::rms () const
{ return std::sqrt ( variance () ) ; }
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::SinhAsinh::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) ; // / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::SinhAsinh::tag () const 
{  
  static const std::string s_name = "SinhAsinh" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_sigma , m_epsilon , m_delta ) ; 
}
// ============================================================================

// ============================================================================
// Johnson-SU
// ============================================================================
/*  constructor with all parameters
 *  @param xi     \f$\xi\f$-parameter       \f$-\inf<\xi<+\inf\f$
 *  @param lambda \f$\lambda\f$-parameter   \f$   0<\lambda<+\inf\f$
 *  @param delta  \f$\delta\f$-parameter    \f$   0<\delta<+\inf\f$
 *  @param gamma  \f$\gamma\f$-parameter    \f$-\inf<\epsilon<+\inf\f$
 */
// ============================================================================
Ostap::Math::JohnsonSU::JohnsonSU  
( const double xi      , // related to location 
  const double lambda  , // related to variance
  const double delta   , // shape 
  const double gamma   ) // shape 
  : m_xi      (            xi) 
  , m_lambda  ( std::abs ( lambda ) ) 
  , m_delta   ( std::abs ( delta  ) ) 
  , m_gamma   (            gamma    ) 
{}
// ============================================================================
// get the mean value
// ============================================================================
double Ostap::Math::JohnsonSU::mean() const 
{
  const double d = 
    std::exp   ( 0.5 / ( m_delta * m_delta ) ) * 
    std::sinh  ( m_gamma / m_delta           ) ;
  //
  return m_xi - m_lambda * d ;
}
// ============================================================================
// get the variance
// ============================================================================
double Ostap::Math::JohnsonSU::variance () const 
{
  const double d1 = std::exp ( 1.0 / ( m_delta * m_delta ) );
  //
  const double d2 = ( d1 - 1 ) * ( d1 * std::cosh ( 2 * m_gamma / m_delta ) + 1 ) ;
  //
  return 0.5 * m_lambda * m_lambda * d2 ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setXi( const double value ) 
{
  if ( s_equal ( value , m_xi ) ) { return false ; }
  m_xi = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setGamma( const double value ) 
{
  if ( s_equal ( value , m_gamma ) ) { return false ; }
  m_gamma = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setLambda ( const double value ) 
{
  const double value_ = std::abs ( value ) ;  
  if ( s_equal ( value_ , m_lambda ) ) { return false ; }
  m_lambda = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setDelta ( const double value ) 
{
  const double value_ = std::abs ( value ) ;  
  if ( s_equal ( value_ , m_delta ) ) { return false ; }
  m_delta = value_ ;
  return true ;
}
// ============================================================================
// evaluate JohnsonSU-distributions
// ============================================================================
double Ostap::Math::JohnsonSU::pdf        ( const double x ) const 
{
  // get z 
  const long double dx  = ( x - m_xi ) / m_lambda ;
  const long double z   = m_gamma + m_delta * std::asinh ( dx ) ;
  //
  const long double res = std::exp ( -0.5 * z * z ) / std::sqrt ( 1 + dx * dx ) ;
  //
  return res * m_delta / ( m_lambda * s_SQRT2PI ) ;
}
// ============================================================================
// evaluate JohnsonSU-distributions
// ============================================================================
double Ostap::Math::JohnsonSU::cdf        ( const double x ) const 
{
  // get z 
  const long double dx  = ( x - m_xi ) / m_lambda ;
  const long double z   = m_gamma + m_delta * std::asinh ( dx ) ;
  //
  return gsl_cdf_ugaussian_P ( z ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::JohnsonSU::integral 
( const double low  ,
  const double high ) const 
{ return  s_equal ( low , high ) ? 0.0 : ( cdf ( high ) - cdf ( low ) ) ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::JohnsonSU::integral () const { return 1 ; } 
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::JohnsonSU::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::JohnsonSU::tag () const 
{ 
  static const std::string s_name = "JohnsonSU" ;
  return Ostap::Utils::hash_combiner ( s_name , m_xi , m_lambda , m_delta , m_gamma ) ; 
}
// ============================================================================


// ============================================================================
/*  Constructor from location and mean 
 *  @param mu location 
 *  @param scale the scale, scale>0
 */
// ============================================================================
Ostap::Math::Slash::Slash
( const double mu    , // location 
  const double scale ) // scale 
  : m_mu    ( mu ) 
  , m_scale ( std::abs ( scale ) ) 
{}
// ============================================================================
bool Ostap::Math::Slash::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Slash::setScale ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( v                                        ,
		  "Parameter 'scale' must be non-zero"     ,
		  "Ostap::Math::Slash::setSigma"           ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //    
  m_scale = v;
  return true ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  const long double s_slash = 0.5L / std::sqrt( 2.0L * M_PI );
  // ==========================================================================
  //  (phi(0)-phi(x))/x^2
  inline long double _slash_pdf_ ( const long double x ) 
  {
    if        ( s_zero ( x ) ) { return s_slash ; }
    else if   ( 0.1 < std::abs ( x ) ) 
    { return  ( 2 * s_slash -  Ostap::Math::gauss_pdf ( x ) ) / ( x * x ) ; }
    //
    const long double z = - 0.5L * x * x ;
    return  s_slash * ( std::expm1 ( z ) / z ) ;
  }
  // ==========================================================================
  //  Phi(x) - (phi(0)-phi(x))/x
  inline long double _slash_cdf_ ( const long double x ) 
  {
    return 
      s_equal ( x , 0 ) ? 0.5 : 
      Ostap::Math::gauss_cdf ( x ) - x * _slash_pdf_ ( x ) ;
  }
}
// ============================================================================
// evaluate slash function
// ============================================================================
double Ostap::Math::Slash::pdf ( const double x ) const
{
  const double y = ( x - m_mu ) / m_scale ;
  return _slash_pdf_ ( y ) / m_scale ;
}
// ============================================================================
// evaluate slash CDF 
// ============================================================================
double Ostap::Math::Slash::cdf ( const double x ) const
{
  const double y = ( x - m_mu ) / m_scale ;
  return _slash_cdf_ ( y ) ;
}
// ============================================================================
// get integral from low to high
// ============================================================================
double Ostap::Math::Slash::integral
( const double low  ,
  const double high ) const 
{ return s_equal ( low ,  high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ===========================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Slash::tag () const 
{ 
  static const std::string s_name = "Slash" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_scale ) ; 
}
// ============================================================================




// ============================================================================
/*  constructor with all arguments 
 *  @param mu  the mean/mode/median of the distribution
 *  @param s   the width-parameteter
 */
// ============================================================================
Ostap::Math::RaisingCosine::RaisingCosine
( const double mu , 
  const double s  ) 
  : m_mu ( mu ) 
  , m_s  ( std::abs ( s ) ) 
{}
// ============================================================================
bool Ostap::Math::RaisingCosine::setS ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_s ) ) { return false ; }
  m_s = v;
  return true ;
}
// ============================================================================
bool Ostap::Math::RaisingCosine::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// evaluate raising cosine distribution
// ============================================================================
double Ostap::Math::RaisingCosine::pdf ( const double x ) const 
{
  return 
    x <= m_mu - m_s ? 0.0 : 
    x >= m_mu + m_s ? 0.0 : 
    ( 1  + std::cos ( M_PI * ( x  - m_mu ) / m_s ) ) / ( 2 * m_s )  ;
}
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::RaisingCosine::variance () const 
{
  static const double s_c1 = ( 1./3 - 2 / ( M_PI * M_PI ) ) ;
  return m_s * m_s * s_c1 ;
}
// ============================================================================
// rms
// ============================================================================
double Ostap::Math::RaisingCosine::rms () const 
{
  static const double s_c2 = std::sqrt ( 1./3 - 2 / ( M_PI * M_PI ) ) ;
  return m_s * s_c2 ;
}
// ============================================================================
// kurtosis
// ============================================================================
double Ostap::Math::RaisingCosine::kurtosis () const 
{
  static const double s_k = 
    1.2 * ( 90. - std::pow ( M_PI , 4 ) ) / std::pow ( M_PI*M_PI  - 6. , 2 ) ;
  return  s_k ;
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::RaisingCosine::cdf      ( const double x ) const 
{
  if      ( x <= m_mu - m_s ) { return 0 ; }
  else if ( x >= m_mu - m_s ) { return 1 ; }
  //
  const double y = ( x - m_mu ) / m_s ;
  return 0.5 * ( 1 + y + std::sin ( y  * M_PI ) / M_PI ) ;
}
// ============================================================================
// evaluate the integral
// ============================================================================
double Ostap::Math::RaisingCosine::integral
( const double low  , 
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( low > high             ) { return -integral ( high , low ) ; }
  else if ( high <  m_mu - m_s     ) { return 0 ; }
  else if ( low  >  m_mu + m_s     ) { return 0 ; }
  //                                            
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// evaluate the integral
// ============================================================================
double Ostap::Math::RaisingCosine::integral () const { return 1 ; } 
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::RaisingCosine::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::RaisingCosine::tag () const 
{ 
  static const std::string s_name = "RasisingCosine" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_s ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor from all parameters 
 *  @param mu  peak location
 *  @param lambdaL ``left''  exponential slope  (lambdaL>0)
 *  @param lambdaR ``right'' exponential slope  (lambdaR>0)
 */
// ============================================================================
Ostap::Math::AsymmetricLaplace::AsymmetricLaplace 
( const double mu      , // location 
  const double lambdaL , // left  exponential slope 
  const double lambdaR ) // right exponential slope 
  : m_mu ( mu ) 
  , m_lambdaL ( std::abs ( lambdaL ) ) 
  , m_lambdaR ( std::abs ( lambdaR ) ) 
{}
// ============================================================================
bool Ostap::Math::AsymmetricLaplace::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::AsymmetricLaplace::setLambdaL ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_lambdaL ) ) { return false ; }
  m_lambdaL = v;
  return true ;
}
// ============================================================================
bool Ostap::Math::AsymmetricLaplace::setLambdaR ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_lambdaR ) ) { return false ; }
  m_lambdaR = v;
  return true ;
}
// ============================================================================
// evaluate asymmetic laplace function
// ============================================================================
double  Ostap::Math::AsymmetricLaplace::pdf ( const double x ) const 
{
  const long double l = 1.0L / ( m_lambdaL + m_lambdaR ) ;
  return 
    x < m_mu ? 
    l * std::exp (   ( x - m_mu ) / m_lambdaL ) :
    l * std::exp ( - ( x - m_mu ) / m_lambdaR ) ;
}
// ============================================================================
// evaluate CDF
// ============================================================================
double  Ostap::Math::AsymmetricLaplace::cdf ( const double x ) const 
{
  const long double l = 1.0L / ( m_lambdaL + m_lambdaR ) ;
  return 
    x < m_mu ? 
    m_lambdaR     * l * std::exp (   ( x - m_mu ) / m_lambdaL ) :
    1 - m_lambdaL * l * std::exp ( - ( x - m_mu ) / m_lambdaR ) ;
}
// ============================================================================
// get integral from low to high
// ============================================================================
double Ostap::Math::AsymmetricLaplace::integral
( const double low  ,
  const double high ) const 
{ return s_equal ( low ,  high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::AsymmetricLaplace::tag () const 
{ 
  static const std::string s_name = "AsymmetricLaplace" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_lambdaL , m_lambdaR ) ; 
}
// ============================================================================
// get mean 
// ============================================================================
double Ostap::Math::AsymmetricLaplace::mean () const
{ return  m_mu + ( m_lambdaR - m_lambdaL ) ; }
// ============================================================================
// get mean 
// ============================================================================
double Ostap::Math::AsymmetricLaplace::variance () const
{ 
  const double l  = lambdaL () ; 
  const double r  = lambdaR () ;
  const double l2 = l * l ;
  const double r2 = r * r ;  
  return  ( l2 * l2 + r2 * r2 ) / ( l * r ) ;  
}
// ============================================================================
// get median
// ============================================================================
double Ostap::Math::AsymmetricLaplace::median () const
{
  if ( s_equal ( m_lambdaL , m_lambdaR ) ) { return m_mu ; }
  //
  const double l2 = m_lambdaL * m_lambdaL ;
  const double r2 = m_lambdaR * m_lambdaR ;
  //
  return ( m_lambdaL >= m_lambdaR ) ? 
    ( m_mu + m_lambdaL * std::log ( 0.5 * ( 1 + r2 / l2 ) ) ) :
    ( m_mu - m_lambdaR * std::log ( 0.5 * ( 1 + l2 / r2 ) ) ) ; 
}
// ============================================================================
// get skewness 
// ============================================================================
double Ostap::Math::AsymmetricLaplace::skewness () const
{
  if ( s_equal ( m_lambdaL , m_lambdaR ) ) { return 0 ; }
  //
  const double K2 = k2()    ;
  const double K4 = K2 * K2 ;
  const double K6 = K2 * K4 ;
  //
  return 2 * ( 1 - K6  ) / std::pow ( 1 + K4 , 1.5 ) ; 
}
// ============================================================================
// get (excess) kurtosis 
// ============================================================================
double Ostap::Math::AsymmetricLaplace::kurtosis () const
{
  //
  const double K2 = k2()    ;
  const double K4 = K2 * K2 ;
  const double K8 = K4 * K4 ;
  //
  return 6 * ( 1 + K8 ) / std::pow ( 1 + K4 , 2 ) ;
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
// double Ostap::Math::AsymmetricLaplace::non_gaussian 
// ( const double xlow  ,
//   const double xhigh ) const
// {
//   if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
//   else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
//   //
//   const double I_CB = integral ( xlow , xhigh ) / integral ();
//   //
//   const double m = mean () ;
//   const double s = rms  () ;
//   //
//   const double I_G  =
//     Ostap::Math::gauss_cdf ( xhigh , m , s ) -
//     Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
//   //
//   return 1 - I_G / I_CB ;
// }
// ============================================================================


// ============================================================================
/*  constructor from all arguments            
 *  @param mean  the mean/mode/location of the peak 
 *  @param q     q-value   (q<3, for q>3, it will be reflected)
 *  @param scale 
 */
// ============================================================================
Ostap::Math::QGaussian::QGaussian 
( const double mean  ,  // mean/mode/location 
  const double scale ,  // scale/sigma
  const double q     )  // q-parameter (shape)  
  : m_mean  ( mean               )
  , m_scale ( std::abs ( scale ) )
  , m_q     ( 1                  )
  , m_cq    ( s_SQRTPI           ) 
{
  setQ ( q ) ; 
}
// ============================================================================
// evaluate  pdf  for q-Gaussian distribution
// ============================================================================
// evaluate PDF of q-Gaussian distribution 
// ============================================================================
double Ostap::Math::QGaussian::pdf ( const  double x ) const 
{
  //
  if ( 1 == m_q || s_equal ( m_q , 1 ) ) { return gauss_pdf ( x , m_mean , m_scale ) ; }
  //
  const double dx  = ( x - m_mean ) / m_scale ;
  //
  static const double s_sq2 = std::sqrt ( 2.0 ) ;
  //
  return Ostap::Math::tsallis_qexp ( -0.5 * dx * dx , m_q ) / 
    ( s_sq2 * m_scale * m_cq ) ;  
}
// ============================================================================
// set mean 
// ============================================================================
bool Ostap::Math::QGaussian::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean ) ) { return false ; }
  m_mean = value ;
  return true ;
}
// ============================================================================
// set q
// ============================================================================
bool Ostap::Math::QGaussian::setQ ( const double value ) 
{
  if ( value > 3 ) { return setQ ( 6  - value ) ; } // ATTENTION! 
  //
  if ( s_equal ( value , m_q ) ) { return false ; }
  //
  m_q  = value ;
  //
  m_cq = s_SQRTPI ; 
  //
  if      ( 1 > m_q ) 
  {
    const long double q  = m_q ;
    const long double g1 = std::lgamma ( 1.0             / ( 1.0L - q ) ) ;
    const long double g2 = std::lgamma ( 0.5 * ( 3 - q ) / ( 1.0L - q ) ) ;
    m_cq *= 2 * std::exp ( g1  - std::log ( 3.0L - q ) - 
                           0.5 * std::log ( 1.0L - q ) - g2 ) ;
  }
  else if ( 1 < m_q )
  {
    const long double q  = m_q ;
    const long double g1 = std::lgamma ( 1.0                / ( q - 1.0L ) ) ;
    const long double g2 = std::lgamma ( 0.5 * ( 3.0L - q ) / ( q - 1.0L ) ) ;
    m_cq *=     std::exp ( g2 - 0.5 * std::log ( q - 1.0L ) - g1 ) ;
    //
  }
  //
  return true ;
}
// ============================================================================
// set scale
// ============================================================================
bool Ostap::Math::QGaussian::setScale ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( v                                        ,
		  "Parameter 'scale' must be non-zero"     ,
		  "Ostap::Math::QGaussian::setScale"       ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //    
  m_scale = v ;
  return true ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::QGaussian::integral
( const double low  , 
  const double high ) const 
{
  ///
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) {  return -integral ( high , low ) ; }
  ///
  if ( 1 == m_q || s_equal ( m_q , 1 ) ) 
  {
    return 
      gauss_cdf  ( high , m_mean , m_scale ) - 
      gauss_cdf  (  low , m_mean , m_scale ) ;
  }
  //
  if ( m_q  > 1 ) 
  {
    if ( low < m_mean && m_mean < high ) 
    {
      const double dx1 = m_mean -  low    ;
      const double dx2 = high   -  m_mean ;
      return 
        dx1 < dx2 ? 
        2 * integral ( low    , m_mean ) + integral ( 2 * m_mean - low ,              high ) :
        2 * integral ( m_mean , high   ) + integral (              low , 2 * m_mean - high ) ;
    }
    else if ( high - low > 3 * m_scale ) 
    {
      const double mid = 0.5 * ( low  + high ) ;
      return integral ( low , mid ) + integral ( mid , high ) ;
    } 
    //
  }
  
  double xlow  = low  ;
  double xhigh = high ;
  
  if ( m_q < 1 ) 
  {
    static const double s_sq2 = std::sqrt ( 2.0 ) ;
    const double win  = s_sq2 * m_scale / std::sqrt ( 1.0L - m_q ) ;
    const double xmin = m_mean - win ;
    const double xmax = m_mean + win ;
    if ( high <= xmin || low >= xmax ) { return 0 ; } // RETURN
    xlow  = std::max ( xmin , xlow  ) ;
    xhigh = std::min ( xmax , xhigh ) ;
  }
  //  are we already in the tail? 
  const bool in_tail = 
    std::min ( std::abs ( xhigh - m_mean ) , std::abs ( m_mean - xlow ) )  > 8 * m_scale ;   
  //
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<QGaussian> s_integrator {} ;
  static char s_message[] = "Integral(QGaussian)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::QGaussian::tag () const 
{ 
  static const std::string s_name = "QGaussian" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_q , m_scale ) ; 
}
// ============================================================================





// ============================================================================
/*  constructor from all arguments            
 *  @param mean  the mean/mode/location of the peak 
 *  @param kappa  kappa-value 
 *  @param scale 
 */
// ============================================================================
Ostap::Math::KGaussian::KGaussian 
( const double mean  ,  // mean/mode/location 
  const double scale ,  // scale/sigma
  const double kappa )  // kappa-parameter (shape)  
  : m_mean  ( mean               )
  , m_scale ( std::abs ( scale ) )
  , m_k     ( 100   )
  , m_kappa ( kappa )
{
  setKappa ( kappa ) ; 
}
// ============================================================================
// evaluate  pdf  for k-Gaussian distribution
// ============================================================================
double Ostap::Math::KGaussian::pdf ( const  double x ) const 
{
  //
  if ( 0 == m_k || s_zero ( m_k ) ) { return gauss_pdf ( x , m_mean , m_scale ) ; }
  //
  const double dx  = ( x - m_mean ) / m_scale ;
  //
  return m_Zk / m_scale * Ostap::Math::kaniadakis_kexp ( -0.5 * dx * dx , m_k ) ;
}
// ============================================================================
// set mean 
// ============================================================================
bool Ostap::Math::KGaussian::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean ) ) { return false ; }
  m_mean = value ;
  return true ;
}
// ============================================================================
// set kappa 
// ============================================================================
bool Ostap::Math::KGaussian::setKappa ( const double value ) 
{
  //
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_kappa ) && std::abs ( m_k ) <= 1 ) { return false ; }
  //
  m_kappa = avalue                ;
  m_k     = std::tanh ( m_kappa ) ;
  //
  if ( s_zero ( m_k ) ) { m_Zk = s_SQRT2PIi ; }
  else 
  {
    m_Zk = std::sqrt ( m_k / M_PI ) 
      * ( 1 + 0.5 * m_k ) 
      * std::exp ( std::lgamma ( 0.5/m_k + 0.25 ) - std::lgamma ( 0.5/m_k - 0.25 ) ) ;  
  }
  //
  return true ;
}
// ============================================================================
// set scale
// ============================================================================
bool Ostap::Math::KGaussian::setScale ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( v                                        ,
		  "Parameter 'scale' must be non-zero"     ,
		  "Ostap::Math::KGaussian::setScale"       ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //    
  m_scale = v ;
  return true ;
}
// ============================================================================
// get variance 
// ============================================================================
double Ostap::Math::KGaussian::variance   () const 
{
  if ( 0 == m_k || s_zero ( m_k ) ) { return m_scale * m_scale ; }
  
  const double f1 = std::exp ( std::lgamma ( 0.5/m_k + 0.25 ) - 
                               std::lgamma ( 0.5/m_k - 0.25 ) ) ;
  const double f2 = 4 * m_k * ( 2 + m_k ) / ( ( 2 - m_k ) * ( 4 - 9 * m_k * m_k ) ) ;
  return 2 * m_scale * m_scale * f2 * f1 * f1 ;  
}
// ============================================================================
// get RMS
// ============================================================================
double Ostap::Math::KGaussian::rms () const 
{ return std::sqrt ( variance () ) ; }  
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::KGaussian::integral
( const double low  , 
  const double high ) const 
{
  ///
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) {  return -integral ( high , low ) ; }
  ///
  if ( 0 == m_k || s_zero ( m_k ) ) 
  {
    return 
      gauss_cdf  ( high , m_mean , m_scale ) - 
      gauss_cdf  (  low , m_mean , m_scale ) ;
  }
  //
  // split into some "reasonable" intervals 
  if ( low <  m_mean  && m_mean < high ) 
  { return integral ( low , m_mean ) + integral ( m_mean , high ) ; }
  //
  { 
    const double x1 = m_mean + 3 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 3 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }    
  }
  //
  { 
    const double x1 = m_mean + 5 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 5 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }    
  }
  //
  { 
    const double x1 = m_mean + 10 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 10 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }    
  }
  //
  { 
    const double x1 = m_mean + 15 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 15 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    
  }
  //
  const double x_low  = m_mean - 15 * m_scale  ;
  const double x_high = m_mean + 15 * m_scale  ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<KGaussian> s_integrator {} ;
  static char s_message[] = "Integral(KGaussian)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::KGaussian::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::KGaussian::tag () const 
{ 
  static const std::string s_name = "KGaussian" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_kappa , m_scale ) ; 
}
// ============================================================================



namespace 
{
  // ==========================================================================
  const double z_SMALL = 1.e-6 ;
  // ==========================================================================
  inline double _knu_ ( const double z , const double nu ) 
  {
    const double zh  = 0.5 * z ;
    const double zh2 = zh  * zh ;
    const double gn  = std::tgamma ( nu ) ;
    return gn * std::pow ( zh , -nu ) * ( 1 + zh2 /(1-nu) + 0.5 * zh2 * zh2 / ( ( 1 - nu ) * ( 2 - nu ) ) ) ;
  }
  // ==========================================================================
  ///  evaluate \f$ K_{\nu}(z) \f$ for small values of z 
  inline double knu  ( const double z , const double nu ) 
  { 
    return 
      z < z_SMALL && !s_zero ( nu ) ?
          0.5 * ( _knu_( z , nu ) + _knu_  ( z , -nu ) ) 
          : Ostap::Math::bessel_Knu ( nu , z ) ;
  }
  // ===========================================================================
  /** evaluate \f$ z^{\nu} K^{*}_{\nu}(z)}\f$ for small values of z,
   *  where \f$ K^*_{\nu}(z)\f$ is a scaled modified Bessel function 
   */
  inline double z_knu_scaled ( const double z , 
                               const double nu ) 
  {
    //
    if ( s_zero ( z ) ) 
    { return ( nu <= 0 ) ? 0.0 : std::pow ( 2 , nu - 1 ) * std::tgamma ( nu  ) ; }
    //
    if      ( z > z_SMALL ) { return std::pow ( z , nu ) * Ostap::Math::bessel_Knu_scaled ( nu , z ) ; }
    //
    if      ( nu >  0.2 ) { return 0.5 * std::pow ( 2       ,  nu ) * std::tgamma (  nu ) ; }
    else if ( nu < -0.2 ) { return 0.5 * std::pow ( 2/(z*z) , -nu ) * std::tgamma ( -nu ) ; }
    //
    if ( s_zero ( nu ) ) { return z * ( -M_EULER - std::log ( 0.5 * z ) ) ; }
    //
    const double zh  = 0.5 * z ;
    const double zh2 = zh  * zh ;
    const double gn1 = std::tgamma (  nu ) ;
    const double gn2 = std::tgamma ( -nu ) ;
    //
    // const double g1  = gn1 * std::pow ( zh , -nu ) * ( 1 + zh2 /(1-nu) + 0.5 * zh2 * zh2 / ( ( 1 - nu ) * ( 2 - nu ) ) ) ;
    // const double g2  = gn2 * std::pow ( zh ,  nu ) * ( 1 + zh2 /(1+nu) + 0.5 * zh2 * zh2 / ( ( 1 + nu ) * ( 2 + nu ) ) ) ;
    //
    const double g1  = gn1 * ( 1 + zh2 /(1-nu) + 0.5 * zh2 * zh2 / ( ( 1 - nu ) * ( 2 - nu ) ) ) ;
    const double g2  = gn2 * ( 1 + zh2 /(1+nu) + 0.5 * zh2 * zh2 / ( ( 1 + nu ) * ( 2 + nu ) ) ) ;
    ///
    return 0 <= nu ? 
      0.5 * ( g1 + std::pow ( zh,  2 * nu ) * g2 ) :
      0.5 * ( g2 + std::pow ( zh, -2 * nu ) * g1 ) ;
  }
  // ==========================================================================
  /// calculate \f$ z K_{nu+1}(z)/K_{nu}(z) \f$ 
  inline double _AL2_ ( const  double nu , const double z ) 
  {
    // 
    // return z * 
    //   Ostap::Math::bessel_Knu_scaled ( nu + 1 , z ) / 
    //   Ostap::Math::bessel_Knu_scaled ( nu     , z ) ;
    //
    if ( z_SMALL  <=  z ) 
    {
      return z * 
        Ostap::Math::bessel_Knu_scaled ( nu + 1 , z ) / 
        Ostap::Math::bessel_Knu_scaled ( nu     , z ) ;
    }
    //
    if ( s_equal ( nu , -1 ) )
    {
      const double zh = 0.5 * z ;
      const double zlog = std::log ( zh ) ;
      return z * z * ( -M_EULER - zlog ) / ( 1 + z * zh * zlog ) ;
    }
    //
    else if ( s_equal ( nu , 0 ) )
    {
      const double zh = 0.5 * z ;
      const double zlog = std::log ( zh ) ;
      return ( 1 + zh * zh * ( 1 + 2 * zlog ) ) / ( -M_EULER + ( 1 - M_EULER ) * zh * zh - zlog ) ;
    }
    //
    else if ( nu < -1.15 ) 
    { 
      return 0.5 * z * z / abs ( nu ) ;   
    }
    else if ( nu < -1 ) 
    {
      return z * 
        Ostap::Math::bessel_Knu_scaled ( nu + 1 , z ) / 
        Ostap::Math::bessel_Knu_scaled ( nu     , z ) ;
    }
    else if ( nu < -0.2 ) 
    {
      const double d  = std::abs  ( nu ) ;
      const double xh = 0.5 * z ;
      return 2 * std::tgamma ( 1 - d ) / std::tgamma ( d ) * std::pow ( xh , 2*d ) ;
    }
    else if ( std::abs ( nu ) <= 0.2 ) 
    {
      return z * knu ( z , nu + 1 ) / knu ( z , nu ) ;
    }
    //
    return 2 * nu ;
  }
  // ==========================================================================
}
// ============================================================================



// ============================================================================
/* constructor from mu, sigma, zeta and kappa 
 *  @param mu    related to location 
 *  @param beta  related to asymmetry
 *  @param sigma related to width 
 *  @param zeta  related to what ?
 */
// ============================================================================
Ostap::Math::Hyperbolic::Hyperbolic
( const double mu     ,   // related to location 
  const double sigma  ,   // related to withs  
  const double zeta   ,   // shape parameter
  const double kappa  )   // related to asymmetry 
  : m_mu    ( mu    ) 
  , m_sigma ( -1    )
  , m_zeta  ( -1    )
  , m_kappa ( kappa )
  , m_AL    ( -1    )  
  , m_N     ( -1    )
{
  setSigma ( sigma ) ;
  setZeta  ( zeta  ) ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setSigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( avalue                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Hyperbolic::setSigma"      ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //              
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setZeta ( const double value ) 
{
  //
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_zeta ) &&  ( 0 < m_AL ) && ( 0 < m_N ) ) { return false ; }
  m_zeta = avalue ;
  //
  m_AL = std::sqrt ( _AL2_ ( 1 , m_zeta ) ) ;
  m_N  = 1 / ( z_knu_scaled ( m_zeta , 1 ) ) ; 
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setKappa ( const double value ) 
{
  if ( s_equal ( value , m_kappa ) ) { return false ; }
  m_kappa = value ;
  return true ;
}
// ============================================================================
/*  set "standard" parameters 
 *  @param mu     mu-parameter, location 
 *  @param beta   beta-parameter, asymmetry 
 *  @param gamma  gamma-parameter, shape 
 *  @param delta  delta-parameter, scale 
 *
 *  \f$ \alpha = \sqrt{\beta^2+\gamma^2} \f$ 
 */ 
// ============================================================================
bool Ostap::Math::Hyperbolic::setStandard
( const double mu     ,
  const double beta   , 
  const double gamma  , 
  const double delta  )
{
  bool modified = !s_equal ( m_mu , mu )  ;
  //
  m_mu     = mu    ;
  //
  const double _zeta   = std::abs ( delta ) * std::abs ( gamma ) ;
  if ( !s_equal ( m_zeta , _zeta ) ) { modified = true ; }
  m_zeta = _zeta ;
  //
  if ( modified ) { m_AL = std::sqrt ( _AL2_ ( 1  , m_zeta ) ) ; }
  //
  const double _sigma = m_AL / std::abs ( gamma ) ;
  if ( !s_equal ( m_sigma , _sigma ) ) { modified = true ; }
  m_sigma = _sigma ;
  //
  if ( modified ) { m_N = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , 1 ) ) ; }
  //
  const double _kappa = beta / m_sigma ;
  if ( !s_equal ( m_kappa , _kappa ) ) { modified = true ; }
  m_kappa = _kappa ;
  //
  return modified ;
}
// ============================================================================
// calculate the mean of the distribution  
// ============================================================================
double Ostap::Math::Hyperbolic::mean () const 
{ return m_mu + m_kappa * m_sigma ; }
// ============================================================================
// get the actual mode of the distribution
// ============================================================================
double Ostap::Math::Hyperbolic::mode () const 
{ return m_mu + m_kappa * m_sigma * m_zeta / ( m_AL * m_AL ) ; }
// ============================================================================
// get the variance/dispersion 
// ============================================================================
double Ostap::Math::Hyperbolic::variance () const 
{
  //
  const double s2 = sigma2 () ;
  const double k2 = kappa2 () ;
  const double z2 = zeta2  () ;
  //
  return s2 + k2 * s2 * ( _AL2_ ( 1 + 1 , m_zeta ) /  ( m_AL * m_AL ) - 1.0 ) ;
}
// ============================================================================
// evaluate  pdf  for the Hyperbolic distribution
// ============================================================================
double Ostap::Math::Hyperbolic::pdf ( const double x ) const 
{
  //
  const double dx =  ( x - m_mu ) / m_sigma ;
  //
  const double a2 = m_AL * m_AL            ;
  const double ka = m_kappa * m_kappa + a2 ;
  //
  const double q  = 
    - std::sqrt ( ka * ( m_zeta * m_zeta / a2 + dx * dx ) )  
    + m_kappa * dx    
    + m_zeta          ;  // from normalization 
  //
  const double aa = 0.5 * a2 / ( m_sigma * std::sqrt ( ka ) ) ;
  //
  return m_N * std::exp ( q ) * aa ;
}
// ============================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::Hyperbolic::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0        ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low  ) ; } // RETURN
  //
  const double m1    = mode () ;
  const double m2    = mean () ;
  const double mmin  = std::min ( m1 , m2 ) ;
  const double mmax  = std::max ( m1 , m2 ) ;
  const double mlow  = mmin - 5 * m_sigma ;
  const double mhigh = mmax + 5 * m_sigma ;
  //
  const double mc [] = { mmin - 3.0 * m_sigma , 
                         mmax + 3.0 * m_sigma , 
                         mlow , mhigh         } ;
  //
  for ( const double c : mc ) 
  { if ( low < c  && c < high ) { return integral ( low , c ) + integral ( c , high ) ; } }
  //
  // in tails 
  const bool in_tail = ( high <= mlow ) || ( low >= mhigh ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Hyperbolic> s_integrator {} ;
  static char s_message[] = "Integral(Hyperbolic)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     ,  
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size () ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::Hyperbolic::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Hyperbolic::tag () const 
{ 
  static const std::string s_name = "Hyperbolic" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_sigma , m_zeta , m_kappa  ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor from mu, sigma, zeta and kappa 
 *  @param mu     related to location 
 *  @param sigma  related to width 
 *  @param zeta   related to kurtosis 
 *  @param kappa  related to asymmetry
 *  @param lambda shape-related    
 */
// ============================================================================
Ostap::Math::GenHyperbolic::GenHyperbolic
( const double mu     ,  // related to location 
  const double sigma  ,  // related to width 
  const double zeta   ,  // related to kurtosis
  const double kappa  ,  // related to asymmetry 
  const double lambda ) 
  : m_mu     ( mu    ) 
  , m_sigma  ( std::abs ( sigma ) )
  , m_zeta   ( zeta  ) 
  , m_kappa  ( kappa ) 
  , m_lambda ( lambda ) 
  , m_AL     ( -1 ) 
  , m_N      ( -1 ) 
{
  setSigma  ( sigma  ) ;  
  setLambda ( lambda ) ;  
  setZeta   ( zeta   ) ;  
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setSigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( avalue                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::GenHyperbolic::setSigma"   ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //              
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setKappa ( const double value ) 
{
  if ( s_equal ( value , m_kappa ) ) { return false ; }
  m_kappa = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setZeta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_zeta ) &&  ( 0 < m_AL ) && ( 0 < m_N ) ) { return false ; }
  m_zeta = avalue ;
  //
  m_AL = std::sqrt ( _AL2_ ( m_lambda , m_zeta ) ) ;
  m_N  = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , m_lambda ) ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setLambda ( const double value ) 
{
  if ( s_equal ( value , m_lambda ) && ( 0 < m_AL ) && ( 0 < m_N ) ) { return false ; }
  m_lambda = value ;
  //
  m_AL = std::sqrt ( _AL2_ ( m_lambda , m_zeta ) ) ;
  m_N  = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , m_lambda ) ) ;
  //
  return true ;
}
// ============================================================================
/*  set "standard" parameters 
 *  @param mu     mu-parameter, location 
 *  @param beta   beta-parameter, asymmetry 
 *  @param gamma  alpha-parameter, shape 
 *  @param delta  delta-parameter, scale 
 *  @param lambda lambda-parameter, shape   
 *
 *  \f$ \alpha = \sqrt{\beta^2+\gamma^2} \f$ 
 *  \f[ \begin{array}{ll} 
 *    \delta \ge 0, left| \beta \right|<  \alpha &  if~\lambda > 0 \\ 
 *    \delta >   0, left| \beta \right|<  \alpha &  if~\lambda = 0 \\ 
 *    \delta >   0, left| \beta \right|\le\alpha &  if~\lambda < 0 
 *  \end{array}\f] 
 */ 
// ============================================================================
bool Ostap::Math::GenHyperbolic::setStandard
( const double mu     ,
  const double beta   , 
  const double gamma  , 
  const double delta  ,
  const double lambda ) 
{
  bool modified = !s_equal ( m_mu , mu ) || !s_equal ( m_lambda , lambda ) ;
  //
  m_mu     = mu     ;
  m_lambda = lambda ;
  //
  const double _zeta   = std::abs ( delta ) * std::abs ( gamma ) ;
  if ( !s_equal ( m_zeta , _zeta ) ) { modified = true ; }
  m_zeta = _zeta ;
  //
  if ( modified ) { m_AL = std::sqrt ( _AL2_ ( m_lambda , m_zeta ) ) ; }
  //
  const double _sigma = m_AL / std::abs ( gamma ) ;
  if ( !s_equal ( m_sigma , _sigma ) ) { modified = true ; }
  m_sigma = _sigma ;
  //
  if ( modified ) { m_N = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , m_lambda ) ) ; }
  //
  const double _kappa = beta / m_sigma ;
  if ( !s_equal ( m_kappa , _kappa ) ) { modified = true ; }
  m_kappa = _kappa ;
  //
  return modified ;
}
// ============================================================================
// evaluate  pdf  for Generalised Hyperbolic distribution
// ============================================================================
double Ostap::Math::GenHyperbolic::pdf ( const  double x ) const 
{  
  const double dx   =  ( x - m_mu ) / m_sigma ;
  //
  const double k2   = m_kappa * m_kappa ;
  const double k2pA = k2 + m_AL * m_AL  ;
  const double z_A  = m_zeta    / m_AL  ;
  //
  const double arg2 =  k2pA * ( dx * dx + z_A * z_A ) ;
  const double arg  =  std::sqrt ( arg2 ) ;
  //
  // NB: we use scaled bessel function here!
  const double kfun = Ostap::Math::bessel_Knu_scaled ( m_lambda - 0.5 , arg ) ;
  //
  const double f   = 
    + std::log ( kfun )             // scaled bessel function 
    - arg                           // "unscale" it 
    + m_zeta                        // from normalzation 
    + m_kappa * dx                  // asymmetry factor  
    + ( m_lambda - 0.5 ) * std::log ( arg * m_sigma * m_sigma / k2pA ) ;
  //
  return m_N * std::exp ( f ) * std::pow ( gamma2() , m_lambda ) ;
}
// ============================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::GenHyperbolic::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0        ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low  ) ; } // RETURN
  //
  const double m1    = mean () ;
  const double mlow  = m1 - 5 * m_sigma ;
  const double mhigh = m1 + 5 * m_sigma ;
  //
  const double mc [] = { m1 - 3.0 * m_sigma , 
                         m1 + 3.0 * m_sigma , 
                         mlow , mhigh       } ;
  //
  for ( const double c : mc ) 
  { if ( low < c  && c < high ) { return integral ( low , c ) + integral ( c , high ) ; } }
  //
  // in tails 
  const bool in_tail = ( high <= mlow ) || ( low >= mhigh ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenHyperbolic> s_integrator {} ;
  static char s_message[] = "Integral(GenHyperbolic)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     ,  
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size () ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::GenHyperbolic::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ==============================================================================
// get mean value 
// ==============================================================================
double Ostap::Math::GenHyperbolic::mean        () const 
{ return m_mu + m_kappa * m_sigma ; }
// ==============================================================================
// get variance 
// ==============================================================================
double Ostap::Math::GenHyperbolic::variance    () const 
{
  //
  const double s2 = sigma2 () ;
  const double k2 = kappa2 () ;
  const double z2 = zeta2  () ;
  //
  return s2 + k2 * s2 * ( _AL2_ ( m_lambda + 1 , m_zeta ) /  ( m_AL * m_AL ) - 1.0 ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::GenHyperbolic::tag () const 
{ 
  static const std::string s_name = "GHD" ;
  return Ostap::Utils::hash_combiner ( s_name   , 
                             m_mu     , 
                             m_sigma  , 
                             m_kappa  , 
                             m_zeta   , 
                             m_lambda ) ; 
}
// ============================================================================



// ============================================================================
/*  constructor with full parameters 
 *  @param mu peak location 
 *  @param sigma sigma for Gaussian Core 
 *  @param kL    left tail parameter 
 *  @param kR    right tail parameter 
 */
// ============================================================================
Ostap::Math::Das::Das 
( const double mu     , // location parameter 
  const double sigma  , // width parameter 
  const double kL     , // left tails 
  const double kR     ) // right tail 
  : m_core  ( mu , sigma )
  , m_left  ( kL )
  , m_right ( kR ) 
{}
// ============================================================================
// evaluate  pdf 
// ============================================================================
double Ostap::Math::Das::pdf ( const  double x ) const 
{
  const double _xL = xL () ;
  if ( x <= _xL ) { return m_left  ( x , _xL , m_core ( _xL ) , m_core.dFoF ( _xL ) ) ;  }
  const double _xR = xR () ;
  if ( x >= _xR ) { return m_right ( x , _xR , m_core( _xR ) , m_core.dFoF ( _xR ) ) ;  }
  //
  return m_core ( x ) ; 
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::Das::integral () const 
{
  //
  const double _xL = xL () ;
  const double _xR = xR () ;
  return m_core .integral ( _xL , _xR ) 
       + m_left .integral ( _xL , _xL , m_core ( _xL ) , m_core.dFoF ( _xL ) ) 
       + m_right.integral ( _xR , _xR , m_core ( _xR ) , m_core.dFoF ( _xR ) ) ; 
}
// =============================================================================
// get the integral 
// =============================================================================
double Ostap::Math::Das::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) { return - integral ( high , low ) ; }
  //
  const double _xL = xL () ;
  const double _xR = xR () ;
  //
  if ( low < _xL && _xL < high ) { return integral ( low , _xL ) + integral ( _xL , high ) ; }
  if ( low < _xR && _xR < high ) { return integral ( low , _xR ) + integral ( _xR , high ) ; }
  //
  if ( high <= _xL ) { return m_left .integral ( low , high , _xL , m_core ( _xL ) , m_core.dFoF ( _xL ) ) ; }
  if ( low  >= _xR ) { return m_right.integral ( low , high , _xR , m_core ( _xR ) , m_core.dFoF ( _xR ) ) ; }
  //
  return m_core.integral ( low , high ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Das::tag () const 
{ 
  static const std::string s_name = "Das" ;
  return Ostap::Utils::hash_combiner ( s_name  , 
                                       m_core .tag  () ,  
                                       m_left .tag  () , 
                                       m_right.tag  () ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor with full parameters 
 *  @param mu peak location 
 *  @param sigma sigma for Gaussian Core 
 *  @param kL    left tail parameter 
 *  @param kR    right tail parameter 
 */
// ============================================================================
Ostap::Math::ADas::ADas 
( const double mu     , // location parameter 
  const double sigmaL , // leftt width parameter 
  const double sigmaR , // right width parameter
  const double kL     , // left tails 
  const double kR     ) // right tail 
  : m_core  ( mu , sigmaL , sigmaR )
  , m_left  ( kL )
  , m_right ( kR ) 
{}
// ============================================================================
// evaluate  pdf 
// ============================================================================
double Ostap::Math::ADas::pdf ( const  double x ) const 
{
  const double _xL = xL () ;
  if ( x <= _xL ) { return m_left  ( x , _xL , m_core ( _xL ) , m_core.dFoF ( _xL ) ) ;  }
  const double _xR = xR () ;
  if ( x >= _xR ) { return m_right ( x , _xR , m_core ( _xR ) , m_core.dFoF ( _xR ) ) ;  }
  //
  return m_core ( x ) ; 
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::ADas::integral () const 
{
  //
  const double _xL = xL () ;
  const double _xR = xR () ;
  return m_core .integral ( _xL , _xR ) 
       + m_left .integral ( _xL , _xL , m_core ( _xL ) , m_core.dFoF ( _xL ) ) 
       + m_right.integral ( _xR , _xR , m_core ( _xR ) , m_core.dFoF ( _xR ) ) ; 
}
// =============================================================================
// get the integral 
// =============================================================================
double Ostap::Math::ADas::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) { return - integral ( high , low ) ; }
  //
  const double _xL = xL () ;
  const double _xR = xR () ;
  //
  if ( low < _xL && _xL < high ) { return integral ( low , _xL ) + integral ( _xL , high ) ; }
  if ( low < _xR && _xR < high ) { return integral ( low , _xR ) + integral ( _xR , high ) ; }
  //
  if ( high <= _xL ) { return m_left .integral ( low , high , _xL , m_core ( _xL ) , m_core.dFoF ( _xL ) ) ; }
  if ( low  >= _xR ) { return m_right.integral ( low , high , _xR , m_core ( _xR ) , m_core.dFoF ( _xR ) ) ; }
  //
  return m_core.integral ( low , high ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::ADas::tag () const 
{ 
  static const std::string s_name = "ADas" ;
  return Ostap::Utils::hash_combiner ( s_name  , 
                                       m_core .tag  () ,  
                                       m_left .tag  () , 
                                       m_right.tag  () ) ; 
}
// ============================================================================





// ============================================================================
/*  constructor with full parameters 
 *  @param mu    related to location 
 *  @param sigma related to RSM/scale/width 
 *  @param xi    related to asymmetry/skewness
 *  @param r     shape parameter 
 *  @param zeta  shape parameter     
 */
// ============================================================================
Ostap::Math::SkewGenT::SkewGenT  
( const double mu     ,    // location parameter 
  const double sigma  ,    // width parameter 
  const double psi    ,    // asymmetry/skewness parameter 
  const double r      ,    // shape parameter 
  const double zeta   )    // shape parameter 
  : m_mu     ( mu    ) 
  , m_sigma  ( sigma ) 
  , m_psi    ( psi   ) 
  , m_r      ( std::abs ( r    ) ) 
  , m_zeta   ( std::abs ( zeta ) )   
  , m_lambda ( -100  ) 
  , m_b1     ( -100  ) 
  , m_b2     ( -100  ) 
  , m_b3     ( -100  )
{
  setMu      ( mu    ) ;
  setSigma   ( sigma ) ;
  setPsi     ( psi   ) ;
  setR       ( r     ) ;
  setZeta    ( zeta  ) ;
}
// ======================================================================
// set mu-parameter
// ======================================================================
bool Ostap::Math::SkewGenT::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ======================================================================
// set sigma-parameter
// ======================================================================
bool Ostap::Math::SkewGenT::setSigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( avalue                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::SkewGenT::setSigma"        ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //              
  m_sigma = avalue ;
  return true ;
}
// ======================================================================
// set psi-parameter
// ======================================================================
bool Ostap::Math::SkewGenT::setPsi
( const double value ) 
{
  if ( s_equal ( value , m_psi )
       && -1<= m_lambda 
       &&      m_lambda <= 1 ) { return false ; }
  //
  m_psi    = value ;
  m_lambda = std::tanh ( value ) ; 
  //
  return true ;
}
// ============================================================================
/*  calculate helper math constants
 *  \f[ \left( \begin{array}{l} 
 *    b_1 \\ b_2 \\ b3
 *   \end{array}\right) = 
 *   \left( \begin{array}{l}
 *    \Beta( \frac{1}{p} , q ) \\ 
 *   \frac{ \Beta( \frac{2}{p} , q -\frac{1}{p})} { \Beta( \frac{1}{p} , q )} \ \
 *   \frac{ \Beta( \frac{3}{p} , q -\frac{3}{p})} { \Beta( \frac{1}{p} , q )} 
 *   \end{array}\right) \f]
 */
// ============================================================================
void Ostap::Math::SkewGenT::calc_b 
( double& b1 ,    // 1/B( 1/p, q ) 
  double& b2 ,    // B(2/p,q-1/p) / B ( 1/p.q) 
  double& b3 )    // B(3/p,q-2/p) / B ( 1/p.q) 
{
  //
  const double qq   = q () ;
  const double lnb1 = Ostap::Math::lnbeta ( m_r , qq ) ;
  b1 = std::exp ( - lnb1 ) ;
  b2 = std::exp ( Ostap::Math::lnbeta ( 2 * m_r , qq -     m_r ) - lnb1 ) ;
  b3 = std::exp ( Ostap::Math::lnbeta ( 3 * m_r , qq - 2 * m_r ) - lnb1 ) ;
  //
}
// ======================================================================
// set r-parameter
// ======================================================================
bool Ostap::Math::SkewGenT::setR
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_r ) 
       && -100 != m_b1  
       && -100 != m_b2  
       && -100 != m_b3 ) { return false ; }
  //
  m_r = avalue  ;
  //
  calc_b ( m_b1 , m_b2 , m_b3 ) ;
  //
  return true ;
}
// ======================================================================
// set zeta-parameter
// ======================================================================
bool Ostap::Math::SkewGenT::setZeta
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_zeta ) 
       && -100 != m_b1  
       && -100 != m_b2  
       && -100 != m_b3 ) { return false ; }
  //
  m_zeta = avalue  ;
  //
  calc_b ( m_b1 , m_b2 , m_b3 ) ;
  //
  return true ;
}
// ============================================================================
// helper scale parameter 
// ============================================================================
double Ostap::Math::SkewGenT::v_scale () const 
{
  return 1 / std::sqrt ( ( 3 * m_lambda * m_lambda + 1 ) * m_b3 
                         - 4 * m_lambda * m_lambda * m_b2 * m_b2 ) ;                       
}
// ============================================================================
/* helper bias parameter 
 *  \f$ m^{\prime} = 2 \sigma \lambda b_2 \f$
 */
// ============================================================================
double Ostap::Math::SkewGenT::m_bias  () const 
{ return 2 * m_sigma * m_lambda * m_b2 ; }
// ============================================================================
// evaluate the pdf 
// ============================================================================
double Ostap::Math::SkewGenT::pdf ( const double x ) const 
{ 
  //
  const double qq = q () ;
  const double pp = p () ;
  const double v  = v_scale () ;
  const double m  = m_bias  () * v ;
  //
  const double dx = ( x - m_mu + m ) / ( v * m_sigma ) ;
  const double t  = std::abs ( dx ) / ( m_lambda * std::copysign ( 1.0 , dx ) + 1 ) ;
  const double tp = std::pow ( t , pp ) ;
  //
  return m_b1 / ( 2 * m_sigma * v * m_r * std::pow ( tp + 1 , m_r + qq ) ) ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::SkewGenT::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::SkewGenT::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  // split into reasonable sub intervals
  //
  const double v  = v_scale ()      ;  
  const double m  = m_bias  () * v  ;
  
  const double mm = m_mu - m ;
  //
  if ( low <  mm && mm < high ) { return integral ( low , mm ) + integral ( mm , high ) ; }
  //
  {
    const double x1 = mm + 3 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm - 3 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = mm  + 5 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm  - 5 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }  
  //
  {
    const double x1 = mm + 10 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm - 10 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = mm + 15 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm -  15 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  const double x1     = mm - 15 * m_sigma  ;
  const double x2     = mm + 15 * m_sigma  ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<SkewGenT> s_integrator {} ;
  static char s_message[] = "Integral(SkewGenT)" ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high             ,    // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,           // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::SkewGenT::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::SkewGenT::skewness   () const 
{
  if ( s_zero ( m_lambda ) || s_zero ( m_psi ) ) { return 0 ; }
  //
  const double qq = q ()  ;
  const double vs = v_scale () * m_sigma ;
  const double l2 = m_lambda * m_lambda  ;
  //
  const double b4 = std::exp ( Ostap::Math::lnbeta ( 4 * m_r , qq - 3 * m_r ) - 
                               Ostap::Math::lnbeta (     m_r , qq           ) ) ;
  //
  return 
    m_lambda * std::pow ( vs , 3 ) * ( 8 * l2  * std::pow ( m_b3 , 3   ) - 
                                       3 * ( 3 * l2 + 1 ) * m_b2 * m_b3  +
                                       2 * (     l2 + 1 ) * b4 ) ;  
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::SkewGenT::tag () const 
{ 
  static const std::string s_name = "SkewGenT" ;
  return Ostap::Utils::hash_combiner ( s_name   , 
                             m_mu     , 
                             m_sigma  , 
                             m_psi    , 
                             m_r      , 
                             m_zeta   ) ;
}
// ============================================================================






// ============================================================================
/*  constructor with full parameters 
 *  @param mu    related to location 
 *  @param sigma related to RSM/scale/width 
 *  @param xi    related to asymmetry/skewness
 *  @param p     shape parameter 
 */
// ============================================================================
Ostap::Math::SkewGenError::SkewGenError  
( const double mu     ,    // location parameter 
  const double sigma  ,    // width parameter 
  const double xi     ,    // asymmetry/skewness parameter 
  const double p      )    // shape parameter 
  : m_mu     ( mu    ) 
  , m_sigma  ( sigma ) 
  , m_xi     ( xi    ) 
  , m_p      ( std::abs ( p    ) ) 
  , m_lambda ( -100  ) 
  , m_b0     ( -100  ) 
  , m_b1     ( -100  ) 
  , m_b2     ( -100  ) 
{
  setMu      ( mu    ) ;
  setSigma   ( sigma ) ;
  setXi      ( xi    ) ;
  setP       ( p     ) ;
}
// ======================================================================
// set mu-parameter
// ======================================================================
bool Ostap::Math::SkewGenError::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ======================================================================
// set sigma-parameter
// ======================================================================
bool Ostap::Math::SkewGenError::setSigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( avalue                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::SkewGenError::setSigma"    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //  
  m_sigma = avalue ;
  return true ;
}
// ======================================================================
// set xi-parameter
// ======================================================================
bool Ostap::Math::SkewGenError::setXi
( const double value ) 
{
  if ( s_equal ( value , m_xi )
       && -1<= m_lambda 
       &&      m_lambda <= 1 ) { return false ; }
  //
  m_xi     = value ;
  m_lambda = std::tanh ( value ) ; 
  //
  return true ;
}
// ============================================================================
/*  calculate helper math constants
 *  \f[ \left( \begin{array}{l} 
 *    b_0 \\ b_1 \\ b_2 
 *   \end{array}\right) = 
 *   \left( \begin{array}{l}
 *   \frac{1}{\Gamma(1/p) } \\ 
 *   \frac{\Gamma(3/p)    }{\Gamma^3(1/p)} \\ 
 *   2^{2/p}\frac{\Gamma(1/2+1/p)}{\Gamma(1/p)  } 
 *   \end{array}\right) \f]
 */
// ============================================================================
void Ostap::Math::SkewGenError::calc_b 
( double& b9 ,    // 1/Gamma(1/p) 
  double& b1 ,    // Gamma(3/p)/Gamma(1/p) 
  double& b2 )    // 2^{2/p} Gamma(1/2+1/p)/Gamma(1/p)
{
  //
  const long double ip  = 1.0L / m_p ;
  const long double lg1 = std::lgamma ( ip ) ;
  //
  m_b0 = Ostap::Math::igamma ( ip ) ;
  m_b1 = std::exp ( std::lgamma ( 3 * ip ) - 3 * lg1 ) ;
  m_b2 = std::exp ( 2 * ip * s_ln2 + std::lgamma ( 0.5L + ip ) - lg1 ) ;
  //
}
// ============================================================================
// set p-parameter
// ============================================================================
bool Ostap::Math::SkewGenError::setP
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_p ) 
       && -100 != m_b0  
       && -100 != m_b1  
       && -100 != m_b2 ) { return false ; }
  //
  m_p = avalue  ;
  //
  calc_b ( m_b0 , m_b1 , m_b2 ) ;
  //
  return true ;
}
// ============================================================================
// helper scale parameter 
// ============================================================================
double Ostap::Math::SkewGenError::v_scale () const 
{
  const double l2 = m_lambda * m_lambda ;
  return std::sqrt ( M_PI / ( M_PI * ( 1.0 + 3.0 * l2 ) * m_b1 - l2 * m_b2 * m_b2 ) ) ;
}
// ============================================================================
/* helper bias parameter 
 *  \f$ m^{\prime} = \sigma \lambda b_2 \f$
 */
// ============================================================================
double Ostap::Math::SkewGenError::m_bias  () const 
{ return m_sigma * m_lambda * m_b2 * s_SQRTPIi ; }
// ============================================================================
// evaluate the pdf 
// ============================================================================
double Ostap::Math::SkewGenError::pdf ( const double x ) const 
{ 
  //
  const double vp = v_scale () ;
  const double mp = m_bias  () * vp ;
  //
  const double dx = ( x - m_mu + mp ) / ( vp * m_sigma * m_b0 ) ;
  const double t  = std::abs ( dx ) / ( m_lambda * std::copysign ( 1.0 , dx ) + 1 ) ;
  const double tp = std::pow ( t , m_p ) ;
  //
  return m_p * std::exp ( -tp )  / ( 2 * m_sigma * vp ) ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::SkewGenError::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::SkewGenError::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  // split into reasonable sub intervals
  //
  const double v  = v_scale ()      ;  
  const double m  = m_bias  () * v  ;
  
  const double mm = m_mu - m ;
  //
  if ( low <  mm && mm < high ) { return integral ( low , mm ) + integral ( mm , high ) ; }
  //
  {
    const double x1 = mm + 3 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm - 3 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = mm  + 5 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm  - 5 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }  
  //
  {
    const double x1 = mm + 10 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm - 10 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = mm + 15 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = mm -  15 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  const double x1     = mm - 15 * m_sigma  ;
  const double x2     = mm + 15 * m_sigma  ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<SkewGenError> s_integrator {} ;
  static char s_message[] = "Integral(SkewGenError)" ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high             ,                  // low & high edges
      workspace ( m_workspace ) ,                  // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()        ,                  // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::SkewGenError::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::SkewGenError::tag () const 
{ 
  static const std::string s_name = "SkewGenError" ;
  return Ostap::Utils::hash_combiner ( s_name   , 
                                       m_mu     , 
                                       m_sigma  , 
                                       m_xi     , 
                                       m_p      ) ;
}
// ============================================================================





// ===========================================================================
// constructor with location and scale parmaeters 
// ===========================================================================
Ostap::Math::Hat::Hat
( const double mu       , 
  const double varsigma ) 
  : m_mu       ( mu )
  , m_varsigma ( std::abs ( varsigma ) ) 
{}
// ============================================================================
// set mu
// ============================================================================
bool Ostap::Math::Hat::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// set varsigma
// ============================================================================
bool Ostap::Math::Hat::setVarsigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_varsigma ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Hat::evaluate 
( const double x ) const
{
  static const double s_norm = 1/0.443993816168079313833061405603 ;
  const double z = ( x - m_mu ) / m_varsigma ;
  return 1 <= std::abs ( z ) ? 0.0 : Ostap::Math::hat ( z )  * s_norm / m_varsigma ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::Hat::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::Hat::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  const double mn = ( low  - m_mu ) / m_varsigma ;
  const double mx = ( high - m_mu ) / m_varsigma ;
  //
  if      ( mx <= -1            ) { return 0 ; }
  else if ( mn >=  1            ) { return 0 ; }
  else if ( mn <= -1 && mx >= 1 ) { return 1 ; }
  //
  const double xmn = std::max ( low  , m_mu - m_varsigma ) ;
  const double xmx = std::min ( high , m_mu + m_varsigma ) ;
  //
   // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Hat> s_integrator {} ;
  static char s_message[] = "Integral(Hat)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     ,  
      xmn    , xmx              , // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size ()       , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
}
// ============================================================================
// get the variance of the distribution 
// ============================================================================
double Ostap::Math::Hat::variance () const 
{ return m_varsigma * m_varsigma * 0.15811363626379668 ; }
// ============================================================================
// get the RMS of the distribution 
// ============================================================================
double Ostap::Math::Hat::rms      () const 
{ return m_varsigma * 0.3976350541184676 ; }
// ============================================================================
// get the (excess) kurtosis 
// ============================================================================
double Ostap::Math::Hat::kurtosis () const 
{ return -0.8807206646393597 ; }
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::Hat::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Hat::tag () const 
{ 
  static const std::string s_name = "Hat" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma ) ;
}
// ============================================================================


// ===========================================================================
// constructor with location and scale parmaeters 
// ===========================================================================
Ostap::Math::Up::Up
( const double mu       , 
  const double varsigma ) 
  : m_mu       ( mu )
  , m_varsigma ( std::abs ( varsigma ) ) 
{}
// ============================================================================
// set mu
// ============================================================================
bool Ostap::Math::Up::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// set varsigma
// ============================================================================
bool Ostap::Math::Up::setVarsigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_varsigma ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Up::evaluate 
( const double x ) const
{
  const double z = ( x - m_mu ) / m_varsigma ;
  return 1 <= std::abs ( z ) ? 0.0 : eval ( z ) / m_varsigma ;
}
// ===========================================================================
// evaluate the "standard" up function 
// ===========================================================================
double Ostap::Math::Up::eval ( const double z )  const 
{
  //
  auto fourrier = []( unsigned int k ) -> double  
  { return k == 0 || 1 == k % 2 ?  Ostap::Math::up_F ( M_PI * k ) : 0.0 ; } ;
  static const std::array<double,120> s_fourrier { make_array ( fourrier , std::make_index_sequence<120>() ) } ;
  //
  return 1 <= std::abs ( z ) ? 0.0 : 
    std::max ( 0.0L , Ostap::Math::Clenshaw::cosine_sum 
               ( s_fourrier.begin () , s_fourrier.end () , z * M_PI ) ) ; 
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::Up::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::Up::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  const double mn = ( low  - m_mu ) / m_varsigma ;
  const double mx = ( high - m_mu ) / m_varsigma ;
  //
  if      ( mx <= -1            ) { return 0 ; }
  else if ( mn >=  1            ) { return 0 ; }
  else if ( mn <= -1 && mx >= 1 ) { return 1 ; }
  //
  const double xmn = std::max ( low  , xmin () ) ;
  const double xmx = std::min ( high , xmax () ) ;
  //
   // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Up> s_integrator {} ;
  static char s_message[] = "Integral(Up)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     ,  
      xmn    , xmx              , // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size ()       , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::Up::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the variance of the distribution 
// ============================================================================
double Ostap::Math::Up::variance () const 
{ return m_varsigma * m_varsigma / 9 ;}
// ============================================================================
// get the RMS of the distribution 
// ============================================================================
double Ostap::Math::Up::rms      () const 
{ return m_varsigma / 3 ; }
// ============================================================================
// get the (excess) kurtosis 
// ============================================================================
double Ostap::Math::Up::kurtosis () const 
{
  static const double s_kurtosis = 19.0 * 9 * 9 / ( std::pow ( 3 , 3 ) * 5 * 5 ) - 3 ;
  return s_kurtosis ;
}
// ============================================================================
// get the value of the derivative 
// ============================================================================
double Ostap::Math::Up::derivative ( const double x ) const 
{
  const double z = ( x - m_mu ) / m_varsigma ;
  return 1 <= std::abs ( z )  ? 0.0 : 
    2 * ( eval ( 2 * z + 1 )  - eval ( 2 * z - 1 ) ) / m_varsigma ;
}  
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Up::tag () const 
{ 
  static const std::string s_name = "Up" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma ) ;
}
// ============================================================================


// ============================================================================
namespace 
{  
  // ==========================================================================
  typedef std::array<double,120>           RESULT ;
  typedef std::map<unsigned short,RESULT>  MAP    ;
  typedef SyncedCache<MAP>                 CACHE  ;
  // =========================================================================
  CACHE s_FupN_cache {} ;
  // =========================================================================
}
// ===========================================================================
// constructor with location and scale parmaeters 
// ===========================================================================
Ostap::Math::FupN::FupN
( const unsigned short N        , 
  const double         mu       , 
  const double         varsigma ) 
  : m_N        ( N  )
  , m_mu       ( mu )
  , m_varsigma ( std::abs ( varsigma ) ) 
{
  CACHE::Lock lock { s_FupN_cache.mutex() } ;
  auto it = s_FupN_cache->find  ( m_N  ) ;
  if ( s_FupN_cache->end() == it ) 
  {
    auto fourrier = [this]( unsigned int k ) -> double
      { return Ostap::Math::fupN_F ( this->m_N , M_PI * k / ( this->m_N + 1 ) ) ; } ;
    //
    const RESULT res { make_array( fourrier , std::make_index_sequence<120>() ) } ;
    s_FupN_cache->insert ( std::make_pair ( m_N , res ) ) ;  
  }
}
// ============================================================================
// set mu
// ============================================================================
bool Ostap::Math::FupN::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// set varsigma
// ============================================================================
bool Ostap::Math::FupN::setVarsigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_varsigma ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::FupN::evaluate 
( const double x ) const
{
  const double z = ( x - m_mu ) / m_varsigma ;
  return 0.5 * ( m_N + 2 )  <= std::abs ( z ) ? 0.0 : eval ( z ) / m_varsigma ;
}
// ===========================================================================
// evaluate the "standard" fupN function 
// ===========================================================================
double Ostap::Math::FupN::eval ( const double z )  const 
{
  //
  auto it = s_FupN_cache->find  ( m_N  ) ;
  Ostap::Assert ( s_FupN_cache->end() != it , 
                  "Cache does not exist!"   , 
                  "Ostap::Math::FupN"       ) ;
  //
  return 0.5 * ( m_N + 2 ) <= std::abs ( z ) ? 0.0 : 
    std::max ( 0.0L , Ostap::Math::Clenshaw::cosine_sum 
               ( it->second.begin () , 
                 it->second.end   () ,  M_PI * z / ( m_N + 1 ) ) / ( m_N + 1 ) ) ; 
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::FupN::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::FupN::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  const double mn = ( low  - m_mu ) / m_varsigma ;
  const double mx = ( high - m_mu ) / m_varsigma ;
  //
  const double nn = 0.5 * ( m_N + 2 ) ;
  //
  if      ( mx <= -1               ) { return 0 ; }
  else if ( mn >=  1               ) { return 0 ; }
  else if ( mn <= -nn  && mx >= nn ) { return 1 ; }
  //
  const double xmn = std::max ( low  , xmin () ) ;
  const double xmx = std::min ( high , xmax () ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<FupN> s_integrator {} ;
  static char s_message[] = "Integral(FupN)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     ,  
      xmn    , xmx              , // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size ()       , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::FupN::tag () const 
{ 
  static const std::string s_name = "FupN" ;
  return Ostap::Utils::hash_combiner ( s_name , m_N , m_mu , m_varsigma ) ;
}
// ============================================================================

                  
// ===========================================================================
// Meixner distribution
// =========================================================================== 
Ostap::Math::Meixner::Meixner
( const double mu    , // location
  const double sigma , // sigma 
  const double psi   , // b = 2 * atan ( psi )
  const double shape ) // d    
  : m_mu    (  0 ) 
  , m_sigma (  1 ) 
  , m_psi   (  0 )
  , m_shape (  1 )
  , m_a     (  1 )
  , m_b     (  0 )
  , m_C     ( -1 ) 
{
  setMu    ( mu    ) ; 
  setSigma ( sigma ) ;
  setPsi   ( psi   ) ; 
  setShape ( shape ) ; 
  //
  m_C = 2 * m_shape * std::log ( 2 * std::cos ( 0.5 * m_b ) ) 
    - std::lgamma ( 2 * m_shape ) - std::log ( 2 * M_PI ) ;
  //
  m_a = m_sigma * std::sqrt ( ( std::cos ( m_b ) + 1 ) / m_shape ) ;
}
// ============================================================================
// set Mu
// ============================================================================
bool  Ostap::Math::Meixner::setMu
( const double value ) 
{
  if ( s_equal( value , m_mu) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// set Sigma 
// ============================================================================
bool Ostap::Math::Meixner::setSigma 
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  //
  Ostap::Assert ( avalue                                   ,
		  "Parameter 'sigma' must be non-zero"     ,
		  "Ostap::Math::Meixner::setSigma"         ,
		  INVALID_PARAMETER , __FILE__ , __LINE__  ) ;
  //             
  m_sigma = avalue ;
  //
  m_a  = m_sigma * std::sqrt ( ( std::cos ( m_b ) + 1 ) / m_shape ) ;
  return true ;
} 
// ============================================================================
// set Shape 
// ============================================================================
bool Ostap::Math::Meixner::setShape 
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_shape ) ) { return false ; }
  m_shape = avalue ;
  //
  m_C = 2 * m_shape * std::log ( 2 * std::cos ( 0.5 * m_b ) )
    - std::lgamma ( 2 * m_shape ) - std::log ( 2 * M_PI ) ;
  //
  m_a = m_sigma * std::sqrt ( ( std::cos ( m_b ) + 1 ) / m_shape ) ;
  return true ;
} 
// ============================================================================
// set psi 
// ============================================================================
bool Ostap::Math::Meixner::setPsi  
( const double value ) 
{
  if ( s_equal ( value , m_psi ) ) { return false ; }
  m_psi = value ;
  if  ( s_zero ( m_psi ) ) { m_psi = 0 ; }
  //
  m_b = m_psi ? 2 * std::atan ( m_psi ) : 0.0 ; 
  //
  m_C = 2 * m_shape * std::log ( 2 * std::cos ( 0.5 * m_b ) ) 
    - std::lgamma ( 2 * m_shape ) - std::log ( 2 * M_PI ) ;
  //
  m_a = m_sigma * std::sqrt ( ( std::cos ( m_b ) + 1 ) / m_shape ) ;
  return true ;
} 
// ============================================================================
// evaluate Meixner function
// ============================================================================
double Ostap::Math::Meixner::evaluate   ( const double x ) const 
{
  const double z  = ( x - m_mu ) /  m_a ;
  const std::complex<double> v { m_shape , z } ; 
  const double r = std::real ( Ostap::Math::lgamma ( v ) ) ;
  //
  const double f = m_C + m_b * z + 2 * r ; 
  return std::exp ( f ) / m_a  ;
}
// ============================================================================
// kappa 
// ============================================================================
double Ostap::Math::Meixner::kappa () const 
{ return m_b / M_PI ; }
// ============================================================================
// mean values  
// ============================================================================
double Ostap::Math::Meixner::mean() const
{ return m_psi ? m_mu + m_a * m_shape + std::tan ( 0.5 * m_b ) : m_mu ; }
// ============================================================================
//  skewness 
// =============================================================================
double Ostap::Math::Meixner::skewness  () const 
{ return m_psi ? std::sin ( m_b ) / std::sqrt ( m_shape * ( std::cos ( m_b ) + 1 ) ) : 0.0 ; }
// ============================================================================
//  (excess) kurtosis
// ============================================================================
double Ostap::Math::Meixner::kurtosis() const
{ return ( 2 - std::cos ( m_b ) ) / m_shape  ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::Meixner::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::Meixner::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  if  ( low < m_mu && m_mu < high )
    { return integral ( low , m_mu   ) + integral ( m_mu   , high ) ; }
  //
  const double x_mean = mean () ;
  //
  if ( m_psi && low < x_mean && x_mean < high )
    { return integral ( low , x_mean ) + integral ( x_mean , high ) ; }
  //
  for  ( unsigned int  j = 1 ; j <= 5 ; j += 2  )
    {
      const double x1 = std::max ( m_mu , x_mean ) + j * m_sigma ; 
      if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
      const double x2 = std::min ( m_mu , x_mean ) - j * m_sigma  ; 
      if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Meixner> s_integrator {} ;
  static char s_message[] = "Integral(Meixner)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     ,  
      low    , high             , // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size ()       , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
}
// ============================================================================
/*  quantify the effect of the tails, difference from Gaussian
 *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
 * where 
 * - \f$ I_{CB} \f$ is integral over Gaussian function 
 * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
 * - Gaussian is centered at mean-valuw with sigma = RMS 
 */
// ============================================================================
double Ostap::Math::Meixner::non_gaussian 
( const double xlow  ,
  const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; } // convention
  else if ( xhigh < xlow             ) { return -non_gaussian ( xhigh , xlow ) ; } 
  //
  const double I_CB = integral ( xlow , xhigh ) / integral ();
  //
  const double m = mean () ;
  const double s = rms  () ;
  //
  const double I_G  =
    Ostap::Math::gauss_cdf ( xhigh , m , s ) -
    Ostap::Math::gauss_cdf ( xlow  , m , s ) ;
  //
  return 1 - I_G / I_CB ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Meixner::tag () const 
{ 
  static const std::string s_name = "Mexner" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_sigma , m_psi , m_shape ) ;
}
// ============================================================================
// Asymptotic 
// ===========================================================================
double Ostap::Math::Meixner::rho         () const { return 2 * d() - 1 ; }
double Ostap::Math::Meixner::sigma_plus  () const { return ( M_PI + m_b ) / m_a ; }
double Ostap::Math::Meixner::sigma_minus () const { return ( M_PI - m_b ) / m_a ; }


// ============================================================================
//                                                                      The END 
// ============================================================================
