// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Peaks.h"
#include "Ostap/MoreMath.h"
// ============================================================================
//  Local
// ============================================================================
#include "local_math.h"
#include "local_gsl.h"
#include "local_hash.h"
#include "gauss.h"      
#include "Integrator1D.h"      
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

  double x_sinh ( const double x , double precision = s_PRECISION )
  {
    //
    if      ( s_equal   ( x , 0 )    ) { return 1 ; }
    else if ( std::fabs ( x ) < 0.1  )
    {
      double result = 1.0  ;
      double delta  = x    ;
      //
      precision = std::fabs ( precision ) ;
      precision = std::min  ( precision , std::fabs ( s_PRECISION_TAIL ) ) ;
      unsigned int n = 1 ;
      //
      do
      {
        delta  *= x * x / ( ( n + 1 )  * ( n + 2 ) ) ;
        result += delta ;
        n      += 2 ;
      }
      while ( std::fabs ( delta ) > 0.1 * precision && n < 10000 ) ;
      //
      return result ;
    }
    //
    if ( 100 < std::fabs ( x )  ) { return  s_INFINITY ; }
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
  // ==========================================================================
  // Studen-t'
  // ==========================================================================
  inline double student_cdf (  const double t , const double nu ) 
  {
    const double xt    = nu / ( t * t + nu ) ;
    const double value = 0.5 * gsl_sf_beta_inc ( 0.5 * nu , 0.5 , xt ) ;
    return t >= 0 ? 1 - value : value ;
  }
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
  , m_sigmaL  ( std::fabs ( sigmaL ) )
  , m_sigmaR  ( std::fabs ( sigmaR ) )
//
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::BifurcatedGauss::~BifurcatedGauss (){}
// ============================================================================
// evaluate Bifurcated Gaussian
// ============================================================================
double Ostap::Math::BifurcatedGauss::evaluate ( const double x ) const
{
  const double dx   = x - m_peak ;
  const double norm = s_SQRTPIHALF * ( sigmaL() + sigmaR() ) ;
  //
  return
    dx < 0 ?
    my_exp ( -0.5 * dx * dx / ( sigmaL () * sigmaL () ) ) / norm :
    my_exp ( -0.5 * dx * dx / ( sigmaR () * sigmaR () ) ) / norm ;
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
    const double sf    = s_SQRT2i  / sigma  ;
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
double Ostap::Math::BifurcatedGauss::integral ( const double low  ,
                                                const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }                         // RETURN
  if (           low > high   ) { return - integral ( high , low ) ; } // RETURN
  //
  // left half-gaussian
  if       ( high <= m_peak )
  {
    const double sigma = sigmaL () ;
    const double sf    = s_SQRT2i  / sigma  ;
    const double nf    = sigma    / ( sigmaL() + sigmaR() ) ;
    const double a     = ( low  - m_peak ) * sf ;
    const double b     = ( high - m_peak ) * sf ;
    return  ( std::erf ( b ) -  std::erf ( a ) ) * nf ; // RETURN
  }
  // right half-gaussian
  else if ( low >= m_peak )
  {
    const double sigma = sigmaR () ;
    const double sf    = s_SQRT2i  / sigma  ;
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
{ return std::hash_combine ( m_peak , m_sigmaL , m_sigmaR ) ; }
// ============================================================================

// ============================================================================
bool Ostap::Math::BifurcatedGauss::setSigmaL ( const double value )
{
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( m_sigmaL , value_ ) ) { return false ; }
  m_sigmaL = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setSigmaR ( const double value )
{
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( m_sigmaR , value_ ) ) { return false ; }
  m_sigmaR = value_ ;
  //
  return true ;
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
{ return std::hash_combine ( m_peak , m_sigma , m_fraction , m_scale ) ; }
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
  , m_sigma ( std::fabs ( sigma ) )
    //
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Gauss::~Gauss (){}
// ============================================================================
// evaluate Bifurcated Gaussian
// ============================================================================
double Ostap::Math::Gauss::evaluate ( const double x ) const
{
  const double dx   = ( x - m_peak ) / m_sigma ;
  const double norm = s_SQRTPIHALF * m_sigma  ;
  //
  return my_exp ( -0.5 * dx * dx ) / norm ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Gauss::integral () const { return 1 ; }
// ============================================================================
//  get CDF 
// ============================================================================
double Ostap::Math::Gauss::cdf ( const double x ) const 
{
  const double dx =  s_SQRT2i * ( x - m_peak ) / m_sigma ;
  return 0.5 * ( 1 + std::erf ( dx ) ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Gauss::integral ( const double low  ,
                                      const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }                         // RETURN
  //
  const double c = s_SQRT2i / m_sigma ;
  const double l = c * ( low  - m_peak ) ;
  const double h = c * ( high - m_peak ) ;
  //
  return 0.5 * ( std::erf ( h ) - std::erf ( l ) ) ;
}
// ============================================================================
bool Ostap::Math::Gauss::setSigma ( const double value )
{
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( m_sigma , value_ ) ) { return false ; }
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
{ return std::hash_combine ( m_peak , m_sigma ) ; }
// ============================================================================





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
// destructor 
// ============================================================================
Ostap::Math::GenGaussV1::~GenGaussV1(){}
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
{ return std::hash_combine ( m_mu , m_alpha , m_beta ) ; }
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
// destructor
// ============================================================================
Ostap::Math::GenGaussV2::~GenGaussV2(){}
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
{ return std::hash_combine ( m_xi , m_alpha , m_kappa ) ; }
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
// desctructor
// ============================================================================
Ostap::Math::SkewGauss::~SkewGauss(){}
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
  return m_omega * m_omega * ( 1 - 2 * delta * delta , M_PI ) ;
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
{ return std::hash_combine ( m_xi , m_omega , m_alpha ) ; }
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
// destructor
// ============================================================================
Ostap::Math::Bukin::~Bukin () {}
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
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
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
  m_B2    = 1/x_log ( beta )  ;
  m_B2   *= m_B2              ;
  m_B2   *= ab*ab             ;
  //
  //
  const double delta = xi + xi2sqrt - 1 ;
  const double tail  =
    0.5 * s_Bukin * xi2sqrt * ( 1 + xi + xi2sqrt ) / ( xi + xi2sqrt ) / x_log ( delta ) ;
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
  const double value_ = std::fabs ( value ) ;
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
  const double A    = x_log ( m_A * dx ) ;
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
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag () , 
      &F     ,  
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
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
  std::tie ( ierror1 , result1 , error1 ) = s_integrator.gaqil_integrate_with_cache
    ( tag () , 
      &F     , 
      m_x1   ,                       // low edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION_TAIL    ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message1          , 
      __FILE__ , __LINE__ ) ;
  //
  int    ierror2  =  0 ;
  double result2  =  1 ;
  double error2   = -1 ;
  std::tie ( ierror2 , result2 , error2 ) = s_integrator.gaqiu_integrate_with_cache
    ( tag () , 
      &F     , 
      m_x2   ,                       // high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION_TAIL    ,          // relative precision
      m_workspace.size()              ,          // size of workspace
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
{ return std::hash_combine ( m_peak , m_sigma , m_xi , m_rho_L , m_rho_R ) ; }
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
( const double m0         ,
  const double sigma      ,
  const double tau        )
  : m_m0        ( m0                   )
  , m_sigma     ( std::fabs ( sigma )  )
  , m_tau       ( std::tanh ( tau   )  )
//
  , m_lambda    ( 0.0   )
//
  , m_integral  ( -1000 )
  , m_workspace ()
{
  //
  m_lambda = x_sinh ( m_tau * s_Novosibirsk ) ;
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Novosibirsk::~Novosibirsk() {}
// ============================================================================
// set parameter m0
// ============================================================================
bool Ostap::Math::Novosibirsk::setM0    ( const double value )
{
  //
  if ( s_equal ( m_m0 ,  value ) ) { return false ; }
  //
  m_m0 = value ;
  //
  return true ;
}
// ============================================================================
// set parameter sigma
// ============================================================================
bool Ostap::Math::Novosibirsk::setSigma ( const double value )
{
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma    = value_ ;
  m_integral = -1000 ;
  //
  return true ;
}
// ============================================================================
// set parameter tau
// ============================================================================
bool Ostap::Math::Novosibirsk::setTau ( const double value )
{
  //
  const double value_ = std::tanh ( value )   ;
  if ( s_equal ( value_ , m_tau ) ) { return false ; }
  //
  m_tau      = value_ ;
  m_integral = -1000 ;
  //
  m_lambda   = x_sinh ( m_tau * s_Novosibirsk ) ;
  //
  return true ;
}
// ============================================================================
// evaluate Novosibirsk's function
// ============================================================================
double Ostap::Math::Novosibirsk::pdf  ( const double x ) const
{
  //
  const double dx     = ( x - m_m0 ) / m_sigma ;
  //
  const double arg    = m_lambda * dx * m_tau ;
  //
  if ( arg <= -1 || s_equal ( arg , -1 ) ) { return 0 ; } // RETURN
  //
  const double l      =  x_log ( arg ) * m_lambda * dx ;
  //
  const double result = l * l  + m_tau * m_tau ;
  //
  return  my_exp ( -0.5 * result ) ;
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
  const double x1     = m_m0 - 10 * m_sigma  ;
  const double x2     = m_m0 + 10 * m_sigma  ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  if      ( low < x_low  && x_low < high )
  {
    return
      integral (   low , x_low  ) +
      integral ( x_low ,   high ) ;
  }
  else if ( low <  x_high && x_high < high )
  {
    return
      integral (   low  , x_high  ) +
      integral ( x_high ,   high  ) ;
  }
  //
  // split, if the interval is too large
  //
  const double width = std::max ( std::abs  ( m_sigma )  , 0.0 ) ;
  if ( 0 < width &&  3 * width < high - low  )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Novosibirsk> s_integrator {} ;
  static char s_message[] = "Integral(Novosibirsk)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag () , 
      &F     , 
      low    , high             ,    // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL :
      s_PRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
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
  if ( m_integral <= 0 )
  {
    Novosibirsk* novosibirsk = const_cast<Novosibirsk*> ( this ) ;
    novosibirsk -> integrate() ;
  }
  //
  return m_integral ;
}
// ============================================================================
// calculate  the integral
// =========================================================================
void Ostap::Math::Novosibirsk::integrate()
{
  //
  const double x1     = m_m0 - 10 * m_sigma ;
  const double x2     = m_m0 + 10 * m_sigma ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
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
  std::tie ( ierror1 , result1 , error1 ) = s_integrator.gaqil_integrate_with_cache
    ( tag () , 
      &F     , 
      x_low  ,                       // low edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION_TAIL    ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message1          , 
      __FILE__ , __LINE__ ) ;
  //
  int    ierror2  =  0 ;
  double result2  =  1 ;
  double error2   = -1 ;
  std::tie ( ierror2 , result2 , error2 ) = s_integrator.gaqiu_integrate_with_cache 
    ( tag () , 
      &F     , 
      x_high ,                       // high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION_TAIL    ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message2          , 
      __FILE__ , __LINE__ ) ;
  //
  m_integral = result1 + result2 + integral ( x_low ,  x_high ) ;
  //
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Novosibirsk::tag () const 
{ return std::hash_combine ( m_m0 , m_sigma , m_tau ) ; }
// ============================================================================


// ============================================================================
// Crystal Ball & Co
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
  : m_m0         ( m0 )
  , m_sigma      (  1 )
  , m_alpha      (  2 )
  , m_n          (  2 )
//
  , m_A  ( -1000 ) 
  , m_B  ( -1000 ) 
  , m_C  ( -1000 ) 
{
  //
  setM0     ( m0     ) ;
  setAlpha  ( alpha  ) ;
  setSigma  ( sigma  ) ;
  setN      ( n      ) ;
  //
  m_A = my_exp ( -0.5 * m_alpha * m_alpha ) ;
  m_B =  0.5 * ( 1 + std::erf ( - m_alpha * s_SQRT2i ) ) ;
  if   ( !s_equal ( m_n , 0 ) && !s_equal ( m_alpha , 0 )  ) 
  { m_C  = ( m_n + 1 )  / aa ()  / m_n  * s_SQRT2PIi ; }
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::CrystalBall::~CrystalBall (){}
// ============================================================================
bool  Ostap::Math::CrystalBall::setM0 ( const double value )
{
  //
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBall::setSigma ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBall::setAlpha  ( const double value )
{
  //
  if ( s_equal ( value , m_alpha ) ) { return false ; }
  //
  m_alpha    = value  ;
  //
  m_A        = my_exp ( -0.5 * alpha () * alpha ()) ;
  //
  // 
  if   ( s_equal ( n () , 0 ) || s_equal ( m_alpha , 0 )  ) { m_C = -1000 ; }
  else { m_C  = np1 () / aa () / n()  * s_SQRT2PIi ; }
  //
  m_B  = 0.5 * ( 1 + std::erf ( - m_alpha * s_SQRT2i ) ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBall::setN      ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n     ) ) { return false ; }
  //
  m_n        = value_ ;
  if ( s_equal ( m_n , 0 ) ) { m_n = 0 ; }
  //
  if   ( s_equal ( n () , 0 ) || s_equal ( m_alpha , 0 )  ) { m_C = -1000 ; }
  else { m_C  = np1()  / aa () / n() * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBall::pdf ( const double x ) const
{
  //
  const double dx    = ( x - m_m0 ) / m_sigma ;
  //
  // the tail
  //
  if  ( dx < -m_alpha )
  {
    const double frac = np1 () / ( np1 () - aa() * ( m_alpha + dx ) )  ;
    return std::pow ( frac , np1 () ) * m_A * s_SQRT2PIi / sigma() ;
  }
  //
  // the peak
  //
  return my_exp ( -0.5 * dx * dx ) * s_SQRT2PIi / sigma() ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBall::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  const double x0 = m_m0 - m_alpha * m_sigma ;
  //
  // split into proper subintervals
  //
  if      ( low < x0 && x0 < high )
  { return integral ( low , x0 ) + integral ( x0 , high ) ; }
  //
  // Z = (x-x0)/sigma 
  //
  const double zlow  = ( low  - m_m0 ) / sigma() ;
  const double zhigh = ( high - m_m0 ) / sigma() ;
  //
  // peak
  //
  if ( x0 <= low  )
  { return s_SQRT2PIi * details::gaussian_int ( 0.5   , 0 , zlow  , zhigh ) ; }
  //
  // tail
  //
  const double A =   np1 () ;
  const double B =   np1 () ;
  const double C = - aa  () ;
  //
  const double result = s_SQRT2PIi * m_A * 
    tail_integral ( A , B , C , np1 () , zlow + alpha() , zhigh + alpha() ) ;
  //
  return result ;
}
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::CrystalBall::integral () const
{
  /// the regular case 
  if ( 0 < m_C ) { return m_C + m_B ; }
  //
  /// trunkate it! 
  const double left = ( 0 < m_alpha ) ?  (-m_alpha-s_TRUNC) : -s_TRUNC ;
  // 
  return m_B + integral ( m0 () + left     * sigma() , 
                          m0 () - alpha () * sigma() ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBall::tag () const 
{ return std::hash_combine ( m_m0 , m_sigma , m_alpha , m_n ) ; }
// ============================================================================


// ============================================================================
// Needham function
// ============================================================================
/* constructor from all parameters
 *  @param m0     m0       parameter
 *  @param sigma  sigma    parameter
 *  @param a0     a0       parameter
 *  @param a1     a1       parameter
 *  @param a2     a2       parameter
 */
// ============================================================================
Ostap::Math::Needham::Needham
( const double m0    ,
  const double sigma ,
  const double a0    ,
  const double a1    ,
  const double a2    )
/// @see Ostap::Math:CrystalBall
  : m_cb  ( m0 , sigma , 1 , 0 ) // Ostap::Math:CrystalBall
  , m_a0  ( std::abs ( a0 )  )
  , m_a1  (            a1    )
  , m_a2  (            a2    )
{
  m_cb.setAlpha ( alpha () ) ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Needham::~Needham(){}
// ============================================================================
bool Ostap::Math::Needham::setA0 ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_a0 ) ) { return false ; }
  m_a0      = value_ ;
  return m_cb.setAlpha ( alpha () ) ;
}
// ============================================================================
bool Ostap::Math::Needham::setA1 ( const double value )
{
  if ( s_equal ( value , m_a1 ) ) { return false ; }
  m_a1 = value ;
  return m_cb.setAlpha ( alpha () ) ;
}
// ============================================================================
bool Ostap::Math::Needham::setA2 ( const double value )
{
  if ( s_equal ( value , m_a2 ) ) { return false ; }
  m_a2 = value ;
  return m_cb.setAlpha ( alpha () ) ;
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
{ return std::hash_combine ( m_cb.tag() ,  m_a0 , m_a1 , m_a2 ) ; }
// ============================================================================


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
  : m_cb         ( m0 , sigma , alpha , n ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::CrystalBallRightSide::~CrystalBallRightSide (){}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallRightSide::pdf ( const double x ) const
{
  const double y = 2 * m0 ()  - x ;
  //
  return  m_cb.pdf ( y ) ;  
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallRightSide::integral
( const double low ,
  const double high ) const
{ return m_cb.integral ( 2 * m0 () - high  , 2 * m0 () - low ) ; }
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::CrystalBallRightSide::integral () const 
{ return m_cb.integral () ; }
// =========================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallRightSide::tag () const 
{ return std::hash_combine ( m_cb.tag() ,  -1 ) ; }
// ============================================================================


// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSided::CrystalBallDoubleSided
( const double m0      ,
  const double sigma   ,
  const double alpha_L ,
  const double n_L     ,
  const double alpha_R ,
  const double n_R     )
  : m_m0         (  m0 )
  , m_sigma      (   1 )
  , m_alpha_L    (   2 )
  , m_n_L        (   2 )
  , m_alpha_R    (   2 )
  , m_n_R        (   2 )
//
  , m_AL         ( -1000 ) 
  , m_AR         ( -1000 )
  , m_B          ( -1000 ) 
  , m_TL         ( -1000 ) 
  , m_TR         ( -1000 ) 
{
  //
  setM0       ( m0      ) ;
  setSigma    ( sigma   ) ;
  setAlpha_L  ( alpha_L ) ;
  setAlpha_R  ( alpha_R ) ;
  setN_L      ( n_L     ) ;
  setN_R      ( n_R     ) ;
  //
  m_AL = my_exp ( -0.5 * m_alpha_L * m_alpha_L ) ;
  m_AR = my_exp ( -0.5 * m_alpha_R * m_alpha_R ) ;
  m_B  = 0.5 *  ( std::erf (  m_alpha_R * s_SQRT2i ) - 
                  std::erf ( -m_alpha_L * s_SQRT2i ) ) ;
  //
  if   ( !s_equal ( m_n_L , 0 ) && !s_equal ( m_alpha_L , 0 )  ) 
  { m_TL  = ( m_n_L + 1 )  / std::abs ( m_alpha_L )  / m_n_L  * s_SQRT2PIi ; }
  if   ( !s_equal ( m_n_R , 0 ) && !s_equal ( m_alpha_R , 0 )  ) 
  { m_TR  = ( m_n_R + 1 )  / std::abs ( m_alpha_R )  / m_n_R  * s_SQRT2PIi ; }
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::CrystalBallDoubleSided::~CrystalBallDoubleSided(){}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setM0 ( const double value )
{
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setSigma ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setAlpha_L ( const double value )
{
  if ( s_equal ( value , m_alpha_L ) ) { return false ; }
  //
  m_alpha_L  = value  ;
  m_AL       = my_exp ( -0.5 * m_alpha_L * m_alpha_L ) ;
  m_B        = 0.5 *  ( std::erf (  m_alpha_R * s_SQRT2i ) - 
                        std::erf ( -m_alpha_L * s_SQRT2i ) ) ;
  //
  if   ( s_equal ( m_n_L , 0 ) || s_equal  ( m_alpha_L , 0 )  ) {  m_TL = -1000 ; }
  else { m_TL  = ( m_n_L + 1 )  / std::abs ( m_alpha_L )  / m_n_L  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setAlpha_R ( const double value )
{
  if ( s_equal ( value , m_alpha_R ) ) { return false ; }
  //
  m_alpha_R  = value  ;
  m_AR       = my_exp ( -0.5 * m_alpha_R * m_alpha_R ) ;
  m_B        = 0.5 *  ( std::erf (  m_alpha_R * s_SQRT2i ) - 
                        std::erf ( -m_alpha_L * s_SQRT2i ) ) ;
  //
  if   ( s_equal ( m_n_R , 0 ) || s_equal ( m_alpha_R , 0 )  ) { m_TR = -1000 ; }
  else { m_TR  = ( m_n_R + 1 )  / std::abs ( m_alpha_R )  / m_n_R  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setN_L     ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n_L    ) ) { return false ; }
  //
  m_n_L      = value_ ;
  if ( s_equal ( m_n_L , 0 ) ) { m_n_L = 0 ; }
  //
  if   ( s_equal ( m_n_L , 0 ) || s_equal ( m_alpha_L , 0 )  ) {  m_TL = -1000 ; }
  else { m_TL  = ( m_n_L + 1 )  / std::abs ( m_alpha_L )  / m_n_L  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setN_R     ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n_R    ) ) { return false ; }
  //
  m_n_R      = value_ ;
  if ( s_equal ( m_n_R , 0 ) ) { m_n_R = 1 ; }
  //
  if   ( s_equal ( m_n_R , 0 ) || s_equal ( m_alpha_R , 0 )  ) { m_TR = -1000 ; }
  else { m_TR  = ( m_n_R + 1 )  / std::abs ( m_alpha_R )  / m_n_R  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::pdf ( const double x ) const
{
  //
  const double dx   = ( x - m_m0 ) / m_sigma ;
  //
  // the left tail
  if      ( dx  < -m_alpha_L )  // left tail
  {
    const double np1  = n_L() + 1 ;
    const double frac = np1 / ( np1 - std::abs ( m_alpha_L ) * ( m_alpha_L + dx ) )  ;
    return std::pow ( frac , np1 ) * m_AL * s_SQRT2PIi / sigma() ;
  }
  // the right tail
  else if  ( dx >  m_alpha_R )  // right tail
  {
    const double np1  = n_R () + 1 ;
    const double frac = np1 / ( np1 - std::abs ( m_alpha_R ) * ( m_alpha_R - dx ) )  ;
    return std::pow ( frac , np1 ) * m_AR * s_SQRT2PIi / sigma() ;
  }
  //
  // the peak
  //
  return my_exp ( -0.5 * dx * dx ) * s_SQRT2PIi / sigma() ; 
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
  const double x_low  = m_m0 - m_alpha_L * m_sigma ;
  const double x_high = m_m0 + m_alpha_R * m_sigma ;
  //
  // split into proper subintervals
  //
  if ( low < x_low  && x_low  < high )
  { return integral ( low , x_low  ) + integral ( x_low  , high ) ; }
  if ( low < x_high && x_high < high )
  { return integral ( low , x_high ) + integral ( x_high , high ) ; }
  //
  // Z = (x-x0)/sigma 
  //
  const double zlow  = ( low  - m_m0 ) / sigma() ;
  const double zhigh = ( high - m_m0 ) / sigma() ;
  //
  // the peak
  //
  if ( x_low <= low && high <= x_high )
  { return  s_SQRT2PIi * details::gaussian_int ( 0.5   , 0 , zlow  , zhigh ) ; }
  //
  // left tail 
  //
  if ( high <= x_low ) 
  {
    const double np1 = n_L () + 1 ;
    //
    const double A   = np1 ;
    const double B   = np1 ;
    const double C   = - std::abs ( alpha_L () ) ;
    //
    return s_SQRT2PIi * m_AL * 
      tail_integral ( A , B , C , np1 , zlow + alpha_L() , zhigh + alpha_L() ) ;
  }
  //
  // right tail 
  // 
  if ( low  >= x_high ) 
  {
    //
    const double np1 = n_R () + 1 ;
    //
    const double A   = np1 ;
    const double B   = np1 ;
    const double C   = std::abs ( alpha_R () ) ;
    //
    return s_SQRT2PIi * m_AR * 
      tail_integral ( A , B , C , np1 , zlow - alpha_R() , zhigh - alpha_R() ) ;
  }
  //
  return 0 ;
}
// ============================================================================
// get the (truncated)  integral
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::integral () const 
{
  //
  if      ( 0 < m_TL && 0 <= m_TR ) { return m_TL + m_TR + m_B ; }
  else if ( 0 < m_TR ) 
  {
    /// truncate it! 
    const double left = ( 0 < alpha_L() ) ?  (-alpha_L()-s_TRUNC) : -s_TRUNC ;
    // 
    return m_TR + m_B + integral ( m0 () + left       * sigma() , 
                                   m0 () - alpha_L () * sigma() ) ;
    
  }
  else if ( 0 < m_TL ) 
  {
    /// truncate it! 
    const double right = ( 0 < alpha_R() ) ?  ( alpha_R () + s_TRUNC) : + s_TRUNC ;
    // 
    return m_TL + m_B + integral ( m0 () + alpha_R () * sigma() , 
                                   m0 () + right      * sigma() ) ;
    
  }
  //
  /// truncate both
  const double left  = ( 0 < alpha_L() ) ?  (-alpha_L () - s_TRUNC ) : -s_TRUNC ;
  /// truncate it! 
  const double right = ( 0 < alpha_R() ) ?  ( alpha_R () + s_TRUNC ) : + s_TRUNC ;
  //
  return integral ( m0 () - left  * sigma () , m0 () + right * sigma () ) ;

}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallDoubleSided::tag () const 
{ return std::hash_combine ( m_m0      , m_sigma , 
                             m_alpha_L , m_n_L   , 
                             m_alpha_R , m_n_R   ) ; }
// ============================================================================




// ============================================================================
// apolonios 
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 *  @param b     b-parameter 
 */
// ============================================================================
Ostap::Math::Apolonios::Apolonios
( const double m0    ,
  const double sigma ,
  const double alpha ,
  const double n     ,
  const double bp    )
  : m_m0         ( m0 )
  , m_sigma      (  1 )
  , m_alpha      (  2 )
  , m_n          (  2 )
  , m_b          (  2 )
//
  , m_A  ( -1000 ) 
//
  , m_workspace () 
{
  //
  setM0     ( m0    ) ;
  setAlpha  ( alpha ) ;
  setSigma  ( sigma ) ;
  setN      ( n     ) ;
  setB      ( bp    ) ;
  //
  m_A = my_exp ( -b () * a1 () ) ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Apolonios::~Apolonios(){}
// ============================================================================
bool  Ostap::Math::Apolonios::setM0 ( const double value )
{
  //
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apolonios::setSigma ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apolonios::setAlpha  ( const double value )
{
  //
  if ( s_equal ( value , m_alpha ) ) { return false ; }
  //
  m_alpha    = value  ;
  //
  m_A = my_exp ( -b() * a1 () ) ; 
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apolonios::setN      ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n     ) ) { return false ; }
  //
  m_n        = value_ ;
  if ( s_equal ( m_n , 0 ) ) { m_n = 0 ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apolonios::setB  ( const double value )
{
  //
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_b ) ) { return false ; }
  //
  m_b    = value_  ;
  //
  if ( s_equal ( m_b , 0 ) ) { m_b = 0 ; }
  if ( s_equal ( m_b , 1 ) ) { m_b = 1 ; }
  //
  m_A = my_exp ( -b () * a1 () ) ;
  //
  return true ;
}
// ============================================================================
//  evaluate Apolonios' function
// ============================================================================
double Ostap::Math::Apolonios::pdf ( const double x ) const
{
  //
  const double dx    = ( x - m_m0 ) / m_sigma ;
  //
  // the tail
  //
  if  ( dx < -m_alpha )
  {
    const double frac = np1 () / ( np1 () - ( m_alpha + dx ) * aa () )  ;
    return std::pow ( frac , np1 () ) * m_A * s_SQRT2PIi / sigma() ;
  }
  //
  // the peak
  //
  return my_exp ( -b() * std::sqrt ( 1 + dx*dx ) ) * s_SQRT2PIi / sigma() ;  
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Apolonios::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  const double x0 = m_m0 - m_alpha * m_sigma ;
  //
  // split into proper subintervals
  //
  if      ( low < x0 && x0 < high )
  { return integral ( low , x0 ) + integral ( x0 , high ) ; }
  //
  // Z = (x-x0)/sigma 
  //
  const double zlow  = ( low  - m_m0 ) / sigma() ;
  const double zhigh = ( high - m_m0 ) / sigma() ;
  //
  // peak
  //
  if ( x0 <= low  )
  {
    //
    // use GSL to evaluate the integral
    //
    static const Ostap::Math::GSL::Integrator1D<Apolonios> s_integrator {} ;
    static char s_message[] = "Integral(Apolonios)" ;
    //
    const auto F = s_integrator.make_function ( this ) ;
    int    ierror   =  0 ;
    double result   =  1 ;
    double error    = -1 ;
    std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
      ( tag () , 
        &F     , 
        low    , high  ,               // low & high edges
        workspace ( m_workspace ) ,    // workspace
        s_PRECISION         ,          // absolute precision
        s_PRECISION         ,          // relative precision
        m_workspace.size()              ,          // size of workspace
        s_message           , 
        __FILE__ , __LINE__ ) ;
    //  
    return result ;
  }
  //
  // tail
  //
  const double A = np1 () ;
  const double B = np1 () ;
  const double C = - std::abs ( alpha () * b () ) / a1 ()  ;
  //
  const double result = s_SQRT2PIi * m_A * 
    tail_integral ( A , B , C , np1 () , zlow + alpha() , zhigh + alpha() ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Apolonios::tag () const 
{ return std::hash_combine ( m_m0      , m_sigma , 
                             m_alpha   , m_n     , m_b ) ; }
// ============================================================================

// ============================================================================
// apolonios2 
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 *  @param b     b-parameter 
 */
// ============================================================================
Ostap::Math::Apolonios2::Apolonios2
( const double m0     ,
  const double sigmaL ,
  const double sigmaR ,
  const double beta   )
  : m_m0         (  0 )
  , m_sigmaL     (  1 )
  , m_sigmaR     (  1 )
  , m_beta       (  1 )
  , m_workspace () 
{
  //
  setM0     ( m0     ) ;
  setSigmaL ( sigmaL ) ;
  setSigmaR ( sigmaR ) ;
  setBeta   ( beta   ) ;
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Apolonios2::~Apolonios2(){}
// ============================================================================
bool  Ostap::Math::Apolonios2::setM0 ( const double value )
{
  //
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apolonios2::setSigmaL ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigmaL ) ) { return false ; }
  //
  m_sigmaL    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apolonios2::setSigmaR ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigmaR ) ) { return false ; }
  //
  m_sigmaR    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apolonios2::setBeta ( const double value )
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
//  evaluate Apolonios' function
// ============================================================================
double Ostap::Math::Apolonios2::pdf ( const double x ) const
{
  //
  const double dx = 
    ( x < m_m0 ) ?
    ( x - m_m0 ) / m_sigmaL :
    ( x - m_m0 ) / m_sigmaR ;
  //
  // the peak
  //
  return my_exp ( beta() * ( beta()  - std::sqrt ( b2 () + dx * dx ) ) ) * s_SQRT2PIi / sigma()  ;  
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Apolonios2::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  const double xR = m_m0 + 4.0 * m_sigmaR ;
  if ( low < xR && xR < high ) 
  { return integral ( low , xR ) + integral ( xR , high ) ; }
  //
  const double xL = m_m0 - 4.0 * m_sigmaL ;
  if ( low < xL && xL < high ) 
  { return integral ( low , xL ) + integral ( xL , high ) ; }
  //
  const double in_tail = ( low >= xR || high <= xL ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Apolonios2> s_integrator {} ;
  static char s_message[] = "Integral(Apolonios2)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag () , 
      &F     , 
      low , high  ,                  // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Apolonios2::tag () const 
{ return std::hash_combine ( m_m0 , m_sigmaL , m_sigmaR , m_beta ) ; }
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
// destructor
// ============================================================================
Ostap::Math::Atlas::~Atlas(){}
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
double Ostap::Math::Atlas::integral ( const double low  ,
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
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
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
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Atlas::tag () const 
{ return std::hash_combine ( m_mean , m_sigma ) ; }
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
// destructor
// ============================================================================
Ostap::Math::Sech::~Sech(){}
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
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Sech::tag () const 
{ return std::hash_combine ( m_mean , m_sigma ) ; }
// ============================================================================


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
// destructor
// ============================================================================
Ostap::Math::Logistic::~Logistic(){}
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
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Logistic::tag () const 
{ return std::hash_combine ( m_mean , m_sigma ) ; }
// ============================================================================




// ============================================================================
// Student-T 
// ============================================================================
/*  constructor from mass, resolution and "n"-parameter 
 *  @param M     mass 
 *  @param sigma width parameter
 *  @param N     n-parameter  ( actually  n=1+|N| ) 
 */
// ============================================================================
Ostap::Math::StudentT::StudentT 
( const double mass  , 
  const double sigma ,
  const double n     ) 
//
  : m_M    (      std::abs ( mass  ) )
  , m_s    (      std::abs ( sigma ) )
  , m_n    ( -1 )
  , m_norm ( -1 ) 
{
  setN ( n ) ;  
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::StudentT::~StudentT (){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setM ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_M ) ) { return false ; }
  //
  m_M = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setSigma ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_s ) ) { return false ; }
  //
  m_s = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setN ( const double x )
{
  //
  const double v = 1 + std::abs ( x ) ;
  //
  if ( m_norm < 0 ) 
  {
    m_norm  = gsl_sf_gamma ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
    m_norm /= std::sqrt    ( M_PI * v ) ;
  }
  //
  if ( s_equal ( v , m_n ) ) { return false ; }
  //
  m_n = v ;
  //
  m_norm  = gsl_sf_gamma ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
  m_norm /= std::sqrt    ( M_PI * v ) ;
  //
  return true ;
}
// ==========================================================================
double Ostap::Math::StudentT::pdf ( const double x ) const
{
  //
  const double y = ( x - M () ) / sigma () ;
  //
  const double f = std::pow (  1 + y * y / nu () ,  -0.5 * ( nu () + 1 ) ) ;
  //
  return m_norm * f / sigma () ; // sigma comes from dx = dy * sigma 
}
// ============================================================================
double Ostap::Math::StudentT::cdf ( const double y ) const
{
  //
  const double  t    = ( y - M () ) / sigma () ;
  return student_cdf ( t , nu() ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::StudentT::integral() const { return 1 ; }
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::StudentT::integral
( const double low  , 
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::StudentT::tag () const 
{ return std::hash_combine ( m_M , m_s , m_n ) ; }
// ============================================================================

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
  , m_normL ( -1 ) 
  , m_normR ( -1 ) 
{
  setNL ( nL ) ;  
  setNR ( nR ) ;  
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::BifurcatedStudentT::~BifurcatedStudentT (){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setM ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_M ) ) { return false ; }
  //
  m_M = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setSigmaL ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_sL ) ) { return false ; }
  //
  m_sL = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setSigmaR ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_sR ) ) { return false ; }
  //
  m_sR = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setNL ( const double x )
{
  //
  const double v = 1 + std::abs ( x ) ;
  //
  if ( m_normL < 0 ) 
  {
    m_normL = gsl_sf_gamma  ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
    m_normL /= std::sqrt    ( M_PI * v ) ;
  }
  //
  if ( s_equal ( v , m_nL ) ) { return false ; }
  //
  m_nL =  v ;
  //
  m_normL = gsl_sf_gamma  ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
  m_normL /= std::sqrt    ( M_PI * v ) ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setNR ( const double x )
{
  //
  const double v = 1 + std::abs ( x ) ;
  //
  if ( m_normR < 0 ) 
  {
    m_normR  = std::tgamma ( 0.5 * ( v + 1 ) ) / std::tgamma ( 0.5 * v ) ;  
    m_normR /= std::sqrt   ( M_PI * v ) ;
  }
  //
  if ( s_equal ( v , m_nR ) ) { return false ; }
  //
  m_nR = v ;
  //
  m_normR  = std::tgamma ( 0.5 * ( v + 1 ) ) / std::tgamma ( 0.5 * v ) ;  
  m_normR /= std::sqrt   ( M_PI * v ) ;
  //
  return true ;
}
// ==========================================================================
double Ostap::Math::BifurcatedStudentT::pdf ( const double x ) const
{
  //
  const double y = ( x <= M() ) ? 
    ( x - M () ) / sigmaL () : ( x - M () ) / sigmaR () ;
  //
  const double f = ( x <= M() ) ? 
    std::pow (  1 + y * y / nuL () ,  -0.5 * ( nuL () + 1 ) ) :
    std::pow (  1 + y * y / nuR () ,  -0.5 * ( nuR () + 1 ) ) ;
  //
  const double n_1 = m_normL       / sigmaL () ;
  const double n_2 = m_normR       / sigmaR () ;
  const double n_t = 2 * n_1 * n_2 / ( n_1 + n_2 ) ;
  //
  return n_t * f ; 
}
// ============================================================================
double Ostap::Math::BifurcatedStudentT::cdf ( const double y ) const
{
  //
  const double n_1 = m_normL / sigmaL () ;
  const double n_2 = m_normR / sigmaR () ;
  //
  if ( y <= M() ) 
  {
    const double  t    = ( y - M () ) / sigmaL () ;
    return     2 * n_2 / ( n_1 + n_2 ) * student_cdf (  t , nuL () ) ;  
  }
  //
  const   double  t    = ( y - M () ) / sigmaR () ;
  return   1 - 2 * n_1 / ( n_1 + n_2 ) * student_cdf ( -t , nuR () ) ;  
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::BifurcatedStudentT::integral() const { return 1 ; }
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::BifurcatedStudentT::integral
( const double low  , 
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::BifurcatedStudentT::tag () const 
{ return std::hash_combine ( m_M , m_sL , m_sR  , m_nL , m_nR ) ; }
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
Ostap::Math::SinhAsinh::~SinhAsinh(){}
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
double Ostap::Math::SinhAsinh::integral ( const double low  ,
                                          const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::SinhAsinh::tag () const 
{ return std::hash_combine ( m_mu , m_sigma , m_epsilon , m_delta ) ; }
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
( const double xi      ,   // related to location 
  const double lambda  ,   // related to variance
  const double delta   ,   // shape 
  const double gamma   )  // shape 
  : m_xi      (            xi) 
  , m_lambda  ( std::abs ( lambda ) ) 
  , m_delta   ( std::abs ( delta  ) ) 
  , m_gamma   (            gamma    ) 
{}
// ============================================================================
// Destructor
// ============================================================================
Ostap::Math::JohnsonSU::~JohnsonSU(){}
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
  const double d2 = ( d1 - 1 ) * ( d1 * std::cosh ( 2  *m_gamma / m_delta ) + 1 ) ;
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
// get the tag 
// ============================================================================
std::size_t Ostap::Math::JohnsonSU::tag () const 
{ return std::hash_combine ( m_xi , m_lambda , m_delta , m_gamma ) ; }
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
// destructor
// ============================================================================
Ostap::Math::Slash::~Slash(){}
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
{ return std::hash_combine ( m_mu , m_scale ) ; }
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
// get the tag 
// ============================================================================
std::size_t Ostap::Math::RaisingCosine::tag () const 
{ return std::hash_combine ( m_mu , m_s ) ; }
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
// destructor 
// ============================================================================
Ostap::Math::AsymmetricLaplace::~AsymmetricLaplace(){}
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
{ return std::hash_combine ( m_mu , m_lambdaL , m_lambdaR ) ; }
// ============================================================================


// ============================================================================
/*  constructor from all arguments            
 *  @param mean  the mean/mode/location of the peak 
 *  @param q     q-value   (q<3, for q>3, it will be reflected)
 *  @param scale 
 */
// ============================================================================
Ostap::Math::QGaussian::QGaussian 
( const double mean  ,   // mean/mode/location 
  const double q     ,   //  q-parameter 
  const double scale )   // scale/sigma
  : m_mean  ( mean               )
  , m_q     ( 1                  ) 
  , m_scale ( std::abs ( scale ) )
{
  setQ ( q ) ; 
}
// ============================================================================
// evaluate  pdf  for q-Gaussian distribution
// ============================================================================
namespace 
{
  // 
  inline double q_exp ( const long double x , const long double q ) 
  { return std::pow ( 1 + ( 1 - q ) * x , 1 / ( 1 - q ) ) ; }
  //
}
// ============================================================================
// evaluate PDF of q-Gaussian distribution 
// ============================================================================
double Ostap::Math::QGaussian::pdf ( const  double x ) const 
{
  //
  if ( s_equal ( m_q , 1 ) ) { return gauss_pdf ( x , m_mean , m_scale ) ; }
  //
  const long double dx  =   ( x - m_mean ) / m_scale ;
  //
  static const double s_sq2 = std::sqrt ( 2.0 ) ;
  //
  return
    ( m_q < 1 && s_sq2 <= std::abs ( dx ) * std::sqrt ( 1 - m_q ) ) ? 0.0 :  
    q_exp ( - 0.5 * dx * dx , m_q ) / ( s_sq2 * m_scale * m_cq ) ;  
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
  static const double s_sqpi = std::sqrt ( M_PI ) ;
  m_cq = s_sqpi ; 
  if      ( 1 > m_q ) 
  {
    const long double q  =  m_q ;
    const long double g1 = std::lgamma ( 1.0             / ( 1 - q ) ) ;
    const long double g2 = std::lgamma ( 0.5 * ( 3 - q ) / ( 1 - q ) ) ;
    //
    m_cq = 2 * s_sqpi * std::exp ( g1  - std::log ( 3 - q ) - 
                                   0.5 * std::log ( 1 - q ) - g2 ) ;
  }
  else if ( 1 < m_q )
  {
    const long double q  =  m_q ;
    const long double g1 = std::lgamma ( 1.0             / ( q - 1 ) ) ;
    const long double g2 = std::lgamma ( 0.5 * ( 3 - q ) / ( q - 1 ) ) ;
    //
    m_cq = s_sqpi * std::exp ( g2 - 0.5 * std::log ( q - 1 ) - g1 ) ;
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
  m_scale = v ;
  return true ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::QGaussian::integral ( const double low  , 
                                          const double high ) const 
{
  ///
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) {  return -integral ( high , low ) ; }
  ///
  if ( s_equal ( m_q , 1 ) ) 
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
    const double win  = s_sq2 * m_scale / ( 1 - m_q ) ;
    const double xmin = m_mean - win ;
    const double xmax = m_mean + win ;
    if ( high <= xmin || low >= xmax ) { return 0 ; } // RETURN
    xlow  = std::max ( xmin , xlow  ) ;
    xhigh = std::min ( xmax , xhigh ) ;
  }
  //  are we already in the tail? 
  const bool in_tail = 
    std::min ( std::abs ( xhigh - m_mean ) , std::abs ( m_mean - xlow ) )  > 5 * m_scale ;   
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
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
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
{ return std::hash_combine ( m_mean , m_q , m_scale ) ; }
// ============================================================================
 

// ============================================================================
//                                                                      The END 
// ============================================================================
