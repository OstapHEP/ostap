// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <map>
#include <limits>
#include <complex>
#include <algorithm>
#include <numeric>
#include <array>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_log.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_psi.h"
#include "gsl/gsl_sf_zeta.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TMath.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Math.h"
#include "Ostap/Power.h"
#include "Ostap/StatusCode.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/MoreMath.h"
#include "Ostap/PhaseSpace.h"
#include "Ostap/Models.h"
// ============================================================================
//  Local 
// ============================================================================
#include "local_math.h"
#include "local_hash.h"
#include "local_gsl.h"
#include "status_codes.h"
#include "Integrator1D.h"
// ============================================================================
/** @file
 *  Implementation file for functions from the file Ostap/Models.h
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-04-19
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /// join two vectors together
  template <class SCALAR>
  inline
  std::vector<SCALAR>
  join
  ( const std::vector<SCALAR>& a ,
    const std::vector<SCALAR>& b )
  {
    // destination 
    std::vector<double> ab ( a.size() + b.size() ) ;
    std::copy ( a.begin () , a.end () , ab.begin ()            ) ;
    std::copy ( b.begin () , b.end () , ab.begin () + a.size() ) ;
    return ab ;    
  }
  // ==========================================================================
} // end of anonymous namespace
// ============================================================================
/*  constructor  from all parameters 
 *  @param mu location, bias parameter 
 *  @param beta scale parameter 
 */
// ============================================================================
Ostap::Math::Gumbel::Gumbel
( const double mu   , 
  const double beta )
  : ShiftAndScale ( beta , mu   , "beta" , "mu"  , typeid ( *this ) , false )
{}
// ============================================================================
// minimal x 
// ============================================================================
double Ostap::Math::Gumbel::xmin () const { return s_NEGINF ; } 
// ============================================================================
// maximal x 
// ============================================================================
double Ostap::Math::Gumbel::xmax () const { return s_POSINF ; } 
// ============================================================================
double Ostap::Math::Gumbel::median () const 
{ return x ( -s_ln_ln2 ) ; }  
// ============================================================================
double Ostap::Math::Gumbel::mean () const 
{ return x ( s_Euler ) ;  }
// ============================================================================
double Ostap::Math::Gumbel::variance () const 
{
  const double ss  = scale () ;
  return ss * ss * s_pi2 / 6 ; 
}
// ============================================================================
double Ostap::Math::Gumbel::rms () const 
{ return std::sqrt ( variance () ) ; }
// ============================================================================
double Ostap::Math::Gumbel::skewness () const 
{  
  static const double s_skew  =  12 * std::sqrt ( 6.0L ) * Ostap::Math::zeta ( 3 ) / s_pi3  ;
  return scale_sign () * s_skew ;
}
// ============================================================================
// get a value for the function      
// ============================================================================
double Ostap::Math::Gumbel::pdf  ( const double x ) const 
{
  const double y = t ( x ) ; 
  const double s = std::abs  ( scale () ) ;
  return std::exp ( - y - std::exp ( -y ) ) / s ;
}
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::Gumbel::cdf ( const double x ) const 
{
  const double y    = t ( x ) ; 
  const double rcdf = std::exp ( - std::exp ( -y ) ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Gumbel::integral
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  return cdf ( high ) - cdf ( low ) ; 
}
// ============================================================================
/*  quantile function
 *  @parameter p  probability \f$ 0 < p < 1 \f$
 */
// ============================================================================
double Ostap::Math::Gumbel::quantile
( const double p ) const
{
  return
    p < 0             ?  s_QUIETNAN :
    p > 1             ?  s_QUIETNAN : 
    s_zero  ( p     ) ?  s_NEGHUGE  : 
    s_equal ( p , 1 ) ?  s_POSHUGE  :
    x ( - std::log ( - std::log ( p ) ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Gumbel::tag () const 
{ 
  static const std::string s_name = "Gumbel" ;
  return Ostap::Utils::hash_combiner ( s_name , ShiftAndScale::tag () ) ; 
}
// ============================================================================
 
// ============================================================================
/* constructor form scale & shape parameters
 *  param k      \f$k\f$ parameter (shape)
 *  param theta  \f$\theta\f$ parameter (scale)
 */
// ============================================================================
Ostap::Math::GammaDist::GammaDist 
( const double k     ,   // shape parameter  
  const double theta ,   // scale parameter
  const double shift )   // shidt parameter
  : ShiftAndScale ( theta , shift , "theta" , "shift" , typeid ( *this ) , false ) 
  , m_k           ( k                       , "k"     , typeid ( *this ) )
{
  // evaluate auxillary parameter 
  m_lgk = Ostap::Math::lgamma ( m_k.value() ) ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::GammaDist::setK ( const double x )
{
  if ( !m_k.setValue ( x ) ) { return false ; }
  // evaluate auxillary parameter 
  m_lgk = Ostap::Math::lgamma ( m_k.value() ) ;
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::GammaDist::setLogK ( const double x )
{
  if ( !m_k.setLogValue ( x ) ) { return false ; }
  // evaluate auxillary parameter 
  m_lgk = Ostap::Math::lgamma ( m_k.value() ) ;
  return true ;
}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::GammaDist::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::GammaDist::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 0 ) ; }
// ============================================================================
// raw moments (unscaled & unshifted)
// ============================================================================
double Ostap::Math::GammaDist::raw_moment
( const unsigned short r ) const
{ return rising_factorial  ( k () , r ) ; }
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::GammaDist::mean       () const
{
  //
  const double kv     = k     () ; 
  const double tv     = theta () ;
  //
  return x ( kv * tv ) ; 
}
// ============================================================================
// dispersion value 
// ============================================================================
double Ostap::Math::GammaDist::variance () const
{
  //
  const double kv     = k     () ; 
  const double ss     = scale () ;
  //
  return kv * ss * ss ; 
}
// ============================================================================
double Ostap::Math::GammaDist::rms     () const
{ return std::sqrt ( dispersion ()  ) ; }
// ============================================================================
double Ostap::Math::GammaDist::skewness () const
{
  const double kv = k () ; 
  const double r  = 2.0 / std::sqrt ( kv ) ;
  return scale_sign () * r ;
}
// ============================================================================
// calculate gamma distribution shape
// ============================================================================
double Ostap::Math::GammaDist::evaluate ( const double x ) const
{
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ;}
  //
  const double kv     = k     () ; 
  if ( 1 < kv && s_zero ( y ) ) { return 0 ; }
  //
  const double logr = -y + ( kv - 1 ) * std::log ( y  ) - m_lgk ;
  //
  const  double s = std::abs ( scale () ) ;
  return std::exp ( logr ) / s ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::GammaDist::integral () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::GammaDist::integral
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low  , high )  ) { return  0 ; }
  else if (           low  > high    ) { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the cdf 
// ============================================================================
double Ostap::Math::GammaDist::cdf 
( const double x ) const 
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ; 
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ; 
}
// ============================================================================
// "raw" cdf (for positive scale)
// ============================================================================
double Ostap::Math::GammaDist::raw_cdf
( const double y ) const
{
  if  ( y <= 0 ) { return 0 ; }
  const double kv = k     () ;  
  return gsl_sf_gamma_inc_P ( kv , y ) ;
}
// ============================================================================
// calculate the mode
// ============================================================================
double Ostap::Math::GammaDist::mode () const
{
  const double kv = k     () ; 
  const double tv = theta () ;
  const double r  =  1 <= kv ? ( kv - 1 ) * tv : 0.0 ;  
  //
  return x ( r ) ; 
}
// ============================================================================
// calculate the quantile   (0<p<1) 
// ============================================================================
double Ostap::Math::GammaDist::quantile ( const double p ) const 
{
  const double kv = k     () ; 
  const double tv = theta () ;

  return
    p <  0            ? s_QUIETNAN :
    p >  1            ? s_QUIETNAN :
    s_zero  ( p     ) ? x ( 0 )    :
    s_equal ( p , 1 ) ? s_POSHUGE  :    
    x ( gsl_cdf_gamma_Pinv ( p , kv , tv ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GammaDist::tag () const 
{ 
  static const std::string s_name = "GammaDist" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       m_k           .tag () , 
				       ShiftAndScale::tag () ) ;
}
// ============================================================================
/*  effective \f$ \chi^2 \f$-parameters
 *  If   \f$ Q  \sim \chi^2(\nu)\f$  and c is a positive constant,
 *  than \f$ cQ \sim \Gamma (k = \nu/2, \theta = 2c) \f$
 *  @attention for shift!=0  NaN is returned
 */
// ============================================================================
double Ostap::Math::GammaDist::nu () const
{
  const double kv = k     () ;  
  const double sv = shift () ;
  //
  return !sv ? 2 * kv : s_QUIETNAN ;
}
// ============================================================================
/*  effective \f$ \chi^2 \f$-parameters
 *  If   \f$ Q  \sim \chi^2(\nu)\f$  and c is a positive constant,
 *  than \f$ cQ \sim \Gamma (k = \nu/2, \theta = 2c) \f$
 *  @attention for shift!=0  NaN is returned
 */
// ============================================================================
double Ostap::Math::GammaDist::c  () const
{
  const double tv = theta () ;  
  const double sv = shift () ;
  //
  return !sv ? 0.5 * tv : s_QUIETNAN  ;
}
// ============================================================================

// ============================================================================
// Generalized Gamma distribtion
// ============================================================================
/*  constructor
 *  param k     \f$k\f$ parameter      (shape)
 *  param theta \f$\theta\f$ parameter (scale)
 *  param p     \f$p\f$ parameter 
 *  param low   bias       
 */
// ============================================================================
Ostap::Math::GenGammaDist::GenGammaDist
( const double k     , 
  const double theta , 
  const double p     , // 1 corresponds to gamma distribution 
  const double shift )
  : ShiftAndScale ( theta , shift , "theta" , "shift" , typeid ( *this ) , false ) 
  , m_k           ( k                       , "p"     , typeid ( *this ) )
  , m_p           ( p                       , "p"     , typeid ( *this ) )
{
  // evaluate auxillary parameter 
  m_lgkp = Ostap::Math::lgamma ( m_k.value() / m_p.value () ) ;
}
// ============================================================================
bool Ostap::Math::GenGammaDist::setK     ( const double value ) 
{
  if ( !m_k.setValue ( value ) ) { return false ; }
  // evaluate auxillary parameter 
  m_lgkp = Ostap::Math::lgamma ( m_k.value() / m_p.value () ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGammaDist::setP    ( const double value ) 
{
  if ( !m_p.setValue ( value ) ) { return false ; }
  // evaluate auxillary parameter 
  m_lgkp = Ostap::Math::lgamma ( m_k.value() / m_p.value () ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGammaDist::setLogK     ( const double value ) 
{
  if ( !m_k.setLogValue ( value ) ) { return false ; }
  // evaluate auxillary parameter 
  m_lgkp = Ostap::Math::lgamma ( m_k.value() / m_p.value () ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGammaDist::setLogP    ( const double value ) 
{
  if ( !m_p.setLogValue ( value ) ) { return false ; }
  // evaluate auxillary parameter 
  m_lgkp = Ostap::Math::lgamma ( m_k.value() / m_p.value () ) ;
  return true ;
}
// ============================================================================
double Ostap::Math::GenGammaDist::evaluate ( const double x ) const 
{
  //
  const double y  = t ( x ) ;
  if ( y <  0 ) { return 0 ; }
  //
  const double kv = k () ; 
  if ( 1 < kv && s_zero ( y ) ) { return 0 ; }
  //
  const double pv = p () ;
  //
  const double yp = std::pow ( y , pv ) ;
  //
  const double lr = ( kv - 1 ) * std::log ( y ) - yp - m_lgkp ;
  //
  const double sv = std::abs ( scale ()  ) ;
  return pv * std::exp ( lr ) / sv ; 
}
// ============================================================================
double Ostap::Math::GenGammaDist::cdf ( const double x ) const 
{
  const double y = t ( x ) ;
  const double r = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? r : 1 - r ;
}
// ============================================================================
double Ostap::Math::GenGammaDist::raw_cdf ( const double y ) const 
{
  if ( y <= 0 ) { return 0 ; }
  //
  const double kv = k () ;
  const double pv = p () ;
  //
  const double yp = std::pow ( y , pv ) ;
  return gsl_sf_gamma_inc_P ( kv / pv , yp ) ;  
}
// ============================================================================
double Ostap::Math::GenGammaDist::integral () const { return 1 ; }
// ============================================================================
double Ostap::Math::GenGammaDist::integral
( const double low  , 
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high < low             ) { return -integral ( high , low ) ; }  
  return cdf ( high ) - cdf ( low ) ;  
}
// ============================================================================
// raw moments (unscaled & unshifted)
// ============================================================================
double Ostap::Math::GenGammaDist::raw_moment
( const unsigned short r ) const
{
  const double kv = k () ;
  const double pv = p () ;
  return std::exp ( Ostap::Math::lgamma ( ( kv + r ) / pv ) -
		    Ostap::Math::lgamma (   kv       / pv ) ) ;
}  
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::GenGammaDist::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::GenGammaDist::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 0 ) ; }
// ============================================================================
// mode 
// ============================================================================
double Ostap::Math::GenGammaDist::mode () const
{
  const double kv = k () ;
  const double pv = p () ;
  //
  const double r  =  1 < kv ? std::pow ( ( kv - 1 ) / pv , 1/pv ) : 0.0 ;
  return x ( r ) ;
}
// ============================================================================
// mean 
// ============================================================================
double Ostap::Math::GenGammaDist::mean () const
{ return x ( raw_moment ( 1 ) ) ; }
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::GenGammaDist::variance () const
{
  const double ss = scale () ;
  return ss * ss * Ostap::Math::variance ( raw_moment ( 2 ) ,
					   raw_moment ( 1 ) ) ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::GenGammaDist::rms () const
{ return std::sqrt ( variance () ) ; }
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::GenGammaDist::skewness () const
{
  return scale_sign () * Ostap::Math::skewness ( raw_moment ( 3 ) ,
						 raw_moment ( 2 ) ,
						 raw_moment ( 1 ) ) ;
}
// ============================================================================
// kurtosis 
// ============================================================================
double Ostap::Math::GenGammaDist::kurtosis  () const
{
  return Ostap::Math::kurtosis ( raw_moment ( 4 ) ,
				 raw_moment ( 3 ) ,
				 raw_moment ( 2 ) ,
				 raw_moment ( 1 ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenGammaDist::tag () const 
{ 
  static const std::string s_name = "GenGammaDist" ;
  return Ostap::Utils::hash_combiner ( s_name ,
				       ShiftAndScale::tag () , 
				       m_p           .tag () ,
				       m_k           .tag () ) ;
}
// ============================================================================


// ============================================================================
// Amoroso 
// ============================================================================
/*  constructor
 *  param a     a-parameter 
 *  param theta \f$\theta\f$-parameter  
 *  param alpha \f$\alpha\f$-parameter (>0)
 *  param beta  \f$\beta\f$-parameter 
 *  Note that   \f$\alpha\beta\f$ is equal to k-parameter 
 */
// ============================================================================
Ostap::Math::Amoroso::Amoroso 
( const double theta , 
  const double alpha , 
  const double beta  ,
  const double a     )
  : ShiftAndScale ( theta , a , "theta" , "a"     , typeid ( *this ) , false ) 
  , m_alpha       ( alpha               , "alpha" , typeid ( *this ) ) 
  , m_beta        ( beta                , "beta"  , typeid ( *this ) ) 
{
  m_iga = Ostap::Math::igamma ( m_alpha.value () ) ; 
}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::Amoroso::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::Amoroso::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 0 ) ; }
// ============================================================================
// Set parameter alpha
// ============================================================================
bool Ostap::Math::Amoroso::setAlpha  ( const double value )
{
  if ( !m_alpha.setValue ( value ) ) { return false ; }
  m_iga = Ostap::Math::igamma ( m_alpha.value () ) ;
  return true ;
}
// ============================================================================
// evaluate Amoroso distribtion
// ============================================================================
double Ostap::Math::Amoroso::evaluate ( const double x ) const 
{
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const double av = alpha () ;
  const double bv = beta  () ;
  //
  const double ln1 = ( av * bv - 1 ) * std::log ( y ) ;
  const double ln2 = std::pow ( y , bv ) ;  
  //
  return m_iga * std::abs ( bv ) * std::exp ( ln1 - ln2 ) / abs_scale() ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::Amoroso::cdf ( const double x ) const 
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1  - rcdf ;
}
// ============================================================================
// "raw" cdf (for positive scale)
// ============================================================================
double Ostap::Math::Amoroso::raw_cdf ( const double y ) const 
{
  if ( y <= 0 ) { return 0 ; }
  //
  const double av = alpha () ;
  const double bv = beta  () ;
  const double tv = theta () ;
  //
  const double xt = std::pow ( y , bv ) ;
  //
  return 0 < bv ? 1 - gsl_sf_gamma_inc_Q ( av , xt ) : gsl_sf_gamma_inc_Q ( av , xt ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::integral () const { return 1 ; }
// ============================================================================
double Ostap::Math::Amoroso::integral
( const double low  , 
  const double high ) const 
{
  if ( s_equal ( low ,high ) ) { return 0 ; }
  return  cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::mode () const 
{
  const double av = alpha () ;
  const double bv = beta  () ;
  const double m  = av * bv < 1 ? 0 : std::pow ( av - 1 / bv , 1 / bv ) ;
  return x ( m ) ;
}
// ============================================================================
// raw moments (unscaled & unshifted)
// ============================================================================
double Ostap::Math::Amoroso::raw_moment
( const unsigned short k ) const
{
  if ( !k ) { return 1 ; }
  //
  const double av = alpha () ;
  const double bv = beta  () ;
  const double vv = av + k / bv ;
  //
  if ( vv <= 0 || s_zero ( vv ) ) { return s_QUIETNAN ; }
  //
  return Ostap::Math::pochhammer ( av , k / bv ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::mean () const 
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  return x ( m1 ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::variance () const 
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  //
  const double ss = scale () ;
  return ss * ss * Ostap::Math::variance ( m2 , m1 ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::rms () const 
{
  const double vv = variance () ;
  if ( !std::isfinite ( vv ) ) { return vv ; } 
  return std::sqrt ( vv ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::skewness () const 
{ 
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m3 = raw_moment ( 3 ) ;
  if ( !std::isfinite ( m3 ) ) { return m3 ; }
  //
  return scale_sign() * Ostap::Math::skewness ( m3 , m2 , m1 ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::kurtosis  () const 
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m3 = raw_moment ( 3 ) ;
  if ( !std::isfinite ( m3 ) ) { return m3 ; }
  const double m4 = raw_moment ( 4 ) ;
  if ( !std::isfinite ( m4 ) ) { return m4 ; }
  //
  return Ostap::Math::kurtosis ( m4 , m3 , m2 , m1 ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Amoroso::tag () const 
{ 
  static const std::string s_name = "Amoroso" ;
  return Ostap::Utils::hash_combiner ( s_name ,
				       ShiftAndScale::tag () ,
				       m_alpha       .tag () ,
				       m_beta        .tag () ) ;
}
// ============================================================================

// ============================================================================
/* constructor from scale & shape parameters
 *  param k      \f$k\f$ parameter (shape)
 *  param theta  \f$\theta\f$ parameter (scale)
 */
// ============================================================================
Ostap::Math::LogGammaDist::LogGammaDist 
( const double k     ,   // shape parameter  
  const double theta ,   // scale parameter
  const double shift )   // shift parameter
  : m_gamma ( k , theta , shift ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::LogGammaDist::~LogGammaDist (){}
// ============================================================================
// calculate log-gamma distribution shape
// ============================================================================
double Ostap::Math::LogGammaDist::evaluate ( const double x ) const
{
  const double z = my_exp ( x ) ;
  return m_gamma ( z ) * z ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::LogGammaDist::integral () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::LogGammaDist::integral
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low  , high ) ) { return 0 ; }
  else if (           low  > high   ) { return -1 * integral ( high , low  ) ; }
  //
  const double z_low  = my_exp ( low  ) ;
  const double z_high = my_exp ( high ) ;
  //
  return m_gamma.integral ( z_low , z_high ) ;
}
// ============================================================================
// calculate the quantile   (0<p<1) 
// ============================================================================
double Ostap::Math::LogGammaDist::quantile ( const double p ) const 
{
  if      ( p < 0 ) { return s_QUIETNAN ; }
  else if ( p > 1 ) { return s_QUIETNAN ; }
  //
  return my_log ( gsl_cdf_gamma_Pinv ( p , k() , theta() ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::LogGammaDist::tag () const 
{ 
  static const std::string s_name = "LogGammaDist" ;
  return Ostap::Utils::hash_combiner ( s_name , m_gamma.tag () ) ;
}
// ============================================================================

// ============================================================================
/* constructor form scale & shape parameters
 *  param k      \f$k\f$ parameter (shape)
 *  param theta  \f$\theta\f$ parameter (scale)
 */
// ============================================================================
Ostap::Math::Log10GammaDist::Log10GammaDist 
( const double k     ,   // shape parameter  
  const double theta ,   // scale parameter
  const double shift )   // shift parameter
  : Ostap::Math::LogGammaDist( k , theta , shift ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Log10GammaDist::~Log10GammaDist (){}
// ============================================================================
// calculate log-gamma distribution shape
// ============================================================================
double Ostap::Math::Log10GammaDist::evaluate ( const double x ) const
{ return LogGammaDist::evaluate ( x * s_ln10 ) * s_ln10 ; }
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::Log10GammaDist::integral () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Log10GammaDist::integral
( const double low  ,
  const double high ) const 
{ return LogGammaDist::integral ( low  * s_ln10 , high * s_ln10 ) ; }
// ============================================================================
// calculate the quantile   (0<p<1) 
// ============================================================================
double Ostap::Math::Log10GammaDist::quantile ( const double p ) const 
{
  if      ( p < 0 ) { return s_QUIETNAN ; }
  else if ( p > 1 ) { return s_QUIETNAN ; }
  return LogGammaDist::quantile ( p ) / s_ln10 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Log10GammaDist::tag () const 
{ 
  static const std::string s_name = "Log10GammaDist" ;
  return Ostap::Utils::hash_combiner ( s_name , m_gamma.tag () ) ; 
}
// ============================================================================



// ============================================================================
// Log-Gamma
// ============================================================================
/*  constructor from scale & shape parameters
 *  param nu      \f$\nu\f$ parameter      (location)
 *  param lambda  \f$\lambda\f$ parameter  
 *  param alpha   \f$\alpha\f$ parameter    (>0)
 */
// ============================================================================
Ostap::Math::LogGamma::LogGamma
( const double nu     , 
  const double lambda , 
  const double alpha  ) 
  : m_nu     ( nu     ) 
  , m_lambda ( lambda ) 
  , m_alpha  ( std::abs ( alpha ) ) 
{
}
// ============================================================================
bool Ostap::Math::LogGamma::setNu   ( const double value ) 
{
  if ( s_equal ( value , m_nu ) ) { return false ; }
  m_nu  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::LogGamma::setLambda ( const double value ) 
{
  if ( s_equal ( value , m_lambda ) ) { return false ; }
  m_lambda = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::LogGamma::setAlpha ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  return true ;
}
// ============================================================================
// calculate log-gamma shape
// ============================================================================
double Ostap::Math::LogGamma::pdf ( const double x ) const 
{
  //
  const double xc  = x  -  nu    () ;
  const double xt  = xc / lambda () ;
  //
  const double arg = alpha() * xt - my_exp ( xt ) ;
  //
  double result  = arg ;
  result        -= gsl_sf_log      ( std::abs ( lambda () ) ) ;
  result        -= gsl_sf_lngamma  (            alpha  ()   ) ;
  //
  return my_exp ( result ) ;
}
// ============================================================================
double Ostap::Math::LogGamma::cdf ( const double x ) const 
{
  //
  const double xc  = x  -  nu    () ;
  const double xt  = xc / lambda () ;
  //
  const double ext = my_exp ( xt ) ;
  //
  return 
    lambda () > 0 ? 
    1 - gsl_sf_gamma_inc_Q ( alpha() , ext ) : gsl_sf_gamma_inc_Q ( alpha() , ext ) ;
}
// ============================================================================
double Ostap::Math::LogGamma::integral ( const double low  , 
                                         const double high ) const
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::LogGamma::integral () const { return 1 ; }
// ============================================================================
double Ostap::Math::LogGamma::mode     () const 
{ return nu() - lambda() * gsl_sf_log ( alpha () ) ; }
// ============================================================================
double Ostap::Math::LogGamma::mean     () const 
{ return nu() + lambda() * gsl_sf_psi ( alpha () ) ; }
// ============================================================================
double Ostap::Math::LogGamma::sigma    () const 
{ return std::sqrt ( variance () ) ; }
// ============================================================================
double Ostap::Math::LogGamma::variance () const 
{ return lambda() * lambda() * gsl_sf_psi_1 ( alpha () ) ; }
// ============================================================================
double Ostap::Math::LogGamma::skewness () const 
{ 
  const double l_p2 = gsl_sf_psi_n ( 2 , alpha () ) ; 
  const double l_p1 = gsl_sf_psi_1 (     alpha () ) ; 
  return 
    lambda() > 0 ?
    l_p2 / std::pow ( l_p1 , 1.5 ) : -1 * l_p2 / std::pow ( l_p1 , 1.5 ) ;
}
// ============================================================================
double Ostap::Math::LogGamma::kurtosis () const 
{ 
  const double l_p3 = gsl_sf_psi_n ( 3 , alpha () ) ; 
  const double l_p1 = gsl_sf_psi_1 (     alpha () ) ; 
  return l_p3 / ( l_p1 * l_p1) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::LogGamma::tag () const 
{ 
  static const std::string s_name = "LogGamma" ;
  return Ostap::Utils::hash_combiner ( s_name , m_nu , m_lambda , m_alpha ) ; 
}
// ============================================================================

// ============================================================================
// Beta
// ============================================================================
/*  constructor with all parameters
 *  @param alpha \f$\alpha\f$-parameter
 *  @param beta  \f$\beta\f$-parameter
 *  @param scale  scale-parameter
 *  @param shift  shift-parameter
 */
// ============================================================================
Ostap::Math::Beta::Beta
( const double p     ,
  const double q     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) )
  , PQ            ( p     , q     , "p"     , "q"     , typeid ( *this ) )
{}
// ============================================================================
// evaluate beta-distribution
// ============================================================================
double Ostap::Math::Beta::evaluate  ( const double x ) const
{
  const double tt = t ( x ) ;
  if ( tt < 0 || 1 < tt ) { return 0 ; }
  //
  const double pv = p () ;
  const double qv = q () ;
  //
  if ( 1 < pv && s_equal ( 1 + tt , 1 ) ) { return 0 ; }
  if ( 1 < qv && s_equal ( 0 + tt , 1 ) ) { return 0 ; }
  //
  // density diverges
  if ( pv < 1 && s_zero (     tt ) ) { return std::numeric_limits<double>::max () ; }
  if ( qv < 1 && s_zero ( 1 - tt ) ) { return std::numeric_limits<double>::max () ; }
  //
  const double pp = std::pow ( 0.0L + tt , pv - 1.0L ) ;
  const double qq = std::pow ( 1.0L - tt , qv - 1.0L ) ;
  //
  const double s  = std::abs ( scale () ) ;
  return pp * qq * inv_Beta () / s ; 
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::Beta::cdf ( const double x ) const 
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  const bool   ps   = Ostap::Math::is_positive ( scale() ) ;
  return ps ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::Beta::raw_cdf ( const double y ) const 
{
  if ( y <=  0 ) { return 0 ; }
  return beta_cdf ( y , p ()  , q ()  ) ;
}
// ============================================================================
// get integral 
// ============================================================================
double Ostap::Math::Beta::integral () const { return 1 ; }
// ============================================================================
// get integral 
// ============================================================================
double Ostap::Math::Beta::integral
( const double low  ,
  const double high ) const
{
  if      ( s_equal ( low  , high ) ) { return  0 ; }
  else if (           high < low    ) { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ===========================================================================
// raw moments (unscaled & unshifted)
// ===========================================================================
double Ostap::Math::Beta::raw_moment
( const unsigned short n ) const
{
  if ( !n ) { return 1 ; }
  //
  const double pv = p () ;
  const double qv = q () ;
  //  
  double result = 1 ;
  for ( unsigned short i = 0 ; i < n ; ++ i )
  { result *= ( pv + i ) / ( pv + qv + i ) ; }
  //
  return result ;
}  
// ===========================================================================
// mean
// ===========================================================================
double Ostap::Math::Beta::mean () const
{
  //
  const double pv = p () ;
  const double qv = q () ;
  //  
  const double  value = pv / ( pv + qv ) ;
  return x ( value ) ;
}
// ===========================================================================
// mode 
// ===========================================================================
double Ostap::Math::Beta::mode () const
{
  //
  const double pv = p () ;
  const double qv = q () ;
  //  
  const double  value =
    pv <  1            ? 0   :
    qv <  1            ? 1   :
    s_equal ( pv , 1 ) && s_equal ( qv , 1 ) ? 0.5 :
    ( pv - 1 ) / ( pv + qv - 2 ) ;
  //
  return x ( value ) ;
}
// ===========================================================================
// median 
// ===========================================================================
double Ostap::Math::Beta::median () const
{
  const double pv = p () ;
  const double qv = q () ;
  //  
  return s_equal ( pv  , qv  ) ? x ( 0.5 ) : quantile ( 0.5 ) ;
}
// ===========================================================================
// variance 
// ===========================================================================
double Ostap::Math::Beta::variance () const
{
  const double pv = p     () ;
  const double qv = q     () ;
  //
  const double sc = scale () ;
  //  
  const double value = pv * qv
    / ( std::pow ( pv + qv  , 2 ) * ( pv + qv  + 1 ) ) ;
  //
  return value  * sc * sc ;
}
// ===========================================================================
// RMS  
// ===========================================================================
double Ostap::Math::Beta::rms   () const
{ return std::sqrt ( variance ()  ) ; }
// ===========================================================================
// skewness 
// ===========================================================================
double Ostap::Math::Beta::skewness () const
{
  //
  const double pv = p     () ;
  const double qv = q     () ;
  //
  const double r = 2 * ( qv - pv )   
    * std::sqrt ( ( pv + qv + 1 ) / ( pv * qv ) )
    / ( pv + qv + 2 ) ;
  //
  return scale_sign ()  * r ;
}
// ===========================================================================
// (excess) kurtosis
// ===========================================================================
double Ostap::Math::Beta::kurtosis () const
{
  //
  const double pv = p     () ;
  const double qv = q     () ;
  //  
  const double pq  = pv * qv ;
  const double pq1 = pv + qv + 1 ;
  const double pq2 = pv + qv + 2 ;
  const double pq3 = pv + qv + 3 ;
  //
  return 6 * ( std::pow ( pv - pq , 2 ) * pq1 - pq * pq2 ) / ( pq * pq2 * pq3 ) ;
}
// ===========================================================================
// quantile \f$ 0 \le p \le 1 \f$
// ===========================================================================
double Ostap::Math::Beta::quantile ( const double z  ) const
{
  //
  if      ( z <  0 ) { return s_QUIETNAN ; }  
  else if ( z >  1 ) { return s_QUIETNAN ; }  
  else if ( z <= 0 ) { return x ( 0 )    ; } 
  else if ( z >= 1 ) { return x ( 1 )    ; } 
  //
  const double tt = beta_quantile ( z , p () , q ()  ) ;
  return x ( tt ) ;
}
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Beta::tag () const 
{ 
  static const std::string s_name = "Beta" ;
  return Ostap::Utils::hash_combiner ( s_name                ,
				       ShiftAndScale::tag () ,
				       PQ::tag ()            ) ;
}
// ============================================================================

// ============================================================================
// Beta' 
// ============================================================================
/*  constructor with all parameters 
 *  @param p     p-parameter 
 *  @param q      q-parameter 
 *  @param scale  scale-parameter
 *  @param shift  shift-parameter
 */
// ============================================================================
Ostap::Math::BetaPrime::BetaPrime
( const double p     ,
  const double q     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , PQ            ( p     , q     , "p"     , "q"     , typeid ( *this ) )
{}
// ============================================================================
// evaluate beta'-distributions 
// ============================================================================
double Ostap::Math::BetaPrime::evaluate ( const double x ) const 
{
  //
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const double s = std::abs ( scale () ) ;
  //
  return inv_Beta ()
    * pow_ratio_a1 ( y , p () )
    * pow_ratio_a2 ( y , q () ) / ( s * y ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BetaPrime::cdf ( const double x ) const
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ; 
}
// ============================================================================
// "raw" cdf (for positive scale)
// ============================================================================
double Ostap::Math::BetaPrime::raw_cdf ( const double y ) const 
{
  if ( y <= 0 ) { return 0 ; }
  const double z = y / ( 1 + y ) ;  
  return beta_inc ( p () , q () , z ) ;  
}
// ============================================================================
double Ostap::Math::BetaPrime::integral
( const double low  , 
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high < low             ) { return -integral ( high , low ) ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::integral () const { return 1 ; } 
// ============================================================================
double Ostap::Math::BetaPrime::mean () const 
{
  const double qv = q () ;
  if ( qv < 1 || s_equal ( qv , 1 ) ) { return s_QUIETNAN  ; } 
  //
  const double pv = p () ;
  const double r = pv / ( qv - 1 ) ;
  return x ( r ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::mode () const 
{
  //
  const double pv = p () ;
  const double qv = q () ;
  //
  if ( pv < 1 ) { return 0 ; }
  //
  const double r = ( pv - 1 ) / ( qv + 1 ) ;
  return x ( r ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::variance () const 
{
  const double qv = q () ;
  if ( qv < 2 || s_equal ( qv , 2 ) ) { return s_QUIETNAN  ; } 
  //
  const double pv = p () ;
  //
  const double q2 = qv - 2 ;
  const double q1 = qv - 1 ;
  const double r  = pv * ( pv + qv - 1 ) / ( q2 * q1 * q1 ) ;
  
  const double s = std::abs ( scale () ) ;
  //
  return r * s * s ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::rms () const 
{
  const double qv = q () ;
  if ( qv < 2 || s_equal ( qv , 2 ) ) { return s_QUIETNAN  ; } 
  //
  return std::sqrt ( variance () ) ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::skewness  () const 
{
  const double qv = q () ;
  if ( qv < 3 || s_equal ( qv , 3 ) ) { return s_QUIETNAN  ; } 
  //
  const double pv = p () ;
  //
  const double q3 = qv - 3 ;
  const double q2 = qv - 2 ;
  //
  const double r = 2 * ( 2 * pv + qv - 1 ) / q3 * std::sqrt( q2 / ( pv * ( pv + qv - 1 ) ) ) ;
  //
  return scale_sign () * r ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::kurtosis  () const 
{
  const double qv = q () ;
  if ( qv < 4 || s_equal ( qv , 4 ) ) { return s_QUIETNAN  ; } 
  //
  const double pv = p () ;
  //
  const double q1 = qv - 1 ;
  const double q2 = qv - 2 ;
  const double q3 = qv - 3 ;
  const double q4 = qv - 4 ;
  //
  const double c0 = ( pv + qv - 1 ) ; 
  const double c1 = c0 * pv * ( 3 * qv - 11 ) ;
  const double c2 = q1 * q1 * q2 ;
  const double c3 = c0 * pv * q3 * q4 ;
  //
  return  c1 * c2  / c3 ;
}
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BetaPrime::tag () const 
{ 
  static const std::string s_name = "BetaPrime" ;
  return Ostap::Utils::hash_combiner ( s_name                ,
				       ShiftAndScale::tag () ,
				       PQ::tag ()            ) ;
}
// ============================================================================

// ============================================================================
/*  constructor from all parameters 
 *  @param d1  d1-parameter d1>0
 *  @param d2  d2-parameter d2>0
 *  @param scale  scale-parameter
 *  @param shift  shift-parameter
 */
// ============================================================================
Ostap::Math::FDistribution::FDistribution
( const double d1    ,
  const double d2    ,
  const double scale ,
  const double shift )
  : m_betap ( 0.5 * d1 , 0.5 * d2 )
  , m_scale ( scale )
  , m_shift ( shift )
{}
// ============================================================================
bool Ostap::Math::FDistribution::setScale ( const double value ) 
{
  if ( s_equal ( value , m_scale  ) ) { return false ; }
  m_scale  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::FDistribution::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift  ) ) { return false ; }
  m_shift  = value ;
  return true ;
}
// ============================================================================
// evaluate F-distribution
// ============================================================================
double Ostap::Math::FDistribution::pdf ( const double x ) const 
{
  //
  const double s = d1 () / d2 () ;
  const double z = ( x - m_shift ) / m_scale ;
  //
  return m_betap.pdf ( z * s ) * s / m_scale ;
}
// ============================================================================
double Ostap::Math::FDistribution::cdf ( const double x ) const 
{
  //
  const double z = ( x - m_shift ) / m_scale ;
  const double s = d1 () / d2 () ;
  return m_betap.cdf ( z * s ) ;
}
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::FDistribution::tag () const 
{ 
  static const std::string s_name = "FDistribution"  ;
  return Ostap::Utils::hash_combiner ( s_name         ,
				       m_betap.tag () ,
				       m_scale        ,
				       m_shift        ) ; 
}
// ============================================================================
double Ostap::Math::FDistribution::mean () const 
{
  const double s = m_scale / ( d1 () / d2 () ) ;
  return m_shift + s * m_betap.mean () ;  
}
// ============================================================================
double Ostap::Math::FDistribution::mode () const 
{
  const double s = m_scale / ( d1 () / d2 () ) ;
  return m_shift + s * m_betap.mode() ;
}
// ============================================================================
double Ostap::Math::FDistribution::variance () const 
{
  const double s = m_scale / ( d1 () / d2 () ) ;
  return s * s * m_betap.variance () ;  
}
// ===========================================================================
double Ostap::Math::FDistribution::sigma () const 
{
  if ( m_betap.beta() <= 2 || s_equal ( m_betap.beta() , 2 ) ) { return -1.e+9 ; }  
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::FDistribution::integral () const { return 1 ; }  
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::FDistribution::integral
( const double low  , 
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================

// ============================================================================
// Generalized Beta' 
// ============================================================================
/*  constructor with all parameters
 *  @param a shape a-parameter
 *  @param p shape p-parameter
 *  @param q shape q-parameter
 *  @param scale   scale-parameter
 *  @param shift   shift-parameter
 */
// ============================================================================
Ostap::Math::GenBetaPrime::GenBetaPrime 
( const double a     , 
  const double p     ,
  const double q     , 
  const double scale , 
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , PQ            ( p     , q     , "p"     , "q"     , typeid ( *this )         )    
  , m_a  ( a                                , "a"     , typeid ( *this ) , false )
{}
// ============================================================================
// evaluate Beta'-distribution
// ============================================================================
double Ostap::Math::GenBetaPrime::evaluate ( const double x ) const 
{
  //
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  //
  const double ap = av * pv ;
  //
  if ( s_zero ( y ) && 1 < ap  ) { return 0 ; }
  //
  const double s = std::abs ( scale () ) ;
  //
  return pv * inv_Beta () 
    * std::pow ( pow_ratio_a1 ( y , av ) , pv )
    * std::pow ( pow_ratio_a2 ( y , av ) , qv ) / ( s * y ) ;
}
// ============================================================================
// xmin: can be infinite or NaN
// ============================================================================
double Ostap::Math::GenBetaPrime::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? shift () : s_NEGINF ;   }
// ============================================================================
// xmax: can be infinite or NaN
// ============================================================================
double Ostap::Math::GenBetaPrime::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : shift () ;   }
// ============================================================================
double Ostap::Math::GenBetaPrime::integral () const { return 1 ; }
// ============================================================================
// "raw" cdf (for positive scale)   
// ============================================================================
double Ostap::Math::GenBetaPrime::raw_cdf ( const double x ) const
{
  const double y = t ( x ) ;
  if ( y <= 0 ) { return 0 ; }
  //
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  //
  const double z = pow_ratio_a1 ( y , av ) ;
  //    
  return beta_inc ( pv , qv , z ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::GenBetaPrime::cdf ( const double x ) const
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale() ) ? rcdf : 1 - rcdf ;
}

// ============================================================================
double Ostap::Math::GenBetaPrime::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::GenBetaPrime::mean () const 
{
  const double av = a () ;
  const double qv = q () ;
  //
  const double aq = av * qv ;
  //
  if ( aq <= 1 || s_equal ( aq , 1 ) ) { return s_QUIETNAN ; }
  //
  const double pv   = p () ;
  const double inva = 1 / av ;
  //
  const double logm =
    Ostap::Math::lgamma ( pv + inva ) - Ostap::Math::lgamma ( pv )  + 
    Ostap::Math::lgamma ( qv - inva ) - Ostap::Math::lgamma ( qv ) ;
  //
  const double r = std::exp ( logm ) ;
  return x ( r ) ;
}
// ============================================================================
double Ostap::Math::GenBetaPrime::mode () const 
{
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  //
  const double ap = av * pv ; 
  if ( ap <= 1 ) { return x ( 0 ) ; }
  //
  const double aq = av * qv ; 
  const double  z = std::pow ( ( ap - 1 ) / ( aq + 1 ) , 1.0 / av ) ;
  //
  return x ( z ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBetaPrime::tag () const 
{ 
  static const std::string s_name = "GenBetaPrime" ;
  return Ostap::Utils::hash_combiner ( s_name                ,
				       m_a .tag ()           ,
				       ShiftAndScale::tag () ,
				       PQ::tag ()            ) ;
}
// ============================================================================

// ============================================================================
// the most general form
// ============================================================================
Ostap::Math::GenBeta::GenBeta
( const double a     , // shape  
  const double r     , // r = beta2/beta1 
  const double p     , // shape 
  const double q     , // shape 
  const double scale , 
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , PQ            ( p     , q     , "p"     , "q"     , typeid ( *this )         )    
  , m_a           ( a     ,                   "a"     , typeid ( *this ) , false ) 
  , m_r           ( r     ,                   "r"     , typeid ( *this )         )
{}
// ============================================================================
// parameter c: beta2/beta1 = r = ( ( 1 - c ) / c ) ^ { 1 / a } 
// ============================================================================
double Ostap::Math::GenBeta::c () const
{
  const double av = a () ;
  const double rv = r () ;
  return pow_ratio_a2 ( rv , av ) ;
}
// ============================================================================
// evaluate GenBeta-distribution
// ============================================================================
double Ostap::Math::GenBeta::evaluate ( const double x ) const 
{
  //
  const double y = t ( x ) ;
  if      ( y < 0        ) { return 0 ; }
  //
  const bool fr = finite_range () ;
  if      (  fr && 1 < y ) { return 0 ; }
  else if ( !fr && y < 1 ) { return 0 ; }
  //
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  const double rv = r () ;
  //
  const double s  = std::abs ( scale () ) ;
  //
  if ( fr )
  {
    const double ap = av * pv ; 
    const double yr = y / rv  ; 
    return std::abs ( av / s ) * inv_Beta ()      
      * std::pow ( pow_ratio_a2 ( 1 / rv , av ) , -pv )       
      * std::pow ( 0.0L + y                     , ap - 1.0L )
      * std::pow ( 1.0L - std::pow ( y   , av ) , qv - 1.0L )
      * std::pow ( pow_ratio_a2 ( yr     , av ) , pv + qv   ) ;      
  }
  //
  const double za  = std::pow ( 1 / y , -av ) ;
  const double yr  = y / rv  ; 
  return std::abs ( av / s ) * inv_Beta ()
    * std::pow ( pow_ratio_a2 ( 1 / rv , av ) , -pv )
    * std::pow ( 0.0L + za , pv     )
    * std::pow ( 1.0L - za , qv - 1 )
    * std::pow ( pow_ratio_a2 ( yr     , av ) , pv + qv   ) / y  ;
}
// ============================================================================
// xmin: can be infinite or NaN
// ============================================================================
double Ostap::Math::GenBeta::xmin () const
{
  const bool fr = finite_range () ;
  const bool ps = Ostap::Math::is_positive ( scale () ) ;
  //
  if      (  fr &&  ps ) { return x ( 0 ) ; } 
  else if (  fr && !ps ) { return x ( 1 ) ; }  
  else if ( !fr &&  ps ) { return x ( 1 ) ; }
  //
  return s_NEGINF ;
}
// ============================================================================
// xmax: can be infinite or NaN
// ============================================================================
double Ostap::Math::GenBeta::xmax () const
{
  const bool fr = finite_range () ;
  const bool ps = Ostap::Math::is_positive ( scale () ) ;
  //
  if      (  fr &&  ps ) { return x ( 1 ) ; } 
  else if (  fr && !ps ) { return x ( 0 ) ; }  
  else if ( !fr && !ps ) { return x ( 1 ) ; }  
  //
  return s_POSINF ; 
}
// ============================================================================
// evaluate the integral 
// ============================================================================
double Ostap::Math::GenBeta::integral () const { return 1 ; }
// ============================================================================
// evaluate integral 
// ============================================================================
double Ostap::Math::GenBeta::integral
( const double low  ,
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return - integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the CDF
// ============================================================================
double Ostap::Math::GenBeta::cdf 
( const double x    ) const
{
  const bool   ps  = Ostap::Math::is_positive ( scale () ) ;
  const double y   = t ( x ) ;
  return ps ? raw_cdf ( y ) : 1.0 - raw_cdf ( y ) ;
}
// ============================================================================
// get the raw CDF for positive scale 
// ============================================================================
double Ostap::Math::GenBeta::raw_cdf 
( const double y ) const
{
  //
  if      ( y <= 0        ) { return 0 ; }
  //
  const bool fr = finite_range () ;
  if      (  fr && 1 <= y ) { return 1 ; }
  else if ( !fr && y <= 1 ) { return 0 ; }
  //
  const double av = a ()   ;
  const double rv = r ()   ;
  const double yr = y / rv ;
  //
  const double z = pow_ratio_a1 ( yr , av ) / pow_ratio_a2 ( rv , av ) ;
  //
  return z <= 0 ? 0 : z >= 1 ? 1 : beta_inc ( p () , q () , z ) ;
} 
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBeta::tag () const 
{ 
  static const std::string s_name = "GenBeta" ;
  return Ostap::Utils::hash_combiner ( s_name                ,
				       ShiftAndScale::tag () ,
				       PQ           ::tag () , 
				       m_a           .tag () ,
				       m_r           .tag () ) ;
}
// ============================================================================

// ============================================================================
// Generalized beta of first kind 
// ============================================================================
Ostap::Math::GenBeta1::GenBeta1
( const double a     , // shape 
  const double p     , // shape 
  const double q     , // shape 
  const double scale , // scale 
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , PQ            ( p     , q     , "p"     , "q"     , typeid ( *this )         )    
  , m_a           ( a     ,                   "a"     , typeid ( *this ) , false ) 
{}
// ============================================================================
// xmin: can be infinite or NaN
// ============================================================================
double Ostap::Math::GenBeta1::xmin () const
{
  const bool fr = finite_range () ;
  const bool ps = Ostap::Math::is_positive ( scale () ) ;
  //
  if      (  fr &&  ps ) { return x ( 0 ) ; } 
  else if (  fr && !ps ) { return x ( 1 ) ; }  
  else if ( !fr &&  ps ) { return x ( 1 ) ; }
  //
  return s_NEGINF ;
}
// ============================================================================
// xmax: can be infinite or NaN
// ============================================================================
double Ostap::Math::GenBeta1::xmax () const
{
  const bool fr = finite_range () ;
  const bool ps = Ostap::Math::is_positive ( scale () ) ;
  //
  if      (  fr &&  ps ) { return x ( 1 ) ; } 
  else if (  fr && !ps ) { return x ( 0 ) ; }  
  else if ( !fr && !ps ) { return x ( 1 ) ; }  
  //
  return s_POSINF ; 
}
// ============================================================================
// evaluate GenBeta1-distribution
// ============================================================================
double Ostap::Math::GenBeta1::evaluate
( const double x ) const 
{
  //
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const bool fr = finite_range () ;
  if      (  fr && 1 < y ) { return 0 ; }
  else if ( !fr && y < 1 ) { return 0 ; }
  //
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  //
  const double s  = std::abs ( scale () ) ;
  //
  const double ya = std::pow ( y , av ) ;
  return std::abs ( av / s ) * inv_Beta ()
    * std::pow ( 0.0L + ya , pv + 0.0L )
    * std::pow ( 0.1L - ya , qv - 1.0L ) / y ;          
}
// ============================================================================
// mean-value @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta1::mean () const
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  return x ( m1 ) ;
}
// ============================================================================
// variance @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta1::variance () const
{
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  //
  const double ss = scale () ; 
  return ss * ss * Ostap::Math::variance ( m2 , m1 ) ;
}
// ============================================================================
// rms @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta1::rms () const
{
  const double vv = variance () ;
  if ( !std::isfinite ( vv ) ) { return vv ; }
  return std::sqrt ( vv ) ;
}
// ============================================================================
// skewness @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta1::skewness () const
{
  const double m3 = raw_moment ( 3 ) ;
  if ( !std::isfinite ( m3 ) ) { return m3 ; }
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  //
  return scale_sign () * Ostap::Math::skewness ( m3 , m2 , m1 ) ;
}
// ============================================================================
// (excess) kurtosis @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta1::kurtosis () const
{
  const double m4 = raw_moment ( 4 ) ;
  if ( !std::isfinite ( m4 ) ) { return m4 ; }
  const double m3 = raw_moment ( 3 ) ;
  if ( !std::isfinite ( m3 ) ) { return m3 ; }
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  //
  return Ostap::Math::kurtosis ( m4 , m3 , m2 , m1 ) ;
}
// ============================================================================
// raw moment  of unscaled/unshifted distribution
// ============================================================================
double Ostap::Math::GenBeta1::raw_moment
( const unsigned short k ) const
{
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  //
  const double vv = pv + k / av ;
  if ( vv <= 0 || s_zero ( vv ) ) { return s_QUIETNAN ; } 
  //
  return Ostap::Math::beta  ( vv , qv  ) * inv_Beta () ; 
}
// ============================================================================
// evaluate integral 
// ============================================================================
double Ostap::Math::GenBeta1::integral () const { return 1 ; }
// ============================================================================
// evaluate integral 
// ============================================================================
double Ostap::Math::GenBeta1::integral
( const double low  ,
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return - integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the CDF
// ============================================================================
double Ostap::Math::GenBeta1::cdf 
( const double x    ) const
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the CDF 
// ============================================================================
double Ostap::Math::GenBeta1::raw_cdf 
( const double y  ) const
{
  //
  if      ( y <= 0        ) { return 0 ; }
  //
  const bool fr = finite_range () ;
  if      (  fr && 1 <= y ) { return 1 ; }
  else if ( !fr && y <= 1 ) { return 0 ; }
  //
  const double av = a ()   ;
  //
  const double z = std::pow ( y , av ) ;
  return z <= 0 ? 0 : z >= 1 ? 1 : beta_inc ( p () , q () , z ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBeta1::tag () const 
{ 
  static const std::string s_name = "GenBeta1" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       ShiftAndScale::tag () ,
				       PQ           ::tag () , 
				       m_a           .tag () ) ;
}
// ============================================================================

// ============================================================================
// Generalized beta of second kind 
// ============================================================================
Ostap::Math::GenBeta2::GenBeta2
( const double a     , // shape 
  const double p     , // shape 
  const double q     , // shape 
  const double scale , // scale 
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , PQ            ( p     , q     , "p"     , "q"     , typeid ( *this )         )    
  , m_a           ( a     ,                   "a"     , typeid ( *this ) , false ) 
{}
// ============================================================================
// xmin
// ============================================================================
double Ostap::Math::GenBeta2::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::GenBeta2::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 1 ) ; }
// ============================================================================
// evaluate GenBeta2-distribution
// ============================================================================
double Ostap::Math::GenBeta2::evaluate
( const double x ) const 
{
  const double y = t ( x ) ;
  if ( y < 0  ) { return 0 ; }
  //
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  //
  const double ap = av * pv ;  
  //
  if ( 1 < ap && s_zero ( y ) ) { return 0 ; }
  //
  const double s = std::abs ( scale () ) ;
  //
  return std::abs ( av / s ) * inv_Beta () 
    * std::pow ( pow_ratio_a1 ( y , av ) , pv + 0.0L )
    * std::pow ( pow_ratio_a2 ( y , av ) , qv + 0.0L ) / y ;
}
// ============================================================================
// raw moment  of unscaled/unshifted distribution
// ============================================================================
double Ostap::Math::GenBeta2::raw_moment
( const unsigned short k ) const
{
  const double av = a () ;
  const double pv = p () ;
  const double qv = q () ;
  //
  const double pa = pv + k / av ;
  const double qa = qv - k / av ;
  //
  return
    pa <= 0 || s_zero ( pa      ) ? s_QUIETNAN : 
    qa <= 0 || s_zero ( qa      ) ? s_QUIETNAN :
    Ostap::Math::beta ( pa , qa ) * inv_Beta () ; 
}
// ============================================================================
// mean-value @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta2::mean () const
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  return x ( m1 ) ;
}
// ============================================================================
// variance @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta2::variance () const
{
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  //
  const double ss = scale () ; 
  return ss * ss * Ostap::Math::variance ( m2 , m1 ) ;
}
// ============================================================================
// rms @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta2::rms () const
{
  const double vv = variance () ;
  if ( !std::isfinite ( vv ) ) { return vv ; }
  return std::sqrt ( vv ) ;
}
// ============================================================================
// skewness @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta2::skewness () const
{
  const double m3 = raw_moment ( 3 ) ;
  if ( !std::isfinite ( m3 ) ) { return m3 ; }
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  //
  return scale_sign () * Ostap::Math::skewness ( m3 , m2 , m1 ) ;
}
// ============================================================================
// (excess) kurtosis @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta2::kurtosis () const
{
  const double m4 = raw_moment ( 4 ) ;
  if ( !std::isfinite ( m4 ) ) { return m4 ; }
  const double m3 = raw_moment ( 3 ) ;
  if ( !std::isfinite ( m3 ) ) { return m3 ; }
  const double m2 = raw_moment ( 2 ) ;
  if ( !std::isfinite ( m2 ) ) { return m2 ; }
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  //
  return Ostap::Math::kurtosis ( m4 , m3 , m2 , m1 ) ;
}
// ============================================================================
// evaluate integral 
// ============================================================================
double Ostap::Math::GenBeta2::integral () const { return 1 ; }
// ============================================================================
// evaluate integral 
// ============================================================================
double Ostap::Math::GenBeta2::integral
( const double low  ,
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return - integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the CDF
// ============================================================================
double Ostap::Math::GenBeta2::cdf 
( const double x    ) const
{
  const bool   ps   = Ostap::Math::is_positive ( scale () ) ;
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return ps ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the raw CDF for positive scale 
// ============================================================================
double Ostap::Math::GenBeta2::raw_cdf 
( const double y ) const
{
  if ( y <= 0 ) { return 0 ; }
  //
  const double z = pow_ratio_a1 ( y , a()  ) ; 
  return z <= 0 ? 0 : z >= 1 ? 1 : beta_inc ( p () , q () , z ) ;
} 
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBeta2::tag () const 
{ 
  static const std::string s_name = "GenBeta2" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       ShiftAndScale::tag () ,
				       PQ           ::tag () , 
				       m_a           .tag () ) ;
}
// ============================================================================

// ============================================================================
// Landau
// ============================================================================
/*  constructor with all parameters 
 */
// ============================================================================
Ostap::Math::Landau::Landau
( const double scale , 
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
{}
// ============================================================================
// evaluate Landau-distributions 
// ============================================================================
double Ostap::Math::Landau::evaluate ( const double x ) const 
{
  const double y = t ( x ) ;
  return gsl_ran_landau_pdf ( y ) / std::abs ( scale () ) ;
}
// ============================================================================
namespace 
{
/* Not needed yet */
/* This function is a translation from the original Fortran of the
 * CERN library routine DISLAN, the integral from -inf to x of the
 * Landau p.d.f.
 */
//static
//double
//gsl_ran_landau_dislan(const double x)
double _dislan(const double x)
  {
  static double P1[5] =
    {
      0.2514091491E0, -0.6250580444E-1,
      0.1458381230E-1, -0.2108817737E-2,
      0.7411247290E-3
    };

  static double P2[4] =
    {
      0.2868328584E0, 0.3564363231E0,
      0.1523518695E0, 0.2251304883E-1
    };

  static double P3[4] =
    {
      0.2868329066E0, 0.3003828436E0,
      0.9950951941E-1, 0.8733827185E-2
    };

  static double P4[4] =
    {
      0.1000351630E1, 0.4503592498E1,
      0.1085883880E2, 0.7536052269E1
    };

  static double P5[4] =
    {
      0.1000006517E1, 0.4909414111E2,
      0.8505544753E2, 0.1532153455E3
    };

  static double P6[4] =
    {
      0.1000000983E1, 0.1329868456E3,
      0.9162149244E3, -0.9605054274E3
    };

  static double Q1[5] =
    {
      1.0, -0.5571175625E-2,
      0.6225310236E-1, -0.3137378427E-2,
      0.1931496439E-2
    };

  static double Q2[4] =
    {
      1.0, 0.6191136137E0,
      0.1720721448E0, 0.2278594771E-1
    };

  static double Q3[4] =
    {
      1.0, 0.4237190502E0,
      0.1095631512E0, 0.8693851567E-2
    };

  static double Q4[4] =
    {
      1.0, 0.5539969678E1,
      0.1933581111E2, 0.2721321508E2
    };

  static double Q5[4] =
    {
      1.0, 0.5009928881E2,
      0.1399819104E3, 0.4200002909E3
    };

  static double Q6[4] =
    {
      1.0, 0.1339887843E3,
      0.1055990413E4, 0.5532224619E3
    };

  static double A1[3] =
    {
      -0.4583333333E0, 0.6675347222E0, -0.1641741416E1
    };

  static double A2[3] =
    {
      1.0, -0.4227843351E0, -0.2043403138E1
    };

  double U, V, DISLAN;

  V = x;
  if (V < -5.5)
    {
      U = exp(V + 1);
      DISLAN = 0.3989422803 * exp( -1 / U) * sqrt(U) *
               (1 + (A1[0] + (A1[1] + A1[2] * U) * U) * U);
    }
  else if (V < -1)
    {
      U = exp( -V - 1);
      DISLAN = (exp( -U) / sqrt(U)) *
               (P1[0] + (P1[1] + (P1[2] + (P1[3] + P1[4] * V) * V) * V) * V) /
               (Q1[0] + (Q1[1] + (Q1[2] + (Q1[3] + Q1[4] * V) * V) * V) * V);
    }
  else if (V < 1)
    {
      DISLAN = (P2[0] + (P2[1] + (P2[2] + P2[3] * V) * V) * V) /
               (Q2[0] + (Q2[1] + (Q2[2] + Q2[3] * V) * V) * V);
    }
  else if (V < 4)
    {
      DISLAN = (P3[0] + (P3[1] + (P3[2] + P3[3] * V) * V) * V) /
               (Q3[0] + (Q3[1] + (Q3[2] + Q3[3] * V) * V) * V);
    }
  else if (V < 12)
    {
      U = 1 / V;
      DISLAN = (P4[0] + (P4[1] + (P4[2] + P4[3] * U) * U) * U) /
               (Q4[0] + (Q4[1] + (Q4[2] + Q4[3] * U) * U) * U);
    }
  else if (V < 50)
    {
      U = 1 / V;
      DISLAN = (P5[0] + (P5[1] + (P5[2] + P5[3] * U) * U) * U) /
               (Q5[0] + (Q5[1] + (Q5[2] + Q5[3] * U) * U) * U);
    }
  else if (V < 300)
    {
      U = 1 / V;
      DISLAN = (P6[0] + (P6[1] + (P6[2] + P6[3] * U) * U) * U) /
               (Q6[0] + (Q6[1] + (Q6[2] + Q6[3] * U) * U) * U);
    }
  else
    {
      U = 1 / (V - V * log(V) / (V + 1));
      DISLAN = 1 - (A2[0] + (A2[1] + A2[2] * U) * U) * U;
    }

  return DISLAN;
  }
}
// ============================================================================
// evaluate Landau-CDF 
// ============================================================================
double Ostap::Math::Landau::cdf ( const double x ) const 
{
  const double y    = t ( x )       ;
  const double rcdf = _dislan ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
double Ostap::Math::Landau::integral
( const double low  , 
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0 ; }  
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Landau::tag () const 
{ 
  static const std::string s_name = "Landau" ;
  return Ostap::Utils::hash_combiner ( s_name , ShiftAndScale::tag () ) ;
}
// ============================================================================



// ============================================================================
// Weibull
// ============================================================================
/*  constructor from all parameters    
 *  @param scale the scale parameter "lambda"  >0 
 *  @param shape the shape parameter "k"       >0
 *  @param shift the shift parameter "x0"
 */
// ============================================================================
Ostap::Math::Weibull::Weibull
(  const double scale ,  
   const double shape ,  
   const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_k           ( shape                   , "k"     , typeid ( *this ) )
{}    
// ============================================================================
// evaluate Weibull-distributions
// ============================================================================
double Ostap::Math::Weibull::evaluate ( const double x ) const 
{
  const double y = t ( x ) ;
  if ( y <= 0 ) { return 0 ; }
  //
  const double kv = k () ;
  const double ss = std::abs ( scale () ) ;
  //
  return ( kv / ss )
    * std::pow (   y            , kv - 1 )
    * std::exp ( - std::pow ( y , kv     ) ) ;  
}
// ============================================================================
double Ostap::Math::Weibull::cdf ( const double x ) const 
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
double Ostap::Math::Weibull::raw_cdf ( const double y ) const 
{
  if ( y <= 0 ) {  return 0 ; }
  return 1 - std::exp ( - std::pow ( y  , k () ) ) ;
}
// ============================================================================
double Ostap::Math::Weibull::integral  () const { return 1 ; } 
// ============================================================================
double Ostap::Math::Weibull::integral 
( const double low  , const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high < low             ) { return -integral ( high, low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;  
}
// ============================================================================
// raw moment  of unscaled/unshifted distribution
// ============================================================================
double Ostap::Math::Weibull::raw_moment
( const unsigned short n ) const
{ return Ostap::Math::gamma ( 1.0 + n / k () ) ; }
// ============================================================================
// the mean 
// ============================================================================
double Ostap::Math::Weibull::mean () const 
{
  const double m1 = raw_moment ( 1 ) ;
  return x ( m1 ) ;
}
// ============================================================================
// the mode 
// ============================================================================
double Ostap::Math::Weibull::mode () const 
{
  const double kv = k () ; 
  const double r = ( 1 < kv ) ?  std::pow ( ( kv - 1 ) / kv , 1 / kv ) : 0.0 ;
  return x ( r ) ;
}
// ============================================================================
// the median 
// ============================================================================
double Ostap::Math::Weibull::median () const 
{
  const double kv = k () ;
  const double r = std::pow ( s_ln2 , 1/ kv ) ;
  return x ( r ) ;
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::Weibull::variance () const 
{
  const double ss = scale () ;
  return ss * ss * Ostap::Math::variance ( raw_moment ( 2 ) ,
					   raw_moment ( 1 ) ) ;
}
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::Weibull::skewness () const 
{
  return scale_sign () * Ostap::Math::skewness ( raw_moment ( 3 ) ,
						 raw_moment ( 2 ) ,
						 raw_moment ( 1 ) ) ;
}
// ============================================================================
// kurtosis
// ============================================================================
double Ostap::Math::Weibull::kurtosis () const 
{
  return Ostap::Math::kurtosis ( raw_moment ( 4 ) ,
				 raw_moment ( 3 ) ,
				 raw_moment ( 2 ) ,
				 raw_moment ( 1 ) ) ;
}
// ============================================================================
// rms  
// ============================================================================
double Ostap::Math::Weibull::rms () const { return  std::sqrt ( variance () ) ; }
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Weibull::tag () const 
{ 
  static const std::string s_name = "Weibull" ;
  return Ostap::Utils::hash_combiner ( s_name                ,
				       ShiftAndScale::tag () ,
				       m_k.tag()             ) ;
}
// ============================================================================


// namespace 
// {
//   //
//   inline double shash ( const double x   , 
//                         const double eps , 
//                         const double dlt ) 
//   {
//     const double y = eps + dlt * std::asinh ( x ) ;
//     return 
//       (     GSL_LOG_DBL_MAX < y ) ?    s_INFINITY :
//       ( -1* GSL_LOG_DBL_MAX > y ) ? -1*s_INFINITY : std::sinh ( y ) ;
//   }
//   //
// }
// ============================================================================



// ============================================================================
// contructor with all parameters
// ============================================================================
Ostap::Math::Rice::Rice
( const double nu       , 
  const double varsigma ,
  const double shift    ) 
  : m_nu         ( std::abs ( nu       ) ) 
  , m_varsigma   ( std::abs ( varsigma ) ) 
  , m_shift      ( shift ) 
  , m_workspace ()
{}
// ============================================================================
// update nu
// ============================================================================
bool Ostap::Math::Rice::setNu ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_nu , avalue ) ) { return false ; }
  m_nu = avalue ;
  return true ;
}
// ============================================================================
// update varsigma
// ============================================================================
bool Ostap::Math::Rice::setVarsigma ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_varsigma , avalue ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
// update shift 
// ============================================================================
bool Ostap::Math::Rice::setShift ( const double value )
{
  if ( s_equal ( m_shift , value ) ) { return false ; }
  m_shift = value ;
  return true ;
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::Rice::mean       () const 
{
  return m_shift + m_varsigma * s_sqrt_pi_2 * 
    Ostap::Math::laguerre_q ( 0.5 , -0.5 * std::pow ( m_nu / m_varsigma , 2 ) ) ;
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::Rice::variance   () const 
{
  const double s2 = m_varsigma * m_varsigma ;
  const double n2 = m_nu       * m_nu       ;
  const double l  = Ostap::Math::laguerre_q ( 0.5 , - 0.5 * n2 / s2 ) ;
  return 2 * s2 + n2 - s_pi_2 * s2 * l * l ;  
}
// ============================================================================
/// evaluate the function
// ============================================================================
double Ostap::Math::Rice::evaluate ( const double x ) const 
{
  if ( x <= m_shift ) { return 0 ; }
  const double s2 = m_varsigma * m_varsigma ;
  const double n2 = m_nu       * m_nu       ;
  const double dx = x - m_shift ;
  return 
    ( dx / s2 ) * std::exp ( -0.5 *  (dx * dx + n2 ) / s2 ) * 
    Ostap::Math::bessel_In ( 0 , dx * m_nu / s2 ) ; 
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Rice::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Rice::integral   
( const double low  ,
  const double high ) const
{
  // 
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low )        { return - integral ( high , low ) ; }
  //
  if ( high <= m_shift ) { return 0 ; }
  //
  const double xmin = std::max ( low , m_shift ) ;
  //
  const double x0 = m_nu + m_shift ;
  //
  // split at x0 
  if ( xmin < x0 && x0 < high ) 
  { return integral ( xmin , x0 ) + integral ( x0 , high ) ; }
  //
  // split at x1 
  const double x1 = x0 - 3 * m_varsigma ;
  if ( xmin < x1 && x1 < high ) 
  { return integral ( xmin , x1 ) + integral ( x1 , high ) ; }
  //
  // split at x2 
  const double x2 = x0 + 4 * m_varsigma ;
  if ( xmin < x2 && x2 < high ) 
  { return integral ( xmin , x2 ) + integral ( x2 , high ) ; }
  //
  // split at x3 
  const double x3 = x0 + 10 * m_varsigma ;
  if ( xmin < x3 && x3 < high ) 
  { return integral ( xmin , x3 ) + integral ( x3 , high ) ; }
  //
  //
  const bool in_tail = ( x0 + 10 * m_varsigma ) <= low ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Rice> s_integrator ;
  static const char s_message[] = "Integral(Rice)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,                             // low & high edges
      workspace ( m_workspace ) ,                  // workspace
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
std::size_t Ostap::Math::Rice::tag () const 
{ 
  static const std::string s_name = "Rice" ;
  return Ostap::Utils::hash_combiner ( s_name , m_nu , m_varsigma , m_shift ) ;
}
  

// ============================================================================
// Constructor 
// ============================================================================
Ostap::Math::GenInvGauss::GenInvGauss 
( const double theta , 
  const double eta   ,
  const double p     ,  
  const double shift )
  : ShiftAndScale ( eta   , shift , "eta" , "shift" , typeid ( *this ) , false )
  , m_theta       ( theta                 , "theta" , typeid ( *this ) )
  , m_p           ( p                     , "p"     , typeid ( *this ) )
{
  m_iKp_theta = 1.0 / Ostap::Math::bessel_Knu ( m_p.value()  , m_theta.value()  ) ;
}
// ============================================================================
/// evaluate the function 
// ============================================================================
double Ostap::Math::GenInvGauss::evaluate ( const double x ) const 
{
  const double y = t ( x ) ;
  if ( y <= 0 ) { return 0 ; }
  //
  const double pv = p     () ;
  const double tv = theta () ;
  //
  const double ss = std::abs ( scale () ) ;
  return 0.5 * m_iKp_theta / ss
    * std::pow ( y , pv - 1 )  
    * std::exp ( -0.5 * tv * ( y + 1 / y ) )  ; 
}
// ============================================================================
// set parameter theta 
// ============================================================================
bool Ostap::Math::GenInvGauss::setTheta  ( const double value ) 
{
  if ( !m_theta.setValue ( value ) ) { return false ; }
  m_iKp_theta = 1.0 / Ostap::Math::bessel_Knu ( p () , theta ()) ;  
  return true ;
}
// ============================================================================
// set parameter P
// ============================================================================
bool Ostap::Math::GenInvGauss::setP  ( const double value ) 
{
  if ( !m_p.setValue ( value ) ) { return false ; }
  m_iKp_theta = 1.0 / Ostap::Math::bessel_Knu ( p () , theta ()) ;  
  return true ;
}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::GenInvGauss::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::GenInvGauss::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 0 ) ; }
// ============================================================================
// mean value
// ============================================================================
double Ostap::Math::GenInvGauss::mean     () const 
{
  const double r = Ostap::Math::bessel_Knu ( p () + 1 , theta () ) * m_iKp_theta ;
  return x ( r ) ;
}
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::GenInvGauss::variance () const 
{
  const double pv = p     () ;
  const double tv = theta () ;
  //
  const double k1 = Ostap::Math::bessel_Knu ( pv + 1 , tv ) ;
  const double k2 = Ostap::Math::bessel_Knu ( pv + 2 , tv ) ;
  //
  const double ss = scale () ;
  //
  return ss * ss * ( k2 * m_iKp_theta - std::pow ( k1 * m_iKp_theta , 2 ) ) ;
}
// ============================================================================
// RMS
// ============================================================================
double Ostap::Math::GenInvGauss::rms () const 
{ return std::sqrt ( variance () ) ; }
// ============================================================================
// mode  
// ============================================================================
double Ostap::Math::GenInvGauss::mode     () const 
{
  const double pv = p     () ;
  const double tv = theta () ;  
  const double r  = ( ( pv - 1 ) + std::hypot ( pv - 1 , tv ) ) / tv ; 
  return x ( r ) ;
}  
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::GenInvGauss::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::GenInvGauss::integral   
( const double low  ,
  const double high ) const
{
  // 
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  const double xmn = xmin () ;
  if ( std::isfinite ( xmn ) )
  {
    if      ( high <= xmn ) { return 0 ; }
    else if ( low  <  xmn ) { return integral ( xmn , high ) ; }
  }
  //
  const double xmx = xmax () ;
  if ( std::isfinite ( xmx ) )
  {
    if      ( low  >= xmx ) { return 0 ; }
    else if ( high >  xmx ) { return  integral ( low , xmx ) ; }
  }
  //
  // (1) split at mean 
  const double vmean = mean () ;
  if ( low < vmean && vmean < high )
  { return integral ( low , vmean ) + integral ( vmean , high ) ; }
  //
  // (2) split at mode
  const double vmode = mode () ;
  if ( low < vmode && vmode < high )
  { return integral ( low , vmode ) + integral ( vmode , high ) ; }
  //
  // (3) more splits
  const double sigma = rms () ;
  // split points 
  static const std::array<double,3> s_splits {{ 2.0 , 5.0 , 10.0 }}  ; 
  for ( const auto split  : s_splits )
  {
    //
    const double xmid1 = vmean + split * sigma ;
    if ( low < xmid1 && xmid1 < high ) 
    { return integral ( low , xmid1 ) + integral ( xmid1 , high ) ; }
    //
    const double xmid2 = vmean - split * sigma ;
    if ( low < xmid2 && xmid2 < high ) 
    { return integral ( low , xmid2 ) + integral ( xmid2 , high ) ; }
  }
  //
  const bool in_tail =
    ( high <= vmean - s_splits.back() * sigma ) ||
    ( low  >= vmean + s_splits.back() * sigma ) ;    
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenInvGauss> s_integrator ;
  static const char s_message[] = "Integral(GenInvGauss)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,                             // low & high edges
      workspace ( m_workspace ) ,                  // workspace
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
std::size_t Ostap::Math::GenInvGauss::tag () const 
{ 
  static const std::string s_name = "GenInvGauss" ;
  return Ostap::Utils::hash_combiner ( s_name                ,
				       ShiftAndScale::tag () ,				       
				       m_theta       .tag () , 
				       m_p           .tag () ) ; 
}
  
// ============================================================================
namespace 
{
  // ==========================================================================
  inline 
  std::vector<double> index_knots 
  ( const unsigned short n ) 
  {
    std::vector<double> result ( n + 1 ) ;
    std::iota ( result.begin() , result.end() , 0 ) ;
    return result ;
  }
  // ==========================================================================
}
// ============================================================================
/*  constructor from n-parameter
 *  @param n  n-parameter (shape) 
 */
// ============================================================================
Ostap::Math::IrwinHall::IrwinHall
( const unsigned short n ) 
  : m_bspline ( index_knots ( n ) , n - 1 ) 
{
  //
  static const std::string s_E {"Shape parameter n must be positive!"} ;
  static const std::string s_R {"Ostap::Math::IrwinHall"} ;
  Ostap::Assert ( 1<= n , s_E , s_R , 360 ) ;
  //
  m_bspline.setPar ( n - 1 , 1 ) ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double 
Ostap::Math::IrwinHall::evaluate ( const double x )  const 
{
  return
    x <= xmin() ? 0.0 :
    x >= xmax() ? 0.0 : m_bspline ( x ) ;
}
// ============================================================================
// get rms 
// ============================================================================
double Ostap::Math::IrwinHall::rms () const { return std::sqrt ( variance () ) ; }
// ============================================================================
// get kurtosis 
// ============================================================================
double Ostap::Math::IrwinHall::kurtosis () const 
{ return 3 - 6.0 / ( 5 * n () )  ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::IrwinHall::integral ()  const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::IrwinHall::integral 
( const double low  ,
  const double high ) const 
{
  if ( high < low ) { return - integral ( high , low ) ; } 
  return m_bspline.integral ( low , high ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::IrwinHall::tag () const 
{ 
  static const std::string s_name = "IrwinHall" ;
  return Ostap::Utils::hash_combiner ( s_name , n () ) ;
}




// ============================================================================
/*  constructor from n-parameter
 *  @param n  n-parameter (shape) 
 */
// ============================================================================
Ostap::Math::Bates::Bates
( const unsigned short n ) 
  : m_ih ( n ) 
{
  //
  static const std::string s_E {"Shape parameter n must be positive!"} ;
  static const std::string s_R {"Ostap::Math::Bates"} ;
  Ostap::Assert ( 1<= n , s_E , s_R , 360 ) ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double 
Ostap::Math::Bates::evaluate ( const double x )  const 
{
  const std::size_t nn = n() ;
  return
    x <= 0 ? 0.0 :
    x >= 1 ? 0.0 : m_ih ( x * nn ) * nn ;
}
// ============================================================================
// get variance 
// ============================================================================
double Ostap::Math::Bates::variance () const { return 1.0 / ( 12.0 * n() ) ; }
// ============================================================================
// get rms 
// ============================================================================
double Ostap::Math::Bates::rms      () const { return std::sqrt ( variance () ) ; }
// ============================================================================
// get kurtosis 
// ============================================================================
double Ostap::Math::Bates::kurtosis () const { return m_ih.kurtosis ()  ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Bates::integral ()  const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Bates::integral 
( const double low  ,
  const double high ) const 
{
  if ( high < low ) { return - integral ( high , low ) ; }
  if ( high <= 0  ) { return 0 ; }
  if ( low  >= 1  ) { return 0 ; }
  const std::size_t nn = n() ;
  return m_ih.integral ( low * nn , high * nn ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Bates::tag () const 
{ 
  static const std::string s_name = "Bates" ;
  return Ostap::Utils::hash_combiner ( s_name , n () ) ;
}
// ============================================================================
/*  constructor from n-parameter
 *  @param n  n-parameter (shape) 
 */
// ============================================================================
Ostap::Math::BatesShape::BatesShape 
( const double         mu    ,  
  const double         sigma , 
  const unsigned short n     ) 
  : m_sq12n ( std::sqrt ( 12.0 * n ) )
  , m_mu    (            mu      )
  , m_sigma ( std::abs ( sigma ) )
  , m_bates (  n ) 
{}
// ============================================================================
bool Ostap::Math::BatesShape::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) {  return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BatesShape::setSigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) {  return false ; }
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double 
Ostap::Math::BatesShape::evaluate ( const double x )  const 
{
  const double y = t ( x ) ;
  return 
    y <= 0 ? 0.0 : 
    y >= 1 ? 0.0 :
    m_bates ( y ) / ( m_sq12n * m_sigma ) ;
}
// ============================================================================
// minimal x 
// ============================================================================
double Ostap::Math::BatesShape::xmin () const { return x ( 0 ) ; }
// ============================================================================
// maximal x 
// ============================================================================
double Ostap::Math::BatesShape::xmax () const { return x ( 1 ) ; }
// ============================================================================
// get kurtosis 
// ============================================================================
double Ostap::Math::BatesShape::kurtosis   () const 
{ return m_bates.kurtosis ()  ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::BatesShape::integral ()  const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::BatesShape::integral 
( const double low  ,
  const double high ) const 
{
  if ( high < low ) { return - integral ( high , low ) ; }
  const double ylow  = t ( low  ) ;
  const double yhigh = t ( high ) ;
  return m_bates.integral ( ylow , yhigh ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BatesShape::tag () const 
{ 
  static const std::string s_name = "BatesShape" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_sigma , n () ) ;
}

// ============================================================================
// Generalized Pareto Distribution
// ============================================================================
Ostap::Math::GenPareto::GenPareto
( const double mu    , 
  const double scale ,
  const double shape ) 
  : m_mu    ( mu ) 
  , m_scale ( std::abs ( scale ) ) 
  , m_shape ( shape  )
{}
// ============================================================================
// evaluate function 
// ============================================================================
double Ostap::Math::GenPareto::evaluate 
( const double x ) const 
{
  if ( x < m_mu ) { return 0 ; }
  //
  const double z =  ( x - m_mu ) / m_scale ;
  if ( s_zero ( m_shape ) ) { return std::exp ( -z ) / m_scale ; }
  //
  if ( m_shape < 0 && x >= m_mu - m_scale / m_shape ) { return 0 ; }
  //
  return std::pow ( 1 + m_shape * z , -1 -1 / m_shape ) / m_scale ;
}
// =============================================================================
// set mu 
// =============================================================================
bool Ostap::Math::GenPareto::setMu ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set scale parameter
// =============================================================================
bool Ostap::Math::GenPareto::setScale ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_scale , avalue ) ) { return false ; }
  m_scale = avalue ;
  return true ;    
}
// =============================================================================
// set shape parameter
// =============================================================================
bool Ostap::Math::GenPareto::setShape ( const double value ) 
{
  if ( s_equal ( m_shape , value ) ) { return false ; }
  m_shape = value ;
  return true ;    
}
// =============================================================================
// xmax: can be infinite 
// =============================================================================
double Ostap::Math::GenPareto::xmax  () const
{ return m_shape < 0 ? m_mu - m_scale/m_shape : s_POSINF ; }
// =============================================================================
// get the integral 
// =============================================================================
double Ostap::Math::GenPareto::integral  () const { return 1 ; }
// =============================================================================
// get cdf 
// =============================================================================
double Ostap::Math::GenPareto::cdf ( const double x ) const
{
  if ( x <= m_mu ) { return 0 ; }
  //
  const double z =  ( x - m_mu ) / m_scale ;
  if ( s_zero ( m_shape ) ) { return 1.0 - std::exp ( z ) ; }
  //
  const double ishape = -1/m_shape ;
  if ( m_shape < 0 && z >= ishape  ) { return 1 ; }
  //
  return 1 - std::pow ( 1 + m_shape * z , ishape ) ;
}
// =============================================================================
// get the integral between low and high
// =============================================================================
double Ostap::Math::GenPareto::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low ,high ) ) { return 0                          ; }
  else if ( high <  low           ) { return - integral ( high , low  ) ; }
  else if ( high <= m_mu          ) { return 0 ; }
  else if ( low  <  m_mu          ) { return   integral ( m_mu , high ) ; }
  //
  const double zl = ( low  - m_mu ) / m_scale ;
  double       zh = ( high - m_mu ) / m_scale ;
  //
  if ( s_zero ( m_shape ) ) 
  { return std::exp ( -zl ) - std::exp ( -zh ) ; }
  //
  const double ishape = -1/m_shape ;
  if ( m_shape < 0 ) 
  {
    if ( zl >= ishape ) { return 0 ; }
    zh = std::min ( zh , ishape ) ;
  }
  //
  return 
    std::pow ( 1 + m_shape * zl , ishape ) - 
    std::pow ( 1 + m_shape * zh , ishape ) ;
}
// ============================================================================
// mean value, defined only for shape<1 
// ============================================================================
double Ostap::Math::GenPareto::mean       () const 
{
  return 
    m_shape < 1 ? m_mu + m_scale / ( 1 - m_shape ) : std::numeric_limits<double>::quiet_NaN() ;
}
// ============================================================================
// median value
// ============================================================================
double Ostap::Math::GenPareto::median () const 
{ return m_mu + m_scale * ( std::pow ( 2 , m_shape ) - 1 ) / m_shape ; }
// ============================================================================
// variance , defined only for shape<1/2 
// ============================================================================
double Ostap::Math::GenPareto::variance () const 
{
  if ( 0.5 <= m_shape ) { return std::numeric_limits<double>::quiet_NaN() ; }
  const double a = m_scale / ( 1 - m_shape ) ;
  return a * a / ( 1 - 2 * m_shape ) ;
}
// ============================================================================
// rms, defined only for shape<1/2 
// ============================================================================
double Ostap::Math::GenPareto::rms () const 
{
  if ( 0.5 <= m_shape ) { return std::numeric_limits<double>::quiet_NaN() ; }
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// skewness, defined only for shape<1/3 
// ============================================================================
double Ostap::Math::GenPareto::skewness() const 
{
  if ( 1 <= 3 * m_shape ) { return std::numeric_limits<double>::quiet_NaN() ; }
  return 2 * ( 1 + m_shape ) * std::sqrt ( 1 - 2 * m_shape ) / ( 1 - 2 * m_shape ) ;
}
// ============================================================================
// (ex)kurtosis, defined only for shape<1/4 
// ============================================================================
double Ostap::Math::GenPareto::kurtosis() const 
{
  if ( 1 <= 4 * m_shape ) { return std::numeric_limits<double>::quiet_NaN() ; }
  return 
    3 * ( 1 - 2 * m_shape ) * ( 2 * m_shape * m_shape + m_shape + 3 ) /
    ( ( 1 - 3 * m_shape ) * ( 1 - 4 * m_shape ) ) ;  
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenPareto::tag () const 
{ 
  static const std::string s_name = "GenPareto" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_scale , m_shape ) ;
}





// ============================================================================
// Generalized Pareto Distribution
// ============================================================================
Ostap::Math::ExGenPareto::ExGenPareto
( const double mu    , 
  const double scale ,
  const double shape ) 
  : m_mu        ( mu ) 
  , m_scale     ( std::abs ( scale ) ) 
  , m_shape     ( shape * 1000 + 1.e+5 ) // gargabe  
  , m_alog_xi   ( 0 ) 
  , m_bias_mean ( 0 ) 
  , m_bias_var  ( 0 ) 
{
  setShape ( shape ) ;  
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::ExGenPareto::evaluate 
( const double x ) const 
{
  //
  const double z = ( x - m_mu ) / m_scale ;
  //
  if ( s_zero ( m_shape ) ) 
  { return GSL_LOG_DBL_MAX < z ? 0.0 : std::exp ( z - std::exp ( z ) ) ; }
  //
  if ( m_shape < 0 && z >= - m_alog_xi ) { return 0 ; }
  //
  return 
    GSL_LOG_DBL_MAX <= z ? 0.0 : 
    std::exp ( z - ( 1 / m_shape + 1 ) * std::log ( 1 + m_shape * std::exp ( z ) ) ) ;  
}
// =============================================================================
// set mu 
// =============================================================================
bool Ostap::Math::ExGenPareto::setMu ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set scale parameter
// =============================================================================
bool Ostap::Math::ExGenPareto::setScale ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_scale , avalue ) ) { return false ; }
  m_scale = avalue ;
  return true ;    
}
// =============================================================================
// set shape parameter
// =============================================================================
bool Ostap::Math::ExGenPareto::setShape ( const double value ) 
{
  if ( s_equal ( m_shape , value ) ) { return false ; }
  m_shape     = value ;
  //
  static const double s_psi1 = Ostap::Math::digamma  ( 1.0 ) ;
  static const double s_psi2 = Ostap::Math::trigamma ( 1.0 ) ;
  //
  if ( s_zero ( value ) )
  {
    m_alog_xi   = 0 ;
    m_bias_mean = -m_alog_xi + s_psi1 ;  
    m_bias_var  =              s_psi2 ;  
  }
  else if ( value < 0 ) 
  { 
    m_alog_xi   = std::log ( std::abs ( value ) ) ;
    m_bias_mean = -m_alog_xi + s_psi1 - Ostap::Math::digamma  ( 1 - 1/m_shape ) ;  
    m_bias_var  =              s_psi2 - Ostap::Math::trigamma ( 1 - 1/m_shape ) ;  
  }
  else if ( value > 0 ) 
  {
    m_alog_xi   = std::log ( value ) ;
    m_bias_mean = -m_alog_xi + s_psi1 - Ostap::Math::digamma  ( 1/m_shape     ) ;  
    m_bias_var  =              s_psi2 + Ostap::Math::trigamma ( 1/m_shape     ) ;  
  }
  //
  return true ;    
}
// =============================================================================
// get the integral 
// =============================================================================
double Ostap::Math::ExGenPareto::integral  () const { return 1 ; }
// =============================================================================
// get the integral between low and high
// =============================================================================
double Ostap::Math::ExGenPareto::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low ,high ) ) { return 0                          ; }
  else if ( high <  low           ) { return - integral ( high , low  ) ; }
  //
  if ( m_shape < 0 && !s_zero ( m_shape ) )
  {
    const double zl = ( low  - m_mu ) / m_scale ;
    if ( zl >= -m_alog_xi ) { return 0 ; }                         // RETURN!!
    const double xmax = m_mu - m_alog_xi * m_scale ;
    if ( xmax < high      ) { return integral ( low , xmax ) ; }   // RETURN
  }
  //
  const double m0    = mean () ;
  const double sigma = rms  () ;
  //
  // split points 
  static const std::array<int,9> s_split = {{ -8 , -5 , -3 , -1 , 0 , 1 , 3 , 5 , 8 }} ; 
  for( const auto p : s_split )
  {
    const double xmid = m0 + p * sigma ;
    if ( low < xmid && xmid < high ) 
    { return integral ( low , xmid ) + integral ( xmid , high ) ; }
  }
  //
  // in tail?
  const bool in_tail = high <= ( m0 - 7 * sigma ) || low  >= ( m0 + 7 * sigma ) ;
  //  
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<ExGenPareto> s_integrator {} ;
  static char s_message[] = "Integral(ExGenPareto)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high  ,              // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,            // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// mode 
// ============================================================================
double Ostap::Math::ExGenPareto::mode () const 
{ return -1.0 < m_shape ? m_mu : m_mu - m_scale * m_alog_xi ; }
// ============================================================================
// mean value
// ============================================================================
double Ostap::Math::ExGenPareto::mean       () const 
{ return m_mu + m_scale * m_bias_mean ; }
// ============================================================================
// variance
// ============================================================================
double Ostap::Math::ExGenPareto::variance () const 
{ return m_scale * m_scale * m_bias_var ; }
// ============================================================================
// rms
// ============================================================================
double Ostap::Math::ExGenPareto::rms () const 
{ return m_scale * std::sqrt ( m_bias_var ) ; }
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::ExGenPareto::tag () const 
{ 
  static const std::string s_name = "ExGenPareto" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_scale , m_shape ) ;
}
// ============================================================================

// ============================================================================
Ostap::Math::Benini::Benini
( const std::size_t          n     ,   // number of shape parameters: n >= 1
  const double               scale ,   // scale parameter 
  const double               shift )  // shift parameter
  : ShiftAndScale ( scale , shift  , "beta" , "mu"  , typeid ( *this )  , false )    
  , m_pars  ( !n ? 1 : n , 0.0 ) 
  , m_pars2 ( !n ? 1 : n , 0.0 )
{
  m_pars  [ 0 ] = 1 ;
  m_pars2 [ 0 ] = 1 ;
}
// ============================================================================

// ============================================================================
// Benini distribution 
// ============================================================================
Ostap::Math::Benini::Benini
( const std::vector<double>& pars  ,   // shape parameters               
  const double               scale ,   // scale parameter 
  const double               shift )   // shift parameter
  : ShiftAndScale ( scale , shift  , "scale" , "shift" , typeid ( *this )  , false )    
  , m_pars  ( pars.size() ) 
  , m_pars2 ( pars.size() ) 
{
  setPars ( pars ) ;
  //
  Ostap::Assert ( !m_pars.empty() && 0 < std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ,
                  "Invalid parameters!"  , 
                  "Ostap::Math::Benini"  ,
		  INVALID_PARAMETERS , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// one shape parameters: alpha
// ============================================================================
Ostap::Math::Benini::Benini
( const double               alpha     , 
  const double               scale     ,      
  const double               shift     )
  : Benini ( std::vector<double>(1,alpha)  , scale , shift )
{}
// ============================================================================
// two shape parameters: alpha and beta 
// ============================================================================
Ostap::Math::Benini::Benini
( const double               alpha     , 
  const double               beta      ,	
  const double               scale     ,      
  const double               shift     )
  : Benini ( {{ alpha , beta }} , scale , shift )
{}
// ======================================================================
// three shape parameters: alpha, beta & gamma 
// ======================================================================
Ostap::Math::Benini::Benini
( const double               alpha     , 
  const double               beta      , 	
  const double               gamma     , 	
  const double               scale     ,  // scale parameter 
  const double               shift     )  // shift parameter
  : Benini ( {{ alpha , beta , gamma }} , scale , shift )
{}
// ======================================================================
// four shape parameters: alpha, beta, gamma & delta 
// ======================================================================
Ostap::Math::Benini::Benini
( const double               alpha     , 
  const double               beta      , 	
  const double               gamma     , 	
  const double               delta     , 	
  const double               scale     ,  // scale parameter 
  const double               shift     )  // shift parameter
  : Benini ( {{ alpha , beta , gamma , delta }} , scale , shift )
{}
// ============================================================================
bool Ostap::Math::Benini::setPar  
( const unsigned short i     , 
  const double         value ) 
{
  if ( m_pars.size() <= i ) { return false ; }
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_pars [ i ] , avalue ) ) { return false ; }
  //
  m_pars  [ i ] =              avalue ;
  m_pars2 [ i ] = ( i + 1 ) *  avalue ;
  //
  Ostap::Assert ( !s_zero ( avalue ) ||
		  0 < std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ,
                  "Invalid parameters/sum of parameters"   , 
                  "Ostap::Math::Benini"                    ,
		  INVALID_PARAMETERS , __FILE__ , __LINE__ ) ;
  //
  return true ; 
}
// ============================================================================
bool Ostap::Math::Benini::setPars  
( const std::vector<double>& values ) 
{
  bool changed = false ;
  const std::size_t n = std::min ( m_pars.size() , values.size() ) ;
  for ( unsigned short i = 0 ; i < n ; ++i ) 
  {
    const double avalue = std::abs ( values [ i ] ) ;
    if ( s_equal ( m_pars [ i ] , avalue ) ) { continue ; }
    //
    m_pars  [ i ] =              avalue ;
    m_pars2 [ i ] = ( i + 1 ) *  avalue ;
    //
    changed = true ;
  }
  //
  Ostap::Assert ( 0 < std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ,
                  "Invalid sum of parameters"   , 
                  "Ostap::Math::Benini"         ,
		  INVALID_PARAMETERS , __FILE__ , __LINE__ ) ;
  //
  return changed ;
}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::Benini::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 1 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::Benini::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 1 ) ; }
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Benini::evaluate 
( const double x ) const 
{
  //
  const double y    = t ( x ) ;
  if ( y <= 1 ) { return 0 ; }
  //
  const double logy = std::log ( y ) ;
  //
  const double p1 = Ostap::Math::Clenshaw::monomial_sum ( m_pars .rbegin() , m_pars .rend  () , logy ).first ;
  const double p2 = Ostap::Math::Clenshaw::monomial_sum ( m_pars2.rbegin() , m_pars2.rend  () , logy).first ;
  //
  return std::exp ( - logy * p1 ) * p2 / ( y * scale_abs () ) ;
}
// =============================================================================
// evaluate the cdf 
// ============================================================================
double Ostap::Math::Benini::cdf
( const double x ) const
{
  const double y    = t ( x )       ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// evaluate the cdf 
// ============================================================================
double Ostap::Math::Benini::raw_cdf
( const double y ) const 
{
  if ( y <= 1 ) { return 0 ; }
  const double logy = std::log ( y ) ;
  // evaluate polynomial 
  const double p1 = Ostap::Math::Clenshaw::monomial_sum ( m_pars.rbegin() , m_pars.rend  () , logy ).first ;
  return 1.0L - std::exp ( - logy * p1 ) ;
}
// =============================================================================
// get the integral 
// =============================================================================
double Ostap::Math::Benini::integral  () const { return 1 ; }
// =============================================================================
// get the integral between low and high
// =============================================================================
double Ostap::Math::Benini::integral 
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low ,high ) ) { return 0 ; }
  else if ( high < low            ) { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Benini::tag () const 
{ 
  static const std::string s_name = "Benini" ;
  return Ostap::Utils::hash_combiner ( s_name  , 
				       ShiftAndScale::tag () , 
                                       Ostap::Utils::hash_range ( m_pars ) ) ;
}
// ============================================================================



// ============================================================================
// Generalized extreme value distribution 
// ============================================================================  
Ostap::Math::GEV::GEV
( const double mu    , 
  const double scale ,
  const double shape ) 
  : m_mu ( mu ) 
  , m_scale ( std::abs ( scale ) )
  , m_shape ( shape ) 
{}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::GEV::evaluate 
( const double x ) const 
{
  //
  const double z = ( x - m_mu ) / m_scale ;
  //
  if ( s_zero ( m_shape ) ) 
  {
    return -z < GSL_LOG_DBL_MAX ? std::exp ( -z - std::exp ( -z ) ) : 0.0 ;  
  }
  else if ( 0 < m_shape && z <= -1 / m_shape ) { return 0.0 ; }
  else if ( 0 > m_shape && z >= -1 / m_shape ) { return 0.0 ; }
  //
  const double tx = std::pow ( 1 + m_shape * z , -1/m_shape ) ;
  //
  return std::pow ( tx , m_shape + 1 ) * std::exp ( -tx ) / m_scale ;
}
// =============================================================================
// set mu 
// =============================================================================
bool Ostap::Math::GEV::setMu ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set scale parameter
// =============================================================================
bool Ostap::Math::GEV::setScale ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_scale , avalue ) ) { return false ; }
  m_scale = avalue ;
  return true ;    
}
// =============================================================================
// set shape parameter
// =============================================================================
bool Ostap::Math::GEV::setShape ( const double value ) 
{
  if ( s_equal ( m_shape , value ) ) { return false ; }
  m_shape     = value ;
  return true ;    
}
// ===========================================================================
// get the CDF 
// ===========================================================================
double Ostap::Math::GEV::cdf ( const double x ) const 
{
  //
  const double z = ( x - m_mu ) / m_scale ;
  //
  if ( s_zero ( m_shape ) ) 
  { return -z < GSL_LOG_DBL_MAX ? std::exp ( - std::exp ( -z ) ) : 0.0 ; }
  else if ( 0 < m_shape && z <= -1 / m_shape ) { return 0.0 ; }
  else if ( 0 > m_shape && z >= -1 / m_shape ) { return 1.0 ; }
  //
  const double tx = std::pow ( 1 + m_shape * z , -1/m_shape ) ;
  //
  return std::exp ( -tx ) ;
}
// =============================================================================
// get the integral 
// =============================================================================
double Ostap::Math::GEV::integral  () const { return 1 ; }
// =============================================================================
// get the integral between low and high
// =============================================================================
double Ostap::Math::GEV::integral 
( const double low  ,
  const double high ) const 
{ return s_equal ( low , high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// mode of distribution 
// ============================================================================
double Ostap::Math::GEV::mode       () const 
{
  return s_zero ( m_shape ) ? m_mu :
    m_mu + m_scale * ( std::pow ( 1 + m_shape , - m_shape ) - 1 ) / m_shape ;
}
// ============================================================================
// median of distribution 
// ============================================================================
double Ostap::Math::GEV::median   () const 
{
  static const double s_ln2   =            std::log ( 2.0 )   ;
  static const double s_lnln2 = std::log ( std::log ( 2.0 ) ) ;
  return s_zero ( m_shape ) ? m_mu - m_scale * s_lnln2 :
    m_mu + m_scale * ( std::pow ( s_ln2 , - m_shape ) - 1 ) / m_shape ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GEV::tag () const 
{ 
  static const std::string s_name = "GEV" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_scale , m_shape ) ;
}
// ============================================================================


// ============================================================================
// constructor from parameters
// ============================================================================
Ostap::Math::FisherZ::FisherZ 
( const double mu     , 
  const double d1     , 
  const double d2     ,
  const double scale  )
  : m_mu    ( mu )
  , m_scale ( std::abs ( scale ) ) 
  , m_d1    ( std::abs ( d1    ) )
  , m_d2    ( std::abs ( d2    ) )
  , m_workspace ()
  , m_C     ( -1 )
{
  Ostap::Assert ( 0 < m_scale , "Invalid scale parameter!" , "Ostap::Math::FisherZ"     ) ;
  Ostap::Assert ( 0 < m_d1    , "Invalid d1-parameter!"    , "Ostap::Math::FisherZ"     ) ;
  Ostap::Assert ( 0 < m_d2    , "Invalid d2-parameter!"    , "Ostap::Math::FisherZ"     ) ;
  //
  m_C = norm() ;
}
// =============================================================================
double Ostap::Math::FisherZ::norm() const
{
  double C  = 0.5 * m_d1 * std::log ( m_d1 ) + 0.5 * m_d2 * std::log ( m_d2 ) ;
  C        -= Ostap::Math::lnbeta ( m_d1 / 2 , m_d2 / 2 ) ;
  return  2 * std::exp ( C ) ;
}
// =============================================================================
// set location parameter
// =============================================================================
bool Ostap::Math::FisherZ::setMu ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu     = value ;
  return true ;    
}
// =============================================================================
// set scale parameter
// =============================================================================
bool Ostap::Math::FisherZ::setScale ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_scale , avalue ) ) { return false ; }
  m_scale = avalue ;
  return true ;    
}
// =============================================================================
// set d1 parameter
// =============================================================================
bool Ostap::Math::FisherZ::setD1 ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_d1 , avalue ) ) { return false ; }
  m_d1 = avalue ;
  m_C  = norm  () ; 
  return true ;    
}
// =============================================================================
// set d1 parameter
// =============================================================================
bool Ostap::Math::FisherZ::setD2 ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_d2 , avalue ) ) { return false ; }
  m_d2 = avalue ;
  m_C  = norm () ;  
  return true ;    
}
// ============================================================================
// get the value of FizherZ distribution
// ============================================================================
double Ostap::Math::FisherZ::evaluate ( const double x ) const
{
  const double z = ( x - m_mu ) / m_scale ;
  //
  const double logres = z <=0 ?
    +m_d1 * z - 0.5 * ( m_d1 + m_d2 ) * std::log ( m_d1 * std::exp (  2 * z ) + m_d2 ) :
    -m_d2 * z - 0.5 * ( m_d1 + m_d2 ) * std::log ( m_d2 * std::exp ( -2 * z ) + m_d1 ) ;
  //
  return  m_C * std::exp ( logres ) / m_scale ;
}
// =================================================================================
/// get the integral 
// =================================================================================
double Ostap::Math::FisherZ::integral () const { return 1 ; }
// =================================================================================
// get the integral 
// =================================================================================
double Ostap::Math::FisherZ::integral
( const double xlow , 
  const double xhigh ) const
{
  if      ( s_equal ( xlow  , xhigh ) ) { return 0 ; }
  else if (           xhigh < xlow    ) { return - integral ( xhigh , xlow ) ; }
  //
  if ( xlow < m_mu && m_mu < xhigh )
    { return integral ( xlow , m_mu  ) + integral ( m_mu , xhigh ) ; }
  //
  const double ss = std::sqrt ( 0.5 * ( 1/m_d1 + 1/m_d2 ) ) ;
  // 
  const double x_l10 = m_mu - 10 * m_scale * ss ;
  const double x_h10 = m_mu + 10 * m_scale * ss ;
  //
  if ( xlow < x_l10 && x_l10 < xhigh )
    { return integral ( xlow , x_l10 ) + integral ( x_l10 , xhigh ) ; }
  if ( xlow < x_h10 && x_h10 < xhigh )
    { return integral ( xlow , x_h10 ) + integral ( x_h10 , xhigh ) ; }
  //
  const double x_l5  = m_mu -  5 * m_scale * ss ;
  const double x_h5  = m_mu +  5 * m_scale * ss ;
  //
  if ( xlow < x_l5 && x_l5  < xhigh )
    { return integral ( xlow , x_l5  ) + integral ( x_l5  , xhigh ) ; }
  if ( xlow < x_h5 && x_h5  < xhigh )
    { return integral ( xlow , x_h5  ) + integral ( x_h5  , xhigh ) ; }
  //
  const double x_l3  = m_mu -  3 * m_scale * ss ;
  const double x_h3  = m_mu +  3 * m_scale * ss ;
  //
  if ( xlow < x_l3 && x_l3  < xhigh )
    { return integral ( xlow , x_l3  ) + integral ( x_l3  , xhigh ) ; }
  if ( xlow < x_h3 && x_h3  < xhigh )
    { return integral ( xlow , x_h3  ) + integral ( x_h3  , xhigh ) ; }
  //
  // we are in tails? 
  const bool in_tail = ( xhigh <= x_l10 ) || ( x_h10 <= xlow ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<FisherZ> s_integrator ;
  static const char s_message[] = "Integral(FisherZ)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xlow    , xhigh ,                                       // low & high edges
      workspace ( m_workspace ) ,                             // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION ,            // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION ,            // relative precision
      m_workspace.size()              ,                       // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::FisherZ::tag () const 
{ 
  static const std::string s_name = "FisherZ" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_scale , m_d1 , m_d2 ) ;
}
// ============================================================================

// ============================================================================
Ostap::Math::BirnbaumSaunders::BirnbaumSaunders
( const double mu    , // location 
  const double beta  , // scale
  const double gamma ) // shape 
  : m_mu ( mu )
  , m_beta  ( std::abs ( beta  ) )
  , m_gamma ( std::abs ( gamma ) )
{}
// =============================================================================
// set scale parameter
// =============================================================================
bool Ostap::Math::BirnbaumSaunders::setMu ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set scale parameter
// =============================================================================
bool Ostap::Math::BirnbaumSaunders::setBeta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_beta , avalue ) ) { return false ; }
  m_beta = avalue ;
  return true ;    
}
// =============================================================================
// set shape parameter
// =============================================================================
bool Ostap::Math::BirnbaumSaunders::setGamma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_gamma , avalue ) ) { return false ; }
  m_gamma = avalue ;
  return true ;    
}
// =============================================================================
// get mean value 
// =============================================================================
double Ostap::Math::BirnbaumSaunders::mean () const
{ return m_mu + m_beta * ( 1 + 0.4 * m_gamma * m_gamma ) ; }
// =============================================================================
// get variance  
// =============================================================================
double Ostap::Math::BirnbaumSaunders::variance () const
{
  const double a2 = m_gamma * m_gamma ;
  return a2 * m_beta * m_beta * ( 1 + 1.25 * a2 ) ;
}
// =============================================================================
// get RMS
// =============================================================================
double Ostap::Math::BirnbaumSaunders::rms  () const
{ return std::sqrt ( variance () ) ;  }
// =============================================================================
// get skewness 
// =============================================================================
double Ostap::Math::BirnbaumSaunders::skewness() const
{
  const double a2 = m_gamma * m_gamma ;
  return 4 * m_gamma * ( 11 * a2 + 6 ) / std::pow ( 5 * a2 + 4 ,   1.5 ) ;
}
// =============================================================================
// get kurtosis 
// =============================================================================
double Ostap::Math::BirnbaumSaunders::kurtosis() const
{
  const double a2 = m_gamma * m_gamma ;
  return 6 * a2 * ( 93 * a2 + 40 ) / std::pow ( 5 * a2 + 4 , 2 ) ;
}
// ============================================================================
// get the value
// ===========================================================================
double Ostap::Math::BirnbaumSaunders::evaluate ( const double x ) const
{
  if ( x <= m_mu ) { return 0 ; }
  //
  const double z2 = ( x - m_mu ) / m_beta ;
  const double z  = std::sqrt ( z2 ) ;
  const double zi = 1/z ;   
  //
  const double qq = Ostap::Math::gauss_pdf ( ( z - zi ) / m_gamma ) ;
  if ( s_zero ( qq )  )  { return 0 ; }
  //
  return qq *  ( z + zi ) / ( 2 * m_gamma * ( x - m_mu ) ) ;
}
// ============================================================================
// get integral 
// ============================================================================
double Ostap::Math::BirnbaumSaunders::integral() const { return  1 ; }
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::BirnbaumSaunders::cdf ( const double x ) const
{
  if ( x <= m_mu ) { return 0 ; }
  const double z2 = ( x - m_mu ) / m_beta ;
  const double z  = std::sqrt ( z2 ) ;
  return Ostap::Math::gauss_cdf (  ( z + 1/z ) / m_gamma ) ;
}
// ============================================================================
/// get the intergral between xmin and xmax 
// ============================================================================
double Ostap::Math::BirnbaumSaunders::integral
( const double xlow , 
  const double xhigh ) const
{
  if ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xhigh < xlow ) { return - integral ( xhigh , xlow ) ; }
  //
  if ( xhigh <= m_mu ) { return 0 ; }
  return cdf ( xhigh ) - cdf ( xlow ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BirnbaumSaunders::tag () const 
{ 
  static const std::string s_name = "BitnbaumSaunders" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_beta , m_gamma ) ;
}
// ============================================================================

// ============================================================================
// Constructor
// ============================================================================
Ostap::Math::Frechet::Frechet
( const double alpha , 
  const double scale ,  
  const double shift )
  : ShiftAndScale ( scale , shift   , "scale" , "shift" , typeid ( *this ) , false ) 
  , m_alpha       ( alpha                     , "alpha" , typeid ( *this ) )
{}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::Frechet::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::Frechet::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 0 ) ; }
// ============================================================================
// evaluate Frechet distribution
// ============================================================================
double Ostap::Math::Frechet::evaluate
( const double x ) const
{
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const double av = alpha     () ;
  const double s  = abs_scale () ;
  //
  const double dd = ( av + 1 ) * std::log ( y ) ;
  //
  if ( 1 < y ) { return av / s * std::exp ( - std::pow ( 1.0 / y , av ) - dd ) ; }
  return av / s * std::exp ( - std::pow ( y , - av ) - dd ) ;
}
// ============================================================================
// raw moments (unscaled & unshifted)
// ============================================================================
double Ostap::Math::Frechet::raw_moment ( const unsigned short r ) const
{
  const double av = alpha () ;
  return Ostap::Math::gamma ( 1.0 - r / av ) ;
}
// ============================================================================
// mean value
// ============================================================================
double Ostap::Math::Frechet::mean () const
{
  const double av = alpha () ;
  if ( av <= 1 || s_equal ( av , 1 ) ) { return s_QUIETNAN ; } 
  const double m1 = raw_moment ( 1 ) ;
  return x ( m1 ) ;
}
// ============================================================================
// median value
// ============================================================================
double Ostap::Math::Frechet::median () const
{
  const double av = alpha () ;
  const double m  = std::pow ( s_ln2 , - 1 / av ) ;
  return x ( m ) ;
}
// ============================================================================
// mode value
// ============================================================================
double Ostap::Math::Frechet::mode  () const
{
  const double av = alpha () ;
  const double m  = std::pow ( av / ( 1 + av ) , 1 / av ) ;
  return x ( m )  ;
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::Frechet::variance  () const
{
  const double av = alpha () ;
  if ( av <= 2 || s_equal ( av , 2 ) ) { return s_QUIETNAN ; } 
  //
  const double ss = scale () ;
  return ss * ss * Ostap::Math::variance
    ( raw_moment ( 2 ) ,
      raw_moment ( 1 ) ) ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::Frechet::rms  () const
{
  const double av = alpha () ;
  if ( av <= 2 || s_equal ( av , 2 ) ) { return s_QUIETNAN ; } 
  //
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// skewnes 
// ============================================================================
double Ostap::Math::Frechet::skewness() const
{
  const double av = alpha () ;
  if ( av <= 3 || s_equal ( av , 3 ) ) { return s_QUIETNAN ; } 
  //
  return scale_sign () * Ostap::Math::skewness
    ( raw_moment ( 3 ) ,
      raw_moment ( 2 ) ,
      raw_moment ( 1 ) ) ;
}
// ============================================================================
// (excess) kurtosis
// ============================================================================
double Ostap::Math::Frechet::kurtosis () const
{
  const double av = alpha () ;
  if ( av <= 4 || s_equal ( av , 4 ) ) { return s_QUIETNAN ; } 
  //
  return Ostap::Math::kurtosis
    ( raw_moment ( 4 ) ,
      raw_moment ( 3 ) ,
      raw_moment ( 2 ) ,
      raw_moment ( 1 ) ) ;
}
// ============================================================================
// evaluate Frechet CDF 
// ============================================================================
double Ostap::Math::Frechet::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return  0 < scale_sign () ? rcdf : 1 - rcdf ;
}
// ============================================================================
// evaluate Frechet CDF 
// ============================================================================
double Ostap::Math::Frechet::raw_cdf 
( const double y ) const
{
  if ( y <= 0 ) { return  0 ; }
  // 
  const double av = alpha () ;
  if ( 1 <= y  ) { return std::exp ( - std::pow ( 1 / y , av ) ) ; }
  return std::exp ( - std::pow ( y , -av ) ) ;  
}
// ============================================================================
// evaluate Frechet integral 
// ============================================================================
double Ostap::Math::Frechet::integral () const { return 1 ; } 
// ============================================================================
// evaluate Frechet integral 
// ============================================================================
double Ostap::Math::Frechet::integral
( const double a ,
  const double b ) const
{
  if ( s_equal ( a , b ) ) { return 0 ; }
  return cdf ( b ) - cdf ( a ) ;
}
// ============================================================================
/*  get quantile 
 *  @param p probability \f$ 0 \le p < 1 \f$
 */
// ============================================================================
double Ostap::Math::Frechet::quantile
( const double p ) const
{
  const double av = alpha () ;
  return
    s_zero  ( p     ) ?  x ( 0 )    :                               
    s_equal ( p , 1 ) ?  s_POSINF   : 
    p <  0            ?  s_QUIETNAN :
    p >  1            ?  s_QUIETNAN :
    x ( std::pow ( -std::log ( p ) , -1 / av ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Frechet::tag () const 
{ 
  static const std::string s_name = "Frechet" ;
  return Ostap::Utils::hash_combiner ( s_name ,
				       ShiftAndScale::tag () , 
				       m_alpha       .tag () ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor from all parameters
 *  @param p shape parameter \f$ 0 < p \f$
 *  @param a shape parameter \f$ 0 < a \f$
 *  @param b scale parameter \f$ 0 < b \f$
 *  @param shift shift parameter        
 */
// ============================================================================
Ostap::Math::Dagum::Dagum
( const double p     ,
  const double a     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift   , "scale" , "shift" , typeid ( *this ) , false ) 
  , m_a           ( a                         , "a"     , typeid ( *this ) )
  , m_p           ( p                         , "p"     , typeid ( *this ) )
{}
// ============================================================================
// evaluate Dagum distribution
// ============================================================================
double Ostap::Math::Dagum::evaluate
( const double x ) const
{
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const double av = a () ;
  const double pv = p () ;
  const double s  = std::abs ( scale () ) ;
  //
  return ( av * pv / s )
    * pow_ratio_a2 ( y , av )
    * std::pow ( pow_ratio_a1 ( y , av ) , pv ) / y ;
}  
// ============================================================================
// evaluate Dagum CDF 
// ============================================================================
double Ostap::Math::Dagum::cdf
( const double x    ) const
{  
  const double y = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ; 
}
// ============================================================================
// evaluate Dagum CDF 
// ============================================================================
double Ostap::Math::Dagum::raw_cdf
( const double y    ) const
{
  if ( y <= 0 ) { return 0 ; }
  const double av = a () ;
  const double pv = p () ;
  return std::pow ( pow_ratio_a2 ( y , -av ) , pv ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::Dagum::integral
( const double low  ,
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return - integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get quantile \f$ 0 < p < 1 \f$
// ============================================================================
double Ostap::Math::Dagum::quantile
( const double u ) const
{
  return
    s_zero  ( u     ) ?  x ( 0 )    : 
    s_equal ( u , 1 ) ?  s_POSHUGE  :    
    u < 0             ?  s_QUIETNAN :
    u > 1             ?  s_QUIETNAN :
    x ( std::pow ( std::pow ( u , -1/ p () ) -1 , -1 / a () ) ) ;
}
// ============================================================================
// mean 
// ============================================================================
double Ostap::Math::Dagum::mean() const
{
  const double av = a () ;  
  if ( av <= 1 || s_equal ( av , 1 ) ) { return s_QUIETNAN ; }
  //
  const double lg1 = Ostap::Math::lgamma ( 1   - 1 / m_a ) ;
  const double lg2 = Ostap::Math::lgamma ( m_p + 1 / m_a ) ;
  const double lgp = Ostap::Math::lgamma ( m_p ) ;
  //
  return x ( std::exp ( lg1 + lg2 - lgp ) ) ;  
}
// ============================================================================
// mode
// ============================================================================
double Ostap::Math::Dagum::mode () const
{
  const double av = a () ;  
  const double pv = p () ;  
  const double ap = av * pv  ;
  if ( ap <= 1 ) { return x ( 0 ) ; }
  //
  return x ( std::pow ( ( ap - 1 ) / ( av + 1 ) , 1 / av  ) ) ;
}
// ============================================================================
// median
// ============================================================================
double Ostap::Math::Dagum::median () const
{
  const double av = a () ;  
  const double pv = p () ;  
  const double q  = 1 / pv ;
  if ( 1 < q )
  {
    const double q2 = std::pow ( 2.0 , -q ) ;
    return x ( std::pow ( q2 / ( 1.0L - q2 ) , 1/ av  ) ) ;
  }
  const double q2 = std::pow ( 2.0 , q ) ;
  return x ( std::pow ( 1 / ( q2 - 1.0L )  , 1 /av ) ) ;
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::Dagum::variance () const
{
  const double av = a () ;  
  if ( av <= 2 || s_equal ( av , 2 ) ) { return s_QUIETNAN ; }
  //
  const double pv = p () ;  
  //
  const double lg1 = Ostap::Math::lgamma ( 1   - 2 / av ) ;
  const double lg2 = Ostap::Math::lgamma ( pv  + 2 / av ) ;
  const double lg3 = Ostap::Math::lgamma ( 1   - 1 / av ) ;
  const double lg4 = Ostap::Math::lgamma ( pv  + 1 / av ) ;
  const double lgp = Ostap::Math::lgamma ( pv  ) ;
  //
  const double t1 = std::exp ( lg1 + lg2 - lgp ) ;
  const double t2 = std::exp ( lg3 + lg4 - lgp ) ;
  //
  const double ss = scale () ;
  return ss * ss * ( t1 - t2 * t2 ) ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::Dagum::rms  () const
{
  const double av = a () ;  
  if ( av <= 2 || s_equal ( av , 2 ) ) { return s_QUIETNAN ; }
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Dagum::tag () const 
{ 
  static const std::string s_name = "Dagum" ;
  return Ostap::Utils::hash_combiner ( s_name                ,
				       ShiftAndScale::tag () , 
				       m_a           .tag () ,
				       m_p           .tag () ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::BenktanderI::BenktanderI
( const double a     ,
  const double r     , 
  const double scale ,
  const double shift )
  : m_a  ( a                       , "a"     , typeid ( *this ) )
  , m_r  ( r                       , "r"     , typeid ( *this ) )
  , m_ss ( scale , shift , "scale" , "shift" , typeid ( *this ) )
  , m_p  ( 1.0 / std::hypot ( r , 1.0 ) )    
{}
// ============================================================================
// Set R parameter  
// ============================================================================
bool Ostap::Math::BenktanderI::setR 
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( !m_r.setValue ( avalue ) ) { return false ; }
  m_p = 1.0 / std::hypot ( r () , 1.0 ) ;
  return true ;
}
// ============================================================================
// original parameter "b" 
// ============================================================================
double Ostap::Math::BenktanderI::b () const
{
  const double av = a () ;
  const double pv = p () ;
  return pv * av * ( av + 1 ) / 2 ;
}
// ============================================================================
// evaluate Benktander Type I distribution
// ============================================================================
double Ostap::Math::BenktanderI::evaluate
( const double x ) const
{
  const double z = t ( x ) ;
  if ( z < 1 ) { return 0 ; }
  //
  const double bp = b () ;
  const double ap = a () ;
  //
  const double lz = std::log ( z ) ;
  //
  const double t1 = 1 +      2 * bp * lz / ap ;
  const double t2 = 1 + ap + 2 * bp * lz      ;
  const double t3 = 2 + ap +     bp * lz      ;
  //
  const double rpdf =  ( t1 * t2 - 2 * bp / ap ) * std::pow ( z , -t3 ) ;
  return rpdf / std::abs ( scale () ) ;
}
// ============================================================================
// evaluate integral for Benktander Type I distribution
// ============================================================================
double Ostap::Math::BenktanderI::integral () const { return 1 ; }
// ============================================================================
// evaluate integral for Benktander Type I distribution
// ============================================================================
double Ostap::Math::BenktanderI::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// evaluate Benktander Type I CDF 
// ============================================================================
double Ostap::Math::BenktanderI::cdf 
( const double x ) const
{
  const double z = t ( x ) ;
  if ( z <= 1 ) { return 0 ; }
  //
  const double lz = std::log ( z ) ;
  //
  const double bp = b () ;
  const double ap = a () ;
  //
  const double t1 = 1 +      2 * bp * lz / ap ;
  const double t2 = 1 + ap +     bp * lz      ;
  //
  return 1 - t1 * std::pow ( x , -t2 ) ;
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::BenktanderI::mean     () const
{
  const double zmean = 1 + 1 / a () ;
  return x ( zmean ) ;
}
// ============================================================================
// variance value 
// ============================================================================
double Ostap::Math::BenktanderI::variance () const
{
  const double ap  = a () ;
  const double bp  = b () ;
  //
  const double sqb = std::sqrt ( bp ) ;
  const double t1  = ( ap - 1 ) / ( 2 * sqb ) ;
  //
  const double t2  = -sqb + ap * s_sqrt_pi * Ostap::Math::erfcx ( t1 ) ;
  //
  const double ss  = scale () ;
  return ss * ss * t2 / ( ap * ap * sqb ) ;
}
// ============================================================================
// RMS value 
// ============================================================================
double Ostap::Math::BenktanderI::rms () const
{ return std::sqrt ( variance () ) ; }
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BenktanderI::tag () const 
{ 
  static const std::string s_name = "BenktanderI" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_a .tag () ,
				       m_r .tag () ,
				       m_ss.tag () ) ;
}
// ============================================================================

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::BenktanderII::BenktanderII
( const double a     ,
  const double r     , 
  const double scale ,
  const double shift )
  : m_a  ( a                       , "a"     , typeid ( *this ) )
  , m_r  ( r                       , "r"     , typeid ( *this ) )
  , m_ss ( scale , shift , "scale" , "shift" , typeid ( *this ) )
  , m_b  ( 1.0 / std::hypot ( r , 1.0 ) )    
{}
// ============================================================================
// Set R parameter  
// ============================================================================
bool Ostap::Math::BenktanderII::setR 
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( !m_r.setValue ( avalue ) ) { return false ; }
  m_b = 1.0 / std::hypot ( r () , 1.0 ) ;
  return true ;
}
// ============================================================================
// evaluate Benktander Type II distribution
// ============================================================================
double Ostap::Math::BenktanderII::evaluate
( const double x ) const
{
  const double z = t ( x ) ;
  if ( z < 1 ) { return 0 ; }
  //
  const double ap = a () ;
  const double bp = b () ;
  //
  const double zb = std::pow ( z , bp ) ;
  //
  const double zpdf =
    std::exp ( ap * ( 1 - zb ) / bp ) *
    zb * ( bp * zb - bp + 1 ) / ( z * z ) ;
  //
  return zpdf / std::abs ( scale () ) ;
}
// ============================================================================
// evaluate integral for Benktander Type II distribution
// ============================================================================
double Ostap::Math::BenktanderII::integral () const { return 1 ; }
// ============================================================================
// evaluate integral for Benktander Type II distribution
// ============================================================================
double Ostap::Math::BenktanderII::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// evaluate Benktander Type II CDF 
// ============================================================================
double Ostap::Math::BenktanderII::cdf 
( const double x ) const
{
  const double z = t ( x ) ;
  if ( z <= 1 ) { return 0 ; }
  //
  const double ap = a () ;
  const double bp = b () ;
  //  
  const double zb   = std::pow ( z , bp  ) ;
  const double zcdf = 1 - zb * std::exp ( ap  * ( 1 - zb ) / bp ) / z ; 
  //
  return zcdf ;
}
// ============================================================================
// mode 
// ============================================================================
double Ostap::Math::BenktanderII::mode () const
{ return x ( 1 ) ;  }
// ============================================================================
// mean
// ============================================================================
double Ostap::Math::BenktanderII::mean () const
{ return x ( 1 + 1 / a() ) ;  }
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::BenktanderII::variance () const
{
  const double av = a  () ;
  const double bv = b  () ;
  const double ab = av / bv ;
  //
  // something wrong here...
  const double r  =  ( -1 + 2 * ab * std::exp ( ab ) * En ( 1.0 - 1.0 / bv , ab ) ) / ( av * av ) ;
  const double ss = scale () ;
  //
  return  ss * ss * r ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::BenktanderII::rms () const
{ return std::sqrt ( variance () ) ; }
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BenktanderII::tag () const 
{ 
  static const std::string s_name = "BenktanderII" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_a .tag () ,
				       m_r .tag () ,
				       m_ss.tag () ) ;
}
// ============================================================================

// ============================================================================
// constructor from two parametersO
// ============================================================================
Ostap::Math::MPERT::MPERT
( const double xmin , 
  const double xmax ) 
  : MPERT ( xmin , xmax , 0.5 * ( xmin + xmax ) ) 
{}
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::MPERT::MPERT
( const double xmin  , 
  const double xmax  ,
  const double xi    , 
  const double gamma ) 
  : m_xmin  ( std::min ( xmin , xmax ) ) 
  , m_xmax  ( std::max ( xmin , xmax ) ) 
  , m_xi    ( xi  ) 
  , m_gamma ( std::abs ( gamma ) )
{
  Ostap::Assert ( m_xmin < m_xmax , 
                  "Invalid setting of xmin/xmax"      ,
                  "Ostap::Math::MPERT" ) ;
  //  
  if      ( m_xi   < m_xmin ) { m_xi = m_xmin ; }
  else if ( m_xmax < m_xi   ) { m_xi = m_xmax ; }
  //
  setXi    ( m_xi  ) ;
  setGamma ( gamma ) ;
}
// ============================================================================
/// set mode 
// ============================================================================
bool Ostap::Math::MPERT::setXi  ( const double value ) 
{
  const double v = 
    value  < m_xmin ? m_xmin :
    m_xmax < value  ? m_xmax : value ;
  if ( s_equal ( v , m_xi ) && 0 < m_alpha1 && 0 < m_alpha2 && 0 < m_N ) 
  { return false ; }
  //
  m_xi     = v ;
  m_alpha1 = 1 + m_gamma * ( m_xi   - m_xmin ) / ( m_xmax - m_xmin ) ; 
  m_alpha2 = 1 + m_gamma * ( m_xmax - m_xi   ) / ( m_xmax - m_xmin ) ; 
  //
  m_N = 1.0 / ( Ostap::Math::beta ( m_alpha1 , m_alpha2 ) *
                std::pow ( m_xmax - m_xmin , m_alpha1 + m_alpha2 - 1 ) ) ;
  //
  return true ;
}
// ============================================================================
/// set gamma/shape
// ============================================================================
bool Ostap::Math::MPERT::setGamma ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_gamma ) && 0 < m_alpha1 && 0 < m_alpha2 && 0 < m_N ) 
  { return false ; }
  //
  m_gamma  = v ;
  m_alpha1 = 1 + m_gamma * ( m_xi   - m_xmin ) / ( m_xmax - m_xmin ) ; 
  m_alpha2 = 1 + m_gamma * ( m_xmax - m_xi   ) / ( m_xmax - m_xmin ) ; 
  //
  m_N = 1.0 / ( Ostap::Math::beta ( m_alpha1 , m_alpha2 ) *
                std::pow ( m_xmax - m_xmin , m_alpha1 + m_alpha2 - 1 ) ) ;
  //
  return true ;
} 
// ============================================================================
// get the value of MPERT distribution 
// ============================================================================
double Ostap::Math::MPERT::evaluate ( const double x ) const 
{
  if ( x < m_xmin || m_xmax < x ) { return 0 ; }
  return m_N * 
           std::pow ( x - m_xmin , m_alpha1 - 1 ) * 
           std::pow ( m_xmax - x , m_alpha2 - 1 ) ;
}
// ============================================================================
// get the variance 
// ============================================================================
double Ostap::Math::MPERT::variance () const 
{
  const double m = mu () ;
  return ( m - m_xmin ) * ( m_xmax - m ) / ( m_gamma + 3 ) ;
}
// ============================================================================
// get the skewness 
// ============================================================================
double Ostap::Math::MPERT::skewness () const 
{
  const double m = mu () ;
  const double a = 0.25 * ( m_xmin + m_xmax + 2 * m ) ;
  const double b = ( m - m_xmin ) * ( m_xmax - m    ) / 7 ;
  return a / std::sqrt ( b ) ;
}
// ============================================================================
// get the (excess) kurtosis 
// ============================================================================
double Ostap::Math::MPERT::kurtosis() const 
{
  const double a1 = m_alpha1 ;
  const double a2 = m_alpha2 ;
  //
  return 3 * ( a1 + a2 + 1 ) 
    * ( 2 * ( a1 + a2 ) * ( a1 + a2 ) + a1 * a2 * ( a1 + a2 - 6 ) )
    / ( a1 * a2 * ( a1 + a2 + 2 ) * ( a1 + a2 + 3 ) ) - 3 ;
}
// ============================================================================
// get the itegral 
// ============================================================================
double Ostap::Math::MPERT::integral () const { return 1 ; }
// ============================================================================
// get the CDF value 
// ============================================================================
double Ostap::Math::MPERT::cdf ( const double x ) const 
{
  //
  if      ( x       <= m_xmin ) { return 0 ; }
  else if ( m_xmax  <= x      ) { return 1 ; }
  //
  const double z = ( x - m_xmin ) / ( m_xmax - m_xmin ) ;
  return Ostap::Math::beta_inc ( m_alpha1 , m_alpha2 , z ) ;
}  
// ============================================================================
// get the integral between xlow and xhigh 
// ============================================================================
double Ostap::Math::MPERT::integral
( const double xlow  , 
  const double xhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh )  ) { return 0 ; }
  else if ( xhigh  < xlow             ) { return - integral ( xhigh , xlow ) ; }
  else if ( xhigh  <= m_xmin          ) { return 0 ; }
  else if ( m_xmax <= xlow            ) { return 0 ; }
  //
  const double xmn = std::max ( m_xmin , xlow  ) ;
  const double xmx = std::min ( m_xmax , xhigh ) ;
  //
  if ( s_equal ( xmn , xmx ) ) { return 0 ; }
  //
  return cdf ( xmx ) - cdf ( xmn ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::MPERT::tag () const 
{ 
  static const std::string s_name = "MPERT" ;
  return Ostap::Utils::hash_combiner ( s_name , m_xmin , m_xmax , m_xi , m_gamma ) ;
}
// ============================================================================


// ============================================================================
// constructor 
// ============================================================================
Ostap::Math::LogNormal::LogNormal
( const double shape ,
  const double scale ,
  const double shift )
  : m_shape ( std::abs ( shape ) )
  , m_scale ( std::abs ( scale ) )
  , m_shift (            shift   )
{
  Ostap::Assert ( m_shape , 
                  "Shape parameter must be positive!"     ,
                  "Ostap::Math::LogNormal"                ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //  
  Ostap::Assert ( m_scale , 
                  "Scale parameter must be positive!"     ,
                  "Ostap::Math::LogNormal"                ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //  
  m_mu = std::log ( m_scale ) ; 
}
// ============================================================================
// set shape 
// ============================================================================
bool Ostap::Math::LogNormal::setShape ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_shape ) ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "Shape parameter must be positive!"     ,
                  "Ostap::Math::LogNormal"                ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_shape = avalue ; 
  return true ;
} 
// ============================================================================
// set scale 
// ============================================================================
bool Ostap::Math::LogNormal::setScale ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "Scale parameter must be positive!"     ,
                  "Ostap::Math::LogNormal"                ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_scale = avalue ; 
  m_mu    = std::log ( m_scale ) ; 
  return true ;
} 
// ============================================================================
// set Mu 
// ============================================================================
bool Ostap::Math::LogNormal::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  //
  m_mu    = value ;
  m_scale = std::exp ( m_mu ) ; 

  return true ;
} 
// ============================================================================
// set shift 
// ============================================================================
bool Ostap::Math::LogNormal::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift  ) ) { return false ; }
  //
  m_shift = value ; 
  return true ;
} 
// ============================================================================
// evaluate log-normal function
// ============================================================================
double Ostap::Math::LogNormal::evaluate
( const double x ) const
{
  if ( x  <= m_shift || s_equal ( x , m_shift ) ) { return 0 ; }
  const double dx = x - m_shift ; 
  if ( dx <= 0       || s_zero  ( dx )          ) { return 0 ; }
  const double lz = std::log ( dx / m_scale ) / m_shape  ;   
  return Ostap::Math::gauss_pdf ( lz ) / ( dx * m_shape ) ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::LogNormal::integral () const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::LogNormal::integral
( const double low  ,
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low ) ; } 
  else if ( high <= m_shift        ) { return  0 ; }
  else if ( low  <  m_shift        ) { return  integral ( m_shift , high ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF  
// ============================================================================
double Ostap::Math::LogNormal::cdf 
( const double x    ) const
{ 
  if ( x  <= m_shift || s_equal ( x , m_shift ) ) { return 0 ; }
  const double dx = x - m_shift ;
  if ( dx <= 0       || s_zero  ( dx )          ) { return 0 ; }
  //
  const double z = std::log ( dx / m_scale      ) / m_shape ;
  return Ostap::Math::gauss_cdf ( z ) ;
}
// ============================================================================
// quantile  function \f$ 0 < p < 1 \f$ 
// ============================================================================
double Ostap::Math::LogNormal::quantile 
( const double p    ) const
{
  return
    p < 0             ?  s_QUIETNAN :
    p > 1             ?  s_QUIETNAN : 
    s_zero  ( p     ) ?  m_shift    : 
    s_equal ( p , 1 ) ?  s_POSHUGE  :
    m_shift + m_scale * std::exp ( m_shape * Ostap::Math::probit ( p )  ) ;
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::LogNormal::mean     () const
{ return m_shift + m_scale * std::exp ( +0.5 * m_shape * m_shape ) ; }
// ============================================================================
// median 
// ============================================================================
double Ostap::Math::LogNormal::median   () const
{ return m_shift + m_scale ; }
// ============================================================================
// mode value 
// ============================================================================
double Ostap::Math::LogNormal::mode      () const
{ return m_shift + m_scale * std::exp ( -1.0 * m_shape * m_shape ) ; }
// ============================================================================
// variance  value 
// ============================================================================
double Ostap::Math::LogNormal::variance () const
{
  const double es2 = std::exp ( m_shape * m_shape ) ;
  return ( es2 - 1 ) * m_scale * m_scale * es2 ;
}
// ============================================================================
// variance  value 
// ============================================================================
double Ostap::Math::LogNormal::rms  () const { return std::sqrt ( variance () ) ; }
// ============================================================================
// skewness value 
// ============================================================================
double Ostap::Math::LogNormal::skewness () const
{
  const double es2 = std::exp ( m_shape * m_shape ) ;
  return ( es2 + 1 ) * std::sqrt ( es2 - 1 ) ;
}
// ============================================================================
// kurtosis value 
// ============================================================================
double Ostap::Math::LogNormal::kurtosis () const
{
  const double es2 = std::exp ( m_shape * m_shape ) ;  
  return
    std::pow ( es2 , 4 )     +
    std::pow ( es2 , 3 ) * 2 +
    std::pow ( es2 , 2 ) * 3 - 6 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::LogNormal::tag () const 
{ 
  static const std::string s_name = "LogNormal" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_shape ,
				       m_scale ,
				       m_shift ) ;
}
// ============================================================================


// ============================================================================
// ============================================================================
Ostap::Math::ExpoLog::ExpoLog
( const double beta  ,   // scale 
  const double psi   ,   // related to p 
  const double shift )   // shift
  : m_beta  (  std::abs ( beta ) )
  , m_psi   ( psi   )
  , m_shift ( shift )
{
  Ostap::Assert ( m_beta , 
                  "Beta parameter must be positive!"      ,
                  "Ostap::Math::ExpoLog"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_p = 0.5 * ( 1.0L + std::tanh ( 1.0L * m_psi ) ) ;
  if ( s_equal ( m_p , 1 ) ) { m_p    = 1 ; m_logp = 0   ; } 
  else                       { m_logp = std::log ( m_p ) ; }
  //
  Ostap::Assert ( 0 <= m_p && m_p <= 1 , 
                  "P-parameter must be netwen 0 and 1"    ,
                  "Ostap::Math::ExpoLog"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( m_logp <= 0 , 
                  "LogP parameter must be non-positive!"  ,
                  "Ostap::Math::ExpoLog"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  
}
// ============================================================================
// set beta 
// ============================================================================
bool Ostap::Math::ExpoLog::setBeta( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_beta ) ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "Beta parameter must be positive!"      ,
                  "Ostap::Math::ExpoLog"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_beta = avalue ; 
  return true ;
} 
// ============================================================================
// set psi
// ============================================================================
bool Ostap::Math::ExpoLog::setPsi ( const double value ) 
{
  if ( s_equal ( value , m_psi ) && 0 <= m_p && m_p <= 1 && m_logp <= 0 ) { return false ; }
  //
  m_psi  = value ;
  m_p    = 0.5 * ( 1.0L + std::tanh ( 1.0L * m_psi ) ) ;
  if ( s_equal ( m_p , 1 ) ) { m_p = 1 ; m_logp = 0      ; } 
  else                       { m_logp = std::log ( m_p ) ; }
  //
  Ostap::Assert ( 0 <= m_p && m_p <= 1 , 
                  "P-parameter must be netwen 0 and 1"    ,
                  "Ostap::Math::ExpoLog"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( m_logp <= 0 , 
                  "LogP parameter must be non-positive!"  ,
                  "Ostap::Math::ExpoLog"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  return true ;
} 
// ============================================================================
// set shift 
// ============================================================================
bool Ostap::Math::ExpoLog::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift  ) ) { return false ; }
  //
  m_shift = value ; 
  return true ;
} 
// ============================================================================
// evaluate expo-log function
// ============================================================================
double Ostap::Math::ExpoLog::evaluate    ( const double x ) const
{
  if ( x < m_shift ) { return 0 ; }
  const double z    = - m_beta * ( x - m_shift ) ;
  const double expz = std::exp ( z ) ;
  if ( s_equal ( m_p , 1 ) ) { return m_beta * expz ;  }
  //
  const long double dp = 1.0L - m_p ;
  return -1 / m_logp * m_beta * dp * expz / ( 1.0L - dp * expz ) ; 
}
// ============================================================================
// intergal 
// ============================================================================
double Ostap::Math::ExpoLog::integral () const { return 1 ; } 
// ============================================================================
// intergal 
// ============================================================================
double Ostap::Math::ExpoLog::integral
( const double low  , 
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high  < low            ) { return - integral ( high , low ) ; }
  else if ( high  <= m_shift       ) { return 0 ; }
  else if ( low   <  m_shift       ) { return cdf ( high ) ; }
  //
  return cdf ( high ) - cdf ( low ) ; 
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::ExpoLog::cdf
( const double x ) const
{
  if ( x <= m_shift ) { return 0 ; }
  const double z    = - m_beta * ( x - m_shift ) ;
  const double expz = std::exp ( z ) ;
  //
  if ( s_equal ( m_p , 1 ) ) { return 1.0 - expz ; }
  //
  return 1.0L - std::log ( 1.0 - ( 1.0 - m_p ) * expz ) / m_logp ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::ExpoLog::tag () const 
{ 
  static const std::string s_name = "ExpoLog" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_beta  ,
				       m_psi   ,
				       m_shift ) ;
}
// ============================================================================


// ============================================================================
// constructor 
// ============================================================================
Ostap::Math::Davis::Davis
( const double b  ,
  const double n  ,
  const double mu )
  : m_b  ( std::abs ( b ) )
  , m_n  ( std::abs ( n ) )
  , m_mu ( mu )
  , m_C  ( -1 )
{
  Ostap::Assert ( 0 < m_b  , 
                  "b-parameter must be positive!"         ,
                  "Ostap::Math::Davis"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 0 < m_n  , 
                  "b-parameter must be positive!"         ,
                  "Ostap::Math::Davis"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_z0 = s_equal ( m_n     , 1 ) ? s_Mascheroni : Ostap::Math::zeta ( m_n     ) * 1.0L ;
  m_z1 = s_equal ( m_n - 1 , 1 ) ? s_Mascheroni : Ostap::Math::zeta ( m_n - 1 ) * 1.0L ;
  m_z2 = s_equal ( m_n - 2 , 1 ) ? s_Mascheroni : Ostap::Math::zeta ( m_n - 2 ) * 1.0L ;
  //
  m_C  = Ostap::Math::igamma ( m_n     ) / m_z0  ;
  //
}
// ============================================================================
// set beta 
// ============================================================================
bool Ostap::Math::Davis::setB ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_b ) ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "b-parameter must be positive!"         ,
                  "Ostap::Math::Davis"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_b = avalue ; 
  return true ;
} 
// ============================================================================
// set n
// ============================================================================
bool Ostap::Math::Davis::setN ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_n ) && 0 < m_C ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "n-parameter must be positive!"         ,
                  "Ostap::Math::Davis"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_n  = avalue ;
  //
  m_z0 = s_equal ( m_n     , 1 ) ? s_Mascheroni : Ostap::Math::zeta ( m_n     ) * 1.0L ;
  m_z1 = s_equal ( m_n - 1 , 1 ) ? s_Mascheroni : Ostap::Math::zeta ( m_n - 1 ) * 1.0L ;
  m_z2 = s_equal ( m_n - 2 , 1 ) ? s_Mascheroni : Ostap::Math::zeta ( m_n - 2 ) * 1.0L ;
  //
  m_C  = Ostap::Math::igamma ( m_n ) / m_z0 ;  
  //
  return true ;
} 
// ============================================================================
// set mu
// ============================================================================
bool Ostap::Math::Davis::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ; 
  return true ;
} 
// ============================================================================
// evaluate davis function
// ============================================================================
double Ostap::Math::Davis::evaluate  ( const double x ) const
{
  if ( x <= m_mu || s_equal ( x , m_mu ) ) { return 0 ; } 
  const double z = m_b / ( x - m_mu ) ;
  //
  if ( z <= 1 ) { return m_C / m_b * std::pow ( z , m_n ) / Ostap::Math::expm1_x ( z ) ; }
  //
  const double e2 = std::exp ( -z )  ;
  const double rr = ( m_n + 1 ) * std::log ( z ) + std::log ( e2 / ( 1.0L - e2 ) ) ;  
  return m_C / m_b * std::exp ( rr ) ;
}
// ============================================================================
// intergal 
// ============================================================================
double Ostap::Math::Davis::integral ()    const { return 1 ; }
// ============================================================================
// intergal 
// ============================================================================
double Ostap::Math::Davis::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low     ) { return -integral ( high , low )  ; }
  else if ( high <= m_mu   ) { return  0 ; }
  else if ( low  <  m_mu   ) { return  integral ( m_mu , high ) ; }
  //
  if ( 2 < m_n )
  {
    const double xmean = mean () ;
    if ( low < xmean && xmean < high ) { return integral ( low , xmean ) + integral ( xmean , high ) ; }
  }
  //
  if ( 3 < m_n )
  {
    const double xmean = mean () ;
    const double xrms  = rms  () ;
    const double x1    = xmean - 3 * xrms ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2    = xmean + 3 * xrms ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    const double x3    = xmean - 6 * xrms ;
    if ( low < x3 && x3 < high ) { return integral ( low , x3 ) + integral ( x3 , high ) ; }
    const double x4    = xmean + 6 * xrms ;
    if ( low < x4 && x4 < high ) { return integral ( low , x4 ) + integral ( x4 , high ) ; }
  }
  //  
  const double dx = high - low ;
  // split 
  if ( dx > 10 * m_b )
  {
    const double mid = 0.5 * ( low + high ) ;
    return integral ( low , mid ) + integral ( mid , high ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Davis> s_integrator{}  ;
  static const char s_message[] = "Integral(Davis)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,            // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size()        , // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ; 
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::Davis::cdf 
( const double x ) const
{ return x <= m_mu ? 0.0 : integral ( m_mu , x ) ; }
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::Davis::mean     () const
{
  if ( m_n <= 2 || s_equal ( m_n , 2 ) ) { return s_QUIETNAN ; }
  return m_mu + m_b * m_z1 / ( ( m_n - 1 ) * m_z0 ) ;
}
// ============================================================================
// variance   
// ============================================================================
double Ostap::Math::Davis::variance () const
{
  if ( m_n <= 3 || s_equal ( m_n , 3 ) ) { return s_QUIETNAN ; }
  return
    m_b * m_b * ( - ( m_n - 2 ) * m_z1 * m_z1  + ( m_n - 1 ) * m_z2 * m_z0 )
    / ( ( m_n - 2 ) * ( m_n - 1 ) * ( m_n -1 ) * m_z0 * m_z0 ) ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::Davis::rms () const
{
  if ( m_n <= 3 || s_equal ( m_n , 3 ) ) { return s_QUIETNAN ; }
  //
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Davis::tag () const 
{ 
  static const std::string s_name = "Davis" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_b     ,
				       m_n     ,
				       m_mu    ) ;
}
// ============================================================================

// ======================================================================
Ostap::Math::Kumaraswami::Kumaraswami
( const double alpha , // 0<alpha 
  const double beta  , // 0<beta 
  const double scale , // 0<scale
  const double shift )
  : ShiftAndScale ( scale , shift   , "scale"  , "shift" , typeid ( *this ) , false )
  , m_alpha       ( alpha                      , "alpha" , typeid ( *this ) )
  , m_beta        ( beta                       , "beta"  , typeid ( *this ) )
{}
// ======================================================================
// evaluate Kumaraswami function
// ======================================================================
double Ostap::Math::Kumaraswami::evaluate    ( const double x ) const
{
  const double y = t ( x ) ;
  if ( y < 0 || 1 < y ) { return 0 ; }
  //
  const double av = alpha () ;
  const double bv = beta  () ;
  const double s  = std::abs ( scale () ) ;
  //
  if ( 1 <= av && s_zero ( y ) ) { return 0 ; }
  //
  const long double ya = std::pow ( y * 1.0L , 0.0L + av ) ;
  //
  return std::abs ( av * bv / s ) * ya * std::pow ( 1.0L - ya , bv - 1.0L ) / y ;
}
// ======================================================================
// intergal 
// ======================================================================
double Ostap::Math::Kumaraswami::integral ()    const { return 1 ; } 
// ======================================================================
// intergal 
// ======================================================================
double Ostap::Math::Kumaraswami::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low  , high ) ) { return  0 ; }
  else if (           high < low    ) { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ======================================================================
// CDF 
// ======================================================================
double Ostap::Math::Kumaraswami::cdf 
( const double x ) const
{
  const double y  = t ( x ) ; 
  const bool   ps = Ostap::Math::is_positive ( scale () ) ;
  return ps ? raw_cdf ( y ) : 1 - raw_cdf ( y ) ;
}
// ======================================================================
// "raw" cdf (for positive scale)
// ======================================================================
double Ostap::Math::Kumaraswami::raw_cdf
( const double y ) const
{
  return
    y <= 0 ? 0.0 : 
    y >= 1 ? 1.0 :
    1 - std::pow ( 1 - std::pow ( y , alpha ()  ) , beta () ) ;
}
// ======================================================================
// quantile  function \f$ 0 < p < 1 \f$ 
// ======================================================================
double Ostap::Math::Kumaraswami::quantile 
( const double p    ) const
{
  if ( p < 0 || 1 < p ) { return s_QUIETNAN ; }
  //
  const double av = alpha () ;
  const double bv = beta  () ;
  //
  const double qv =
    s_zero  ( p     ) ? 0.0 :
    s_equal ( p , 1 ) ? 1.0 :
    std::pow ( 1.0L - std::pow ( 1.0L - p , 1.0L / beta () ) , 1.0L / alpha () ) ;
  //
  return x ( qv ) ;
}
// ============================================================================
// xmin
// ============================================================================
double Ostap::Math::Kumaraswami::xmin  () const
{
  const bool ps = Ostap::Math::is_positive ( scale () ) ;
  return ps ? x ( 0 ) : x ( 1 ) ;
}
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::Kumaraswami::xmax () const
{
  const bool ps = Ostap::Math::is_positive ( scale () ) ;
  return ps ? x ( 1 ) : x ( 0 ) ;
}
// ============================================================================
// raw moment for standartized Kumaraswami
// ============================================================================
double Ostap::Math::Kumaraswami::raw_moment
( const unsigned short n ) const
{
  const double av = alpha () ;
  const double bv = beta  () ;
  return bv * Ostap::Math::beta ( 1 + n / av , bv ) ;
}
// ============================================================================
// mean-value 
// ============================================================================
double Ostap::Math::Kumaraswami::mean     () const
{ return x ( raw_moment ( 1 ) ) ; }
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::Kumaraswami::variance  () const
{
  const double ss = scale () ;
  return ss * ss * Ostap::Math::variance ( raw_moment ( 2 ) ,
					   raw_moment ( 1 ) ) ;
}
// ============================================================================
// rms-value 
// ============================================================================
double Ostap::Math::Kumaraswami::rms     () const
{ return std::sqrt ( variance () ) ; }
// ============================================================================
// median-value 
// ============================================================================
double Ostap::Math::Kumaraswami::median() const
{
  const double av = alpha () ;
  const double bv = beta  () ;
  //
  const double m = std::pow ( 1.0L - std::pow ( 2.0L , -1.0L / bv  ) , 1.0L / av ) ; 
  return x ( m ) ;
}
// ============================================================================
// mode-value 
// ============================================================================
double Ostap::Math::Kumaraswami::mode () const
{
  const double av = alpha () ;
  const double bv = beta  () ;
  //
  if      ( av < 1 || bv < 1 ) { return s_QUIETNAN ; }
  else if ( s_equal ( av , 1 ) && s_equal ( bv , 1 ) ) { return x ( 0.5 ) ; }
  //  
  const double m = std::pow ( ( av - 1.0L ) / ( av * bv - 1.0L ) , 1.0L/ av ) ;
  return x ( m ) ; 
}
// ============================================================================
// skewness-value 
// ============================================================================
double Ostap::Math::Kumaraswami::skewness () const
{
  return scale_sign () * Ostap::Math::skewness ( raw_moment ( 3 ) ,
						 raw_moment ( 2 ) ,
						 raw_moment ( 1 ) ) ;
}
// ============================================================================
// kurtosis-value 
// ============================================================================
double Ostap::Math::Kumaraswami::kurtosis() const
{ return Ostap::Math::kurtosis ( raw_moment ( 4 ) ,
				 raw_moment ( 3 ) ,
				 raw_moment ( 2 ) ,
				 raw_moment ( 1 ) ) ; } 
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Kumaraswami::tag () const 
{ 
  static const std::string s_name = "Kumaraswami" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_alpha       .tag () ,
				       m_beta        .tag () ) ; 
}
// ============================================================================

// ============================================================================
// full constructot 
// ============================================================================
Ostap::Math::InverseGamma::InverseGamma
( const double alpha ,
  const double beta  ,
  const double shift )
  : m_alpha ( std::abs ( alpha ) )
  , m_beta  ( std::abs ( beta  ) )
  , m_shift (            shift   )
{
  //
  Ostap::Assert ( m_alpha  , 
                  "Alpha parameter must be positive!"     ,
                  "Ostap::Math::InverseGamma"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( m_beta  , 
                  "Beta parameter must be positive!"      ,
                  "Ostap::Math::InverseGamma"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
// set alpha 
// ============================================================================
bool Ostap::Math::InverseGamma::setAlpha ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_alpha ) ) { return false ; }
  //
  Ostap::Assert ( avalue   , 
                  "Alpha parameter must be positive!"     ,
                  "Ostap::Math::InverseGamma"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_alpha = avalue ; 
  return true ;
} 
// ======================================================================
// set beta 
// ============================================================================
bool Ostap::Math::InverseGamma::setBeta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_beta ) ) { return false ; }
  //
  Ostap::Assert ( avalue   , 
                  "Beta parameter must be positive!"      ,
                  "Ostap::Math::InverseGamma"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_beta = avalue ; 
  return true ;
} 
// ======================================================================
// set shift
// ============================================================================
bool Ostap::Math::InverseGamma::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift ) ) { return false ; }
  m_shift = value ; 
  return true ;
}
// ============================================================================
// evaluate Inverse-Gamma function
// ============================================================================
double Ostap::Math::InverseGamma::evaluate    ( const double x ) const
{
  if ( x <= m_shift ) { return 0 ; }
  const double z = ( x - m_shift ) / m_beta ;
  if ( s_zero ( z ) ) { return 0 ; }
  //
  const double lnr =  -1.0 / z - ( m_alpha + 1 ) * std::log ( z ) ;      
  return std::exp ( lnr ) / m_beta ; 
}
// ============================================================================
// mean
// ============================================================================
double Ostap::Math::InverseGamma::mean     () const
{
  if ( m_alpha <= 1 || s_equal ( m_alpha , 1 ) ) { return s_QUIETNAN ; }
  return m_shift + m_beta / ( m_alpha - 1 ) ;
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::InverseGamma::variance () const
{
  if ( m_alpha <= 2 || s_equal ( m_alpha , 2 ) ) { return s_QUIETNAN ; }
  return m_beta * m_beta / ( ( m_alpha - 1 ) * ( m_alpha - 1 ) * ( m_alpha - 2 ) ) ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::InverseGamma::rms  () const
{
  if ( m_alpha <= 2 || s_equal ( m_alpha , 2 ) ) { return s_QUIETNAN ; }
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::InverseGamma::skewness () const
{
  if ( m_alpha <= 3 || s_equal ( m_alpha , 3 ) ) { return s_QUIETNAN ; }
  return 4 * std::sqrt ( m_alpha - 2 ) / ( m_alpha - 3 ) ;
}
// ============================================================================
// kurtosis
// ============================================================================
double Ostap::Math::InverseGamma::kurtosis () const
{
  if ( m_alpha <= 4 || s_equal ( m_alpha , 4 ) ) { return s_QUIETNAN ; }
  return 6 * ( 5 * m_alpha - 11 ) / ( ( m_alpha - 3 ) * ( m_alpha - 4 ) ) ;
}
// ============================================================================
// intergal 
// ============================================================================
double Ostap::Math::InverseGamma::integral ()    const { return 1 ; } 
// ============================================================================
// intergal 
// ============================================================================
double Ostap::Math::InverseGamma::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low  , high ) ) { return  0 ; }
  else if (           high < low    ) { return -integral ( high    , low  ) ; } 
  else if ( high    <= xmin ()      ) { return  0 ; }
  else if ( low     <  xmin ()      ) { return  integral ( xmin () , high ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;  
}
// ============================================================================
// cdf 
// ============================================================================
double Ostap::Math::InverseGamma::cdf 
( const double x ) const
{
  if ( x <= xmin ()           ) { return 0  ; }
  const double z = ( x - m_shift ) / m_beta ;
  if ( z <= 0 || s_zero ( z ) ) { return 0  ; } 
  //
  return Ostap::Math::gamma_inc_Q ( m_alpha , 1 / z ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::InverseGamma::tag () const 
{ 
  static const std::string s_name = "InverseGamma" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_alpha ,
				       m_beta  ,
				       m_shift ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrI::BurrI
( const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false ) 
{}
// ============================================================================
// evaluate BurrXII function
// ============================================================================
double Ostap::Math::BurrI::evaluate    ( const double x ) const
{
  const double y = t ( x ) ;
  if ( y < 0 || 1 < y ) { return 0 ; }
  return 1 / abs_scale() ;  
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrI::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrI::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrI::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double rcdf = y <= 0 ? 0.0 :  y >= 1 ? 1.0 : y ; 
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// raw-moments of unscaled/unshifted distributions 
// ============================================================================
double Ostap::Math::BurrI::raw_moment
( const unsigned short r ) const { return 1 / ( r + 1 ) ; }
// ============================================================================
// mean              
// ============================================================================
double Ostap::Math::BurrI::mean   () const { return x ( 0.5 ) ; }
// ============================================================================
// mode
// ============================================================================
double Ostap::Math::BurrI::mode   () const { return x ( 0.5 ) ; }
// ============================================================================
// median
// ============================================================================
double Ostap::Math::BurrI::median () const { return x ( 0.5 ) ; }
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::BurrI::variance () const
{
  const double ss = scale () ;
  return ss * ss / 12 ;
}
// ============================================================================
// rms              
// ============================================================================
double Ostap::Math::BurrI::rms () const { return scale () * s_1_sqrt12 ; }
// ============================================================================
// skewness            
// ============================================================================
double Ostap::Math::BurrI::skewness () const { return 0 ; }
// ============================================================================
// kurtosis         
// ============================================================================
double Ostap::Math::BurrI::kurtosis () const
{
  //
  const double m4 = raw_moment ( 4 ) ;
  const double m3 = raw_moment ( 3 ) ;
  const double m2 = raw_moment ( 2 ) ;
  const double m1 = raw_moment ( 1 ) ;
  //
  return Ostap::Math::kurtosis ( m4 , m3 , m2 , m1 ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrI::tag () const 
{ 
  static const std::string s_name = "BurrI" ;
  return Ostap::Utils::hash_combiner ( s_name  , ShiftAndScale::tag () ) ;
}
// ============================================================================


// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrII::BurrII
( const double r     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
{}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::BurrII::xmin () const { return s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::BurrII::xmax () const { return s_POSINF ; }
// ============================================================================
// evaluate BurrII function
// ============================================================================
double Ostap::Math::BurrII::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  //
  const double rv = r () ;
  //
  const double r1 = std::pow ( Ostap::Math::logistic (  y ) , rv ) ;
  const double r2 =            Ostap::Math::logistic ( -y )        ;
  //
  return rv * r1 * r2 / abs_scale() ;  
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrII::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrII::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrII::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  // raw-cdf 
  const double rcdf = std::pow ( Ostap::Math::logistic ( y ) , r () ) ;   
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrII::tag () const 
{ 
  static const std::string s_name = "BurrII" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrIII::BurrIII
( const double r     ,
  const double k     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
  , m_k           ( k                       , "k"     , typeid ( *this ) )
{}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::BurrIII::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::BurrIII::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 0 ) ; }
// ============================================================================
// evaluate Burr-III function
// ============================================================================
double Ostap::Math::BurrIII::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  if ( y <= 0 ) { return 0 ; }
  //
  const double kv = k () ;
  const double rv = r () ;
  //
  const double f1 = pow_ratio_a1 ( y , -kv ) ;
  const double f2 = pow_ratio_a2 ( y , -kv ) ;
  //
  return rv * kv * f1 * std::pow ( f2 , rv ) / ( y * abs_scale() ) ; 
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrIII::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrIII::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrIII::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrIII::raw_cdf 
( const double y ) const
{
  if ( y <= 0 ) { return 0 ; }
  //
  const double kv = k () ;
  const double rv = r () ;
  //
  return std::pow ( pow_ratio_a2 ( y , -kv ) , rv  ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrIII::tag () const 
{ 
  static const std::string s_name = "BurrIII" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () , 
				       m_k           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrIV::BurrIV
( const double r     ,
  const double c     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
  , m_c           ( c                       , "c"     , typeid ( *this ) )
{}
// ============================================================================
// evaluate Burr-IV function
// ============================================================================
double Ostap::Math::BurrIV::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  if      ( y <= 0 ) { return 0 ; }
  else if ( y >= 1 ) { return 0 ; }
  //
  const double cv = c () ;
  const double rv = r () ;
  //
  const double z = ( 1 - y ) / y ;
  //
  const double f1 = pow_ratio_a1 ( z , cv ) ;
  const double f2 = pow_ratio_a2 ( z , cv ) ;
  //
  return rv * cv * f1 * std::pow ( f2 , rv ) / ( y * ( 1 - y ) * abs_scale() ) ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrIV::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrIV::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrIV::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrIV::raw_cdf 
( const double y ) const
{
  if      ( y <= 0 ) { return 0 ; }
  else if ( y >= 1 ) { return 1 ; }
  //
  const double cv = c () ;
  const double rv = r () ;
  //
  const double z = ( 1 - y ) / y ;
  return std::pow ( pow_ratio_a2 ( z , cv ) , rv  ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrIV::tag () const 
{ 
  static const std::string s_name = "BurrIV" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () , 
				       m_c           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrV::BurrV
( const double r     ,
  const double k     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
  , m_k           ( k                       , "k"     , typeid ( *this ) )
{}
// ============================================================================
// evaluate Burr-V function
// ============================================================================
double Ostap::Math::BurrV::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  if      ( y <= -1 ) { return 0 ; }
  else if ( y >=  1 ) { return 0 ; }
  //
  const double rv = r () ;
  const double lk = m_k.logValue () ;
  //
  // exponent argument
  const double ea = lk - std::tan ( s_pi_2 * y ) ;
  // 
  const double sc = Ostap::Math::sec ( s_pi_2 * y ) ;
  //
  if ( ea <= 0 )
  {
    const double f = std::exp ( ea ) ;
    return rv * s_pi_2 * sc * sc
      * std::pow ( 1 / ( f + 1 ) , rv )
      *            f / ( f + 1 )       
      / abs_scale() ;		    
  }
  //
  const double f = std::exp ( -ea ) ;
  return rv * s_pi_2 * sc * sc
    * std::pow ( f / ( f + 1 ) , rv )
    *            1 / ( f + 1 )       
    / abs_scale() ;		    
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrV::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrV::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrV::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrV::raw_cdf 
( const double y ) const
{
  if      ( y <= -1 ) { return 0 ; }
  else if ( y >=  1 ) { return 1 ; }
  //
  const double rv = r () ;
  const double lk = m_k.logValue () ;
  //
  /// exponent argument 
  const double ea = lk - std::tan ( s_pi_2 * y ) ;
  //
  if ( ea <= 0 )
  {
    const double f  = std::exp ( ea ) ;
    return std::pow ( 1 / ( 1 + f ) , rv ) / abs_scale () ;
  }
  //
  const double f = std::exp ( - ea  ) ;
  return std::pow   ( f / ( 1 + f ) , rv ) / abs_scale () ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrV::tag () const 
{ 
  static const std::string s_name = "BurrV" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () , 
				       m_k           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrVI::BurrVI
( const double r     ,
  const double k     ,
  const double c     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
  , m_k           ( k                       , "k"     , typeid ( *this ) )
  , m_c           ( c                       , "c"     , typeid ( *this ) )
{}
// ============================================================================
// minimal x 
// ============================================================================
double Ostap::Math::BurrVI::xmin () const { return s_NEGINF ; } 
// ============================================================================
// maximal x 
// ============================================================================
double Ostap::Math::BurrVI::xmax () const { return s_POSINF ; } 
// ============================================================================
// evaluate Burr-VI function
// ============================================================================
double Ostap::Math::BurrVI::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  //
  const double rv = r () ;
  const double cv = c () ;
  const double lk = m_k.logValue () ;
  //
  const double ea = lk - cv * std::sinh ( y ) ;
  //
  double result = rv * cv / abs_scale () ;
  if ( ea < 0 )
  {
    const double f = std::exp ( ea ) ;
    result *= std::pow ( 1 / ( 1 + f ) , rv ) * f / ( 1 + f ) ;       
  }
  else
  {
    const double f = std::exp ( -ea ) ;
    result *= std::pow ( f / ( 1 + f ) , rv ) * 1 / ( 1 + f ) ;       
  }
  //
  return result * ( result ? std::cosh ( y ) : 1.0 ) ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrVI::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrVI::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrVI::cdf 
( const double x ) const
{
  const double y  = t ( x ) ;
  //
  const double cv = c () ;
  const double rv = r () ;
  const double lk = m_k.logValue() ;
  //
  const double ea = lk - cv * std::sinh ( y ) ;
  //
  double rcdf = 1 ; 
  if ( ea < 0 )
  {
    const double f = std::exp ( ea ) ;
    rcdf *= std::pow ( 1 / ( f + 1 ) , rv ) ;
  }
  else 
  {
    const double f = std::exp ( -ea ) ;
    rcdf *= std::pow ( f / ( f + 1 ) , rv ) ;      
  }
  //
  return Ostap::Math::is_positive ( scale() ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrVI::tag () const 
{ 
  static const std::string s_name = "BurrVI" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () , 
				       m_k           .tag () ,
				       m_c           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrVII::BurrVII
( const double r     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
{}
// ============================================================================
// minimal x 
// ============================================================================
double Ostap::Math::BurrVII::xmin () const { return s_NEGINF ; } 
// ============================================================================
// maximal x 
// ============================================================================
double Ostap::Math::BurrVII::xmax () const { return s_POSINF ; } 
// ============================================================================
// evaluate Burr-VII function
// ============================================================================
double Ostap::Math::BurrVII::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  //
  const double rv = r () ;
  //
  // dz/dy = sh^2
  const double sh = Ostap::Math::sech ( y ) ; 
  if ( !sh ) { return 0 ; }
  //
  if ( 1 < rv  && y < - s_TANH_XXXFLOW     ) { return 0 ; }
  const double z  = 0.5 * ( 1 + std::tanh ( y ) ) ;
  //
  return 0.5 * rv * std::pow ( z , rv - 1 ) * sh * sh / abs_scale () ;
  //
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrVII::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrVII::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrVII::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double z    = 0.5 * ( 1 + std::tanh ( y ) ) ;
  const double rv   = r () ;  
  const double rcdf = std::pow ( z , rv ) ;
  return Ostap::Math::is_positive ( scale() ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrVII::tag () const 
{ 
  static const std::string s_name = "BurrVII" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrVIII::BurrVIII
( const double r     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
{}
// ============================================================================
// minimal x 
// ============================================================================
double Ostap::Math::BurrVIII::xmin () const { return s_NEGINF ; } 
// ============================================================================
// maximal x 
// ============================================================================
double Ostap::Math::BurrVIII::xmax () const { return s_POSINF ; } 
// ============================================================================
// evaluate Burr-VIII function
// ============================================================================
double Ostap::Math::BurrVIII::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  //
  // from dz/dy 
  const double sh = Ostap::Math::sech ( y ) ;  
  if ( !sh ) { return 0 ; }
  //  
  const double z  = y <= 0 ? std::atan ( std::exp ( y ) ) : s_pi_2 - std::atan ( std::exp ( -y ) ) ;
  //
  const double rv = r () ;
  //
  return rv * std::pow ( s_2_pi * z , rv - 1 ) * sh * s_1_pi / abs_scale () ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrVIII::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrVIII::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrVIII::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  //
  const double z    = y <= 0 ? std::atan ( std::exp ( y ) ) : s_pi_2 - std::atan ( std::exp ( -y ) ) ;
  //
  const double rv   = r () ;  
  const double rcdf = std::pow ( s_2_pi * z , rv ) ;
  return Ostap::Math::is_positive ( scale() ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrVIII::tag () const 
{ 
  static const std::string s_name = "BurrVIII" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrIX::BurrIX
( const double r     ,
  const double k     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
  , m_k           ( k                       , "k"     , typeid ( *this ) )
{}
// ============================================================================
// minimal x 
// ============================================================================
double Ostap::Math::BurrIX::xmin () const { return s_NEGINF ; } 
// ============================================================================
// maximal x 
// ============================================================================
double Ostap::Math::BurrIX::xmax () const { return s_POSINF ; } 
// ============================================================================
// evaluate Burr-IX function
// ============================================================================
double Ostap::Math::BurrIX::evaluate    ( const double x ) const
{
  const double y  = t ( x ) ;
  //
  const double rv   = r () ;
  const double kv   = r () ;
  //
  double frac  = 1 ;
  double term  = 1 ;
  if ( y <= 0 )
  {
    const long double ey = std::exp ( 1.0L * y ) ;
    frac = 1.0 / ( kv * ( std::pow ( 1.0L + ey , 1.0L * rv ) - 1.0L ) + 2 ) ;
    term = ey  / ( 1 + ey ) * std::pow ( ( 1 + ey ) , rv ) ; 
  }
  else
  {
    const long double ey  = std::exp ( -1.0L * y ) ;
    const long double rey = std::pow ( ey , rv ) ;
    frac = 1.0 * rey / ( kv * ( std::pow ( ey + 1.0L , rv ) - rey ) + 2 * rey ) ;
    term = 1  / ( 1 + ey ) * std::pow ( ey / ( 1 + ey ) , -rv ) ; 
  }
  //
  return 2.0 * rv * kv * frac * frac * term / abs_scale () ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrIX::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrIX::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrIX::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  //
  const double rv   = r () ;
  const double kv   = r () ;
  //
  double frac  = 1 ;
  if ( y <= 0 )
  {
    const long double ey = std::exp ( 1.0L * y ) ;
    frac = 1.0 / ( kv * ( std::pow ( 1.0L + ey , 1.0L * rv ) - 1.0L ) + 2 ) ;
  }
  else
  {
    const long double ey  = std::exp ( -1.0L * y ) ;
    const long double rey = std::pow ( ey , rv ) ;
    frac = 1.0 * rey / ( kv * ( std::pow ( ey + 1.0L , rv ) - rey ) + 2 * rey ) ;
  }
  //
  const double rcdf = 1 - 2.0 * frac ;
  return Ostap::Math::is_positive ( scale() ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrIX::tag () const 
{ 
  static const std::string s_name = "BurrIX" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () , 
				       m_k           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrX::BurrX
( const double r     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
{}
// ============================================================================
// xmin  
// ============================================================================
double Ostap::Math::BurrX::xmin () const
{ return Ostap::Math::is_positive ( scale () ) ? x ( 0 ) : s_NEGINF ; }
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::BurrX::xmax () const
{ return Ostap::Math::is_positive ( scale () ) ? s_POSINF : x ( 0 ) ; }
// ============================================================================
// evaluate Burr-X function
// ============================================================================
double Ostap::Math::BurrX::evaluate    ( const double x ) const
{
  const double y    = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const double rv   = r () ;
  //
  const double z = std::exp ( -1 * y * y ) ;
  return 2 * rv * y * std::pow ( 1 - z , rv ) * ( z / ( 1 - z ) ) / abs_scale () ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrX::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrX::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrX::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double rv   = r () ;
  //
  const double rcdf = y < 0 ? 0.0 : std::pow ( 1 - std::exp ( -1 * y * y ) , rv ) ;
  return Ostap::Math::is_positive ( scale() ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrX::tag () const 
{ 
  static const std::string s_name = "BurrX" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrXI::BurrXI
( const double r     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
  , m_r           ( r                       , "r"     , typeid ( *this ) )
{}
// ============================================================================
// evaluate Burr-XI function
// ============================================================================
double Ostap::Math::BurrXI::evaluate    ( const double x ) const
{
  const double y    = t ( x ) ;
  if ( y < 0 || 1 < y ) { return 0 ; }
  //
  const double rv   = r () ;
  const double y2pi = s_2pi * y ; 
  const double z    = y - s_1_2pi * std::sin ( y2pi ) ;
  return rv * std::pow ( z , rv - 1 ) * ( 1 - std::cos ( y2pi ) ) / abs_scale () ;
  //
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrXI::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrXI::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrXI::cdf 
( const double x ) const
{
  const double y    = t ( x ) ;
  const double rv   = r () ;
  //
  const double rcdf =
    y <= 0 ? 0.0 :
    y >= 1 ? 1.0 :
    std::pow ( y - s_1_2pi * std::sin ( s_2pi * y ) , rv ) ;
  //
  return Ostap::Math::is_positive ( scale() ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrXI::tag () const 
{ 
  static const std::string s_name = "BurrXI" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () , 
				       m_r           .tag () ) ;
}
// ============================================================================

// ============================================================================
/*  constructor
 *  @param c           shape parameter 
 *  @param k           shape parameter 
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::BurrXII::BurrXII
( const double c     ,
  const double k     ,
  const double scale ,
  const double shift )
  : ShiftAndScale ( scale , shift , "scale" , "shift" , typeid ( *this ) , false ) 
  , m_c           ( c                       , "c"     , typeid ( *this ) ) 
  , m_k           ( k                       , "k"     , typeid ( *this ) ) 
{}
// ============================================================================
// evaluate BurrXII function
// ============================================================================
double Ostap::Math::BurrXII::evaluate    ( const double x ) const
{
  const double y = t ( x ) ;
  if ( y < 0 ) { return 0 ; }
  //
  const double cv  = c () ;
  if ( 1 < cv && s_zero ( y ) ) { return 0 ; }
  //
  const double kv = k () ;
  const double as = abs_scale() ;  
  //
  return cv * cv
    *            pow_ratio_a1 ( y , cv )
    * std::pow ( pow_ratio_a2 ( y , cv ) , kv ) / ( as * y ) ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrXII::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::BurrXII::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrXII::cdf 
( const double x ) const
{
  const double y = t ( x ) ;
  const double rcdf = raw_cdf ( y ) ;
  return Ostap::Math::is_positive ( scale () ) ? rcdf : 1 - rcdf ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::BurrXII::raw_cdf 
( const double y ) const
{
  if ( y <= 0 ) { return 0 ; }
  const double cv = c () ;
  const double kv = k () ;  
  return 1 - std::pow ( pow_ratio_a2 ( y , cv ) , kv ) ;
}
// ============================================================================
// raw-moments of unscaled/unshifted distributions 
// ============================================================================
double Ostap::Math::BurrXII::raw_moment
( const unsigned short r ) const
{
  const double cv = c () ;
  const double kv = k () ;  
  const double ck = ck * cv ;
  if ( ck <= r || s_equal ( ck , r ) ) { return s_QUIETNAN ; }
  //
  return m_k * Ostap::Math::beta ( ( ck - r ) / cv , 1 + 1 / cv ) ;
}
// ============================================================================
// mean              (for ck>1) 
// ============================================================================
double Ostap::Math::BurrXII::mean     () const
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; } 
  return x ( m1 ) ;
}
// ============================================================================
// mode
// ============================================================================
double Ostap::Math::BurrXII::mode () const
{
  const double cv = c () ;
  const double kv = k () ;
  const double ck = cv * kv ;
  const double m  = cv <= 1 ? 0.0 : std::pow ( ( cv - 1.0L ) / ( ck + 1.0L ) , 1 / cv ) ;
  return x ( m ) ;
}
// ============================================================================
// median
// ============================================================================
double Ostap::Math::BurrXII::median () const
{
  const double cv = c () ;
  const double kv = k () ;
  const double m  = std::pow ( std::pow ( 2.0L , 1.0L / kv  ) - 1.0 , 1.0 / cv  ) ;
  return x ( m ) ; 
}
// ============================================================================
// variance              (for ck>2) 
// ============================================================================
double Ostap::Math::BurrXII::variance () const
{
  //
  const double cv = c () ;
  const double kv = k () ;  
  const double ck = ck * cv ;
  if ( ck <= 2 || s_equal ( ck , 2 ) ) { return s_QUIETNAN ; }
  //
  const double ss = scale () ;
  return ss * ss * Ostap::Math::variance ( raw_moment ( 2 ) ,
					   raw_moment ( 1 ) ); 
}
// ============================================================================
// rms              (for ck>2) 
// ============================================================================
double Ostap::Math::BurrXII::rms () const
{
  //
  const double cv = c () ;
  const double kv = k () ;  
  const double ck = ck * cv ;
  if ( ck <= 2 || s_equal ( ck , 2 ) ) { return s_QUIETNAN ; }
  //
  const double vv = variance () ; 
  if ( !std::isfinite ( vv ) ) { return vv ; }
  return std::sqrt ( vv ) ;
}
// ============================================================================
// skewness              (for ck>3) 
// ============================================================================
double Ostap::Math::BurrXII::skewness () const
{
  //
  const double cv = c () ;
  const double kv = k () ;  
  const double ck = ck * cv ;
  if ( ck <= 3 || s_equal ( ck , 3 ) ) { return s_QUIETNAN ; }
  //
  const double m3 = raw_moment ( 3 ) ;
  const double m2 = raw_moment ( 2 ) ;
  const double m1 = raw_moment ( 1 ) ;
  //
  return scale_sign () * Ostap::Math::skewness ( m3 , m2 , m1 ) ;
}
// ============================================================================
// kurtosis              (for ck>4) 
// ============================================================================
double Ostap::Math::BurrXII::kurtosis () const
{
  //
  const double cv = c () ;
  const double kv = k () ;  
  const double ck = ck * cv ;
  if ( ck <= 4 || s_equal ( ck , 4 ) ) { return s_QUIETNAN ; }
  //
  const double m4 = raw_moment ( 4 ) ;
  const double m3 = raw_moment ( 3 ) ;
  const double m2 = raw_moment ( 2 ) ;
  const double m1 = raw_moment ( 1 ) ;
  //
  return Ostap::Math::kurtosis ( m4 , m3 , m2 , m1 ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BurrXII::tag () const 
{ 
  static const std::string s_name = "BurrXII" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       ShiftAndScale::tag () ,
				       m_c           .tag () ,
				       m_k           .tag () ) ;
}

// ============================================================================
//                                                                      The END
// ============================================================================


