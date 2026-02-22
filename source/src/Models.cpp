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
  : m_ss ( beta , mu   , "beta" , "mu"   , *this , false )
{}
// ============================================================================
double Ostap::Math::Gumbel::median () const 
{
  static const double s_lnln2 = std::log ( std::log ( 2.0L ) ) ;
  //
  const double _mu   = mu   () ;
  const double _beta = beta () ;
  //
  return _mu - _beta * s_lnln2 ;
}  
// ============================================================================
double Ostap::Math::Gumbel::mean () const 
{
  static const double s_gamma = M_EULER ;
  //
  const double _mu   = mu   () ;
  const double _beta = beta () ;
  //  
  return _mu  + _beta * s_gamma ;
}
// ============================================================================
double Ostap::Math::Gumbel::variance () const 
{
  static const double s_pisq6 = s_pi2 / 6.0L ;
  //
  const double _mu   = mu   () ;
  const double _beta = beta () ;
  //    
  return _beta * _beta * s_pisq6 ;
}
// ============================================================================
double Ostap::Math::Gumbel::rms () const 
{
  static const double s_pisqr6 = s_pi / std::sqrt ( 6.0L ) ;
  //
  const double _beta = beta () ;
  //      
  return std::abs ( _beta ) * s_pisqr6 ;
}
// ============================================================================
double Ostap::Math::Gumbel::skewness () const 
{
  static const double s_skew  = 
    12 * std::sqrt( 6.0L ) * Ostap::Math::zeta ( 3 ) / s_pi3  ;
  //
  const double _beta = beta () ;
  //      
  return std::copysign ( s_skew , _beta ) ;
}
// ============================================================================
// get a value for the function      
// ============================================================================
double Ostap::Math::Gumbel::pdf  ( const double x ) const 
{
  ///
  const double _mu   = mu   () ;
  const double _beta = beta () ;
  //      
  const long double ibeta = 1.0 / _beta ;
  const long double z     = ( x - _mu  ) * ibeta ;
  //
  return std::abs ( _beta ) * std::exp ( -( z + std::exp ( -z ) ) ) ;
}
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::Gumbel::cdf ( const double x ) const 
{
  //
  const double _mu   = mu   () ;
  const double _beta = beta () ;
  //        
  const long double z     = ( x - _mu ) / _beta  ;
  return 0 < _beta  ? 
    std::exp ( -std::exp ( -z ) ) : 1 - std::exp( -std::exp ( -z ) ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Gumbel::integral
( const double low  ,
  const double high ) const 
{
  if  ( s_equal ( low , high ) ) { return 0 ; }
  //
  const double _mu   = mu   () ;
  const double _beta = beta () ;
  //
  const long double ibeta = 1 / _beta ;
  const long double zmin = ( low  - _mu ) * ibeta ;
  const long double zmax = ( high - _mu ) * ibeta ;
  //
  return 0 < _beta ? 
    std::exp ( -std::exp ( -zmax ) ) - std::exp ( -std::exp ( -zmin ) ) : 
    std::exp ( -std::exp ( -zmin ) ) - std::exp ( -std::exp ( -zmax ) ) ;
}
// ============================================================================
/*  quantile function
 *  @parameter p  probability \f$ 0 < p < 1 \f$
 */
// ============================================================================
double Ostap::Math::Gumbel::quantile
( const double p ) const
{
  const double _mu   = mu   () ;
  const double _beta = beta () ;
  //
  return
    p < 0             ?  s_QUIETNAN :
	p > 1             ?  s_QUIETNAN : 
    s_zero  ( p     ) ?  s_NEGHUGE  : 
    s_equal ( p , 1 ) ?  s_POSHUGE  :
    _mu - _beta * std::log ( - std::log ( p ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Gumbel::tag () const 
{ 
  static const std::string s_name = "Gumbel" ;
  return Ostap::Utils::hash_combiner ( s_name , m_ss.tag () ) ; 
}
// ============================================================================
 
// ============================================================================
// Gram-Charlier type A
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::GramCharlierA::GramCharlierA
( const double mean   ,
  const double sigma  ,
  const double kappa3 ,
  const double kappa4 )
  : m_mean   ( mean )
  , m_sigma  ( std::fabs ( sigma ) )
  , m_kappa3 ( kappa3 )
  , m_kappa4 ( kappa4 )
//
  , m_workspace ()
//
{}
// ============================================================================
namespace 
{ 
  constexpr Ostap::Math::Hermite_<3>  s_h3{} ;
  constexpr Ostap::Math::Hermite_<4>  s_h4{} ;
}
// ============================================================================
// evaluate Gram-Charlier type A approximation
// ============================================================================
double Ostap::Math::GramCharlierA::pdf ( const double x ) const
{
  //
  const double dx = ( x - m_mean ) / m_sigma ;
  //
  const double result_0 = my_exp ( -0.5 * dx * dx ) / m_sigma * s_sqrt_1_2pi ;
  //
  double correction = 1 ;
  //
  correction += m_kappa3 * s_h3 ( dx ) /  6 ;
  //
  correction += m_kappa4 * s_h4 ( dx ) / 24 ;
  //
  return correction * result_0 ;
}
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::GramCharlierA::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::GramCharlierA::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                         0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low   ) ; } // RETURN
  //
  const double x_low  = m_mean - 5 * m_sigma ;
  const double x_high = m_mean + 5 * m_sigma ;
  //
  // split for the reasonable sub intervals:
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
  //
  // split, if the interval is too large
  //
  const double width = std::max ( std::abs  ( m_sigma)  , 0.0 ) ;
  if ( 0 < width &&  3 * width < high - low  )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GramCharlierA> s_integrator ;
  static const char s_message[] = "Integral(GramCharlierA)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,                                        // low & high edges
      workspace ( m_workspace ) ,                             // workspace
      ( high   <= x_low  ) ? s_APRECISION_TAIL :
      ( x_high <=   low  ) ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      ( high   <= x_low  ) ? s_RPRECISION_TAIL :
      ( x_high <=   low  ) ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,                                   // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
bool Ostap::Math::GramCharlierA::setM0  ( const double value )
{
  //
  if ( s_equal ( m_mean , value ) ) { return false ; }
  //
  m_mean  = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::GramCharlierA::setSigma  ( const double value )
{
  //
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( m_sigma , value_ ) ) { return false ; }
  //
  m_sigma  = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::GramCharlierA::setKappa3 ( const double value )
{
  if ( s_equal ( m_kappa3 , value )  ) { return false ; }
  //
  m_kappa3  = value ;
  //
  return false ;
}
// ============================================================================
bool Ostap::Math::GramCharlierA::setKappa4 ( const double value )
{
  if ( s_equal ( m_kappa4 , value )  ) { return false ; }
  //
  m_kappa4  = value ;
  //
  return false ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GramCharlierA::tag () const 
{ 
  static const std::string s_name = "GramCharlier" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_sigma , m_kappa3 , m_kappa4 ) ;
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
  : m_k     ( k                       , "k"     , *this )
  , m_ss    ( theta , shift , "theta" , "shift" , *this ) 
{
  // evaluate auxillary parameter 
  m_lgk = - std::lgamma ( m_k.value() ) ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::GammaDist::setK ( const double x )
{
  if ( !m_k.setValue ( x ) ) { return false ; }
  // evaluate auxillary parameter 
  m_lgk = - std::lgamma ( m_k.value() ) ;
  return true ;
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::GammaDist::mean       () const
{
  //
  const double _k     = k     () ; 
  const double _theta = theta () ;
  const double _shift = shift () ; 
  //
  return _shift + _k * _theta  ;
}
// ============================================================================
// dispersion value 
// ============================================================================
double Ostap::Math::GammaDist::variance () const
{
  //
  const double _k     = k     () ; 
  const double _theta = theta () ;
  //
  return _k * _theta * _theta  ; 
}
// ============================================================================
double Ostap::Math::GammaDist::rms     () const
{ return std::sqrt ( dispersion ()  ) ; }
// ============================================================================
double Ostap::Math::GammaDist::skewness () const
{
  const double _k = k () ; 
  return 2.0 / std::sqrt ( _k )  ;
}
// ============================================================================
// calculate gamma distribution shape
// ============================================================================
double Ostap::Math::GammaDist::evaluate ( const double x ) const
{
  const double _k     = k     () ;  
  const double _theta = theta () ;  
  const double _shift = shift () ;  
  //
  // simple cases 
  if ( x <= _shift ) { return 0 ; }
  const double t  = m_ss.t ( x ) ;  
  if ( t <= 0       ) { return 0 ; }
  //
  const double logr = -t + ( _k - 1 ) * std::log ( t )  + m_lgk ;
  return std::exp ( logr ) / _theta ;
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
  const double _shift = shift() ;
  //
  if      ( s_equal ( low  , high )  ) { return 0 ; }
  else if (           low  > high    ) { return -1 * integral ( high   , low  ) ; }
  else if (           high <= _shift ) { return 0 ; }
  else if (           low  <  _shift ) { return      integral ( _shift , high ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the cdf 
// ============================================================================
double Ostap::Math::GammaDist::cdf 
( const double x ) const 
{
  const double _k     = k     () ;  
  const double _shift = shift () ;
  //
  if ( x <= _shift ) { return 0 ; }
  const double z = m_ss.t ( x ) ;
  //
  return gsl_sf_gamma_inc_P ( _k , z ) ;
}
// ============================================================================
// calculate the mode
// ============================================================================
double Ostap::Math::GammaDist::mode () const
{
  const double _k     = k     () ;  
  const double _theta = theta () ;  
  const double _shift = shift () ;  
  //
  return _shift + ( 1 <= _k  ? ( _k  - 1 ) * _theta : 0.0 ) ;
}
// ============================================================================
// calculate the quantile   (0<p<1) 
// ============================================================================
double Ostap::Math::GammaDist::quantile ( const double p ) const 
{
  const double _k     = k     () ;  
  const double _theta = theta () ;  
  const double _shift = shift () ;
  //
  return
    p <  0            ? s_QUIETNAN :
    p >  1            ? s_QUIETNAN :
    s_zero  ( p     ) ? _shift    :
    s_equal ( p , 1 ) ? s_POSHUGE  :    
    _shift + gsl_cdf_gamma_Pinv ( p , _k , _theta ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GammaDist::tag () const 
{ 
  static const std::string s_name = "GammaDist" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       m_k .tag () , 
				       m_ss.tag () ) ;
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
  const double _k     = k     () ;  
  const double _shift = shift () ;
  //
  return ! _shift ? 2   * _k     : s_QUIETNAN ;
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
  const double _theta = theta () ;  
  const double _shift = shift () ;
  //
  return !_shift ? 0.5 * _theta : s_QUIETNAN  ;
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
  const double low   ) 
  : m_k     ( std::abs ( k     ) ) 
  , m_theta ( std::abs ( theta ) ) 
  , m_p     ( std::abs ( p     ) ) 
  , m_low   ( low ) 
{}
// ============================================================================
bool Ostap::Math::GenGammaDist::setK     ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_k ) ) { return false ; }
  m_k   = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGammaDist::setTheta ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_theta ) ) { return false ; }
  m_theta = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGammaDist::setP    ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_p ) ) { return false ; }
  m_p     = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGammaDist::setLow ( const double value ) 
{
  if ( s_equal ( value , m_low ) ) { return false ; }
  m_low   = value ;
  return true ;
}
// ============================================================================
double Ostap::Math::GenGammaDist::pdf ( const double x ) const 
{
  //
  if ( x <= m_low || s_equal ( x , m_low ) ) { return 0 ; }
  //
  const double xc = ( x - m_low ) / m_p ;  
  const double xt = std::pow ( xc , m_p ) ;  
  //
  double result   = ( k () - 1 ) * gsl_sf_log ( xc ) - xt ;
  result         +=  gsl_sf_log     ( m_p  / m_theta ) ;
  result         -=  gsl_sf_lngamma ( m_k  / m_p     ) ;
  //return gsl_sf_exp ( result ) ;
  return my_exp ( result ) ;
}
// ============================================================================
double Ostap::Math::GenGammaDist::cdf ( const double x ) const 
{
  //
  if ( x <= m_low || s_equal ( x , m_low ) ) { return 0 ; }
  //
  const double xc = ( x - m_low ) / m_theta ;  
  const double xt = std::pow ( xc , m_p  ) ;
  //
  return gsl_sf_gamma_inc_P ( k () / m_p , xt ) ;
}
// ============================================================================
double Ostap::Math::GenGammaDist::integral () const { return 1 ; }
// ============================================================================
double Ostap::Math::GenGammaDist::integral
( const double low  , 
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;  
}
// ============================================================================
double Ostap::Math::GenGammaDist::mode () const
{
  const double d_ = d () ;
  const double a_ = a () ;
  return 1 < d_ ? a_ * std::pow ( ( d_ - 1 ) / m_p , 1 / m_p ) : 0.0 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenGammaDist::tag () const 
{ 
  static const std::string s_name = "GenGammaDist" ;
  return Ostap::Utils::hash_combiner ( s_name , m_k , m_theta , m_p , m_low ) ;
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
  : m_a     (            a       ) 
  , m_theta ( theta              ) 
  , m_alpha ( std::abs ( alpha ) ) 
  , m_beta  (            beta    ) 
{}
// ============================================================================
bool Ostap::Math::Amoroso::setA      ( const double value ) 
{
  if ( s_equal ( value , m_a ) ) { return false ; }
  m_a = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Amoroso::setTheta ( const double value ) 
{
  if ( s_equal ( value , m_theta ) ) { return false ; }
  m_theta = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Amoroso::setAlpha ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Amoroso::setBeta ( const double value ) 
{
  if ( s_equal ( value , m_beta ) ) { return false ; }
  m_beta  = value ;
  return true ;
}
// ============================================================================
// evaluate Amoroso distribtion
// ============================================================================
double Ostap::Math::Amoroso::pdf ( const double x ) const 
{
  //
  if      ( m_theta > 0 && ( x <= m_a || s_equal ( x , m_a ) ) ) { return 0 ; }
  else if ( m_theta < 0 && ( x >= m_a || s_equal ( x , m_a ) ) ) { return 0 ; }
  //
  const double xc = ( x - m_a ) / m_theta ;
  const double xt = std::pow ( xc , m_beta ) ;
  //
  double result   = ( m_alpha * m_beta - 1 ) * gsl_sf_log ( xc )  - xt ; 
  result += gsl_sf_log     ( std::abs ( m_beta  / m_theta ) ) ;
  result -= gsl_sf_lngamma (            m_alpha ) ;
  //
  return my_exp ( result ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::cdf ( const double x ) const 
{
  //
  if      ( theta () > 0 && ( x <= m_a || s_equal ( x , m_a ) ) ) { return 0 ; }
  else if ( theta () < 0 && ( x >= m_a || s_equal ( x , m_a ) ) ) { return 1 ; }
  //
  const double xc = ( x - m_a ) / theta ()    ;
  const double xt = std::pow ( xc , beta() ) ;
  //
  return 
    beta() * theta() > 0 ? 
    1 - gsl_sf_gamma_inc_Q ( alpha() , xt ) :
    gsl_sf_gamma_inc_Q ( alpha() , xt ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::integral () const { return 1 ; }
// ============================================================================
double Ostap::Math::Amoroso::integral ( const double low  , 
                                        const double high ) const 
{
  if ( s_equal ( low ,high ) ) { return 0 ; }
  return  cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::mode () const 
{
  if ( alpha() * beta() <= 1 ) { return a () ; }
  return a () + theta() * std::pow ( alpha() - 1./beta() , 1./beta () ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::mean () const 
{
  const double x = alpha() + 1/beta() ;
  if ( x <= 0 || s_equal ( x , 0 )  ) { return -1.e+9 ; } // RETURN 
  //
  if ( x       < 0.2 * GSL_SF_GAMMA_XMAX && 
       alpha() < 0.2 * GSL_SF_GAMMA_XMAX  ) 
  {
    return a () + theta() * gsl_sf_gamma ( x ) / gsl_sf_gamma ( alpha() ) ;
  }
  //
  double aux = gsl_sf_lngamma ( x       ) ;
  aux -= gsl_sf_lngamma       ( alpha() ) ;
  //
  return a() + theta() * gsl_sf_exp ( aux ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::variance () const 
{
  //
  const double x2 = alpha() + 2/beta() ;
  if ( x2 <= 0 || s_equal ( x2 , 0 )  ) { return -1.e+9 ; } // RETURN 
  //
  const double x1 = alpha() + 1/beta() ;
  if ( x1 <= 0 || s_equal ( x1 , 0 )  ) { return -1.e+9 ; } // RETURN 
  //
  //
  if ( x1      < 0.2 * GSL_SF_GAMMA_XMAX && 
       x2      < 0.2 * GSL_SF_GAMMA_XMAX && 
       alpha() < 0.2 * GSL_SF_GAMMA_XMAX  ) 
  {
    const double ga  = gsl_sf_gamma ( alpha () ) ;
    const double gx1 = gsl_sf_gamma ( x1       ) ;
    const double gx2 = gsl_sf_gamma ( x2       ) ;
    //
    return theta2() * ( gx2 / ga - Ostap::Math::POW ( gx1 / ga , 2 ) ) ;
  }
  //
  const double lnga = gsl_sf_lngamma ( alpha () ) ;
  //
  double aux1  = gsl_sf_lngamma ( x1   ) ;
  aux1        -= lnga ;
  aux1         = gsl_sf_exp     ( aux1 ) ;
  //
  double aux2  = gsl_sf_lngamma ( x2   ) ;
  aux2        -= lnga ;
  aux2         = gsl_sf_exp     ( aux2 ) ;
  //
  return theta2() * ( aux2 - aux1 * aux1 ) ;
}
// ============================================================================
double Ostap::Math::Amoroso::sigma () const 
{
  //
  const double x2 = alpha() + 2/beta() ;
  if ( x2 <= 0 || s_equal ( x2 , 0 )  ) { return -1.e+9 ; } // RETURN 
  //
  const double x1 = alpha() + 1/beta() ;
  if ( x1 <= 0 || s_equal ( x1 , 0 )  ) { return -1.e+9 ; } // RETURN 
  //
  return std::sqrt ( variance() ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Amoroso::tag () const 
{ 
  static const std::string s_name = "Amoroso" ;
  return Ostap::Utils::hash_combiner ( s_name , m_a , m_theta , m_alpha , m_beta ) ; 
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
double Ostap::Math::LogGammaDist::operator() ( const double x ) const
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
double Ostap::Math::Log10GammaDist::operator() ( const double x ) const
{ return LogGammaDist::operator() ( x * s_ln10 ) * s_ln10 ; }
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
( const double loga ,
  const double logb  ,
  const double scale ,
  const double shift )
  : m_ab ( loga  , logb  , "a"     , "b"     , typeid ( *this ) )
  , m_ss ( scale , shift , "scale" , "shift" , typeid ( *this ) )
{}
// ============================================================================
// evaluate beta-distribution
// ============================================================================
double Ostap::Math::Beta::evaluate  ( const double x ) const
{
  const double tt = t ( x ) ;
  if ( tt < 0 || 1 < tt ) { return 0 ; }
  //
  const double _a = a () ;
  const double _b = b () ;
  //
  // density diverges
  if ( _a < 1 && s_zero (     tt ) ) { return std::numeric_limits<double>::max () ; }
  if ( _b < 1 && s_zero ( 1 - tt ) ) { return std::numeric_limits<double>::max () ; }
  //
  if ( 1 < _a && s_equal ( 1 + tt , 1 ) ) { return 0 ; }
  if ( 1 < _b && s_equal ( 0 + tt , 1 ) ) { return 0 ; }
  //
  const double a = std::pow ( 0.0L + tt , _a - 1.0L ) ;
  const double b = std::pow ( 1.0L - tt , _b - 1.0L ) ;
  //
  const double _scale = scale () ;
  //
  return a * b * m_ab.inv_Beta_pq () / _scale ;
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::Beta::cdf ( const double x ) const 
{
  //
  if      ( x <= xmin () ) { return 0 ; } 
  else if ( x >= xmax () ) { return 1 ; }
  //
  const double tt = t ( x ) ;
  //
  if      ( tt <= 0 || s_equal ( 1 + tt , 1 ) ) { return 0 ; }
  else if ( tt >= 1 || s_equal ( 0 + tt , 1 ) ) { return 1 ; }
  //
  const double _a = a () ;
  const double _b = b () ;
  //
  return beta_cdf ( tt , _a , _b ) ;
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
  else if ( high <= xmin ()         ) { return  0 ; }
  else if ( low  >= xmax ()         ) { return  0 ; }
  //s
  return cdf ( high ) - cdf ( low ) ;
}
// ===========================================================================
// mean
// ===========================================================================
double Ostap::Math::Beta::mean () const
{
  //
  const double _a = a () ;
  const double _b = b () ;
  //  
  const double  value = _a / ( _a + _b ) ;
  return x ( value ) ;
}
// ===========================================================================
// mode 
// ===========================================================================
double Ostap::Math::Beta::mode () const
{
  //
  const double _a = a () ;
  const double _b = b () ;
  //  
  const double  value =
    _a <  1            ? 0   :
    _b <  1            ? 1   :
    1 == _a && 1 == _b ? 0.5 :
    ( _a - 1 ) / ( _a + _b - 2 ) ;
  //
  return x ( value ) ;
}
// ===========================================================================
// median 
// ===========================================================================
double Ostap::Math::Beta::median () const
{
  const double _a = a () ;
  const double _b = b () ;
  //  
  return s_equal ( _a  , _b  ) ? x ( 0.5 ) : quantile ( 0.5 ) ;
}
// ===========================================================================
// variance 
// ===========================================================================
double Ostap::Math::Beta::variance () const
{
  const double _a     = a     () ;
  const double _b     = b     () ;
  const double _scale = scale () ;
  //  
  const double value = _a * _b /
  ( std::pow ( _a + _b  , 2 ) * ( _a + _b + 1 ) );
  //
  return value  * _scale * _scale ;
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
  const double _a     = a     () ;
  const double _b     = b     () ;
  //
  return
    2 * ( _b - _a )   
    * std::sqrt ( ( _a + _b + 1 ) / ( _a * _b ) )
    / ( _a + _b + 2 ) ;   
}
// ===========================================================================
// (excess) kurtosis
// ===========================================================================
double Ostap::Math::Beta::kurtosis () const
{
  //
  const double _a     = a     () ;
  const double _b     = b     () ;
  //  
  const double ab  = _a * _b     ;
  const double ab1 = _a + _b + 1 ;
  const double ab2 = ab1 + 1 ;
  const double ab3 = ab2 + 1 ;
  //
  return 6 * ( std::pow ( _a - _b , 2 ) * ab1 - ab * ab2 ) / ( ab * ab2 * ab3 ) ;
}
// ===========================================================================
// quantile \f$ 0 \le p \le 1 \f$
// ===========================================================================
double Ostap::Math::Beta::quantile ( const double p  ) const
{
  //
  if      ( p <= 0 ) { return x ( 0 ) ; } 
  else if ( p >= 1 ) { return x ( 1 ) ; } 
  //
  const double _a     = a     () ;
  const double _b     = b     () ;
  //  
  const double tt = beta_quantile ( p , _a , _b ) ;
  return x ( tt ) ;
}
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Beta::tag () const 
{ 
  static const std::string s_name = "Beta" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       m_ab.tag () ,
				       m_ss.tag () ) ;
}
// ============================================================================

// ============================================================================
// Beta' 
// ============================================================================
/*  constructor with all parameters 
 *  @param alpha \f$\alpha\f$-parameter 
 *  @param beta  \f$\beta\f$-parameter 
 *  @param scale  scale-parameter
 *  @param shift  shift-parameter
 */
// ============================================================================
Ostap::Math::BetaPrime::BetaPrime
( const double loga  ,
  const double logb  ,
  const double scale ,
  const double shift )
  : m_ab ( loga  , logb  , "a"     , "b"     , typeid ( *this )         )
  , m_ss ( scale , shift , "scale" , "shift" , typeid ( *this ) , false )
{}
// ============================================================================
// evaluate beta'-distributions 
// ============================================================================
double Ostap::Math::BetaPrime::evaluate ( const double x ) const 
{
  //
  const double _scale = scale () ;
  const double _shift = shift () ;
  //
  if      ( _scale >= 0 && x <= _shift ) { return 0 ; }
  else if ( _scale <= 0 && x >= _shift ) { return 0 ; }
  else if ( s_equal ( x , _shift )     ) { return 0 ; }
  //
  const double y = ( x - _shift ) / _scale ;
  //
  const double _a = a () ;
  const double _b = b () ;
  //
  return m_ab.inv_Beta_pq ()  / std::abs ( _scale ) 
    * std::pow ( 0.0L + y ,   _a -  1 ) 
    * std::pow ( 1.0L + y , - _a - _b ) ;  
}
// ============================================================================
double Ostap::Math::BetaPrime::cdf ( const double x ) const 
{
  //
  const double _scale = scale () ;
  const double _shift = shift () ;
  //
  const double z = ( x - _shift ) /_scale ;
  //
  if ( z <= 0 || s_equal ( z , 0 ) ) { return 0 ; }
  //
  const double y = z / ( 1 + z ) ;
  //
  const double _a = a () ;
  const double _b = b () ;
  //
  Sentry sentry ;
  //
  return
    0 < _scale ? 
    gsl_sf_beta_inc ( _a , _b , y ) : 1 - gsl_sf_beta_inc ( _a , _b , y ) ; 
}
// ============================================================================
double Ostap::Math::BetaPrime::integral
( const double low  , 
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::integral () const { return 1 ; } 
// ============================================================================
double Ostap::Math::BetaPrime::mean () const 
{
  const double _b = b () ;
  if ( _b <= 1 || s_equal ( _b , 1 ) ) { return s_QUIETNAN ; } 
  //
  const double _a     = a () ;
  const double _scale = scale () ;
  const double _shift = shift () ;
  //
  return _shift + _scale * _a  / ( _b - 1 ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::mode () const 
{
  const double _a     = a () ;
  //
  if ( _a < 1 ) { return 0 ; }
  //
  const double _b     = b () ;
  const double _scale = scale () ;
  const double _shift = shift () ;
  //  
  return _shift + _scale * ( _a - 1 ) / ( _b + 1 ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::variance () const 
{
  const double _b     = b () ;
  if ( _b <= 2 || s_equal ( _b , 2 ) ) { return s_QUIETNAN ; } 
  //
  const double _a     = a () ;
  const double _scale = scale () ;
  //
  return _scale * _scale * _a *  ( _a + _b + 1 ) / ( _b - 2 ) / Ostap::Math::POW ( _b - 1 , 2 ) ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::rms () const 
{
  const double _b     = b () ;
  if ( _b <= 2 || s_equal ( _b , 2 ) ) { return s_QUIETNAN ; } 
  return std::sqrt ( variance () ) ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::skewness  () const 
{
  const double _b     = b () ;
  if ( _b <= 3 || s_equal ( _b , 3 ) ) { return s_QUIETNAN ; } 
  //
  const double _a = a () ;
  //
  return 2 * ( 2 * _a + _b - 1 ) / ( _b - 3 ) * std::sqrt( ( _b - 2 ) / _a / ( _a + _b - 1 ) ) ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::kurtosis  () const 
{
  const double _b     = b () ;
  if ( _b <= 4 || s_equal ( _b , 4 ) ) { return s_QUIETNAN ; }
  //
  const double _a = alpha () ;
  //  
  return 6 * ( _a * ( _a + _b - 1 ) * ( 5 * _b - 11 ) + ( _b - 1 ) * ( _b - 1 ) * ( _b - 2 ) ) /
    ( _a * ( _a + _b - 1 ) * ( _b - 3 ) * ( _b - 4 ) ) ;
}
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BetaPrime::tag () const 
{ 
  static const std::string s_name = "BetaPrime" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       m_ab.tag () ,
				       m_ss.tag () ) ;
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
 *  @param alpha \f$\alpha\f$-parameter 
 *  @param p     p-parameter 
 *  @param q     p-parameter 
 *  @param beta  \f$\beta\f$-parameter 
 *  @param scale  scale-parameter
 *  @param shift  shift-parameter
 */
// ============================================================================
Ostap::Math::GenBetaPrime::GenBetaPrime 
( const double alpha , 
  const double beta  ,
  const double p     , 
  const double q     , 
  const double scale , 
  const double shift )
  : m_alpha ( std::abs ( alpha ) )
  , m_beta  ( std::abs ( beta  ) )
  , m_p     ( std::abs ( p     ) ) 
  , m_q     ( std::abs ( q     ) ) 
  , m_scale ( scale )
  , m_shift ( shift )
  , m_aux       ()
  , m_workspace () 
{
  Ostap::Assert ( 0 < m_alpha && !s_zero ( m_alpha ) ,
		  "Invalid value of alpha (must be positive" ,
		  "Ostap::Math::GenBetaPrime" ) ;
  Ostap::Assert ( 0 < m_beta  && !s_zero ( m_beta  ),
		  "Invalid value of beta  (must be positive" ,
		  "Ostap::Math::GenBetaPrime" ) ;  
  Ostap::Assert ( 0 < m_p     && !s_zero ( m_p  ) ,
		  "Invalid value of p    (must be positive" ,
		  "Ostap::Math::GenBetaPrime" ) ;
  Ostap::Assert ( 0 < m_q     && !s_zero ( m_q ),
		  "Invalid value of q    (must be positive" ,
		  "Ostap::Math::GenBetaPrime" ) ;  
  // 
  m_aux = 1.0 / Ostap::Math::beta ( m_alpha , m_beta ) ;
}
// ============================================================================
bool Ostap::Math::GenBetaPrime::setAlpha ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  m_aux   = 1.0/Ostap::Math::beta ( m_alpha , m_beta ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenBetaPrime::setBeta  ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_beta  ) ) { return false ; }
  m_beta  = value_ ;
  m_aux   = 1.0/Ostap::Math::beta ( m_alpha , m_beta ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenBetaPrime::setScale ( const double value ) 
{
  if ( s_equal ( value , m_scale  ) ) { return false ; }
  m_scale  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenBetaPrime::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift  ) ) { return false ; }
  m_shift  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenBetaPrime::setP  ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_p  ) ) { return false ; }
  m_p  = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenBetaPrime::setQ  ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_q  ) ) { return false ; }
  m_q  = value_ ;
  return true ;
}
// ============================================================================
// evaluate beta'-distributions 
// ============================================================================
double Ostap::Math::GenBetaPrime::pdf ( const double x ) const 
{
  //
  if      ( m_scale >= 0 && x <= m_shift ) { return 0 ; }
  else if ( m_scale <= 0 && x >= m_shift ) { return 0 ; }
  else if ( s_equal ( x , m_shift )      ) { return 0 ; }
  //
  const double y = ( x - m_shift ) / m_scale ;
  const double z = y / m_q ;
  //
  return m_aux * m_p / std::abs ( m_q * m_scale )
    * std::pow ( z                         ,   m_alpha * m_p - 1 )
    * std::pow ( 1 + std::pow ( z , m_p )  , - m_alpha - m_beta  ) ;  
}
// ============================================================================
double Ostap::Math::GenBetaPrime::integral () const { return 1 ; }
// ============================================================================
double Ostap::Math::GenBetaPrime::cdf ( const double x ) const
{ return integral ( m_shift , x ) ; }
// ============================================================================
double Ostap::Math::GenBetaPrime::integral
( const double low  , 
  const double high ) const 
{
  if      ( s_equal ( low , high )             ) { return 0 ; }
  else if ( 0 < m_scale && high    <= m_shift ) { return 0 ; }
  else if ( 0 > m_scale && m_shift <= low     ) { return 0 ; }
  else if ( high < low                         ) { return -integral ( high , low ) ; }
  //
  const double xmin = 0 < m_scale ? std::max ( low , m_shift ) : low   ;
  const double xmax = 0 < m_scale ? high : std::min ( high , m_shift ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenBetaPrime> s_integrator ;
  static const char s_message[] = "Integral(GenBetaPrime)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xmin    , xmax ,                                        // low & high edges
      workspace ( m_workspace ) ,                             // workspace
      s_APRECISION , // absolute precision
      s_RPRECISION , // relative precision
      m_workspace.size()              ,                                   // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
double Ostap::Math::GenBetaPrime::mean () const 
{
  const double bp = m_beta * m_p ;
  if ( bp <= 1 || s_equal ( bp , 1 ) )
    { return std::numeric_limits<double>::quiet_NaN(); } 
  //
  const double ip = 1.0/m_p ;
  const double logm =
    Ostap::Math::lgamma ( m_alpha + ip ) + 
    Ostap::Math::lgamma ( m_beta  - ip ) -
    Ostap::Math::lgamma ( m_alpha      ) -
    Ostap::Math::lgamma ( m_beta       ) ;
  //
  return m_shift + m_scale * m_q * std::exp ( logm ) ;
}
// ============================================================================
double Ostap::Math::GenBetaPrime::mode () const 
{
  if ( m_alpha * m_p <= 1 ) { return m_shift ; }
  const double z = ( m_alpha * m_p - 1 ) / ( m_beta * m_p + 1 ) ;
  return m_shift + m_scale * m_q * std::pow ( z , 1.0/m_p ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBetaPrime::tag () const 
{ 
  static const std::string s_name = "GenBetaPrime" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_alpha ,
				       m_beta  ,
				       m_p     ,
				       m_q     ,
				       m_scale ,
				       m_shift ) ; 
}
// ============================================================================

// ============================================================================
// the most general form
// ============================================================================
Ostap::Math::GenBeta::GenBeta
( const double a     , // shape 
  const double gamma , // c = sin^2 gamma
  const double logp  , // shape 
  const double logq  , // shape 
  const double scale , 
  const double shift )
  : m_a      ( a     ,  "a"                    , typeid ( *this ) , false )
  , m_pq     ( logp  , logq    , "p" , "q"     , typeid ( *this ) )
  , m_ss     ( scale , shift   , "b" , "shift" , typeid ( *this ) )
  , m_c      ( gamma ,  0  , 1 , "c"           , typeid ( *this ) )
    //
  , m_c1     (  false ) 
  , m_ac     ( -1     ) 
{
  //
  m_c1 = s_equal  ( m_c.value() , 1 ) ;
  m_ac = m_c1 ? -1.0 : std::pow ( 1.0L - m_c.value() , 1.0L / m_a.value()  ) ;
  //
}
// ============================================================================
// set parameter A 
// ============================================================================
bool Ostap::Math::GenBeta::setA ( const double value ) 
{
  if ( !m_a.setValue ( value ) ) { return false ; }
  m_ac = m_c1 ? -1.0 : std::pow ( 1.0L - m_c.value () , 1.0L / m_a.value ()  ) ;
  return true ;
}
// ============================================================================
// set parameter C
// ============================================================================
bool Ostap::Math::GenBeta::setC ( const double value ) 
{
  if ( !m_c.setValue ( value ) ) { return false ; }
  //
  m_c1 = s_equal ( m_c.value( ) , 1.0 ) ;
  m_ac = m_c1 ? -1.0 : std::pow ( 1.0L - m_c.value () , 1.0L / m_a.value ()  ) ;
  //
  return true ;
}
// ============================================================================
// set parameter Gamma/C
// ============================================================================
bool Ostap::Math::GenBeta::setGamma ( const double value )
{
  if ( !m_c.setExternal ( value ) ) { return value ; }
  //
  m_c1 = s_equal ( m_c.value( ) , 1.0 ) ;
  m_ac = m_c1  ? -1.0 : std::pow ( 1.0L - m_c.value()  , 1 / m_a.value () ) ;
  //
  return true ; 
}
// ============================================================================
// xmin
// ============================================================================
double Ostap::Math::GenBeta::xmin() const 
{
  const double _a     = a     ()  ;
  const double _shift = shift ()  ;
  if ( 0 < _a || m_c1 ) { return _shift ; }
  const double _b     = b     ()  ;  
  return _shift + _b / m_ac ;
}
// ============================================================================
// xmax
// ============================================================================
double Ostap::Math::GenBeta::xmax() const 
{
  const double _a     = a     ()  ;  
  if ( 0 < _a && !m_c1 )
    {
      const double _b     = b     ()  ;  
      const double _shift = shift ()  ;
      return _shift + _b / m_ac ;
    } 
  return s_POSINF ;
}
// ============================================================================
// evaluate GenBeta-distribution
// ============================================================================
double Ostap::Math::GenBeta::evaluate ( const double x ) const 
{
  //
  const double _b     = b     ()  ;  
  const double _shift = shift ()  ;
  //
  const double z = ( x - _shift ) / _b ;
  if ( z <= 0  || s_zero ( z )  ) { return 0 ; }             // RETURN 
  //
  if      ( x <= xmin ()                   ) { return 0 ; }
  else if ( finite_range() && xmax() <= x  ) { return 0 ; }
  //
  const double _a     = a     ()  ;  
  const double _p     = p     ()  ;  
  const double _q     = q     ()  ;  
  //
  const long double za   = std::pow ( z * 1.0L , m_a * 1.0L ) ;
  if ( m_c1 ) 
  {
    const double lnR = - _p * std::log ( 1.0L + 1.0L / za ) - _q * std::log ( 1.0L + 1.0 * za ) ;
    return std::abs ( m_a ) / _b *  std::exp ( lnR - m_pq.log_Beta_pq () ) / z ;
  }
  //
  const double _c     = c     ()  ;  
  //
  const long double d = 1.0L - ( 1.0L - _c ) * za ;
  if ( d <= 0 || s_zero ( d ) ) { return 0 ; } 
  //
  const long double l1 =      std::log ( d              ) ;
  const long double l2 =      std::log ( 1.0 + _c * za ) ;
  const long double l3 = _a * std::log ( z * 1.0L       ) ;
  //
  const double lnR = - _p * l3 + ( _q - 1 ) *  l1 - ( _q + _p ) * l2 ;  
  return std::abs ( _a ) / _b *  std::exp ( lnR - m_pq.log_Beta_pq () ) / z ;
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
  const double _shift = shift () ;
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return - integral ( high , low ) ; }
  else if ( high <=  _shift        ) { return 0 ; }
  ///
  const double xmn = xmin () ;
  if       ( high <= xmn               ) { return 0 ; }
  else if  ( low  <  xmn && xmn < high ) { return integral ( xmn , high ) ; }
  //
  const double d2 = high - low ;
  if ( finite_range () )
  {
    const double xmx = xmax () ;
    if      ( xmx <= low              ) { return 0 ; }
    else if ( low < xmx && xmx < high ) { return integral ( low , xmx  ) ; }
    //
    static const std::array<double,5> s_split1 { 0.1 , 0.25 , 0.5 , 0.75 , 0.9 } ;
    for ( const double s : s_split1 )
    {
      const double split = split * xmn + ( 1 - split ) * xmx ;
      if ( low < split && split < high )
      { return integral ( low , split ) + integral ( split  , high ) ; } 
    }
  }
  //
  const double width = 2 * b ()  ;
  static const std::array<double,5> s_splits { 3 , 6 , 15 , 30 } ;
  for ( const double s : s_splits )
  {
    const double split = xmn + s * width  ;
    if ( low < split && split  < high ) { return integral ( low , split ) + integral ( split , high ) ; }  
  }
  //
  const bool in_tail = !finite_range() && ( xmn + s_splits.back () * width <= low ) ;  
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenBeta> s_integrator ;
  static const char s_message[] = "Integral(GenBeta)" ;
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
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()        , // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the CDF 
// ============================================================================
double Ostap::Math::GenBeta::cdf 
( const double x    ) const
{
  const double xmn = xmin () ;
  if      ( x <= xmn                       ) { return 0 ; }
  else if ( finite_range() && xmax () <= x ) { return 0 ; }
  return integral ( xmn , x ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBeta::tag () const 
{ 
  static const std::string s_name = "GenBeta" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_a.tag  () ,
				       m_pq.tag () ,
				       m_ss.tag () ,
				       m_c .tag () ) ;
}
// ============================================================================


// ============================================================================
// Generalized beta of first kind 
// ============================================================================
Ostap::Math::GenBeta1::GenBeta1
( const double a     , // shape 
  const double logp  , // shape 
  const double logq  , // shape 
  const double scale , // scale 
  const double shift )
  : m_a      ( a                   , "a"     , typeid ( *this ) , false )
  , m_pq     ( logp  , logq  , "p" , "q"     , typeid ( *this ) )
  , m_ss     ( scale , shift , "b" , "shift" , typeid ( *this ) )
{}
// ============================================================================
// x-max
// ============================================================================
double Ostap::Math::GenBeta1::xmax () const
{
  //
  const double _a      = a     () ;
  const double _shift  = shift () ;
  const double _scale  = scale () ;
  //
  return 0 <_a ? _shift + _scale : s_POSINF ;
}
// ============================================================================
// evaluate GenBeta1-distribution
// ============================================================================
double Ostap::Math::GenBeta1::evaluate
( const double x ) const 
{
  //
  const double _a      = a     () ;
  const double _shift  = shift () ;
  const double _scale  = scale () ;
  //
  if      ( 0 < _a && ( x <= _shift || _shift + _scale <= x ) ) { return 0 ; }
  else if ( 0 > _a &&                  _shift + _scale >= x   ) { return 0 ; }
  //
  const double z = m_ss.t ( x ) ;
  //
  if      ( 0 < _a && ( z <=  0 || 1 <= z ) ) { return 0 ; }
  else if ( 0 > _a &&              z <= 1   ) { return 0 ; }
  //
  const long double za  = std::pow ( z * 1.0L , _a * 1.0L ) ;
  //
  const long double lz1 = std::log ( 0.0L + za ) ;
  const long double lz2 = std::log ( 1.0L - za ) ; 
  //
  const double _p      = p     () ;
  const double _q      = q     () ;
  //
  const long double lr = ( _p - 1.0L ) * lz1 + ( _q - 1.0L ) * lz2 - m_pq.log_Beta_pq () ;
  //
  return std::abs ( m_a ) * std::exp ( lr ) / _scale ;
  //
}
// ============================================================================
// mean-value @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta1::mean () const
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  return shift () + scale () * m1 ;
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
  const double _scale = scale () ; 
  return _scale * _scale * Ostap::Math::variance ( m2 , m1 ) ;
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
  return Ostap::Math::skewness ( m3 , m2 , m1 ) ;
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
  const double _a = a () ;
  const double _p = p () ;
  const double _q = q () ;
  //
  const double _pa = _p + k / _a ;
  //
  return
    _pa <= 0 || s_zero ( _pa      ) ? s_QUIETNAN       : 
    Ostap::Math::beta  ( _pa , _q ) * m_pq.inv_Beta_pq () ; 
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
  //
  const double _shift  = shift () ;
  //
  if       ( s_equal ( low , high ) ) { return 0 ; }
  else if  ( high <  low            ) { return - integral ( high , low ) ; }
  else if  ( high <= _shift         ) { return 0 ; }
  ///
  const double xmn = xmin () ;  
  if       ( high <= xmn ) { return 0                       ; }
  else if  ( low  <  xmn ) { return integral ( xmn , high ) ; }
  //
  /// (1) split at mean if exists and finite
  { //
    const double vmean = mean () ;
    if ( std::isfinite ( vmean ) && low < vmean && vmean < high )
    { return integral ( low , vmean ) + integral ( vmean , high ) ; } 
  }
  //
  const double _a = a () ;
  const double _b = b () ;
  //
  // interval length 
  const double d2 = high - low ;
  if ( finite_range () )
  {
    const double xmx = xmax () ;
    if      ( xmx <= low  ) { return 0 ; }
    else if ( xmx <  high ) { return integral ( low , xmx ) ; }
    //
    static const std::array<double,5> s_split1 { 0.1 , 0.25 , 0.5 , 0.75 , 0.9 } ;
    for ( const double s : s_split1 )
    {
      const double split = split * xmn + ( 1 - split ) * xmx ;
      if ( low < split && split < high )
      { return integral ( low , split ) + integral ( split  , high ) ; } 
    }
  }
  //
  // check the sigma
  const double sigma = rms () ;
  const double width = std::isfinite ( sigma ) ? sigma : 2 * _b ;
  //
  static const std::array<double,3> s_splits { 3 , 6 , 15 } ;
  for ( const double s : s_splits )
  {
    const double split = xmn + s * width  ;
    if ( low < split && split  < high ) { return integral ( low , split ) + integral ( split , high ) ; }  
  }
  //
  const bool in_tail = !finite_range() && ( xmn + s_splits.back () * width <= low ) ;  
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenBeta1> s_integrator ;
  static const char s_message[] = "Integral(GenBeta1)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,                             // low & high edges
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
// get the CDF 
// ============================================================================
double Ostap::Math::GenBeta1::cdf 
( const double x    ) const
{
  const double _shift  = shift () ;
  if      ( x <= _shift                    ) { return 0 ; }
  //
  const double xmn = xmin () ;
  if      ( x <= xmn                       ) { return 0 ; }
  else if ( finite_range() && xmax () <= x ) { return 1 ; }
  //
  return integral ( xmn , x ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBeta1::tag () const 
{ 
  static const std::string s_name = "GenBeta1" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       m_a .tag () ,
				       m_pq.tag () ,
				       m_ss.tag () ) ;
}
// ============================================================================

// ============================================================================
// Generalized beta of second kind 
// ============================================================================
Ostap::Math::GenBeta2::GenBeta2
( const double a     , // shape 
  const double logp  , // shape 
  const double logq  , // shape 
  const double scale , // scale 
  const double shift )
  : m_a      ( a                   , "a"     , typeid ( *this ) , false )
  , m_pq     ( logp  , logq  , "p" , "q"     , typeid ( *this ) )
  , m_ss     ( scale , shift , "b" , "shift" , typeid ( *this ) )
{}
// ============================================================================
// evaluate GenBeta2-distribution
// ============================================================================
double Ostap::Math::GenBeta2::evaluate
( const double x ) const 
{
  if ( x < xmin ()  ) { return 0 ; }
  // 
  const double  z = m_ss.t ( x ) ;
  const double _a = a () ;
  //
  if      ( z <   0                 ) { return 0 ; }
  else if ( 1 <= _a && s_zero ( z ) ) { return 0 ; } 
  //
  const double _p = p () ;
  const double _q = q () ;
  const double _b = b () ;
  //
  /// case 1:   1 <= z**a
  if ( ( 1.0 <= z && 0 < _a ) || ( 1.0 >= z && 0 > _a ) )
  {
    const long double zia = std::pow ( z * 1.0L   , -1.0L * _a ) ;
    const long double c1  = std::pow ( 1.0L + zia , _p + _q    ) ; 
    const long double c2  = std::pow ( zia        ,      _q    ) ; 
    //
    return std::abs ( _a ) * m_pq.inv_Beta_pq () * c1 * c2 / ( z * _b ) ;
  }
  //
  const long double za = std::pow ( z * 1.0L            , _a      ) ;
  const long double c1 = std::pow ( za                  , _p      ) ;
  const long double c2 = std::pow ( 1.0 / ( za + 1.0L ) , _p + _q ) ; 
  //
  return std::abs ( _a ) * m_pq.inv_Beta_pq () * c1 * c2 / ( z * _b ) ;
  //
}
// ============================================================================
// raw moment  of unscaled/unshifted distribution
// ============================================================================
double Ostap::Math::GenBeta2::raw_moment
( const unsigned short k ) const
{
  const double _a = a () ;
  const double _p = p () ;
  const double _q = q () ;
  //
  const double _pa = _p + k / _a ;
  const double _qa = _q - k / _a ;
  //
  return
    _pa <= 0 || s_zero ( _pa       ) ? s_QUIETNAN : 
    _qa <= 0 || s_zero ( _qa       ) ? s_QUIETNAN :
    Ostap::Math::beta  ( _pa , _qa ) * m_pq.inv_Beta_pq () ; 
}
// ============================================================================
// mean-value @attention it can be infinite/NaN!
// ============================================================================
double Ostap::Math::GenBeta2::mean () const
{
  const double m1 = raw_moment ( 1 ) ;
  if ( !std::isfinite ( m1 ) ) { return m1 ; }
  return shift () + scale () * m1 ;
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
  const double _scale = scale () ; 
  return _scale * _scale * Ostap::Math::variance ( m2 , m1 ) ;
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
  return Ostap::Math::skewness ( m3 , m2 , m1 ) ;
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
  //
  if       ( s_equal ( low , high ) ) { return 0 ; }
  else if  ( high <  low            ) { return - integral ( high , low ) ; }
  ///
  const double xmn = xmin () ;  
  if       ( high <= xmn ) { return 0                       ; }
  else if  ( low  <  xmn ) { return integral ( xmn , high ) ; }
  // 
  const double _a = a () ;
  const double _b = b () ;
  //
  // split the range at mean-value (if defined and finite ) 
  //
  { // check mean-value
    const double vmean = mean () ;
    if ( std::isfinite ( vmean ) && low < vmean && vmean < high )
    { return integral ( low , vmean ) + integral ( vmean , high ) ; }	 
  }
  //
  // check the sigma
  const double sigma = rms () ;
  const double width = std::isfinite ( sigma ) ? sigma : 2 * _b ;
  //
  static const std::array<double,3> s_splits { 3 , 6 , 15 } ;
  for ( auto s : s_splits )
  {
    const double split = xmn + s * width  ;
    if ( low < split && split  < high ) { return integral ( low , split ) + integral ( split , high ) ; }  
  }
  //
  const bool in_tail = xmn + s_splits.back () * width <= low ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenBeta2> s_integrator ;
  static const char s_message[] = "Integral(GenBeta2)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,                             // low & high edges
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
// get the CDF 
// ============================================================================
double Ostap::Math::GenBeta2::cdf 
( const double x    ) const
{
  const double xmn = xmin () ;
  return x <= xmn ? 0.0 : integral ( xmn , x ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenBeta2::tag () const 
{ 
  static const std::string s_name = "GenBeta2" ;
  return Ostap::Utils::hash_combiner ( s_name      ,
				       m_a .tag () ,
				       m_pq.tag () ,
				       m_ss.tag () ) ;
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
  : m_scale ( scale )
  , m_shift ( shift )
{}
// ============================================================================
bool Ostap::Math::Landau::setScale ( const double value ) 
{
  if ( s_equal ( value , m_scale  ) ) { return false ; }
  m_scale  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Landau::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift  ) ) { return false ; }
  m_shift  = value ;
  return true ;
}
// ============================================================================
// evaluate Landau-distributions 
// ============================================================================
double Ostap::Math::Landau::pdf ( const double x ) const 
{
  //
  const double y = ( x - m_shift ) / m_scale ;
  //  
  return gsl_ran_landau_pdf ( y ) / m_scale ;
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
  //
  const double y = ( x - m_shift ) / m_scale ;
  //
  return _dislan ( y ) ;
}
// ============================================================================
double Ostap::Math::Landau::integral ( const double low  , 
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
  return Ostap::Utils::hash_combiner ( s_name , m_scale , m_shift ) ;
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
  : m_scale (  std::abs ( scale ) ) 
  , m_shape (  std::abs ( shape ) ) 
  , m_shift (             shift   )
{}
// ============================================================================
// evaluate Weibull-distributions
// ============================================================================
double Ostap::Math::Weibull::pdf ( const double x ) const 
{
  if ( x <= m_shift ) {  return 0 ; }
  const double y = ( x - m_shift ) / m_scale ;
  return ( m_shape / m_scale ) * 
    std::pow (  y , m_shape - 1 ) *
    std::exp ( - std::pow ( y , m_shape ) ) ;  
}
// ============================================================================
bool Ostap::Math::Weibull::setScale  ( const double value ) 
{
  const double v = std::abs ( value ) ;  
  if ( s_equal ( v , m_scale  ) ) { return false ; }
  m_scale  = v ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Weibull::setShape  ( const double value ) 
{
  const double v = std::abs ( value ) ;  
  if ( s_equal ( v , m_shape  ) ) { return false ; }
  m_shape  = v ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Weibull::setShift  ( const double value ) 
{
  if ( s_equal ( value , m_shift  ) ) { return false ; }
  m_shift  = value ;
  return true ;
}
// ============================================================================
double Ostap::Math::Weibull::cdf ( const double x ) const 
{
  if ( x <= m_shift ) {  return 0 ; }
  const double y = ( x - m_shift ) / m_scale ;
  return 1 - std::exp ( - std::pow ( y  , m_shape ) ) ;
}
// ============================================================================
double Ostap::Math::Weibull::integral 
( const double low  , const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  if ( high < low  ) { return -integral ( high, low ) ; }
  //
  return high <= m_shift ? 0.0 : cdf ( high ) - cdf ( low ) ;  
}
// ============================================================================
// the mean 
// ============================================================================
double Ostap::Math::Weibull::mean () const 
{ return m_shift + m_scale * std::tgamma ( 1  + 1/m_shape ) ; }
// ============================================================================
// the mode 
// ============================================================================
double Ostap::Math::Weibull::mode () const 
{  
  return 
    1 < m_shape ? m_shift : 
    m_shift + m_scale * std::pow ( ( m_shape - 1 ) / m_shape , 1/m_shape ) ;
}
// ============================================================================
// the median 
// ============================================================================
double Ostap::Math::Weibull::median () const 
{ return m_shift + m_scale * std::pow ( std::log ( 2.0 ) , 1/m_shape ) ; }
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::Weibull::variance () const 
{ 
  const double g1 = std::tgamma ( 1 + 2. / m_shape ) ;
  const double g2 = std::tgamma ( 1 + 1. / m_shape ) ;  
  return m_scale * m_scale * ( g1  - g2  * g2 ) ;
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
  return Ostap::Utils::hash_combiner ( s_name , m_scale , m_shape , m_shift ) ;
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
  : m_theta ( std::abs ( theta ) ) 
  , m_eta   ( std::abs ( eta   ) ) 
  , m_p     ( p     )
  , m_shift ( shift )
  , m_iKp_theta ( -1 ) 
  , m_workspace () 
{
  m_iKp_theta = 1.0 / Ostap::Math::bessel_Knu ( m_p , m_theta ) ;
}
// ============================================================================
/// evaluate the function 
// ============================================================================
double Ostap::Math::GenInvGauss::evaluate ( const double x ) const 
{
  if ( x <= m_shift ) { return 0 ; }
  //
  const double dx = x - m_shift ;
  const double xp = dx / m_eta  ;
  //
  return 0.5 * m_iKp_theta / m_eta 
    * std::pow ( xp , m_p - 1 )  
    * std::exp ( -0.5 * m_theta * ( xp + 1 / xp ) )  ; 
}
// ============================================================================
// parameter a  
// ============================================================================
double Ostap::Math::GenInvGauss::a () const { return m_theta / m_eta ; }
// ============================================================================
// parameter b  
// ============================================================================
double Ostap::Math::GenInvGauss::b () const { return m_theta * m_eta ; }
// ============================================================================
// set parameter theta 
// ============================================================================
bool Ostap::Math::GenInvGauss::setTheta  ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_theta , avalue ) ) { return false ; }
  m_theta     = avalue ;
  m_iKp_theta = 1.0 / Ostap::Math::bessel_Knu ( m_p , m_theta ) ;  
  return true ;
}
// ============================================================================
// set parameter eta 
// ============================================================================
bool Ostap::Math::GenInvGauss::setEta  ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_eta , avalue ) ) { return false ; }
  m_eta     = avalue ;
  return true ;
}
// ============================================================================
// set parameter P
// ============================================================================
bool Ostap::Math::GenInvGauss::setP  ( const double value ) 
{
  if ( s_equal ( m_p , value ) ) { return false ; }
  m_p         = value ;
  m_iKp_theta = 1.0 / Ostap::Math::bessel_Knu ( m_p , m_theta ) ;
  return true ;
}
// ============================================================================
// set shift parameter 
// ============================================================================
bool Ostap::Math::GenInvGauss::setShift  ( const double value ) 
{
  if ( s_equal ( m_shift , value ) ) { return false ; }
  m_shift = value ;
  return true ;
}
// ============================================================================
// mean value
// ============================================================================
double Ostap::Math::GenInvGauss::mean     () const 
{ return m_shift + m_eta * Ostap::Math::bessel_Knu ( m_p + 1 , m_theta ) * m_iKp_theta ; }
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::GenInvGauss::variance () const 
{
  const double k1 = Ostap::Math::bessel_Knu ( m_p + 1 , m_theta ) ;
  const double k2 = Ostap::Math::bessel_Knu ( m_p + 2 , m_theta ) ;
  return m_eta * m_eta * ( k2 * m_iKp_theta - std::pow ( k1 * m_iKp_theta , 2 ) ) ;
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
{ return m_shift + m_eta * ( ( m_p - 1 ) +  std::hypot ( m_p - 1 , m_theta ) ) / m_theta ; }
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
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low )        { return - integral ( high , low ) ; }
  //
  if ( high <= m_shift ) { return 0 ; }
  //
  const double xmin = std::max ( low , m_shift ) ;
  //
  const double x0 = m_eta + m_shift ;
  //
  // split at m_eta  
  if ( xmin < x0 && x0 < high ) 
  { return integral ( xmin , x0 ) + integral ( x0 , high ) ; }
  //
  const double sigma = rms () ;
  //
  // split at x1 
  const double x1 = x0 - 3 * sigma ;
  if ( xmin < x1 && x1 < high ) 
  { return integral ( xmin , x1 ) + integral ( x1 , high ) ; }
  //
  // split at x2 
  const double x2 = x0 + 4 * sigma ;
  if ( xmin < x2 && x2 < high ) 
  { return integral ( xmin , x2 ) + integral ( x2 , high ) ; }
  //
  // split at x3 
  const double x3 = x0 + 10 * sigma ;
  if ( xmin < x3 && x3 < high ) 
  { return integral ( xmin , x3 ) + integral ( x3 , high ) ; }
  //
  //
  const bool in_tail = ( x0 + 10 * sigma ) <= low ;
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
  return Ostap::Utils::hash_combiner ( s_name , m_theta , m_eta , m_p , m_shift ) ;
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
// Benini distribution 
// ============================================================================
Ostap::Math::Benini::Benini
( const std::vector<double>& pars  ,   // shape parameters               
  const double               scale ,   // scale parameter 
  const double               shift )   // shift parameter 
  : m_pars  ( pars.size() ) 
  , m_pars2 ( pars.size() ) 
  , m_scale ( std::abs ( scale ) ) 
  , m_shift (            shift   ) 
{
  setPars ( pars ) ;
  //
  Ostap::Assert ( 2 <= m_pars.size() , 
                  "Invalid number of parameters/1!" , 
                  "Ostap::Math::Benini" ) ;
  std::transform ( m_pars.begin () ,
		   m_pars.end   () ,
		   m_pars.begin () ,
		   []( const double v ) -> double { return std::abs ( v ) ;} ) ;
  Ostap::Assert ( 0 < std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ,
                  "Invalid number of parameters!" , 
                  "Ostap::Math::Benini" ) ;
}
// ============================================================================
// only two shape parameters: alpha and beta 
// ============================================================================
Ostap::Math::Benini::Benini
( const double               scale ,             // scale parameter 
  const double               shift )             // shift parameter
  : Benini ( std::vector<double>( 2 , 1.0 ) , scale , shift ) 
{}
// ============================================================================
// n parameteters 
// ============================================================================
Ostap::Math::Benini::Benini
( const unsigned short       n     ,             // shape parameters 
  const double               scale ,             // scale parameter 
  const double               shift )             // shift parameter
  : Benini ( std::vector<double> ( n , 1.0 ) , scale , shift ) 
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
  const double               scale     ,              // scale parameter 
  const double               shift     )             // shift parameter
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
  const double               scale     ,              // scale parameter 
  const double               shift     )             // shift parameter
  : Benini ( {{ alpha , beta , gamma , delta }} , scale , shift )
{}    
// ============================================================================
bool Ostap::Math::Benini::_setPar  
( const unsigned short i     , 
  const double         value ) 
{
  if ( m_pars.size() <= i ) { return false ; }
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_pars [ i ] , avalue ) ) { return false ; }
  //
  m_pars  [ i ] = avalue ;
  m_pars2 [ i ] = ( i + 1 ) *  avalue ;
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
  return changed ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Benini::evaluate 
( const double x ) const 
{
  //
  if ( x < m_shift + m_scale  ) { return 0 ; }
  //
  const double lxs = std::log ( ( x - m_shift ) / m_scale ) ;
  //
  const double p1 = Ostap::Math::Clenshaw::monomial_sum ( m_pars .rbegin() , m_pars .rend  () , lxs ).first ;
  const double p2 = Ostap::Math::Clenshaw::monomial_sum ( m_pars2.rbegin() , m_pars2.rend  () , lxs ).first ;
  //
  return std::exp ( - lxs * p1 ) * p2 / x ;
}
// =============================================================================
// set scale parameter
// =============================================================================
bool Ostap::Math::Benini::setScale ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_scale , avalue ) ) { return false ; }
  m_scale     = avalue ;
  return true ;    
}
// =============================================================================
// set shift  parameter
// =============================================================================
bool Ostap::Math::Benini::setShift ( const double value ) 
{
  if ( s_equal ( m_shift , value ) ) { return false ; }
  m_shift     = value ;
  return true ;    
}
// ============================================================================
// evaluate the cdf 
// ============================================================================
double Ostap::Math::Benini::cdf
( const double x ) const 
{
  if  ( x <= m_shift + m_scale ) { return 0 ; }
  //
  const double lxs = std::log ( ( x - m_shift ) / m_scale ) ;
  //
  const double p1 = Ostap::Math::Clenshaw::monomial_sum ( m_pars.rbegin() , m_pars.rend  () , lxs ).first ;
  //
  return 1.0L - std::exp ( -lxs * p1 ) ;
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
{ return 
    std::max ( low , high ) <= m_shift + m_scale ? 0.0 :
    s_equal ( low , high )                       ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Benini::tag () const 
{ 
  static const std::string s_name = "Benini" ;
  return Ostap::Utils::hash_combiner ( s_name  , 
                                       m_scale , 
                                       m_shift , 
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
  : m_alpha ( std::abs ( alpha ) )
  , m_scale ( std::abs ( scale ) )
  , m_shift  ( shift )
{
  Ostap::Assert ( m_alpha ,
		  "Alpha parameter must be non-zero!"     ,
		  "Ostap::Math::Frechet"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( m_scale                                 ,
		  "Scale parameter must be non-zero!"     ,
		  "Ostap::Math::Frechet"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// evaluate Frechet distribution
// ============================================================================
double Ostap::Math::Frechet::evaluate
( const double x ) const
{
  if ( x < m_shift ) { return  0 ; }
  //
  const double delta = ( x - m_shift ) / m_scale ;
  return
    m_alpha * std::pow ( delta , -1 - m_alpha ) *
    std::exp ( - std::pow ( delta , -m_alpha ) ) / m_scale ;
}
// ============================================================================
// Set parameter alpha
// ============================================================================
bool Ostap::Math::Frechet::setAlpha
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_alpha ) ) { return false ; }
  Ostap::Assert ( avalue ,
		  "alpha parameter must be non-zero!"     ,
		  "Ostap::Math::Frechet"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_alpha = avalue ;
  return true ;
}
// ============================================================================
// Set scale 
// ============================================================================
bool Ostap::Math::Frechet::setScale
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( avalue ,
		  "Scale parameter must be non-zero!"     ,
		  "Ostap::Math::Frechet"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_scale = avalue ;
  return true ;
}
// ============================================================================
// Set shift 
// ============================================================================
bool Ostap::Math::Frechet::setShift
( const double value )
{
  if ( s_equal ( value , m_shift ) ) { return false ; }
  m_shift = value ;
  return true ;
}
// ============================================================================
// mean value
// ============================================================================
double Ostap::Math::Frechet::mean () const
{
  return
    m_alpha <= 1 || s_equal ( m_alpha , 1 ) ? s_INFINITY : 
    m_shift + m_scale * Ostap::Math::gamma ( 1.0 - 1.0 / m_alpha ) ;
}
// ============================================================================
// median value
// ============================================================================
double Ostap::Math::Frechet::median () const
{
  static const double s_ln2  { std::log ( 2.0 ) } ;
  return m_shift + m_scale * std::pow ( s_ln2 , - 1 / m_alpha ) ;
}
// ============================================================================
// mode value
// ============================================================================
double Ostap::Math::Frechet::mode  () const
{ return m_shift + m_scale * std::pow ( m_alpha / ( 1 + m_alpha ) , 1/m_alpha ) ; }
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::Frechet::variance  () const
{
  if ( m_alpha <= 2 || s_equal ( m_alpha , 2 ) ) { return s_INFINITY ; }
  //
  const double g1 = Ostap::Math::gamma ( 1.0 - 1.0 / m_alpha ) ;
  const double g2 = Ostap::Math::gamma ( 1.0 - 2.0 / m_alpha ) ;
  //
  return m_scale * m_scale  * ( g2  - g1 * g1 ) ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::Frechet::rms  () const
{
  if ( m_alpha <= 2 || s_equal ( m_alpha , 2 ) ) { return s_INFINITY ; }
  //
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// skewnes 
// ============================================================================
double Ostap::Math::Frechet::skewness() const
{
  if ( m_alpha <= 3 || s_equal ( m_alpha , 3 ) ) { return s_INFINITY ; }
  //
  const double g1 = Ostap::Math::gamma ( 1.0 - 1.0 / m_alpha ) ;
  const double g2 = Ostap::Math::gamma ( 1.0 - 2.0 / m_alpha ) ;
  const double g3 = Ostap::Math::gamma ( 1.0 - 3.0 / m_alpha ) ;
  //
  const double A = g3 - 3 * g2 * g1 + 2 * g1 * g1 * g1 ;
  const double B = g2 -     g1 * g1 ;
  //
  return A / std::sqrt ( B * B * B ) ;
}
// ============================================================================
// (excess) kurtosis
// ============================================================================
double Ostap::Math::Frechet::kurtosis () const
{
  if ( m_alpha <= 4 || s_equal ( m_alpha , 4 ) ) { return s_INFINITY ; }
  //
  const double g1 = Ostap::Math::gamma ( 1.0 - 1.0 / m_alpha ) ;
  const double g2 = Ostap::Math::gamma ( 1.0 - 2.0 / m_alpha ) ;
  const double g3 = Ostap::Math::gamma ( 1.0 - 3.0 / m_alpha ) ;
  const double g4 = Ostap::Math::gamma ( 1.0 - 4.0 / m_alpha ) ;
  //
  const double C = g4 - 4 * g3 * g1 + 3 * g2 * g2 ;
  const double D = g2 -     g1 * g1 ;
  //
  return -6 + C / ( D * D ) ; 
}
// ============================================================================
// evaluate Frechet CDF 
// ============================================================================
double Ostap::Math::Frechet::cdf 
( const double x ) const
{
  if ( x <= m_shift ) { return  0 ; }
  // 
  const double delta = ( x - m_shift ) / m_scale ;
  return std::exp ( - std::pow ( delta , -m_alpha ) ) ;
  
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
{ return cdf ( b ) - cdf ( a ) ; }
// ============================================================================
/*  get quantile 
 *  @param p probability \f$ 0 \le p < 1 \f$
 */
// ============================================================================
double Ostap::Math::Frechet::quantile
( const double p ) const
{
  return
    s_zero  ( p     ) ?  m_shift                                   :              
    s_equal ( p , 1 ) ?  std::numeric_limits<double>::max       () : 
    p <= 0            ?  std::numeric_limits<double>::quiet_NaN () :
    p >= 1            ?  std::numeric_limits<double>::quiet_NaN () :
    //
    m_shift + m_scale * std::pow ( -std::log ( p ) , -1 / m_alpha ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Frechet::tag () const 
{ 
  static const std::string s_name = "Frechet" ;
  return Ostap::Utils::hash_combiner ( s_name , m_alpha , m_scale , m_shift ) ;
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
  const double b     ,
  const double shift )
  : m_p     ( std::abs ( p ) )
  , m_a     ( std::abs ( a ) )
  , m_b     ( std::abs ( b ) )
  , m_shift ( shift )
{
  Ostap::Assert ( 0 < m_p , 
		  "p-parameter must be positive!"         ,
		  "Ostap::Math::Dagum"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 0 < m_a , 
		  "a-parameter must be positive!"         ,
		  "Ostap::Math::Dagum"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 0 < m_b , 
		  "b-parameter must be positive!"         ,
		  "Ostap::Math::Dagum"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// Set shape parameter p 
// ============================================================================
bool Ostap::Math::Dagum::setP
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_p ) ) { return false ; }
  //
  Ostap::Assert ( 0 < avalue ,
		  "p-parameter must be positive!"         ,
		  "Ostap::Math::Dagum"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_p = avalue ;
  return true ;
}
// ============================================================================
// Set shape parameter a 
// ============================================================================
bool Ostap::Math::Dagum::setA
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_a ) ) { return false ; }
  //
  Ostap::Assert ( 0 < avalue ,
		  "a-parameter must be positive!"         ,
		  "Ostap::Math::Dagum"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_a = avalue ;
  return true ;
}
// ============================================================================
// Set scale parameter b
// ============================================================================
bool Ostap::Math::Dagum::setB
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_b ) ) { return false ; }
  //
  Ostap::Assert ( 0 < avalue ,
		  "b-parameter must be positive!"         ,
		  "Ostap::Math::Dagum"                    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_b = avalue ;
  return true ;
}
// ============================================================================
// Set shift 
// ============================================================================
bool Ostap::Math::Dagum::setShift
( const double value )
{
  if ( s_equal ( value , m_shift ) ) { return false ; }
  m_shift = value ;
  return true ;
}
// ============================================================================
// evaluate Dagum distribution
// ============================================================================
double Ostap::Math::Dagum::evaluate
( const double x ) const
{
  if ( x <= m_shift ) { return 0 ; } // RETURN 
  //
  const double d  = ( x - m_shift ) / m_b ;
  if ( d <= 0     ) { return 0 ; } // RETURN
  //
  if ( 1 <= d && 1 <= m_a )
  {
    const double t1  = std::pow ( d , -m_a ) ;
    return m_a * m_p * std::pow ( 1 / ( 1 + t1 ) , m_p ) * t1 / ( 1 + t1 ) / d ; 
  }
  //
  const double t1  = std::pow ( d , m_a  ) ;
  return m_a * m_p * std::pow ( t1 / ( 1 + t1 ) , m_p ) * 1 / ( 1 + t1 ) / d ; 
}  
// ============================================================================
// evaluate Dagum CDF 
// ============================================================================
double Ostap::Math::Dagum::cdf
( const double x    ) const
{
  if ( x <= m_shift ) { return 0 ; }
  //
  const double d = ( x - m_shift ) / m_b ;
  if ( d <= 0      ) { return 0 ; } // RETURN
  //
  if ( 1 <= d && 1 <= m_a )
  {
    const double t1 = std::pow ( d , -m_a ) ;
    return std::pow ( 1 / ( t1 + 1 ) , m_p ) ;
  }
  const double t1 = std::pow ( d , m_a ) ;
  return std::pow ( t1 / ( 1 + t1 ) , m_p ) ;
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
  else if ( high <= m_shift        ) { return 0 ; }
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
    u < 0             ?  s_QUIETNAN :
    u > 1             ?  s_QUIETNAN : 
    s_zero  ( u     ) ?  m_shift    : 
    s_equal ( u , 1 ) ?  s_POSHUGE  :
    m_shift + m_b * std::pow ( std::pow ( u , -1/m_p ) -1 , -1 / m_a ) ;
 }
// ============================================================================
// mean 
// ============================================================================
double Ostap::Math::Dagum::mean() const
{
  if ( m_a <= 1 || s_equal ( m_a , 1 ) ) { return s_QUIETNAN ; }
  //
  const double lg1 = Ostap::Math::lgamma ( 1   - 1 / m_a ) ;
  const double lg2 = Ostap::Math::lgamma ( m_p + 1 / m_a ) ;
  const double lgp = Ostap::Math::lgamma ( m_p ) ;
  //
  return m_shift + m_b * std::exp ( lg1 + lg2 - lgp ) ;  
}
// ============================================================================
// mode
// ============================================================================
double Ostap::Math::Dagum::mode () const
{
  const double ap = m_a * m_p ;
  if ( ap <= 1 ) { return m_shift ; }
  return m_shift + m_b * std::pow ( ( ap - 1 ) / ( m_a + 1 ) , 1/m_a ) ;
}
// ============================================================================
// median
// ============================================================================
double Ostap::Math::Dagum::median () const
{
  const double q = 1 / m_p ;
  if ( 1 < q )
  {
    const double q2 = std::pow ( 2.0 , -q ) ;
    return m_shift + m_b * std::pow ( q2 / ( 1.0L - q2 ) , 1/m_a ) ;
  }
  const double q2 = std::pow ( 2.0 , q ) ;
  return m_shift + m_b * std::pow ( 1 / ( q2 - 1.0L )  , 1/m_a ) ;
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::Dagum::variance () const
{
  if ( m_a <= 2 || s_equal ( m_a , 2 ) ) { return s_QUIETNAN ; }
  //
  const double lg1 = Ostap::Math::lgamma ( 1   - 2 / m_a ) ;
  const double lg2 = Ostap::Math::lgamma ( m_p + 2 / m_a ) ;
  const double lg3 = Ostap::Math::lgamma ( 1   - 1 / m_a ) ;
  const double lg4 = Ostap::Math::lgamma ( m_p + 1 / m_a ) ;
  const double lgp = Ostap::Math::lgamma ( m_p ) ;
  //
  const double t1 = std::exp ( lg1 + lg2 - lgp ) ;
  const double t2 = std::exp ( lg3 + lg4 - lgp ) ;
  //
  return m_b * m_b * ( t1 - t2 * t2 ) ;
}
// ============================================================================
// rms 
// ============================================================================
double Ostap::Math::Dagum::rms  () const
{
  if ( m_a <= 2 || s_equal ( m_a , 2 ) ) { return s_QUIETNAN ; }
  //
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Dagum::tag () const 
{ 
  static const std::string s_name = "Dagum" ;
  return Ostap::Utils::hash_combiner ( s_name , m_p , m_a , m_b , m_shift ) ;
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
  : m_a         ( std::abs ( a )     )
  , m_r         ( std::abs ( r )     )
  , m_scale     ( std::abs ( scale ) )
  , m_shift     ( shift )
  , m_p         ( -1 )
{
  //
  Ostap::Assert ( 0 < m_a ,
		  "a-parameter must be positive!"         ,
		  "Ostap::Math::BenktanderI"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( 0 < m_scale ,
		  "scale-parameter must be positive!"     ,
		  "Ostap::Math::BenktanderI"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  setR ( r ) ;
  //
  Ostap::Assert ( 0 < m_p && m_p <= 1 ,
		  "p-parameter must be positive and <=1"  ,
		  "Ostap::Math::BenktanderI"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
}
// ============================================================================
// Set shape parameter a 
// ============================================================================
bool Ostap::Math::BenktanderI::setA
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_a ) ) { return false ; }
  //
  Ostap::Assert ( 0 < avalue ,
		  "a-parameter must be positive!"         ,
		  "Ostap::Math::BenktanderI"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_a = avalue ;
  return true ;
}    
// ============================================================================
// Set delta parameter  
// ============================================================================
bool Ostap::Math::BenktanderI::setR 
( const double value )
{
  const double avalue = std::abs ( value ) ;  
  if ( s_equal ( avalue , m_r ) && 0 < m_p && m_p <= 1 ) { return false ; }
  //
  m_r = avalue ;
  m_p = 1.0 / std::hypot ( m_r , 1.0 ) ;
  //
  Ostap::Assert ( 0 < m_p && m_p <= 1 ,
		  "p-parameter must be positive and <=1"  ,
		  "Ostap::Math::BenktanderI"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
  //
  return true ;
}
// ============================================================================
// Set shift 
// ============================================================================
bool Ostap::Math::BenktanderI::setShift
( const double value )
{
  if ( s_equal ( value , m_shift ) ) { return false ; }
  m_shift = value ;
  return true ;
}
// ============================================================================
// Set scale
// ============================================================================
bool Ostap::Math::BenktanderI::setScale
( const double value )
{
  const double avalue = std::abs ( value ) ; 
  if ( s_equal ( avalue , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( avalue ,
		  "scale-parameter must be positive!"     ,
		  "Ostap::Math::BenktanderI"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_scale = avalue ;
  return true ;
}
// ============================================================================
// evaluate Benktander Type I distribution
// ============================================================================
double Ostap::Math::BenktanderI::evaluate
( const double x ) const
{
  if ( x < m_shift ) { return 0 ; }
  //
  const double z  =  ( x - m_shift ) / m_scale + 1.0 ;
  const double lz = std::log ( z ) ;
  //
  const double _b   = b () ;
  const double t1 = 1 +       2 * _b * lz / m_a ;
  const double t2 = 1 + m_a + 2 * _b * lz       ;
  const double t3 = 2 + m_a +     _b * lz       ;
  //
  const double rpdf =  ( t1 * t2 - 2 * _b / m_a ) * std::pow ( z , -t3 ) ;
  //
  return rpdf / m_scale ;
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
  else if ( high <= m_shift        ) { return 0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// evaluate Benktander Type I CDF 
// ============================================================================
double Ostap::Math::BenktanderI::cdf 
( const double x ) const
{
  if ( x < m_shift ) { return 0 ; }
  //
  const double z  =  ( x - m_shift ) / m_scale + 1.0 ;  
  const double lz = std::log ( z ) ;
  //
  const double _b = b () ;
  const double t1 = 1 +       2 * _b * lz / m_a ;
  const double t2 = 1 + m_a +     _b * lz       ;
  //
  return 1 - t1 * std::pow ( x , -t2 ) ;
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::BenktanderI::mean     () const
{
  const double zmean = 1 + 1 / m_a ;
  return m_shift + m_scale * ( zmean - 1 ) ;
}
// ============================================================================
// variance value 
// ============================================================================
double Ostap::Math::BenktanderI::variance () const
{
  const double _b  = b () ;
  const double sqb = std::sqrt ( _b ) ;
  const double t1  = ( m_a - 1 ) / ( 2 * sqb ) ;
  //
  const double t2  = -sqb + m_a * s_sqrt_pi * Ostap::Math::erfcx ( t1 ) ;
  return m_scale * m_scale * t2 / ( m_a * m_a * sqb ) ;
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
				       m_a     ,
				       m_r     ,
				       m_scale ,
				       m_shift ) ;
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
  : m_a         ( std::abs ( a )     )
  , m_r         ( std::abs ( r )     )
  , m_scale     ( std::abs ( scale ) )
  , m_shift     ( shift )
  , m_b         ( -1 ) 
{
  
  //
  Ostap::Assert ( 0 < m_a ,
		  "a-parameter must be positive!"         ,
		  "Ostap::Math::BenktanderII"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( 0 < m_scale ,
		  "scale-parameter must be positive!"     ,
		  "Ostap::Math::BenktanderII"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  setR ( r ) ;
  //
  Ostap::Assert ( 0 < m_b  && m_b <= 1 		          , 
		  "b-parameter must be positive!"         ,
		  "Ostap::Math::BenktanderII"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
}
// ============================================================================
// Set shape parameter a 
// ============================================================================
bool Ostap::Math::BenktanderII::setA
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_a ) ) { return false ; }
  //
  Ostap::Assert ( 0 < m_a  ,
		  "a-parameter must be positive!"         ,
		  "Ostap::Math::BenktanderII"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_a = avalue ;
  return true ;
}    
// ============================================================================
// Set R parameter  
// ============================================================================
bool Ostap::Math::BenktanderII::setR 
( const double value )
{
  const double avalue = std::abs ( value ) ;  
  if ( s_equal ( avalue , m_r    ) && 0 < m_b && m_b <= 1 ) { return false ; }
  //
  m_r = avalue ;
  m_b = 1.0 / std::hypot ( m_r , 1.0 ) ;
  //
  Ostap::Assert ( 0 < m_b  && m_b <= 1 ,
		  "b-parameter must be positive!"         ,
		  "Ostap::Math::BenktanderII"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  return true ;
}
// ============================================================================
// Set shift 
// ============================================================================
bool Ostap::Math::BenktanderII::setShift
( const double value )
{
  if ( s_equal ( value , m_shift ) ) { return false ; }
  m_shift = value ;
  return true ;
}
// ============================================================================
// Set scale
// ============================================================================
bool Ostap::Math::BenktanderII::setScale
( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( avalue ,
		  "scale-parameter must be positive!"     ,
		  "Ostap::Math::BenktanderII"             ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_scale = avalue ;
  return true ;
}
// ============================================================================
// evaluate Benktander Type II distribution
// ============================================================================
double Ostap::Math::BenktanderII::evaluate
( const double x ) const
{
  if ( x < m_shift ) { return 0 ; }
  //
  const double z  =  ( x - m_shift ) / m_scale + 1.0 ;
  const double zb = std::pow ( z , m_b ) ;
  //
  const double zpdf =
    std::exp ( m_a * ( 1 - zb ) / m_b ) *
    zb * ( m_a * zb - m_b + 1 ) / ( z * z ) ;
  //
  return zpdf / m_scale ;
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
  else if ( high <= m_shift        ) { return 0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// evaluate Benktander Type II CDF 
// ============================================================================
double Ostap::Math::BenktanderII::cdf 
( const double x ) const
{
  if ( x <= m_shift ) { return 0 ; }
  //
  const double z    =  ( x - m_shift ) / m_scale + 1.0 ;
  //
  const double zb   = std::pow ( z , m_b ) ;
  const double zcdf = 1 - zb * std::exp ( m_a * ( 1 - zb ) / m_b ) / z ; 
  //
  return zcdf ;
}
// ============================================================================
// mode 
// ============================================================================
double Ostap::Math::BenktanderII::mode     () const
{ return m_shift ;  }
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BenktanderII::tag () const 
{ 
  static const std::string s_name = "BenktanderII" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_a     ,
				       m_r     ,
				       m_scale ,
				       m_shift ) ;
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
  : m_alpha ( std::abs ( alpha ) )
  , m_beta  ( std::abs ( beta  ) )
  , m_scale ( std::abs ( scale ) )
  , m_shift (            shift   )
{
  Ostap::Assert ( 0 < m_alpha , 
                  "alpha-parameter must be positive!"     ,
                  "Ostap::Math::Kumaraswami"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 0 < m_beta  , 
                  "beta-parameter must be positive!"      ,
                  "Ostap::Math::Kumaraswami"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 0 < m_scale  , 
                  "Scale parameter must be positive!"     ,
                  "Ostap::Math::Kumaraswami"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ======================================================================
// set a 
// ============================================================================
bool Ostap::Math::Kumaraswami::setAlpha ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_alpha ) ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "alpha-parameter must be positive!"     ,
                  "Ostap::Math::Kumaraswami"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_alpha = avalue ; 
  return true ;
} 
// ======================================================================
// set beta
// ============================================================================
bool Ostap::Math::Kumaraswami::setBeta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_beta ) ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "beta-parameter must be positive!"      ,
                  "Ostap::Math::Kumaraswami"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_beta = avalue ; 
  return true ;
}
// ======================================================================
// set scale
// ============================================================================
bool Ostap::Math::Kumaraswami::setScale ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( avalue , 
                  "Scale parameter must be positive!"         ,
                  "Ostap::Math::Kumaraswami"              ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_scale = avalue ; 
  return true ;
}
// ======================================================================
// set shift
// ============================================================================
bool Ostap::Math::Kumaraswami::setShift ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_shift ) ) { return false ; }
  //
  m_shift = avalue ; 
  return true ;
}
// ======================================================================
// evaluate Kumaraswami function
// ======================================================================
double Ostap::Math::Kumaraswami::evaluate    ( const double x ) const
{
  if ( x <= m_shift || m_scale + m_shift <= x ) { return 0 ; }
  const double z = ( x - m_shift ) / m_scale ;
  //
  const long double za = std::pow ( z * 1.0L , m_alpha * 1.0L ) ;
  return m_alpha * m_beta * za * std::pow ( 1.0L - za , m_beta - 1.0L ) / ( z * m_scale ) ;
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
  else if ( high    <= xmin ()      ) { return  0 ; }
  else if ( low     >= xmax ()      ) { return  0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ======================================================================
// CDF 
// ======================================================================
double Ostap::Math::Kumaraswami::cdf 
( const double x ) const
{
  if      ( x <= xmin () ) { return 0 ; }
  else if ( x >= xmax () ) { return 1 ; }
  //
  const double z = ( x - m_shift ) / m_scale ;
  //
  const long double za = std::pow ( z * 1.0L , m_alpha * 1.0L ) ;
  return 1.0 - std::pow ( 1.0L- za , m_beta * 1.0L ) ;
}
// ======================================================================
// quantile  function \f$ 0 < p < 1 \f$ 
// ======================================================================
double Ostap::Math::Kumaraswami::quantile 
( const double p    ) const
{
  return
    p < 0             ?  s_QUIETNAN :
    p > 1             ?  s_QUIETNAN : 
    s_zero  ( p     ) ?  xmin ()    : 
    s_equal ( p , 1 ) ?  xmax ()    :
    m_shift + m_scale * std::pow ( 1.0L - std::pow ( 1.0L - p , 1.0L/m_beta ) , 1.0L/m_alpha ) ;
}
// ============================================================================
// raw moment for standartized Kumaraswami
// ============================================================================
double Ostap::Math::Kumaraswami::moment
( const unsigned int n ) const
{ return m_beta * Ostap::Math::beta ( 1 + n / m_alpha , m_beta ) ; }
// ============================================================================
// mean-value 
// ============================================================================
double Ostap::Math::Kumaraswami::mean     () const
{ return m_shift + m_scale * moment ( 1 ) ; }
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::Kumaraswami::variance  () const
{
  const double m1 = moment ( 1 ) ;
  const double m2 = moment ( 1 ) ;
  //
  return m_scale * m_scale * ( m2 - m1 * m1 ) ; 
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
  const double m = std::pow ( 1.0L - std::pow ( 2.0L , -1.0L/m_beta ) , 1.0L/m_alpha ) ; 
  return m_shift + m_scale * m ;
}
// ============================================================================
// mode-value 
// ============================================================================
double Ostap::Math::Kumaraswami::mode () const
{
  if      ( m_alpha < 1 || m_beta < 1                         ) { return s_QUIETNAN ; }
  else if ( s_equal ( m_alpha , 1 ) && s_equal ( m_beta , 1 ) ) { return m_shift + 0.5 * m_scale ; }
  //  
  const double m = std::pow ( ( m_alpha - 1.0L ) / ( m_alpha * m_beta - 1.0L ) , 1.0L/m_alpha ) ;
  return m_shift + m_scale * m ;
}
// ============================================================================
// skewness-value 
// ============================================================================
double Ostap::Math::Kumaraswami::skewness () const
{
  const double m1 = moment ( 1 ) ;
  const double m2 = moment ( 2 ) ;
  const double m3 = moment ( 3 ) ;
  //
  const double s2 = m2 - m1 * m1 ;
  //
  return ( m3 - 3 * m2 * m1 + 2 * m1 * m1 * m1 ) / std::pow ( s2 , 1.5 ) ; 
}
// ============================================================================
// kurtosis-value 
// ============================================================================
double Ostap::Math::Kumaraswami::kurtosis() const
{
  const double m1 = moment ( 1 ) ;
  const double m2 = moment ( 2 ) ;
  const double m3 = moment ( 3 ) ;
  const double m4 = moment ( 4 ) ;
  //
  const double s2 = m2 - m1 * m1 ;
  //
  return ( m4 - 4 * m3 * m1 + 6 * m2 * m1 * m1 - 3 * m1 * m1 * m1 * m1  ) / ( s2 * s2 ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Kumaraswami::tag () const 
{ 
  static const std::string s_name = "Kumaraswami" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_alpha ,
				       m_beta  ,
				       m_scale ,
				       m_shift ) ; 
}

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
// ======================================================================
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
 *  @param c           shape parameter 
 *  @param k           shape parameter 
 *  @param scale scale parameter 
 *  @param shift shift parameter
 */
// ============================================================================
Ostap::Math::Burr::Burr
( const double c     ,
  const double k     ,
  const double scale ,
  const double shift )
  : m_c     ( std::abs ( c     ) )
  , m_k     ( std::abs ( k     ) )
  , m_scale ( std::abs ( scale ) )
  , m_shift (            shift   )
{
  Ostap::Assert ( m_c   , 
                  "c-parameter must be positive!"         ,
                  "Ostap::Math::Burr"                     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( m_k   , 
                  "k-parameter must be positive!"         ,
                  "Ostap::Math::Burr"                     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( m_scale , 
                  "scale-parameter must be positive!"     ,
                  "Ostap::Math::Burr"                     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
}
// ======================================================================
// set C 
// ============================================================================
bool Ostap::Math::Burr::setC ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_c ) ) { return false ; }
  //
  Ostap::Assert ( avalue   , 
                  "c-parameter must be positive!"         ,
                  "Ostap::Math::Burr"                     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_c = avalue ; 
  return true ;
}
// ======================================================================
// set K 
// ============================================================================
bool Ostap::Math::Burr::setK ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_k ) ) { return false ; }
  //
  Ostap::Assert ( avalue   , 
                  "k-parameter must be positive!"         ,
                  "Ostap::Math::Burr"                     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_k = avalue ; 
  return true ;
} 
// ======================================================================
// set scale 
// ============================================================================
bool Ostap::Math::Burr::setScale( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_scale ) ) { return false ; }
  //
  Ostap::Assert ( avalue   , 
                  "scale-parameter must be positive!"     ,
                  "Ostap::Math::Burr"                     ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_scale = avalue ; 
  return true ;
} 
// ======================================================================
// set shift
// ============================================================================
bool Ostap::Math::Burr::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift ) ) { return false ; }
  m_shift = value ; 
  return true ;
}
// ============================================================================
// evaluate Burr function
// ============================================================================
double Ostap::Math::Burr::evaluate    ( const double x ) const
{
  if ( x < m_shift ) { return 0 ; }
  double z = ( x - m_shift ) / m_scale ;
  if ( z < 0   ) { return 0 ; }
  ///
  if ( m_c < 1 ) { z = std::max ( z , s_NONZERO ) ; } // NB!!! 
  //
  if ( z <= 1 )
  {
    const double zcm1 = std::pow ( z , m_c - 1 ) ;
    return m_c * m_k * zcm1 / std::pow ( 1 + z * zcm1 , m_k + 1 ) / m_scale ;       
  }
  //
  const double zc = std::pow ( z , -m_c ) ;
  const double lr = ( m_c - 1 ) * std::log ( z ) -
    ( m_k + 1 ) * std::log ( zc / ( 1.0L + zc ) ) ;
  //
  return m_c * m_k * std::exp ( lr ) / m_scale ;    
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::Burr::integral ()    const { return 1 ; }
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::Burr::integral
( const double low  ,
  const double high )  const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if ( high <  low            ) { return -integral ( high    , low  ) ; }
  else if ( high <= m_shift        ) { return  0 ; }
  else if ( low  <  m_shift        ) { return  integral ( m_shift , high ) ; } 
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// CDF 
// ============================================================================
double Ostap::Math::Burr::cdf 
( const double x ) const
{
  if ( x <= m_shift ) { return 0 ; }
  const double z = ( x - m_shift ) / m_scale ;
  if ( z <= 0 || s_zero ( z ) ) { return 0 ; }
  //
  if ( z <= 1 )
  {
    const double zc = std::pow ( z , m_c ) ;
    return 1 - std::pow (1 + zc , - m_k ) ;
  }
  //
  const double zc = std::pow ( z , -m_c ) ;
  return 1 - std::pow ( zc / ( 1 + zc ) , m_k ) ;
}
// ============================================================================
// raw-moments of unscaled/unshifted distributions 
// ============================================================================
double Ostap::Math::Burr::raw_moment
( const unsigned short r ) const
{
  const double ck = m_c * m_k ;
  if ( ck <= r || s_equal ( ck , r ) ) { return s_QUIETNAN ; }
  //
  return m_k * Ostap::Math::beta ( ( ck - r ) / m_c , 1 + 1 / m_c ) ;
}
// ============================================================================
// mean              (for ck>1) 
// ============================================================================
double Ostap::Math::Burr::mean     () const
{
  const double ck = m_c * m_k ;
  if ( ck <= 1 || s_equal ( ck , 1 ) ) { return s_QUIETNAN ; }
  //
  const double m1 = raw_moment ( 1 ) ;
  return m_shift + m_scale * m1 ;
}
// ============================================================================
// mode
// ============================================================================
double Ostap::Math::Burr::mode () const
{
  if ( m_c <= 1 ) { return m_shift ; }
  const double ck = m_c * m_k ;
  return m_shift + m_scale * std::pow ( ( m_c - 1.0L ) / ( ck + 1.0L ) , 1/m_c ) ;
}
// ============================================================================
// median
// ============================================================================
double Ostap::Math::Burr::median () const
{ return m_shift + m_scale * std::pow ( std::pow ( 2.0L , 1.0L/m_k ) - 1.0 , 1.0/m_c ) ; }
// ============================================================================
// variance              (for ck>2) 
// ============================================================================
double Ostap::Math::Burr::variance () const
{
  //
  const double ck = m_c * m_k ;
  if ( ck <= 2 || s_equal ( ck , 2 ) ) { return s_QUIETNAN ; }
  //
  const double m1 = raw_moment ( 1 ) ;
  const double m2 = raw_moment ( 2 ) ;
  //
  return m_scale * m_scale * ( m2 - m1 * m1 ) ;
}
// ============================================================================
// rms              (for ck>2) 
// ============================================================================
double Ostap::Math::Burr::rms () const
{
  //
  const double ck = m_c * m_k ;
  if ( ck <= 2 || s_equal ( ck , 2 ) ) { return s_QUIETNAN ; }
  //
  return std::sqrt ( variance () ) ;
}
// ============================================================================
// skewness              (for ck>3) 
// ============================================================================
double Ostap::Math::Burr::skewness () const
{
  //
  const double ck = m_c * m_k ;
  if ( ck <= 3 || s_equal ( ck , 3 ) ) { return s_QUIETNAN ; }
  //
  const double m1 = raw_moment ( 1 ) ;
  const double m2 = raw_moment ( 2 ) ;
  const double m3 = raw_moment ( 3 ) ;
  //
  return ( m3 - 3 * m1 * m2 + 2 * m1 * m1 * m1 ) /
    std::pow ( m2 - m1 * m1 , 3.0/2 ) ;
}
// ============================================================================
// kurtosis              (for ck>4) 
// ============================================================================
double Ostap::Math::Burr::kurtosis () const
{
  //
  const double ck = m_c * m_k ;
  if ( ck <= 4 || s_equal ( ck , 4 ) ) { return s_QUIETNAN ; }
  //
  const double m1 = raw_moment ( 1 ) ;
  const double m2 = raw_moment ( 2 ) ;
  const double m3 = raw_moment ( 3 ) ;
  const double m4 = raw_moment ( 4 ) ;
  //  
  return ( m4 - 4 * m3 * m1 + 6 * m2 * m1 * m1 - 3 * m1 * m1 * m1 * m1 ) /
    std::pow ( m2 - m1 * m1 , 2 ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Burr::tag () const 
{ 
  static const std::string s_name = "Burr" ;
  return Ostap::Utils::hash_combiner ( s_name  ,
				       m_c     ,
				       m_k     ,
				       m_scale ,
				       m_shift ) ; 
}

// ============================================================================
//                                                                      The END
// ============================================================================


