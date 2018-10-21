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
#include "Ostap/Math.h"
#include "Ostap/Power.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/MoreMath.h"
#include "Ostap/PhaseSpace.h"
#include "Ostap/Models.h"
// ============================================================================
//  Local 
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "local_hash.h"
#include "local_gsl.h"
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
  /** @var x_sqrt2
   *  \f$\sqrt{2}\f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_sqrt2  = s_SQRT2 ;
  // ==========================================================================
} //                                                 end of anonymous namespace
// ============================================================================


// ============================================================================
/*  constructor  from all parameters 
 *  @param mu location, bias parameter 
 *  @param beta scale parameter 
 */
// ============================================================================
Ostap::Math::Gumbel::Gumbel
( const double mu   , 
  const double beta )
  : m_mu   ( mu   ) 
  , m_beta ( beta ) 
{}
// ============================================================================
bool Ostap::Math::Gumbel::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) {  return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Gumbel::setBeta ( const double value ) 
{
  if ( s_equal ( value , m_beta ) ) {  return false ; }
  m_beta = value ;
  return true ;
}
// ============================================================================
double Ostap::Math::Gumbel::median () const 
{
  static const double s_lnln2 = std::log( std::log ( 2.0L ) ) ;
  return mu() - m_beta * s_lnln2 ;
}  
// ============================================================================
double Ostap::Math::Gumbel::mean () const 
{
  static const double s_gamma = M_EULER ;
  return mu() + m_beta * s_gamma ;
}
// ============================================================================
double Ostap::Math::Gumbel::variance () const 
{
  static const double s_pisq6 = M_PI * M_PI / 6.0L ;
  return m_beta * m_beta * s_pisq6 ;
}
// ============================================================================
double Ostap::Math::Gumbel::sigma () const 
{
  static const double s_pisqr6 = M_PI / std::sqrt ( 6.0L ) ;
  return std::abs ( m_beta ) * s_pisqr6 ;
}
// ============================================================================
double Ostap::Math::Gumbel::skewness () const 
{
  static const double s_skew  = 
    12 * std::sqrt( 6.0L ) * gsl_sf_zeta_int ( 3 ) / ( M_PI * M_PI * M_PI ) ;
  return std::copysign ( s_skew , m_beta ) ;
}
// ============================================================================
// get a value for the function      
// ============================================================================
double Ostap::Math::Gumbel::pdf  ( const double x ) const 
{
  const long double ibeta = 1/m_beta ;
  const long double z     = ( x - m_mu ) * ibeta ;
  return std::abs ( m_beta ) * std::exp ( -( z + std::exp ( -z ) ) ) ;
}
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::Gumbel::cdf ( const double x ) const 
{
  const long double z     = ( x - m_mu ) / m_beta ;
  return 0 < m_beta ? 
    std::exp ( -std::exp ( -z ) ) : 1 - std::exp( -std::exp ( -z ) ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Gumbel::integral ( const double low  ,
                                       const double high ) const 
{
  if  ( s_equal ( low , high ) ) { return 0 ; }
  //
  const long double ibeta = 1/m_beta ;
  const long double zmin = ( low  - m_mu ) * ibeta ;
  const long double zmax = ( high - m_mu ) * ibeta ;
  //
  return 0 < m_beta ? 
    std::exp ( -std::exp ( -zmax ) ) - std::exp ( -std::exp ( -zmin ) ) : 
    std::exp ( -std::exp ( -zmin ) ) - std::exp ( -std::exp ( -zmax ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Gumbel::tag () const 
{ return std::hash_combine ( m_mu , m_beta ) ; }
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
// destructor
// ============================================================================
Ostap::Math::GramCharlierA::~GramCharlierA() {}
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
  const double result_0 = my_exp ( -0.5 * dx * dx ) / m_sigma / s_SQRT2PI ;
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
  static const char s_message[] = "Ingegral(GramCharlierA)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag  () , 
      &F      , 
      low     , high ,                                        // low & high edges
      workspace ( m_workspace ) ,                             // workspace
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE              ,                                   // size of workspace
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
{ return std::hash_combine ( m_mean , m_sigma , m_kappa3 , m_kappa4 ) ; }
// ============================================================================




// ======================================================================
/*  constructor from thresholds and number of particles
 *  @param threshold_L the low-mass  threshold
 *  @param threshold_H the high-mass threshold
 *  @param l           how many particles we consider
 *  @param n           total number of particles ( n>l!)
 *  @param N           degree of polynomial 
 */
// ======================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const double         threshold1  ,
  const double         threshold2  ,
  const unsigned short l           ,
  const unsigned short n           , 
  const unsigned short N           )  // degree of polynomial
  : m_phasespace ( threshold1 , threshold2 , l , n ) 
  , m_positive   ( N , 
                   std::min ( std::abs ( threshold1 ) , std::abs ( threshold2 ) ) , 
                   std::max ( std::abs ( threshold1 ) , std::abs ( threshold2 ) ) ) 
  , m_workspace  ()
{}
// =====================================================================
/*  constructor from the phase space and polynomial degree 
 *  @param ps          phase space factor 
 *  @param N           degree of polynomial 
 */
// =====================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const Ostap::Math::PhaseSpaceNL& ps ,
  const unsigned short             N  )  // degree of polynomial
  : m_phasespace ( ps ) 
  , m_positive   ( N  , ps.lowEdge() , ps.highEdge() ) 
  , m_workspace  ()
{}
// ======================================================================
/*  constructor from phase space and polynomial degree 
 *  @param ps          phase space factor 
 *  @param N           degree of polynomial 
 */
// =====================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const Ostap::Math::PhaseSpaceNL& ps    ,
  const unsigned short             N     , 
  const double                     xlow  , 
  const double                     xhigh ) 
  : m_phasespace ( ps ) 
  , m_positive   ( N  , 
                   std::max ( ps. lowEdge() , std::min ( xlow , xhigh ) ) ,
                   std::min ( ps.highEdge() , std::max ( xlow , xhigh ) ) )
  , m_workspace  ()
{}
// =====================================================================
// evaluate N/L-body modulated phase space
// =====================================================================
double Ostap::Math::PhaseSpacePol::evaluate ( const double x ) const 
{
  //
  if      ( x < m_phasespace . lowEdge () ) { return 0 ; }
  else if ( x > m_phasespace .highEdge () ) { return 0 ; }
  else if ( x < m_positive   .   xmin  () ) { return 0 ; }
  else if ( x > m_positive   .   xmax  () ) { return 0 ; }
  //
  return m_positive ( x ) * m_phasespace ( x ) ;
}
// =====================================================================
// destructor 
// =====================================================================
Ostap::Math::PhaseSpacePol::~PhaseSpacePol(){}
// =====================================================================
// get the integral
// =====================================================================
double Ostap::Math::PhaseSpacePol::integral () const 
{
  //
  if      ( m_phasespace.highEdge() <= m_positive.xmin() ) { return 0 ; }
  else if ( m_phasespace. lowEdge() >= m_positive.xmax() ) { return 0 ; }
  //
  const double mn = std::max ( m_phasespace. lowEdge() ,  m_positive.xmin () ) ;
  const double mx = std::min ( m_phasespace.highEdge() ,  m_positive.xmax () ) ;
  //
  return integral ( mn , mx ) ;
}
// =====================================================================
// get the integral between low and high limits
// =====================================================================
double  Ostap::Math::PhaseSpacePol::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if (           low > high   ) { return -1 * integral ( high , low ) ; }
  //
  if      ( high <= m_phasespace .  lowEdge () ) { return 0 ; }
  else if ( high <= m_positive   .     xmin () ) { return 0 ; }
  else if ( low  >= m_phasespace . highEdge () ) { return 0 ; }
  else if ( low  >= m_positive   .     xmax () ) { return 0 ; }
  //
  const double mn    = std::max ( m_phasespace. lowEdge() ,  m_positive.xmin () ) ;
  const double mx    = std::min ( m_phasespace.highEdge() ,  m_positive.xmax () ) ;
  //
  const double xlow  = std::max ( low  , mn ) ;
  const double xhigh = std::min ( high , mx ) ;
  //
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpacePol> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpacePol)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag () ,  
      &F     , 
      xlow   , xhigh      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::PhaseSpacePol::tag () const 
{ return std::hash_combine ( m_phasespace.tag () , m_positive.tag () ) ; }
// ============================================================================




// ============================================================================
/* constructor from threshold and number of particles
 *  @param threshold_L the low-mass  threshold
 *  @param l           how many particles we consider
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const double         threshold_L ,   // low threshold 
  const unsigned short l           ,   // number of particles 
  const unsigned short N           ,   // degree of polynomial
  const double         tau         ,   // the exponent 
  const double         xhigh       )   // high edge 
  : PhaseSpaceLeftExpoPol ( threshold_L , l , N , tau , threshold_L , xhigh ) 
{}
// ============================================================================
/*  constructor from the phase space and polynomial degree
 *  @param ps          phase space factor
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const Ostap::Math::PhaseSpaceLeft& ps ,
  const unsigned short  N     ,   // degree of polynomial
  const double          tau   ,   // the exponent 
  const double          xhigh )   // high edge 
  : PhaseSpaceLeftExpoPol ( ps , N , tau , ps.threshold () , xhigh ) 
{}
// ============================================================================
/* constructor from threshold and number of particles
 *  @param threshold_L the low-mass  threshold
 *  @param l           how many particles we consider
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xlow        the low  edge 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const double         threshold_L ,   // low threshold 
  const unsigned short l           ,   // number of particles 
  const unsigned short N           ,   // degree of polynomial
  const double         tau         ,   // the exponent 
  const double         xlow        ,   // low edge 
  const double         xhigh       )   // high edge 
  : PhaseSpaceLeftExpoPol ( Ostap::Math::PhaseSpaceLeft ( threshold_L , l ) , 
                            N , tau , xlow , xhigh ) 
{}
// ============================================================================
/* constructor from the phase space and polynomial degree
 *  @param ps          phase space factor
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xlow        the low  edge 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const PhaseSpaceLeft& ps    ,
  const unsigned short  N     ,   // degree of polynomial
  const double          tau   ,   // the exponent 
  const double          xlow  ,   // low edge 
  const double          xhigh )  // high edge
  : m_phasespace  ( ps ) 
  , m_positive    ( N  ,
                    std::max ( ps.threshold () , std::min ( xlow , xhigh ) ) , 
                    std::max ( xlow , xhigh ) ) 
  , m_tau         ( std::abs ( tau ) )
  , m_workspace   ( ) 
{
  Ostap::Assert ( m_phasespace.threshold() <= m_positive.xmin () , 
                  "Invalid setting of threshold/xmin/xmax" ,
                  "Ostap::Math::PhaseSpaceLeftPol" ) ;  
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::~PhaseSpaceLeftExpoPol(){}
// ============================================================================
// evaluate the modulated phase space
// ============================================================================
double Ostap::Math::PhaseSpaceLeftExpoPol::evaluate ( const double x ) const 
{
  if  ( x <= xmin () || x >= xmax () ) { return 0 ; }
  const double xc = 0.5 * ( xmin() + xmax() ) ;
  return  
    m_phasespace ( x  ) * 
    // m_phasespace ( xc ) *  
    m_positive   ( x  ) *  
    std::exp     ( -1 * m_tau  * ( x - xc ) ) ;
}
// ============================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::PhaseSpaceLeftExpoPol::setTau ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_tau ) ) { return false ; }
  m_tau = avalue ;
  return true ;
}
// =============================================================================
// get the integral between low and high limits
// =============================================================================
double Ostap::Math::PhaseSpaceLeftExpoPol::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; } // return 
  else if ( high < low             ) { return -1 * integral ( high , low ) ; } // RETURN
  else if ( high <= xmin ()        ) { return 0 ; }
  else if ( low  >= xmax ()        ) { return 0 ; }
  //
  const double xlow  = std::max ( low  , xmin () ) ;
  const double xhigh = std::min ( high , xmax () ) ;
  //
  // if the exponent plays important role, split the interval 
  if ( !s_zero ( m_tau ) ) 
  {
    if  ( 2 < ( xhigh - xlow ) * m_tau )  
    {
      const double xc = 0.5 * ( xhigh + xlow ) ;
      return integral ( xlow , xc ) + integral ( xc , xhigh ) ;
    }
  }
  //
  /// split near-threshold region 
  const double delta =  xmax() - threshold() ;
  const double len   =  xhigh  -  xlow  ;
  const double x1 = threshold () + 0.05 * delta ;
  if ( 0.05 * delta < len && xlow < x1 && x1 < xhigh ) 
  { return integral ( xlow , x1 ) + integral ( x1 , xhigh ) ; }
  const double x2 = threshold () + 0.15 * delta ;
  if ( 0.10 * delta < len && xlow < x2 && x2 < xhigh ) 
  { return integral ( xlow , x2 ) + integral ( x2 , xhigh ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpaceLeftExpoPol> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpaceLeftExpoPol)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag () ,  
      &F     , 
      xlow   , xhigh      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::PhaseSpaceLeftExpoPol::tag () const 
{ return std::hash_combine ( m_phasespace.tag () , m_positive.tag () , m_tau ) ; }
// ============================================================================






// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ExpoPositive::ExpoPositive
( const unsigned short      N    ,
  const double              tau  ,   
  const double              xmin ,
  const double              xmax )
  : m_positive  ( N , xmin , xmax )
  , m_tau       ( tau ) 
{}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ExpoPositive::ExpoPositive
( const std::vector<double>& pars ,
  const double               tau  ,
  const double               xmin ,
  const double               xmax )
  : m_positive  ( pars , xmin , xmax )
  , m_tau       ( tau ) 
{}
// ============================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::ExpoPositive::setTau ( const double value )
{
  if ( s_equal ( value , m_tau ) ) { return false ; }
  m_tau = value ;
  return true ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::ExpoPositive::operator () ( const double x ) const 
{
  //
  if ( x < xmin() || x > xmax() ) { return 0 ; }
  //
  return my_exp ( m_tau * x ) * m_positive ( x ) ;
}
// ============================================================================
double Ostap::Math::ExpoPositive::integral ( const double low  , 
                                             const double high ) const 
{
  return Ostap::Math::integrate ( m_positive.bernstein() , m_tau , low ,high ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::ExpoPositive::tag () const 
{ return std::hash_combine ( m_positive.tag () , m_tau ) ; }
// ============================================================================



// ============================================================================
/* constructor form scale & shape parameters
 *  param k      \f$k\f$ parameter (shape)
 *  param theta  \f$\theta\f$ parameter (scale)
 */
// ============================================================================
Ostap::Math::GammaDist::GammaDist 
( const double k     ,   // shape parameter  
  const double theta )   // scale parameter
  : m_k     ( std::abs ( k     ) )
  , m_theta ( std::abs ( theta ) )
  , m_aux   ( 0 ) 
{
  // evaluate auxillary parameter 
  m_aux = - m_k * std::log ( m_theta ) - std::lgamma ( m_k ) ;
}
// ============================================================================
// destrructor
// ============================================================================
Ostap::Math::GammaDist::~GammaDist (){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::GammaDist::setK ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  //
  if ( s_equal ( v , m_k ) ) { return false ; }
  //
  m_k = v ;
  //
  if ( s_equal ( 1 , m_k ) ) { m_k    = 1 ; }
  //
  // evaluate auxillary parameter 
  m_aux = -m_k * std::log ( m_theta ) - std::lgamma ( m_k ) ;
  //
  return true ;
}
// ============================================================================
double Ostap::Math::GammaDist::sigma    () const
{ return std::sqrt ( dispersion ()  ) ; }
// ============================================================================
double Ostap::Math::GammaDist::skewness () const
{ return 2.0 / std::sqrt ( m_k )      ; }
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::GammaDist::setTheta ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  //
  if ( s_equal ( v , m_theta ) ) { return false ; }
  //
  m_theta = v ;
  //
  // evaluate auxillary parameter 
  m_aux = -m_k * std::log ( m_theta ) - std::lgamma ( m_k ) ;
  //
  return true ;
}
// ============================================================================
// calculate gamma distribution shape
// ============================================================================
double Ostap::Math::GammaDist::pdf ( const double x ) const
{
  // simple cases 
  if ( x <= 0 ) { return 0 ; }
  // 
  double result = m_aux - x / m_theta  + ( m_k - 1 ) * my_log( x ) ;
  //
  return my_exp ( result ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::GammaDist::integral () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::GammaDist::integral ( const double low  ,
                                          const double high ) const 
{
  //
  if      ( s_equal ( low  , high ) ) { return 0 ; }
  else if (           low  > high   ) { return -1 * integral ( high , low  ) ; }
  else if (           high <= 0     ) { return 0 ; }
  else if (           low  < 0      ) { return      integral ( 0    , high ) ; }
  //
  return 
    gsl_sf_gamma_inc_P ( m_k , high / m_theta ) - 
    gsl_sf_gamma_inc_P ( m_k , low  / m_theta ) ;
}
// ============================================================================
// calculatye the quantile   (0<p<1) 
// ============================================================================
double Ostap::Math::GammaDist::quantile ( const double p ) const 
{
  if      ( p <= 0 ) { return          0 ; }
  else if ( p >= 1 ) { return s_INFINITY ; }
  //
  return gsl_cdf_gamma_Pinv ( p , m_k , m_theta ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GammaDist::tag () const 
{ return std::hash_combine ( m_k , m_theta ) ; }
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
// destructor 
// ============================================================================
Ostap::Math::GenGammaDist::~GenGammaDist(){}
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
  const double xc = ( x - m_low ) / theta() ;  
  const double xt = std::pow ( xc , p () ) ;  
  //
  double result   = ( k () - 1 ) * gsl_sf_log ( xc ) - xt ;
  result         +=  gsl_sf_log     ( p ()  / theta  () ) ;
  result         -=  gsl_sf_lngamma ( k ()  / p      () ) ;
  //return gsl_sf_exp ( result ) ;
  return my_exp ( result ) ;
}
// ============================================================================
double Ostap::Math::GenGammaDist::cdf ( const double x ) const 
{
  //
  if ( x <= m_low || s_equal ( x , m_low ) ) { return 0 ; }
  //
  const double xc = ( x - m_low ) / theta() ;  
  const double xt = std::pow ( xc , p () ) ;
  //
  return gsl_sf_gamma_inc_P ( k () / p () , xt ) ;
}
// ============================================================================
double Ostap::Math::GenGammaDist::integral () const { return 1 ; }
// ============================================================================
double Ostap::Math::GenGammaDist::integral ( const double low  , 
                                             const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;  
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenGammaDist::tag () const 
{ return std::hash_combine ( m_k , m_theta , m_p , m_low ) ; }
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
{
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::Amoroso::~Amoroso(){}
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
  if      ( theta () > 0 && ( x <= m_a || s_equal ( x , m_a ) ) ) { return 0 ; }
  else if ( theta () < 0 && ( x >= m_a || s_equal ( x , m_a ) ) ) { return 0 ; }
  //
  const double xc = ( x - m_a ) / theta ()    ;
  const double xt = std::pow ( xc , beta() ) ;
  //
  double result   = ( alpha() * beta() - 1 ) * gsl_sf_log ( xc )  - xt ; 
  result += gsl_sf_log     ( std::abs ( beta  () / theta() ) ) ;
  result -= gsl_sf_lngamma (            alpha ()   ) ;
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
    return theta2() * ( gx2 / ga - Ostap::Math::pow ( gx1 / ga , 2 ) ) ;
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
{ return std::hash_combine ( m_a , m_theta , m_alpha , m_beta ) ; }
// ============================================================================




// ============================================================================
/* constructor from scale & shape parameters
 *  param k      \f$k\f$ parameter (shape)
 *  param theta  \f$\theta\f$ parameter (scale)
 */
// ============================================================================
Ostap::Math::LogGammaDist::LogGammaDist 
( const double k     ,   // shape parameter  
  const double theta )   // scale parameter
  : m_gamma ( k , theta ) 
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
  // 
  const double z = my_exp ( x ) ;
  return m_gamma ( z ) * z ;
  //
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::LogGammaDist::integral () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::LogGammaDist::integral ( const double low  ,
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
  if      ( p <= 0 ) { return -s_INFINITY ; }
  else if ( p >= 1 ) { return  s_INFINITY ; }
  //
  return my_log ( gsl_cdf_gamma_Pinv ( p , k() , theta() ) ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::LogGammaDist::tag () const 
{ return std::hash_combine ( 1 , m_gamma.tag () ) ; }
// ============================================================================


// ============================================================================
/* constructor form scale & shape parameters
 *  param k      \f$k\f$ parameter (shape)
 *  param theta  \f$\theta\f$ parameter (scale)
 */
// ============================================================================
Ostap::Math::Log10GammaDist::Log10GammaDist 
( const double k     ,   // shape parameter  
  const double theta )   // scale parameter
  : Ostap::Math::LogGammaDist( k , theta ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Log10GammaDist::~Log10GammaDist (){}
// ============================================================================
// calculate log-gamma distribution shape
// ============================================================================
double Ostap::Math::Log10GammaDist::operator() ( const double x ) const
{ return LogGammaDist::operator() ( x * s_LN10 ) * s_LN10 ; }
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::Log10GammaDist::integral () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Log10GammaDist::integral ( const double low  ,
                                               const double high ) const 
{ return LogGammaDist::integral ( low  * s_LN10 , high * s_LN10 ) ; }
// ============================================================================
// calculate the quantile   (0<p<1) 
// ============================================================================
double Ostap::Math::Log10GammaDist::quantile ( const double p ) const 
{
  if      ( p <= 0 ) { return -s_INFINITY ; }
  else if ( p >= 1 ) { return  s_INFINITY ; }
  return LogGammaDist::quantile ( p ) / s_LN10 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Log10GammaDist::tag () const 
{ return std::hash_combine ( 10 , m_gamma.tag () ) ; }
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
// destructor 
// ============================================================================
Ostap::Math::LogGamma::~LogGamma(){}
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
{ return std::hash_combine ( m_nu , m_lambda , m_alpha ) ; }
// ============================================================================



// ============================================================================
// Beta' 
// ============================================================================
/*  constructor with all parameters 
 *  @param alpha \f$\alpha\f$-parameter 
 *  @param beta  \f$\beta\f$-parameter 
 */
// ============================================================================
Ostap::Math::BetaPrime::BetaPrime 
( const double alpha , 
  const double beta  , 
  const double scale , 
  const double shift )
  : m_alpha ( std::abs ( alpha ) )
  , m_beta  ( std::abs ( beta  ) )
  , m_scale ( scale )
  , m_shift ( shift )
  , m_aux () 
{
  m_aux = 1/gsl_sf_beta ( m_alpha , m_beta ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::BetaPrime::~BetaPrime (){}
// ============================================================================
bool Ostap::Math::BetaPrime::setAlpha ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  m_aux   = 1/gsl_sf_beta ( m_alpha , m_beta ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BetaPrime::setBeta  ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_beta  ) ) { return false ; }
  m_beta  = value_ ;
  m_aux   = 1/gsl_sf_beta ( m_alpha , m_beta ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BetaPrime::setScale ( const double value ) 
{
  if ( s_equal ( value , m_scale  ) ) { return false ; }
  m_scale  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BetaPrime::setShift ( const double value ) 
{
  if ( s_equal ( value , m_shift  ) ) { return false ; }
  m_shift  = value ;
  return true ;
}
// ============================================================================
// evaluate beta'-distributions 
// ============================================================================
double Ostap::Math::BetaPrime::pdf ( const double x ) const 
{
  //
  if      ( m_scale >= 0 && x <= m_shift ) { return 0 ; }
  else if ( m_scale <= 0 && x >= m_shift ) { return 0 ; }
  else if ( s_equal ( x , m_shift )      ) { return 0 ; }
  //
  const double y = ( x - m_shift ) / m_scale ;
  //
  return m_aux / std::abs ( m_scale  ) 
    * std::pow (     y ,   alpha () - 1       ) 
    * std::pow ( 1 + y , - alpha () - beta () ) ;  
}
// ============================================================================
double Ostap::Math::BetaPrime::cdf ( const double x ) const 
{
  //
  const double z = ( x - m_shift ) / m_scale ;
  //
  if ( z <= 0 || s_equal ( z , 0 ) ) { return 0 ; }
  //
  const double y = z / ( 1 + z ) ;
  //
  Sentry sentry ;
  //
  return
    0 < m_scale ? 
    gsl_sf_beta_inc (  alpha() , beta() , y ) :
    1 - gsl_sf_beta_inc (  alpha() , beta() , y ) ; 
}
// ============================================================================
double Ostap::Math::BetaPrime::integral ( const double low  , 
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
  if ( beta() <= 1 || s_equal ( beta() , 1 ) ) { return -1.e+9 ; }  
  //
  return m_shift + m_scale * alpha() / ( beta() - 1 ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::mode () const 
{
  if ( alpha() < 1 ) { return 0 ; }
  return m_shift + m_scale * ( alpha() - 1 ) / ( beta() + 1 ) ;
}
// ============================================================================
double Ostap::Math::BetaPrime::variance () const 
{
  if ( beta() <= 2 || s_equal ( beta() , 2 ) ) { return -1.e+9 ; }  
  //
  const double a = alpha () ;
  const double b = beta  () ;
  //
  return m_scale * m_scale * a *  ( a + b + 1 ) / ( b - 2 ) / Ostap::Math::pow ( b - 1 , 2 ) ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::sigma () const 
{
  if ( beta() <= 2 || s_equal ( beta() , 2 ) ) { return -1.e+9 ; }  
  return std::sqrt ( variance () ) ;
}
// ===========================================================================
double Ostap::Math::BetaPrime::skewness  () const 
{
  if ( beta() <= 3 || s_equal ( beta() , 3 ) ) { return -1.e+9 ; }  
  //
  const double a = alpha () ;
  const double b = beta  () ;
  //
  return 2 * ( 2 * a + b - 1 ) / ( b - 3 ) * std::sqrt( ( b - 2 ) / a / ( a + b - 1 ) ) ;
}
// ===========================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::BetaPrime::tag () const 
{ return std::hash_combine ( m_alpha , m_beta , m_scale , m_shift ) ; }
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
// destructor 
// ============================================================================
Ostap::Math::Landau::~Landau (){}
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
{ return std::hash_combine ( m_scale , m_shift ) ; }
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
{ return std::hash_combine ( m_scale , m_shape , m_shift ) ; }
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
// Argus
// ============================================================================
/*  constructor with all parameters 
 */
// ============================================================================
Ostap::Math::Argus::Argus
( const double shape , 
  const double high  ,
  const double low   )
  : m_shape ( std::abs ( shape ) ) 
  , m_high  ( std::abs ( high  ) ) 
  , m_low   ( std::abs ( low   ) ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::Argus::~Argus (){}
// ============================================================================
bool Ostap::Math::Argus::setShape ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_shape  ) ) { return false ; }
  m_shape  = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Argus::setLow ( const double value ) 
{
  if ( s_equal ( value , m_low  ) ) { return false ; }
  m_low  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Argus::setHigh ( const double value ) 
{
  if ( s_equal ( value , m_high  ) ) { return false ; }
  m_high  = value ;
  return true ;
}
// ============================================================================
// evaluate Argus-distributions 
// ============================================================================
namespace 
{
  // ==========================================================================
  inline double phi_ ( const double x ) 
  { return gsl_ran_gaussian_pdf  ( x , 1 ) ; }
  inline double Phi_ ( const double x ) 
  { return gsl_cdf_ugaussian_P   ( x     ) ; }
  inline double Psi_ ( const double x ) 
  { return Phi_ ( x ) - x * phi_ ( x ) - 0.5 ; } 
  // ==========================================================================
} // ==========================================================================
// ============================================================================
double Ostap::Math::Argus::pdf ( const double x ) const 
{
  //
  if      ( x >= std::max ( m_high , m_low ) ) { return 0 ; }
  else if ( x <= std::min ( m_high , m_low ) ) { return 0 ; }
  //
  const double y = y_ ( x ) ;
  if ( y <= 0 || y >= 1 ) { return 0 ; }
  //
  double res   = s_SQRT2PIi ;
  res         *= Ostap::Math::pow ( m_shape , 3 ) ;
  res         /= Psi_   ( m_shape ) ;
  res         *= y ;
  //
  const double y2 = 1 - y * y  ;
  res         *= std::sqrt ( y2 ) ;
  res         *= my_exp ( -0.5 * m_shape * m_shape * y2 ) ;
  //
  return     res / std::abs ( m_high - m_low ) ;
}
// ============================================================================
// evaluate Argus-CDF 
// ============================================================================
double Ostap::Math::Argus::cdf ( const double x ) const 
{
  //
  if      ( x > std::max ( m_high , m_low ) ) { return 1 ; }
  else if ( x < std::min ( m_high , m_low ) ) { return 0 ; }
  //
  const double y  = y_ ( x )  ;
  //
  const double y2 = 1 - y * y ;
  //
  const double res =  Psi_ ( m_shape * y2 ) / Psi_( m_shape ) ;
  return m_high > m_low ?  ( 1 - res ) : res ;
}
// ============================================================================
double Ostap::Math::Argus::integral ( const double low  , 
                                      const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Argus::tag () const 
{ return std::hash_combine ( m_shape , m_high , m_low ) ; }
// ============================================================================



// ============================================================================
// constructor from polynom and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid 
( const Ostap::Math::Positive& poly  , 
  const double                 alpha ,
  const double                 x0    ) 
  : m_positive ( poly  )
  , m_alpha    ( alpha )
  , m_x0       ( x0    )
  , m_workspace() 
{}
// ============================================================================
// constructor from polynom and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid 
( const unsigned short             N , 
  const double                 xmin  , 
  const double                 xmax  , 
  const double                 alpha , 
  const double                 x0    ) 
  : m_positive ( N , xmin , xmax )
  , m_alpha    ( alpha )
  , m_x0       ( x0    )
  , m_workspace() 
{}
// ============================================================================
// constructor from polynom and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid 
( const std::vector<double>&   pars  ,
  const double                 xmin  , 
  const double                 xmax  , 
  const double                 alpha , 
  const double                 x0    ) 
  : m_positive ( pars , xmin , xmax )
  , m_alpha    ( alpha )
  , m_x0       ( x0    )
  , m_workspace() 
{}
// ============================================================================
// set new valeu for alpha 
// ============================================================================
bool Ostap::Math::Sigmoid::setAlpha( const double value )
{
  if ( s_equal ( m_alpha, value ) ) { return false ; }
  m_alpha = value ;
  //
  return true ;
}
// ============================================================================
// set new valeu for x0
// ============================================================================
bool Ostap::Math::Sigmoid::setX0 ( const double value )
{
  if ( s_equal ( m_x0, value ) ) { return false ; }
  m_x0 = value ;
  //
  return true ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Sigmoid::operator () ( const double x ) const
{
  return 
    x < xmin () ? 0               :
    x > xmax () ? 0               :
    s_zero  ( m_alpha )    ? 
    0.5 * m_positive ( x ) :
    0.5 * m_positive ( x ) * ( 1 + std::tanh ( m_alpha * ( x - m_x0 ) ) ) ;
}
// ============================================================================
// get the integral between xmin and xmax 
// ============================================================================
double Ostap::Math::Sigmoid::integral   () const 
{ return integral ( m_positive.xmin () , m_positive.xmax() ) ; }
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::Sigmoid::integral  
( const double low  , 
  const double high ) const 
{
  //
  if      ( high < low                ) { return -integral ( high , low ) ; }
  else if ( s_equal ( low , high )    ) { return 0 ; }
  else if ( high < xmin ()            ) { return 0 ; }
  else if ( low  > xmax ()            ) { return 0 ; }
  //
  else if ( s_zero ( m_alpha ) ) { return m_positive.integral ( low , high ) ; }
  //
  // split it, if needed 
  if ( low < m_x0 && m_x0 < high ) 
  { return integral ( low , m_x0 ) + integral ( m_x0 , high ) ; }
  // split further, if needed 
  const double a1 = m_x0 + 3 / m_alpha ;
  if ( low < a1 && a1 < high ) { return integral ( low , a1 ) + integral ( a1 , high ) ; }
  // split further, if needed  
  const double a2 = m_x0 - 3 / m_alpha ;
  if ( low < a2 && a2 < high ) { return integral ( low , a2 ) + integral ( a2 , high ) ; }
  //
  //
  static const Ostap::Math::GSL::Integrator1D<Sigmoid> s_integrator {} ;
  static char s_message[] = "Integral(Sigmoid)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag  () , 
      &F      , 
      low     , high ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Sigmoid::tag () const 
{ return std::hash_combine ( m_positive.tag () , m_alpha , m_x0 ) ; }
// ============================================================================








// ============================================================================
Ostap::Math::TwoExpos::TwoExpos
( const double alpha ,
  const double delta , 
  const double x0    ) 
  : m_alpha ( std::abs ( alpha ) ) 
  , m_delta ( std::abs ( delta ) ) 
  , m_x0    ( x0 ) 
{}
// ============================================================================
// set new value for x0
// ============================================================================
bool Ostap::Math::TwoExpos::setX0 ( const double value )
{
  if ( s_equal ( m_x0, value ) ) { return false ; }
  m_x0 = value ;
  //
  return true ;
}
// ============================================================================
// set new value for alpha
// ============================================================================
bool Ostap::Math::TwoExpos::setAlpha ( const double value )
{
  const double nv = std::abs ( value ) ;
  if ( s_equal ( m_alpha, nv ) ) { return false ; }
  m_alpha = nv ;
  //
  return true ;
}
// ============================================================================
// set new value for delta
// ============================================================================
bool Ostap::Math::TwoExpos::setDelta ( const double value )
{
  const double nv = std::abs ( value ) ;
  if ( s_equal ( m_delta, nv ) ) { return false ; }
  m_delta = nv ;
  //
  return true ;
}
// ============================================================================
// get the value 
// ============================================================================
double Ostap::Math::TwoExpos::operator() ( const double x ) const 
{ return x < m_x0 ? 0 : derivative ( x , 0 ) ; }
// ============================================================================
// get the integral between -inf and +inf
// ============================================================================
double Ostap::Math::TwoExpos::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::TwoExpos::integral   
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( low  > high            ) { return -integral ( high , low  ) ; }
  else if ( high <= m_x0           ) { return 0 ; }
  else if ( low  <  m_x0           ) { return  integral ( m_x0 , high ) ; }
  //
  const double a     = m_alpha            ;
  const double b     = m_alpha + m_delta  ;
  //
  const double xlow  = low  - m_x0 ;
  const double xhigh = high - m_x0 ;  
  //
  const double norm  = 1.0 / m_alpha - 1.0 / ( m_alpha + m_delta ) ;
  return 
    ( ( std::exp ( -b * xhigh ) - std::exp ( -b * xlow ) ) / b -
      ( std::exp ( -a * xhigh ) - std::exp ( -a * xlow ) ) / a ) / norm ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  inline unsigned long long _factorial_ ( const unsigned short N ) 
  {
    return 
      0 == N ?  1 : 
      1 == N ?  1 : 
      2 == N ?  2 : 
      3 == N ?  6 : 
      4 == N ? 24 : N * _factorial_ ( N - 1 ) ;
  }
  // ==========================================================================
  /// get (un-normalized) moment 
  inline long double _moment_ 
  ( const long double    alpha , 
    const long double    delta , 
    const unsigned short N     ) 
  {
    return _factorial_ ( N ) *  
      ( 1 / Ostap::Math::pow ( alpha         , N + 1 ) - 
        1 / Ostap::Math::pow ( alpha + delta , N + 1 ) ) ;  
  }
  // ==========================================================================
}
// ============================================================================
// get normalization constant
// ============================================================================
double Ostap::Math::TwoExpos::norm () const 
{ return 1.L / _moment_ ( m_alpha , m_delta , 0 ) ; } 
// ============================================================================
// mean-value (for -inf,+inf) interval 
// ============================================================================
double Ostap::Math::TwoExpos::mean  () const 
{
  const long double n0 = _moment_ ( m_alpha , m_delta , 0 ) ;
  const long double n1 = _moment_ ( m_alpha , m_delta , 1 ) ;
  //
  return m_x0 + n1 / n0 ;
}
// ============================================================================
// mode 
// ============================================================================
double  Ostap::Math::TwoExpos::mode  () const 
{
  const long double delta = m_delta ;
  return m_x0 + std::log1p ( delta / m_alpha ) / delta ; 
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::TwoExpos::variance () const 
{
  const long double n0 = _moment_ ( m_alpha , m_delta , 0 ) ;
  const long double n1 = _moment_ ( m_alpha , m_delta , 1 ) ;
  const long double n2 = _moment_ ( m_alpha , m_delta , 2 ) ;
  //
  return ( n2 * n0 - n1 * n1 ) / ( n0 * n0 )  ;
}
// ============================================================================
// sigma 
// ============================================================================
double Ostap::Math::TwoExpos::sigma () const { return std::sqrt ( variance() ) ; }
// ============================================================================
// get the derivative at given value 
// ============================================================================
double Ostap::Math::TwoExpos::derivative  ( const double x    ) const 
{ return x < m_x0 ? 0 : derivative ( x , 1 ) ; }
// ============================================================================
// get the second derivative at given value
// ============================================================================
double Ostap::Math::TwoExpos::derivative2 ( const double x    ) const
{ return x < m_x0 ? 0 : derivative ( x , 2 ) ; }
// ============================================================================
// get the Nth derivative at given value
// ============================================================================
double Ostap::Math::TwoExpos::derivative
( const double   x , 
  const unsigned N ) const 
{
  if      ( x <  m_x0 ) { return            0 ; }
  //
  const long double n0 = _moment_ ( m_alpha , m_delta , 0 ) ;
  const long double dx = x - m_x0 ;
  //
  const long double a  = tau1 () ;
  const long double b  = tau2 () ;
  //
  return 
    ( Ostap::Math::pow ( a , N ) *  std::exp ( a * dx ) - 
      Ostap::Math::pow ( b , N ) *  std::exp ( b * dx ) ) / n0 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::TwoExpos::tag () const 
{ return std::hash_combine ( m_alpha , m_delta , m_x0 ) ; }
// ============================================================================


// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const unsigned short N     , 
  const double         alpha , 
  const double         delta , 
  const double         x0    ,
  const double         xmin  , 
  const double         xmax  ) 
  : m_positive ( N , xmin , xmax    )
  , m_2exp     ( alpha , delta , x0 )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const std::vector<double>& pars  ,
  const double               alpha , 
  const double               delta , 
  const double               x0    ,
  const double               xmin  , 
  const double               xmax  ) 
  : m_positive ( pars  , xmin  , xmax )
  , m_2exp     ( alpha , delta , x0   )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const Ostap::Math::Positive& poly  ,  
  const double                 alpha , 
  const double                 delta , 
  const double                 x0    ) 
  : m_positive ( poly                 )
  , m_2exp     ( alpha , delta , x0   )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const Ostap::Math::Positive& poly   , 
  const Ostap::Math::TwoExpos& expos  ) 
  : m_positive ( poly  )
  , m_2exp     ( expos )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const Ostap::Math::TwoExpos& expos  , 
  const Ostap::Math::Positive& poly   )
  : m_positive ( poly  )
  , m_2exp     ( expos )
{}
// ============================================================================
// get the value 
// ============================================================================
double Ostap::Math::TwoExpoPositive::operator() ( const double x ) const 
{
  return 
    x < x0   () ? 0 :  
    x < xmin () ? 0 : 
    x > xmax () ? 0 : m_positive ( x ) * m_2exp ( x ) ;
}
// ============================================================================
// get the integral between xmin and xmax
// ============================================================================ 
double Ostap::Math::TwoExpoPositive::integral () const
{
  const double xlow = std::max ( x0() , xmin () ) ;
  return xlow < xmax() ? integral ( xlow , xmax () ) : 0 ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::TwoExpoPositive::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low, high ) ) { return 0 ; }
  else if ( low > high            ) { return -integral ( high , low ) ; }
  //
  const long double r1 = 
    Ostap::Math::integrate ( m_positive.bernstein() , tau1 () , low , high ) ;
  const long double r2 = 
    Ostap::Math::integrate ( m_positive.bernstein() , tau2 () , low , high ) ;
  //
  return ( r1 - r2 ) / _moment_ ( alpha() , delta () , 0 ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::TwoExpoPositive::tag () const 
{ return std::hash_combine ( m_positive.tag () , m_2exp.tag () ) ; }
// ============================================================================





// ============================================================================
// Tsallis function 
// ============================================================================
/*  constructor from all parameters 
 *  @param mass particle mass  (M>0)
 *  @param n    n-parameter    (N>1)  
 *  @param T    T-parameter    (T>0)
 */
// ============================================================================
Ostap::Math::Tsallis::Tsallis 
( const double mass  , 
  const double n     ,  
  const double T     ) 
  : m_mass ( std::abs ( mass ) ) 
  , m_n    ( std::abs ( n    ) )
  , m_T    ( std::abs ( T    ) )
  , m_workspace() 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Tsallis::~Tsallis(){}
// ============================================================================
// set new value for mass
// ============================================================================
bool Ostap::Math::Tsallis::setMass ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_mass , avalue ) ) { return false ; }
  m_mass = avalue ;
  return true ;
}
// ============================================================================
// set new value for n-parameter
// ============================================================================
bool Ostap::Math::Tsallis::setN ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_n , avalue ) ) { return false ; }
  m_n = avalue ;
  return true ;
}
// ============================================================================
// set new value for T-parameter
// ============================================================================
bool Ostap::Math::Tsallis::setT ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_T , avalue ) ) { return false ; }
  m_T = avalue ;
  return true ;
}
// ============================================================================
//  get Tsallis PDF  
// ============================================================================
double Ostap::Math::Tsallis::pdf ( const double x ) const 
{ return x <= 0 ? 0.0 : x * std::pow ( 1.0 + eTkin ( x ) / ( m_T * m_n ) , -m_n ) ; }
// ============================================================================
//  get Tsallis integrals  
// ============================================================================
double Ostap::Math::Tsallis::integral 
( const double low  , 
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  else if ( high <= xmin ()        ) { return 0 ; }
  //
  const double _low = std::max ( low , xmin () ) ;
  //
  // split too large intervals
  if ( 0 < m_mass ) 
  {
    // split points 
    static const std::array<int,5> s_split = {{ 1 ,  3  , 10 , 20 , 50 }} ;
    for( const auto p : s_split )
    {
      const double middle = m_mass * p ;
      if (  _low < middle && middle < high ) 
      { return integral ( _low , middle ) + integral ( middle , high ) ; }
    }
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Tsallis> s_integrator {} ;
  static char s_message[] = "Integral(Tsallis)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag  () , 
      &F      , 
      low     , high  ,              // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Tsallis::tag () const 
{ return std::hash_combine ( m_mass , m_n , m_T ) ; }
// ============================================================================





// ============================================================================
// QGSM function 
// ============================================================================
/*  constructor from all parameters 
 *  @param mass particle mass  (M>0)
 *  @param n    n-parameter    (N>1)  
 *  @param T    T-parameter    (T>0)
 */
// ============================================================================
Ostap::Math::QGSM::QGSM 
( const double mass , 
  const double b    ) 
  : m_mass ( std::abs ( mass ) ) 
  , m_b    ( std::abs ( b    ) )
  , m_workspace() 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::QGSM::~QGSM(){}
// ============================================================================
// set new value for mass
// ============================================================================
bool Ostap::Math::QGSM::setMass ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_mass , avalue ) ) { return false ; }
  m_mass = avalue ;
  return true ;
}
// ============================================================================
// set new value for b-parameter
// ============================================================================
bool Ostap::Math::QGSM::setB ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_b , avalue ) ) { return false ; }
  m_b = avalue ;
  return true ;
}
// ============================================================================
//  get QGSM PDF  
// ============================================================================
double Ostap::Math::QGSM::pdf ( const double x ) const 
{ return x <= 0 ? 0.0 : x * std::exp ( -m_b * eTkin ( x ) ) ; }
// ============================================================================
//  get QGSM integrals  
// ============================================================================
double Ostap::Math::QGSM::integral 
( const double low  , 
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  else if ( high <= xmin()         ) { return 0 ; }
  //
  const double _low = std::max ( low , xmin() ) ;
  //
  // split too large intervals
  if ( 0 < m_mass ) 
  {
    // split points 
    static const std::array<int,5> s_split = {{ 1 ,  3  , 10 , 20 , 50 }} ;
    for( const auto p : s_split )
    {
      const double middle = m_mass * p ;
      if (  _low < middle && middle < high ) 
      { return integral ( _low , middle ) + integral ( middle , high ) ; }
    }
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<QGSM> s_integrator {} ;
  static char s_message [] = "Integral(QGSM)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag  () , 
      &F      , 
      low     , high      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::QGSM::tag () const 
{ return std::hash_combine ( m_mass , m_b ) ; }
// ============================================================================





// ============================================================================
// The END
// ============================================================================


