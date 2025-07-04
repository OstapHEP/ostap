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
{ 
  static const std::string s_name = "Gumbel" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_beta ) ; 
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
// ======================================================================
// constructor from phase space and polynomial
// ======================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const PhaseSpaceNL&          ps  ,
  const Ostap::Math::Positive& pol ) 
  : m_phasespace ( ps  ) 
  , m_positive   ( pol ) 
  , m_workspace  () 
{
  Ostap::Assert ( m_phasespace.lowEdge () < m_positive.xmax      () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PhaseSpacePol" ) ;
  Ostap::Assert ( m_positive.xmin      () < m_phasespace.highEdge() , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PhaseSpacePol" ) ;                 
}
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
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () ,  
      &F     , 
      xlow   , xhigh      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::PhaseSpacePol::tag () const 
{
  static const std::string s_name = "PhaseSpacePol" ;
  return Ostap::Utils::hash_combiner ( s_name , m_phasespace.tag () , m_positive.tag () ) ; 
}
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
/*  constructor from the phase space and polynomial
 *  @param ps          phase space factor
 *  @param poly        polynomial
 *  @param tau         the exponent 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const PhaseSpaceLeft&        ps  ,   // pjase space 
  const Ostap::Math::Positive& pol ,   // polynomial 
  const double                 tau )  // the exponent 
  : m_phasespace  ( ps  ) 
  , m_positive    ( pol ) 
  , m_tau         ( std::abs ( tau ) )
  , m_workspace   ( ) 
{
  Ostap::Assert ( m_phasespace.threshold() < m_positive.xmax () , 
                  "Invalid setting of threshold/xmin/xmax" ,
                  "Ostap::Math::PhaseSpaceLeftPol" ) ;  
}
// ======================================================================




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
    m_phasespace ( x  ) / 
    m_phasespace ( xc ) * 
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
// ============================================================================
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
    if  ( 3 < ( xhigh - xlow ) * m_tau )  
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
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () ,  
      &F     , 
      xlow   , xhigh      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::PhaseSpaceLeftExpoPol::tag () const 
{
  static const std::string s_name = "PhaseSpaceLeftExpoPol" ;
  return Ostap::Utils::hash_combiner ( s_name , m_phasespace.tag () , m_positive.tag () , m_tau ) ;
}
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
// ======================================================================
// constructor from polynom and exponential 
// ======================================================================
Ostap::Math::ExpoPositive::ExpoPositive
( const Ostap::Math::Positive& pol , 
  const double                 tau ) 
  : m_positive  ( pol )
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
{
  static const std::string s_name = "ExpoPositive" ;
  return Ostap::Utils::hash_combiner ( s_name , m_positive.tag () , m_tau ) ; 
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
{ 
  static const std::string s_name = "GammaDist" ;
  return Ostap::Utils::hash_combiner ( s_name  , m_k , m_theta ) ;
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
{ 
  static const std::string s_name = "LogGamma" ;
  return Ostap::Utils::hash_combiner ( s_name , m_nu , m_lambda , m_alpha ) ; 
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
  Ostap::Assert ( 0 < m_alpha && !s_zero ( m_alpha ) ,
		  "Invalid value of alpha (must be positive" ,
		  "Ostap::Math::BetaPrime" ) ;
  Ostap::Assert ( 0 < m_beta  && !s_zero ( m_beta ) ,
		  "Invalid value of beta  (must be positive" ,
		  "Ostap::Math::BetaPrime" ) ;  
  // 
  m_aux = 1.0 / Ostap::Math::beta ( m_alpha , m_beta ) ;
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
  m_aux   = 1.0/Ostap::Math::beta ( m_alpha , m_beta ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BetaPrime::setBeta  ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_beta  ) ) { return false ; }
  m_beta  = value_ ;
  m_aux   = 1.0/Ostap::Math::beta ( m_alpha , m_beta ) ;
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
  return m_aux / std::abs ( m_scale ) 
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
    gsl_sf_beta_inc (  alpha() , beta() , y ) : 1 - gsl_sf_beta_inc ( alpha() , beta() , y ) ; 
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
  if ( beta() <= 1 || s_equal ( beta() , 1 ) )
    { return std::numeric_limits<double>::quiet_NaN(); } 
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
  if ( beta() <= 2 || s_equal ( beta() , 2 ) )
    { return std::numeric_limits<double>::quiet_NaN(); } 
  //
  const double a = alpha () ;
  const double b = beta  () ;
  //
  return m_scale * m_scale * a *  ( a + b + 1 ) / ( b - 2 ) / Ostap::Math::POW ( b - 1 , 2 ) ;
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
  if ( beta() <= 3 || s_equal ( beta() , 3 ) )
    { return std::numeric_limits<double>::quiet_NaN(); } 
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
{ 
  static const std::string s_name = "BetaPrime" ;
  return Ostap::Utils::hash_combiner ( s_name , m_alpha , m_beta , m_scale , m_shift ) ; 
}
// ============================================================================


// ============================================================================
// generalized Beta' 
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
// destructor 
// ============================================================================
Ostap::Math::GenBetaPrime::~GenBetaPrime (){}
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
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Sigmoid::tag () const 
{
  static const std::string s_name = "Sigmoid" ;
  return Ostap::Utils::hash_combiner ( s_name , m_positive.tag () , m_alpha , m_x0 ) ; 
}
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
      ( 1 / Ostap::Math::POW ( alpha         , N + 1 ) - 
        1 / Ostap::Math::POW ( alpha + delta , N + 1 ) ) ;  
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
    ( Ostap::Math::POW ( a , N ) *  std::exp ( a * dx ) - 
      Ostap::Math::POW ( b , N ) *  std::exp ( b * dx ) ) / n0 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::TwoExpos::tag () const 
{
  static const std::string s_name = "TwoExpos" ;
  return Ostap::Utils::hash_combiner ( s_name , m_alpha , m_delta , m_x0 ) ; 
}
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
{ 
  static const std::string s_name = "TwoExposPositive" ;
  return Ostap::Utils::hash_combiner ( s_name , m_positive.tag () , m_2exp.tag () ) ; 
}
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
double Ostap::Math::Tsallis::evaluate ( const double x ) const 
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
    static const std::array<int,7> s_split = {{ 1 ,  3  , 10 , 20 , 50 , 100 , 1000 }} ;
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
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high  ,              // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
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
{ 
  static const std::string s_name = "Tsallis" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mass , m_n , m_T ) ; 
}
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
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
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
{ 
  static const std::string s_name = "QGSM" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mass , m_b ) ; 
}
// ============================================================================



// ============================================================================
/*  Constrcutor for all parameters
 *  @param mass   mas sof th eparticle 
 *  @param beta   inverse temperature
 */
// ============================================================================
Ostap::Math::Hagedorn::Hagedorn 
( const double mass , 
  const double beta ) 
  : m_mass ( std::abs ( mass ) ) 
  , m_beta ( std::abs ( beta ) ) 
{}
// ============================================================================
// set new value for mass
// ============================================================================
bool Ostap::Math::Hagedorn::setMass ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_mass , avalue ) ) { return false ; }
  m_mass = avalue ;
  return true ;
}
// ============================================================================
// set new value for inverse temporature 
// ============================================================================
bool Ostap::Math::Hagedorn::setBeta ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_beta , avalue ) ) { return false ; }
  m_beta = avalue ;
  return true ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::Hagedorn::evaluate ( const double x ) const 
{
  if ( x <= 0 ) { return 0 ; }
  //
  const double mt  = mT ( x )    ;
  const double arg = m_beta * mt ;
  //
  return 
    // GSL_LOG_DBL_MAX < arg ? 0.0 : 
    300                < arg ? 0.0 : 
    x * mt * Ostap::Math::bessel_Kn ( 1 , arg ) / m_beta ;
}
// ============================================================================
//  get Hagedorn integral  
// ============================================================================
double Ostap::Math::Hagedorn::integral 
( const double low  , 
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  else if ( high <= 0              ) { return 0 ; }
  //
  const double xmin = std::max ( 0.0 , low ) ;
  const double xmax = high ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Hagedorn> s_integrator {} ;
  static char s_message [] = "Integral(HAgedorn)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  ()                   , 
      &F                        , 
      xmin    , xmax            ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION              ,   // absolute precision
      s_RPRECISION              ,   // relative precision
      m_workspace.size()        ,   // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
  //
}
// ============================================================================
// get the mean value 
// ============================================================================
double Ostap::Math::Hagedorn::mean () const 
{
  const double mt = m_mass * m_beta ;
  return 
    s_SQRTPIHALF * std::sqrt ( m_mass / m_beta ) * 
    Ostap::Math::bessel_Knu ( 5.0/2 , mt ) /
    Ostap::Math::bessel_Kn  ( 2     , mt ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Hagedorn::tag () const 
{ 
  static const std::string s_name = "Hagedorn" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mass , m_beta ) ; 
}
// ============================================================================

  


// ============================================================================
/*  Constructor from all parameters
 *  @param right dump direction
 *  @param x0    threshold value 
 *  @param sigma sigma  
 */
// ============================================================================
Ostap::Math::CutOffGauss::CutOffGauss
( const bool   right , 
  const double x0    , 
  const double sigma ) 
  : m_right ( right ) 
  , m_x0    ( x0                 )
  , m_sigma ( std::abs ( sigma ) ) 
{}
// =========================================================================
// update sigma
// ============================================================================
bool Ostap::Math::CutOffGauss::setSigma ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigma , avalue ) ) { return false ; }
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
// update x0
// ============================================================================
bool Ostap::Math::CutOffGauss::setX0 ( const double value )
{
  if ( s_equal ( m_x0 , value ) ) { return false ; }
  m_x0 = value ;
  return true ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::CutOffGauss::operator() ( const double x ) const 
{
  // 
  if      (  m_right && x <= m_x0 ) { return 1 ; }
  else if ( !m_right && x >= m_x0 ) { return 1 ; }
  //
  const double dx = ( x - m_x0 ) / m_sigma ;
  return std::exp ( -0.5 * dx * dx ) ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CutOffGauss::integral
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return 0                        ; }
  else if ( low > high             ) { return -integral ( high , low ) ; }
  //
  if ( low < m_x0 && m_x0 < high ) 
  { return integral ( low , m_x0 ) + integral ( m_x0 , high ) ; }
  //
  if      (  m_right && high <= m_x0 ) { return high - low ; }
  else if ( !m_right && low  >= m_x0 ) { return high - low ; }
  //
  static const double s_norm = std::sqrt ( 2.0 * M_PI ) ;
  return s_norm * m_sigma * 
    ( Ostap::Math::gauss_cdf ( high , m_x0 , m_sigma ) -
      Ostap::Math::gauss_cdf ( low  , m_x0 , m_sigma ) ) ;  
}

// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::CutOffGauss::tag () const 
{  
  static const std::string s_name = "CutOffGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_right , m_x0 , m_sigma ) ; 
}
// ============================================================================




// ============================================================================
/*  Constructor from all parameters
 *  @param right dump direction
 *  @param x0    threshold value 
 *  @param nu    parameter nu 
 *  @param sigma parameter sigma  
 */
// ============================================================================
Ostap::Math::CutOffStudent::CutOffStudent 
( const bool   right , 
  const double x0    , 
  const double nu    ,
  const double sigma ) 
  : m_right  ( right ) 
  , m_x0     ( x0    ) 
  , m_nu     ( -1    ) 
  , m_sigma  ( std::abs ( sigma ) ) 
  , m_C      ( -1    )
{
  setNu ( nu ) ;
}
// =========================================================================
// update sigma
// ============================================================================
bool Ostap::Math::CutOffStudent::setSigma ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigma , avalue ) ) { return false ; }
  m_sigma = avalue ;
  return true ;
}
// =========================================================================
// update nu
// ============================================================================
bool Ostap::Math::CutOffStudent::setNu ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_nu , avalue ) ) { return false ; }
  m_nu = avalue ;
  //
  m_C  = std::exp ( - std::lgamma (  0.5 * ( m_nu + 1 ) )  
                    + std::lgamma (  0.5 * ( m_nu     ) ) 
                    + 0.5 * std::log ( m_nu * M_PI ) ) ;
  return true ;
}
// ============================================================================
// update x0
// ============================================================================
bool Ostap::Math::CutOffStudent::setX0 ( const double value )
{
  if ( s_equal ( m_x0 , value ) ) { return false ; }
  m_x0 = value ;
  return true ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::CutOffStudent::operator() ( const double x ) const 
{
  // 
  if      (  m_right && x <= m_x0 ) { return 1 ; }
  else if ( !m_right && x >= m_x0 ) { return 1 ; }
  //
  const double dx = ( x - m_x0 ) / m_sigma ;
  return std::pow ( 1 + dx * dx / m_nu , -0.5 * ( m_nu + 1 ) ) ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CutOffStudent::integral
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return 0                        ; }
  else if ( low > high             ) { return -integral ( high , low ) ; }
  //
  if ( low < m_x0 && m_x0 < high ) 
  { return integral ( low , m_x0 ) + integral ( m_x0 , high ) ; }
  //
  if      (  m_right && high <= m_x0 ) { return high - low ; }
  else if ( !m_right && low  >= m_x0 ) { return high - low ; }
  //
  const double xl = ( low  - m_x0 ) / m_sigma ;
  const double xh = ( high - m_x0 ) / m_sigma ;
  //
  return m_C * m_sigma * ( Ostap::Math::student_cdf ( xh , m_nu ) - 
                           Ostap::Math::student_cdf ( xl , m_nu ) )  ;
}

// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::CutOffStudent::tag () const 
{ 
  static const std::string s_name = "CutOffStudent" ;
  return Ostap::Utils::hash_combiner ( s_name , m_right , m_x0 , m_nu , m_sigma ) ;
}
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
  return m_shift + m_varsigma * s_SQRTPIHALF * 
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
  return 2 * s2 + n2 - 0.5 * M_PI * s2 * l * l ;  
}
// ============================================================================
/// evaluate the function
// ============================================================================
double Ostap::Math::Rice::operator() ( const double x ) const 
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
  
// =============================================================================
// constructor for all elements 
// =============================================================================
Ostap::Math::Argus::Argus 
( const double mu  , 
  const double  c  , 
  const double chi ) 
  : m_mu   ( mu ) 
  , m_c    ( std::abs ( c   ) )
  , m_chi  ( std::abs ( chi ) )  
  , m_norm ( -1 ) 
{ setChi ( chi ) ; }
// =============================================================================
// set mu parameter
// =============================================================================
bool Ostap::Math::Argus::setMu  ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set c parameter
// =============================================================================
bool Ostap::Math::Argus::setC   ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_c , avalue ) ) { return false ; }
  m_c = avalue ;
  return true ;    
}
// =============================================================================
// set chi parameter
// =============================================================================
bool Ostap::Math::Argus::setChi ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_chi , avalue ) && 0 < m_norm ) { return false ; }
  m_chi  = avalue ;
  m_norm = std::pow ( m_chi , 3 ) / psi ( m_chi ) * s_SQRT2PIi  ;
  return true ;    
}
// ============================================================================
/*  helper function 
 *  \f$ \Psi ( \chi ) = \Phi(\chi )  - \chi \phi  (\chi ) - \frac{1}{2} \f$ 
 */
// ============================================================================
double Ostap::Math::Argus::psi ( const double value ) const 
{ return Ostap::Math::gauss_cdf ( value ) - value * Ostap::Math::gauss_pdf ( value ) - 0.5 ; }
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Argus::evaluate   ( const double x ) const 
{
  if ( x + m_c <= m_mu || m_mu <= x ) { return 0 ; }
  const double dx = ( x + m_c - m_mu ) / m_c ; 
  const double dd = 1 - dx * dx ;
  return m_norm * dx * std::sqrt ( dd ) * std::exp ( -0.5 * m_chi * m_chi * dd ) / m_c ;  
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::Argus::cdf        ( const double x ) const 
{
  //
  if       ( x + m_c <= m_mu ) { return 0 ; }
  else  if ( m_mu <= x       ) { return 1 ; }
  //
  const double dx = ( x + m_c - m_mu ) / m_c ;
  const double dd = std::sqrt ( 1 - dx * dx ) ;
  //
  return 1 - psi ( m_chi * dd ) / psi ( m_chi ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::Argus::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Argus::integral  
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high )             ) { return 0 ; }
  else if ( high < low                         ) { return -integral ( high , low ) ; }
  else if ( high + m_c <=  m_mu                ) { return 0 ; }
  else if ( m_mu       <=  low                 ) { return 0 ; }
  else if ( low  + m_c <= m_mu && m_mu <= high ) { return 1 ; }    
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ===========================================================================
// mean of the distribution 
// ===========================================================================
double Ostap::Math::Argus::mean     () const 
{
  const double c2  = 0.25 * m_chi * m_chi ;
  return ( m_mu - m_c ) + 
    0.5 * m_c * m_chi * s_SQRTPIHALF * std::exp ( - c2 ) * Ostap::Math::bessel_In ( 1 , c2 ) / psi ( m_chi );
}
// ===========================================================================
// mode of the distribution 
// ===========================================================================
double Ostap::Math::Argus::mode     () const 
{
  const double c2  = m_chi * m_chi ;
  return ( m_mu - m_c ) + 
    m_c * s_SQRT2i * std::sqrt ( ( c2 - 2 ) + std::sqrt ( c2 * c2 + 4 ) ) / m_chi ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Argus::tag () const 
{ 
  static const std::string s_name = "Argus" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_c , m_chi ) ;
}






// =============================================================================
// constructor for all elements 
// =============================================================================
Ostap::Math::GenArgus::GenArgus 
( const double mu  , 
  const double  c  , 
  const double chi ,
  const double dp  ) 
  : m_mu   ( mu ) 
  , m_c    ( std::abs ( c   ) )
  , m_chi  ( std::abs ( chi ) )  
  , m_dp   ( std::abs ( dp  ) )  
  , m_norm ( -1 ) 
{ 
  setChi ( chi ) ;
  setDp  ( dp  ) ;
}
// =============================================================================
// set mu parameter
// =============================================================================
bool Ostap::Math::GenArgus::setMu  ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set c parameter
// =============================================================================
bool Ostap::Math::GenArgus::setC   ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_c , avalue ) ) { return false ; }
  m_c = avalue ;
  return true ;    
}
// =============================================================================
// set chi parameter
// =============================================================================
bool Ostap::Math::GenArgus::setChi ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_chi , avalue ) && 0 < m_norm ) { return false ; }
  m_chi  = avalue ;
  //
  const double c2 = m_chi * m_chi ;
  const double p1 = p() + 1 ;
  m_norm = 2 * std::pow ( 0.5 * c2 , p1 ) /
    ( std::tgamma ( p() + 1 ) * ( 1 - Ostap::Math::gamma_inc_Q (   p1 , 0.5 * c2 ) ) ) ;
  //
  return true ;    
}
// =============================================================================
// set dp parameter
// =============================================================================
bool Ostap::Math::GenArgus::setDp ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_dp , avalue ) && 0 < m_norm ) { return false ; }
  m_dp   = avalue ;
  //
  const double c2 = m_chi * m_chi ;
  const double p1 = p() + 1 ;
  m_norm = 2 * std::pow ( 0.5 * c2 , p1 ) /
    ( std::tgamma ( p() + 1 ) * ( 1 - Ostap::Math::gamma_inc_Q (   p1 , 0.5 * c2 ) ) ) ;
  //
  return true ;    
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::GenArgus::evaluate   ( const double x ) const 
{
  if ( x + m_c <= m_mu || m_mu <= x ) { return 0 ; }
  const double dx = ( x + m_c - m_mu ) / m_c ; 
  const double dd = 1 - dx * dx ;
  return m_norm * dx * std::pow ( dd , p () ) * std::exp ( -0.5 * m_chi * m_chi * dd ) / m_c ;  
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::GenArgus::cdf        ( const double x ) const 
{
  //
  if       ( x + m_c <= m_mu ) { return 0 ; }
  else  if ( m_mu <= x       ) { return 1 ; }
  //
  const double dx = ( x + m_c - m_mu ) / m_c ;
  const double dd = 1 - dx * dx  ;
  //
  const double p1 = p() + 1 ;
  const double c2 = 0.5 * m_chi * m_chi ;
  //
  const double a1 = Ostap::Math::gamma_inc_Q ( p1 , c2 * dd ) ;
  const double a2 = Ostap::Math::gamma_inc_Q ( p1 , c2      ) ;
  //
  return ( a1 - a2 ) / ( 1 - a2 ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::GenArgus::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::GenArgus::integral  
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high )             ) { return 0 ; }
  else if ( high < low                         ) { return -integral ( high , low ) ; }
  else if ( high + m_c <=  m_mu                ) { return 0 ; }
  else if ( m_mu       <=  low                 ) { return 0 ; }
  else if ( low  + m_c <= m_mu && m_mu <= high ) { return 1 ; }    
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenArgus::tag () const 
{ 
  static const std::string s_name = "GenArgus" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_c , m_chi , m_dp ) ;
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
/*  constructor fron all parameters
 *  @param a position of the left paraboilic horn
 *  @param delta half-distance fron left to right parabolic horn 
 *  @param phi   linear correction parameter ("efficiency")
 */
// ============================================================================
Ostap::Math::HORNSdini::HORNSdini
( const double a     , 
  const double delta ,
  const double phi   ) 
  : m_a        ( a   ) 
  , m_delta    ( std::abs ( delta ) )
  , m_phi      ( phi ) 
  , m_cos2_phi ( std::pow ( std::cos ( phi + 0.25 * M_PI ) , 2 ) ) 
  , m_sin2_phi ( std::pow ( std::sin ( phi + 0.25 * M_PI ) , 2 ) ) 
{}
// ======================================================================
// evaluate the function 
// ======================================================================
double Ostap::Math::HORNSdini::evaluate ( const double x )  const 
{
  if ( x < xmin () || xmax () < x ) { return 0 ; }
  const double z = ( x - m_a ) / m_delta - 1 ;
  return 1.5 * z * z * ( 1 + z * ( m_cos2_phi - m_sin2_phi ) ) / m_delta ;                        
}
// =============================================================================
// set a parameter
// =============================================================================
bool Ostap::Math::HORNSdini::setA   ( const double value ) 
{
  if ( s_equal ( m_a , value ) ) { return false ; }
  m_a = value ;
  return true ;    
}
// =============================================================================
// set delta parameter
// =============================================================================
bool Ostap::Math::HORNSdini::setDelta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_delta , avalue ) ) { return false ; }
  m_delta = avalue ;
  return true ;    
}
// =============================================================================
// set phi parameter
// =============================================================================
bool Ostap::Math::HORNSdini::setPhi ( const double value ) 
{
  if ( s_equal ( m_phi , value ) ) { return false ; }
  //
  m_phi      = value ;
  m_cos2_phi = std::pow ( std::cos ( m_phi + 0.25 * M_PI ) , 2 ) ;
  m_sin2_phi = std::pow ( std::sin ( m_phi + 0.25 * M_PI ) , 2 ) ;
  //
  return true ;    
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::HORNSdini::integral  () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::HORNSdini::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high )           ) { return 0 ; }
  else if ( high < low                       ) { return - integral ( high , low ) ; }
  else if ( high < xmin ()                   ) { return 0 ; }
  else if ( low  > xmax ()                   ) { return 0 ; }
  else if ( low <= xmin () && high>= xmax () ) { return 1 ; }
  //
  const double xl  = std::max ( low  , xmin () ) ;
  const double xh  = std::min ( high , xmax () ) ;
  //
  const double zl  = ( xl - m_a ) / m_delta - 1  ;
  const double zh  = ( xh - m_a ) / m_delta - 1 ;
  //
  const double zh3 = std::pow ( zh , 3 ) ;
  const double zl3 = std::pow ( zl , 3 ) ;
  //
  return ( ( zh3 - zl3 ) / 3 + 
           ( m_cos2_phi - m_sin2_phi ) * ( zh3 * zh - zl3 * zl ) / 4 ) * 1.5 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::HORNSdini::tag () const 
{ 
  static const std::string s_name = "HORNSdini" ;
  return Ostap::Utils::hash_combiner ( s_name , m_a , m_delta , m_phi ) ;
}


// ============================================================================
/*  constructor fron all parameters
 *  @param a position of the left paraboilic horn
 *  @param delta half-distance fron left to right parabolic horn 
 *  @param phi   linear correction parameter ("efficiency")
 */
// ============================================================================
Ostap::Math::HILLdini::HILLdini
( const double a     , 
  const double delta ,
  const double phi   ) 
  : m_a        ( a   ) 
  , m_delta    ( std::abs ( delta ) )
  , m_phi      ( phi ) 
  , m_cos2_phi ( std::pow ( std::cos ( phi + 0.25 * M_PI ) , 2 ) ) 
  , m_sin2_phi ( std::pow ( std::sin ( phi + 0.25 * M_PI ) , 2 ) ) 
{}
// ======================================================================
// evaluate the function 
// ======================================================================
double Ostap::Math::HILLdini::evaluate ( const double x )  const 
{
  if ( x < xmin () || xmax () < x ) { return 0 ; }
  const double z  = ( x - m_a ) / m_delta - 1 ;
  return 0.75 * ( 1 - z * z ) * ( 1 + z * ( m_cos2_phi - m_sin2_phi ) ) / m_delta ;                        
}
// =============================================================================
// set a parameter
// =============================================================================
bool Ostap::Math::HILLdini::setA   ( const double value ) 
{
  if ( s_equal ( m_a , value ) ) { return false ; }
  m_a = value ;
  return true ;    
}
// =============================================================================
// set delta parameter
// =============================================================================
bool Ostap::Math::HILLdini::setDelta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_delta , avalue ) ) { return false ; }
  m_delta = avalue ;
  return true ;    
}
// =============================================================================
// set phi parameter
// =============================================================================
bool Ostap::Math::HILLdini::setPhi ( const double value ) 
{
  if ( s_equal ( m_phi , value ) ) { return false ; }
  //
  m_phi      = value ;
  m_cos2_phi = std::pow ( std::cos ( m_phi + 0.25 * M_PI ) , 2 ) ;
  m_sin2_phi = std::pow ( std::sin ( m_phi + 0.25 * M_PI ) , 2 ) ;
  //
  return true ;    
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::HILLdini::integral  () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::HILLdini::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high )           ) { return 0 ; }
  else if ( high < low                       ) { return - integral ( high , low ) ; }
  else if ( high < xmin ()                   ) { return 0 ; }
  else if ( low  > xmax ()                   ) { return 0 ; }
  else if ( low <= xmin () && high>= xmax () ) { return 1 ; }
  //
  const double xl  = std::max ( low  , xmin () ) ;
  const double xh  = std::min ( high , xmax () ) ;
  //
  const double zh  = ( xh - m_a ) / m_delta - 1 ; 
  const double zl  = ( xl - m_a ) / m_delta - 1 ;
  //
  const double zh3 = std::pow ( zh , 3 ) ;
  const double zl3 = std::pow ( zl , 3 ) ;
  //
  const double A   = m_cos2_phi - m_sin2_phi ;
  //
  return ( (   zh       - zl       )     
           + ( zh  * zh - zl  * zl ) * A / 2 
           - ( zh3      - zl3      )     / 3
           - ( zh3 * zh - zl3 * zl ) * A / 4 ) * 0.75 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::HILLdini::tag () const 
{ 
  static const std::string s_name = "HILLdini" ;
  return Ostap::Utils::hash_combiner ( s_name , m_a , m_delta , m_phi ) ;
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
  const double z  = std::sqrt ( z ) ;
  const double zi = 1/z ;   
  //
  const double qq = Ostap::Math::gauss_pdf ( ( z - zi ) / m_gamma ) ;
  if ( s_zero ( qq )  )  { return 0 ; }
  //
  return qq *  ( z + zi ) / ( 2 * m_gamma * x ) ;
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
//                                                                      The END
// ============================================================================


