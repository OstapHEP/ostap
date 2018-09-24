// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PhaseSpace.h"
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_sf_gamma.h"
// ============================================================================
// Local
// ============================================================================
#include "local_gsl.h"
#include "local_math.h"
#include "Integrator1D.h"
// ============================================================================
/** @file 
 *  Implementation file for classes from the file Ostap/PhaseSpace.h
 *  @author Vanya Belyaev
 */  
// ============================================================================
// constructor from two masses
// ============================================================================
Ostap::Math::PhaseSpace2::PhaseSpace2
( const double m1 ,
  const double m2 )
  : m_m1 ( std::abs ( m1 ) )
  , m_m2 ( std::abs ( m2 ) )
{}
// ============================================================================
// desctructor
// ============================================================================
Ostap::Math::PhaseSpace2::~PhaseSpace2(){}
// ============================================================================
// evaluate 2-body phase space
// ============================================================================
double Ostap::Math::PhaseSpace2::operator () ( const double x ) const
{ return phasespace ( x , m_m1 , m_m2 ) ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::PhaseSpace2::integral
( const double low  ,
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return  0 ; }                          // RETURN
  else if (           low > high   ) { return -1*integral ( high , low  ) ; } // RETURN
  //
  if ( lowEdge() >= high  ) { return 0 ; }
  //
  const double xlow  = std::max ( lowEdge() , low  ) ;
  const double xhigh = std::max ( lowEdge() , high ) ;
  //
  if ( xlow >= xhigh ) { return 0.0 ; }
  //
  if ( 0 < lowEdge()
       && !s_equal ( std::min ( m_m1 , m_m2 ) , 0 )
       && ( xhigh - xlow ) > 20 * lowEdge() ) 
  {
    return 
      integral ( xlow , 0.5 * ( xhigh + xlow )         ) + 
      integral (        0.5 * ( xhigh + xlow ) , xhigh ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpace2> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpace2)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      xlow , xhigh        ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================a
// get the momentum at center of mass 
// ============================================================================a
double Ostap::Math::PhaseSpace2::q_  ( const double x ) const 
{ return q ( x , m1() , m2() ) ; }
// ============================================================================a
// get the momentum at center of mass 
// ============================================================================a
std::complex<double>
Ostap::Math::PhaseSpace2::q1_ ( const double x ) const 
{ return q1 ( x , m1() , m2() ) ; }
// ============================================================================
/* calculate the phase space for   m -> m1 + m2
 *  \f$ \Phi = \frac{1}{8\pi} \frac{ \lambda^{\frac{1}{2}} 
 *  \left( m^2 , m_1^2, m_2_2 \right) }{ m^2 }\f$,
 *  where \f$\lambda\f$ is a triangle function
 */
// ============================================================================
double Ostap::Math::PhaseSpace2::phasespace
( const double         m  ,
  const double         m1 ,
  const double         m2 ,
  const unsigned short L  )
{
  //
  if ( 0 >= m || 0 > m1 || 0 > m2 ) { return 0 ; } // RETURN
  if ( m < m1 + m2                ) { return 0 ; } // RETURN
  //
  const double msq = m * m ;
  const double lam = triangle ( msq  , m1 * m1 , m2 * m2 ) ;
  //
  static const double s_inv8pi = 1.0 / ( 8 * M_PI ) ;
  //
  return 0 < lam ?
    s_inv8pi * Ostap::Math::pow ( std::sqrt ( lam ) / msq , 2 * L + 1 ) : 0.0 ;
}
// ============================================================================
/*  calculate the triangle function
 *  \f$ \lambda ( a , b, c ) = a^2 + b^2 + c^2 - 2ab - 2bc - 2 ca \f$
 *  @param a parameter a
 *  @param b parameter b
 *  @param c parameter b
 */
// ============================================================================
double
Ostap::Math::PhaseSpace2::triangle
( const double a ,
  const double b ,
  const double c )
{ return a * a + b * b + c * c - 2 * a * b - 2 * b * c - 2 * a * c ; }
// ============================================================================
/*  calculate the particle momentum in rest frame
 *  @param m  the mass
 *  @param m1 the mass of the first particle
 *  @param m2 the mass of the second particle
 *  @return the momentum in rest frame (physical values only)
 */
// ============================================================================
double Ostap::Math::PhaseSpace2::q
( const double m  ,
  const double m1 ,
  const double m2 )
{
  //
  if ( 0 >= m || 0 > m1 || 0 > m2 ) { return 0 ; }
  //
  const double lam = triangle ( m * m  , m1 * m1 , m2 * m2 ) ;
  //
  return 0 < lam ? 0.5  * std::sqrt (  lam ) / m : 0 ;
}
// ============================================================================
/*  calculate the particle momentum in rest frame
 *  @param m the mass
 *  @param m1 the mass of the first particle
 *  @param m2 the mass of the second particle
 *  @return the momentum in rest frame
 *  @return the momentum in rest frame  (imaginary for non-physical branch)
 */
// ============================================================================
std::complex<double>
Ostap::Math::PhaseSpace2::q1
( const double m  ,
  const double m1 ,
  const double m2 )
{
  //
  if ( 0 >= m || 0 > m1 || 0 > m2 ) { return 0 ; }
  //
  const double lam = triangle ( m * m , m1 * m1 , m2 * m2 ) ;
  //
  return
    0 <= lam ?
    std::complex<double> (     0.5  * std::sqrt (  lam ) / m , 0 ) :
    std::complex<double> ( 0 , 0.5  * std::sqrt ( -lam ) / m     ) ;
}
// ============================================================================

// ============================================================================
/*  constructor from three masses
 *  @param m1 the mass of the first  particle
 *  @param m2 the mass of the second particle
 *  @param m3 the mass of the third  particle
 *  @param l1 the angular momentum between 1st and 2nd particle
 *  @param l2 the angular momentum between the pair and 3rd particle
 */
// ============================================================================
Ostap::Math::PhaseSpace3::PhaseSpace3
( const double         m1 ,
  const double         m2 ,
  const double         m3 ,
  const unsigned short l1 ,
  const unsigned short l2 )
  : m_m1  ( std::abs ( m1 ) )
  , m_m2  ( std::abs ( m2 ) )
  , m_m3  ( std::abs ( m3 ) )
  , m_l1  ( l1 )
  , m_l2  ( l2 )
  , m_tmp ( 0  )   
{}
// ============================================================================
// deststructor
// ============================================================================
Ostap::Math::PhaseSpace3::~PhaseSpace3 () {}
// ============================================================================
// evaluate 3-body phase space
// ============================================================================
namespace 
{
  // ==========================================================================
  struct PS32aux 
  {
    PS32aux (  const Ostap::Math::PhaseSpace3* ps ) : m_ps ( ps ) {} ;
    PS32aux () = delete ;
    double operator () ( const double x ) const { return  m_ps->ps2_aux( x ) ; }
    const Ostap::Math::PhaseSpace3* m_ps ;
  } ;
  // ==========================================================================
}
// ============================================================================
double Ostap::Math::PhaseSpace3::operator () ( const double x ) const
{
  //
  if ( x <= lowEdge() ) { return 0 ; }
  //
  /// set the temporary mass
  m_tmp = x ;
  //
  // make integral of ps2_aux from m_m1+m_m2 till x-m_m3
  //
  const double low  = m_m1 + m_m2 ;
  const double high = x    - m_m3 ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PS32aux> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpace3/2aux)" ;
  //
  const PS32aux aux { this } ;
  const auto F = s_integrator.make_function ( &aux ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low  , high         ,          // low & high edges
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
// helper function to get the phase space as
// ============================================================================
double Ostap::Math::PhaseSpace3::ps2_aux
( const double m12 ) const
{
  //
  if ( m_tmp <= lowEdge()    ) { return 0 ; }
  //
  if ( m12   <= m_m1  + m_m2 ) { return 0 ; }
  if ( m12   >= m_tmp - m_m3 ) { return 0 ; }
  //
  // represent 3-body phase space as extention of 2-body phase space
  return  m12 / M_PI *
    Ostap::Math::PhaseSpace2::phasespace ( m12   , m_m1 , m_m2 , m_l1 ) *
    Ostap::Math::PhaseSpace2::phasespace ( m_tmp , m12  , m_m3 , m_l2 ) ;
  //
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::PhaseSpace3::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  if ( lowEdge() >= high  ) { return 0 ; }
  if ( lowEdge() >  low   ) { return integral ( lowEdge() , high ) ; }

  //
  if ( 0 < lowEdge() && 5 * lowEdge() < ( high - low ) )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpace3> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpace3)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
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
// constructor from threshold and number of particles
// ============================================================================
Ostap::Math::PhaseSpaceLeft::PhaseSpaceLeft
( const double         threshold ,
  const unsigned short num       )
  : m_threshold ( std::abs ( threshold ) )
  , m_num       ( num )
{}
// ============================================================================
// constructor from list of masses
// ============================================================================
Ostap::Math::PhaseSpaceLeft::PhaseSpaceLeft
( const std::vector<double>& masses )
  : m_threshold ( 0              )
  , m_num       ( masses.size()  )
{
  //
  for ( std::vector<double>::const_iterator im = masses.begin() ;
        masses.end() != im ; ++im )
  { m_threshold += std::abs ( *im ) ; }
  //
}
// ============================================================================
// desctructor
// ============================================================================
Ostap::Math::PhaseSpaceLeft::~PhaseSpaceLeft(){}
// ============================================================================
// evaluate N-body phase space near left threhsold
// ============================================================================
double Ostap::Math::PhaseSpaceLeft::operator () ( const double x ) const
{
  //
  if ( m_threshold >= x ) { return 0 ; }
  //
  return std::pow ( x - m_threshold , 3 * 0.5 * m_num - 5 * 0.5  ) ;
}
// ============================================================================
double Ostap::Math::PhaseSpaceLeft::integral 
( const double xmin , const double xmax ) const 
{
  //
  if      ( s_equal ( xmin , xmax ) ) { return  0 ; }
  else if (           xmin > xmax   ) { return -1 * integral ( xmax , xmin ) ; }
  else if ( xmax <= m_threshold     ) { return  0 ; }
  //
  const double xlow   = std::max ( xmin , m_threshold ) ;
  const double xhigh  = std::max ( xmax , m_threshold ) ;
  //
  const double n      =  ( 3 * m_num - 5 ) * 0.5 ;
  //
  const double tlow   = xlow  - m_threshold ;
  const double thigh  = xhigh - m_threshold ;
  //
  return ( std::pow ( thigh , n + 1 ) - 
           std::pow ( tlow  , n + 1 ) ) / ( n + 1 ) ;
}
// ============================================================================
// set the new value for threshold
// ============================================================================
bool Ostap::Math::PhaseSpaceLeft::setThreshold ( const double x )
{
  //
  if ( s_equal ( x , m_threshold ) ) { return false ; } // RETURN
  //
  m_threshold = x ;
  //
  return true ;
  //
}
// ============================================================================
// constructor from threshold and number of particles
// ============================================================================
Ostap::Math::PhaseSpaceRight::PhaseSpaceRight
( const double         threshold ,
  const unsigned short l         ,
  const unsigned short n         )
  : m_threshold ( std::abs ( threshold ) )
  , m_N         ( std::max ( l , n ) )
  , m_L         ( std::min ( l , n ) )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::PhaseSpaceRight::~PhaseSpaceRight (){}
// ============================================================================
// evaluate N-body phase space near right threshold
// ============================================================================
double Ostap::Math::PhaseSpaceRight::operator () ( const double x ) const
{
  //
  if ( m_threshold <= x ) { return 0 ; }
  //
  return std::pow ( m_threshold - x , 1.5 * ( m_N - m_L ) - 1  ) ;
}
// ============================================================================
double Ostap::Math::PhaseSpaceRight::integral 
( const double xmin , const double xmax ) const 
{
  //
  if      ( s_equal ( xmin , xmax ) ) { return  0 ; }
  else if (           xmin > xmax   ) { return -1 * integral ( xmax , xmin ) ; }
  else if ( xmin >= m_threshold     ) { return  0 ; }
  //
  const double xlow   = std::min ( xmin , m_threshold ) ;
  const double xhigh  = std::min ( xmax , m_threshold ) ;
  //
  const double n      = 1.5 * ( m_N - m_L ) - 1 ;
  //
  const double thigh  = m_threshold - xlow ;
  //
  const double tlow   = m_threshold - xhigh ;
  //
  return ( std::pow ( thigh , n + 1 ) - 
           std::pow ( tlow  , n + 1 ) ) / ( n + 1 ) ;
}
// ============================================================================
// set the new value for threshold
// ============================================================================
bool Ostap::Math::PhaseSpaceRight::setThreshold ( const double x )
{
  //
  if ( s_equal ( x , m_threshold ) ) { return false ; } // RETURN
  //
  m_threshold = x ;
  //
  return true ;
  //
}

// ============================================================================
// constructor from thresholds and number of particles
// ============================================================================
Ostap::Math::PhaseSpaceNL::PhaseSpaceNL
( const double         threshold1 ,
  const double         threshold2 ,
  const unsigned short l          ,
  const unsigned short n          )
  : m_threshold1 ( std::min ( std::abs ( threshold1 ) , std::abs ( threshold2 ) ) )
  , m_threshold2 ( std::max ( std::abs ( threshold1 ) , std::abs ( threshold2 ) ) )
  , m_N          ( std::max ( l , n ) )
  , m_L          ( std::min ( l , n ) )
  , m_norm       ( 1 )
//
  , m_workspace  ()
//
{
  long double _norm = 1;
  if ( ( 3 * m_N * 0.5 - 3       * 0.5 ) < GSL_SF_GAMMA_XMAX &&
       ( 3 * m_L * 0.5 - 3       * 0.5 ) < GSL_SF_GAMMA_XMAX &&
       ( 3 * m_N * 0.5 - 3 * m_L * 0.5 ) < GSL_SF_GAMMA_XMAX )
  {
    //m_norm  = gsl_sf_gamma   ( 3 * m_N * 0.5 - 3       * 0.5 ) ;
    //m_norm /= gsl_sf_gamma   ( 3 * m_L * 0.5 - 3       * 0.5 ) ;
    //m_norm /= gsl_sf_gamma   ( 3 * m_N * 0.5 - 3 * m_L * 0.5 ) ;
    _norm  = std::tgamma   ( 3.0L * m_N * 0.5 - 3       * 0.5 ) ;
    _norm /= std::tgamma   ( 3.0L * m_L * 0.5 - 3       * 0.5 ) ;
    _norm /= std::tgamma   ( 3.0L * m_N * 0.5 - 3 * m_L * 0.5 ) ;
  }
  else
  {
    //_norm  = gsl_sf_lngamma ( 3 * m_N * 0.5 - 3       * 0.5 ) ;
    //_norm -= gsl_sf_lngamma ( 3 * m_L * 0.5 - 3       * 0.5 ) ;
    //_norm -= gsl_sf_lngamma ( 3 * m_N * 0.5 - 3 * m_L * 0.5 ) ;
    //_norm  = gsl_sf_exp     ( m_norm ) ;
    _norm  = std::lgamma ( 3.0L * m_N * 0.5 - 3       * 0.5 ) ;
    _norm -= std::lgamma ( 3.0L * m_L * 0.5 - 3       * 0.5 ) ;
    _norm -= std::lgamma ( 3.0L * m_N * 0.5 - 3 * m_L * 0.5 ) ;
    _norm  = std::exp    ( _norm ) ;
  }
  m_norm = _norm  ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::PhaseSpaceNL::~PhaseSpaceNL() {}
// ============================================================================
// evaluate N/L-body phase space
// ============================================================================
double Ostap::Math::PhaseSpaceNL::operator () ( const double x ) const
{
  //
  if ( m_threshold1 >= x ) { return 0 ; }
  if ( m_threshold2 <= x ) { return 0 ; }
  //
  const double y = (  x - m_threshold1 ) / ( m_threshold2 - m_threshold1 ) ;
  if ( 0 >= y || 1 <= y )  { return 0 ; }
  //
  return m_norm
    / std::abs ( m_threshold2 - m_threshold1               )
    * std::pow (     y , 3 * 0.5 *   m_L         - 5 * 0.5 )
    * std::pow ( 1 - y , 3 * 0.5 * ( m_N - m_L ) - 1       ) ;
}
// =======================================================================
// set the thresholds
// =======================================================================
bool Ostap::Math::PhaseSpaceNL::setThresholds
( const double mn ,
  const double mx )
{
  const double v1 = std::min ( std::abs ( mn ) ,std::abs ( mx ) ) ;
  const double v2 = std::max ( std::abs ( mn ) ,std::abs ( mx ) ) ;
  //
  if ( s_equal ( v1 , m_threshold1 ) &&
       s_equal ( v2 , m_threshold2 ) ) { return false ; }
  //
  m_threshold1 = v1 ;
  m_threshold2 = v2 ;
  //
  return true ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::PhaseSpaceNL::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  if ( m_threshold2 <= low  ) { return 0 ; }
  if ( m_threshold1 >= high ) { return 0 ; }
  //
  if ( m_threshold1 >  low  ) { return integral ( m_threshold1 ,  high        ) ; }
  if ( m_threshold2 <  high ) { return integral ( low          , m_threshold2 ) ; }
  //
  // split, if the interval is too large
  //
  const double width = 0.2 * std::abs  ( m_threshold2 - m_threshold1 ) ;
  if ( 0 < width &&  width < high - low  )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpaceNL> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpaceNL)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low   , high      ,            // low & high edges
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
// get the integral
// ============================================================================
double  Ostap::Math::PhaseSpaceNL::integral() const
{ return integral ( m_threshold1 , m_threshold2 ) ; }
// ============================================================================



// ============================================================================
/*  constructor from four masses and angular momenta
 *  @param m1 the mass of the first  particle
 *  @param m2 the mass of the second particle
 *  @param m3 the mass of the third  particle
 *  @param m4 the mass of the mother particle (m4>m1+m2+m3)
 *  @param L  the angular momentum between the first pair and
 *  the third particle
 *  @param l  the angular momentum between the first and the second particle
 */
// ============================================================================
Ostap::Math::PhaseSpace23L::PhaseSpace23L
( const double         m1 ,
  const double         m2 ,
  const double         m3 ,
  const double         m  ,
  const unsigned short L  ,
  const unsigned short l  )
  : m_m1   ( std::abs ( m1 ) )
  , m_m2   ( std::abs ( m2 ) )
  , m_m3   ( std::abs ( m3 ) )
  , m_m    ( std::abs ( m  ) )
  , m_l    (            l    )
  , m_L    (            L    )
//
  , m_norm ( -1 )
//
  , m_workspace  ()
//
{
  m_norm = integral() ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::PhaseSpace23L::~PhaseSpace23L() {}
// get the momentum of 1st particle in rest frame of (1,2)
// ============================================================================
double Ostap::Math::PhaseSpace23L::q ( const double x ) const
{ return Ostap::Math::PhaseSpace2::q ( x , m_m1 , m_m2 ) ; }
// ============================================================================
// get the momentum of 3rd particle in rest frame of mother
// ============================================================================
double Ostap::Math::PhaseSpace23L::p ( const double x ) const
{ return Ostap::Math::PhaseSpace2::q ( m_m , x , m_m3 ) ; }
// ============================================================================
// calculate the phase space
// ============================================================================
double Ostap::Math::PhaseSpace23L::operator () ( const double x ) const
{ return ps23L( x ) ; }
// ============================================================================
// calculate the phase space
// ============================================================================
double Ostap::Math::PhaseSpace23L::ps23L ( const double x ) const
{
  //
  if ( lowEdge() >= x || highEdge() <= x ) { return  0 ; }
  //
  // represent 3-body phase space as extention of 2-body phase space
  double ps =  x / M_PI *
    Ostap::Math::PhaseSpace2::phasespace ( x   , m_m1 , m_m2 , m_l  ) *
    Ostap::Math::PhaseSpace2::phasespace ( m_m ,    x , m_m3 , m_L  ) ;
  //
  return 0 < m_norm ? ps / m_norm : ps ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::PhaseSpace23L::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  if ( high <= lowEdge  () ) { return 0 ; }
  if ( low  >= highEdge () ) { return 0 ; }
  //
  if ( low  <  lowEdge  () ) { return integral ( lowEdge() , high        ) ; }
  if ( high >  highEdge () ) { return integral ( low       , highEdge () ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpace23L> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpace23L)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low   , high      ,            // low & high edges
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
// get the integral
// ============================================================================
double  Ostap::Math::PhaseSpace23L::integral () const
{ return integral ( lowEdge () , highEdge() ) ; }
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
