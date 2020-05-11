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
#include "Ostap/Dalitz.h"
#include "Ostap/Kinematics.h"
#include "Ostap/MoreMath.h"
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_sf_gamma.h"
// ============================================================================
// Local
// ============================================================================
#include "local_gsl.h"
#include "local_math.h"
#include "local_hash.h"
#include "Integrator1D.h"
#include "Exception.h"
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
      m_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================a
//  set the first mass
// ============================================================================a
bool Ostap::Math::PhaseSpace2::setM1 ( const double value )
{
  const double a = std::abs ( value ) ;
  if ( s_equal ( a , m_m1 ) ) { return false ; }
  m_m1 = a ;
  return true ;
}
// ============================================================================a
//  set the second mass
// ============================================================================a
bool Ostap::Math::PhaseSpace2::setM2 ( const double value )
{
  const double a = std::abs ( value ) ;
  if ( s_equal ( a , m_m2 ) ) { return false ; }
  m_m2 = a ;
  return true ;
}
// ============================================================================a
/* get the mass for the given momentum
 *  \f$ m = \sqrt{m_1^2+q^2} + \sqrt{m_2^2+q^2}\f$
 */
// ============================================================================a
double Ostap::Math::PhaseSpace2::m_ ( const double q ) const
{
  if ( q <=-0 ) { return 0 ; }
  const  double q2 = q * q ;
  return
    m_m1 == m_m2 ? 2 * std::sqrt ( m_m1 * m_m1 + q2 ) :
    std::sqrt ( m_m1 * m_m1 + q2 ) + std::sqrt ( m_m2 * m_m2 + q2) ;  
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
std::size_t Ostap::Math::PhaseSpace2::tag() const 
{ return std::hash_combine ( m_m1 , m_m2 ) ; }
// ============================================================================
/*  calculate the phase space for   m -> m1 + m2
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
  if ( m <= m1 + m2               ) { return 0 ; } // RETURN
  //
  const double msq = m * m ;
  const double lam = Ostap::Kinematics::triangle ( msq  , m1 * m1 , m2 * m2 ) ;
  //
  static const double s_inv8pi = 1.0 / ( 8 * M_PI ) ;
  //
  return 0 < lam ?
    s_inv8pi * Ostap::Math::POW ( std::sqrt ( lam ) / msq , 2 * L + 1 ) : 0.0 ;
}
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
{ return Ostap::Kinematics::q ( m , m1 , m2 ) ; }
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
  const double lam = Ostap::Kinematics::triangle ( m * m , m1 * m1 , m2 * m2 ) ;
  //
  return
    0 <= lam ?
    std::complex<double> (     0.5  * std::sqrt (  lam ) / m , 0 ) :
    std::complex<double> ( 0 , 0.5  * std::sqrt ( -lam ) / m     ) ;
}
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
// ===============================================================
/*  constructor from three masses
 *  @param l1 the angular momentum between 1st and 2nd particle
 *  @param l2 the angular momentum between the pair and 3rd particle
 */
// ===============================================================
Ostap::Math::PhaseSpace3::PhaseSpace3
( const PhaseSpace3s&  ps3 , 
  const unsigned short l1  ,
  const unsigned short l2  ) 
  : m_m1 ( ps3.m1 () ) 
  , m_m2 ( ps3.m2 () ) 
  , m_m3 ( ps3.m3 () ) 
  , m_l1  ( l1 )
  , m_l2  ( l2 )
  , m_tmp ( 0  )   
{}
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
/*  evaluate 3-body phase space
 *  \f[ R_3 ( M ) = \frac{pi^2}{4M^2}\int_{m2+m3}^{M-m_1} \drac{ds_2}{s_2}
 *   \lambda^{1/2}\left ( s_2 , M^2   , m_1^2\right) 
 *   \lambda^{1/2}\left ( s_2 , m_2^2 , m_3^2\right) 
 *  \f] 
 *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
 *              London, New York, Sydney, Toronto, 1973, Eq. (V.2.17)
 */
double Ostap::Math::PhaseSpace3::evaluate  ( const double x ) const
{
  //
  if ( x <= lowEdge() ) { return 0 ; }
  //
  static const double s_norm = 0.25 *  M_PI * M_PI ;
  //
  const double norm = s_norm / ( x * x ) ;
  //
  // all masses are zero 
  if ( s_zero ( m_m1 ) && s_zero ( m_m2 ) && s_zero ( m_m3 ) ) { return norm ; }
  //
  /// set the temporary mass
  m_tmp = x ;
  //
  // make integral of ps2_aux from m_m1 + m_m2 till  x - m_m3
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
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( std::hash_combine ( tag () , x ) , 
      &F     , 
      low    , high ,                // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      m_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result * norm ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PhaseSpace3::tag ()  const 
{ return std::hash_combine ( m_m1 , m_m2 , m_m3 , m_l1 , m_l2 ) ; }
// ============================================================================
// helper function to get the phase space as
// ============================================================================
double Ostap::Math::PhaseSpace3::ps2_aux
( const double m12 ) const
{
  //
  if ( m_tmp <= lowEdge()    ) { return 0 ; }
  if ( m12   <= m_m1  + m_m2 ) { return 0 ; }
  if ( m12   >= m_tmp - m_m3 ) { return 0 ; }
  //
  // represent 3-body phase space as extention of 2-body phase space
  //
  const double lam1 = Ostap::Kinematics::triangle 
    ( m12   * m12   , m_m1 * m_m1 , m_m2 * m_m2  ) ;
  if ( lam1 <= 0 ) { return 0 ; }
  //
  const double lam2 = Ostap::Kinematics::triangle 
    ( m_tmp * m_tmp ,  m12 * m12  , m_m3 * m_m3  ) ;
  if ( lam2 <= 0 ) { return 0 ; }
  //
  //
  /** True integral is
   *  \f[ \int_{(m_1+m_2)^2}^{(M-m_3)^2} \frac{ds_1}{s_1}
   *  \lambda^{1/2}\left ( M^2 , s_1   , m_3^2\right)
   *  \lambda^{1/2}\left ( s_1 , m_1^2 , m_2^2\right) \f] 
   *  It is rewritten as 
   *  \f[ \int_{m_1+m_2}^{M-m_3} \frac{2m_{12}dm_{12}}{m_{12}^2}
   *  \lambda^{1/2}\left ( M^2      , m_{12}^2   , m_3^2\right)
   *  \lambda^{1/2}\left ( m_{12}^2 , m_1^2      , m_2^2\right) \f] 
   *  then \f$ \lambda^{1/2}\f$ are written as  \f$q\f$:
   *  \f$ q = \frac{\lambda^{1/2}\left(s,m_a^2,m_b^2\right)}{2s}\$ :
   *  \f[ \int_{m_1+m_2}^{M-m_3} \frac{2m_{12}dm_{12}}{m_{12}^2}
   *  2M      q (M\rightarrow m_{12} m_3)  
   *  2m_{12} q ( m_{12}\rightarrow m_{1} m_3) = 
   *  \int_{m_1_m_2}^{M-m_3} 8 M dm_12 
   *  q (M\rightarrow m_{12} m_3)  
   *  q ( m_{12}\rightarrow m_{1} m_3) = \f] 
   *  and as last step all q-s are exponentiated :
   *  \f[ \int_{m_1+m_2}^{M-m_3}  8M d m_{12}
   *  q^{2L_1+1}(M\rightarrow m_{12} m_3)  
   *  q^{2L_2+1}( m_{12}\rightarrow m_{1} m_3) \f]
   */
  // 
  const double q1 = std::sqrt ( lam1 ) / ( 2 * m12   ) ;
  const double q2 = std::sqrt ( lam2 ) / ( 2 * m_tmp ) ;
  //
  return 8 * std::pow ( q1 , 2 * m_l1 + 1 ) * std::pow ( q2 , 2 * m_l2 + 1 ) * m_tmp ;
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
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag () , 
      &F     , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      m_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/** constructor from three masses
 *  @param m1 the mass of the first  particle
 *  @param m2 the mass of the second particle
 *  @param m3 the mass of the third  particle
 */
// ============================================================================
Ostap::Math::PhaseSpace3s::PhaseSpace3s
( const double m1 ,
  const double m2 ,
  const double m3 ) 
  : m_m1  ( std::abs ( m1 ) ) 
  , m_m2  ( std::abs ( m2 ) ) 
  , m_m3  ( std::abs ( m3 ) ) 
{}
// ==============================================================================
// evaluate 3-body phase space
// ==============================================================================
double Ostap::Math::PhaseSpace3s::phasespace
( const double x  , 
  const double m1 , 
  const double m2 , 
  const double m3 )
{ return Ostap::Kinematics::phasespace3 ( x , m1 , m2 , m3 ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PhaseSpace3s::tag ()  const 
{ return std::hash_combine ( m_m1 , m_m2 , m_m3 ) ; }
// ==============================================================================
// evaluate 3-body phase space
// ==============================================================================
double Ostap::Math::PhaseSpace3s::evaluate ( const double x ) const
{ return x <= lowEdge () ? 0.0 : phasespace ( x , m_m1 , m_m2 , m_m3 ) ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::PhaseSpace3s::integral 
( const double low  ,
  const double high ) const 
{
  //
  if ( s_equal  ( low , high ) ) { return 0 ; }
  else if ( low  > high )        { return - integral (   high , low ) ; }
  else if ( high <= lowEdge()  ) { return 0 ; }
  //
  const double xlow  = low  ;
  const double xhigh = high ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpace3s> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpace3s)" ;
  //
  const auto F = s_integrator.make_function ( this) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag () , 
      &F     , 
      xlow   , xhigh ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      m_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
// constructor from threshold and number of particles
// ============================================================================
Ostap::Math::PhaseSpaceLeft::PhaseSpaceLeft
( const double         threshold ,
  const unsigned short num       , 
  const double         scale     ) 
  : m_threshold ( std::abs ( threshold ) )
  , m_num       ( num   )
  , m_scale     ( scale ) 
  , m_ps2       () 
{
  Ostap::Assert ( 2 <= m_num , 
                  "Invalid number of particles" , 
                  "Ostap::Math::PhaseSpaceLeft" ) ;
}
// ============================================================================
// constructor from the list of masses
// ============================================================================
Ostap::Math::PhaseSpaceLeft::PhaseSpaceLeft
( const std::vector<double>& masses ,
  const double               scale  )
  : m_threshold (  0            )
  , m_num       ( masses.size() )
  , m_scale     ( scale ) 
  , m_ps2       () 
{
  Ostap::Assert ( 2 <= m_num                    , 
                  "Invalid number of particles" , 
                  "Ostap::Math::PhaseSpaceLeft" ) ;
  if ( 2 == m_num ) 
  {
    m_ps2.setM1 ( std::abs ( masses[0] ) ) ;
    m_ps2.setM2 ( std::abs ( masses[1] ) ) ;
    m_threshold = m_ps2.lowEdge() ;
  }
  else 
  {
    for ( std::vector<double>::const_iterator im = masses.begin() ;
          masses.end() != im ; ++im )
    { m_threshold += std::abs ( *im ) ; }
  }
}
// ============================================================================
// special case: true 2-body phasespace 
// ============================================================================
Ostap::Math::PhaseSpaceLeft::PhaseSpaceLeft
( const Ostap::Math::PhaseSpace2& ps2   , 
  const double                    scale ) 
  : m_threshold ( ps2.m1()  + ps2.m2()  ) 
  , m_num       ( 0       ) // ATTENTION HERE!! 
  , m_scale     ( scale   ) 
  , m_ps2       ( ps2     ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::PhaseSpaceLeft::~PhaseSpaceLeft(){}
// ============================================================================
// evaluate N-body phase space near left threshold
// ============================================================================
double Ostap::Math::PhaseSpaceLeft::operator () ( const double x ) const
{
  //
  const double  t = threshold ()  ;
  if ( t >= x ) { return 0 ; }
  //
  const double y = t + m_scale * ( x - t ) ;
  //
  if ( 0 == m_num ) { return m_ps2 ( y ) ; }
  //
  return std::pow ( ( y - t ) / y , 3 * 0.5 * m_num - 5 * 0.5  ) ;
}
// ============================================================================
double Ostap::Math::PhaseSpaceLeft::integral 
( const double xmin , const double xmax ) const 
{
  //
  const double t = threshold () ;
  //
  if      ( s_equal ( xmin , xmax ) ) { return  0 ; }
  else if (           xmin > xmax   ) { return -1 * integral ( xmax , xmin ) ; }
  else if ( xmax <= t               ) { return  0 ; }
  //
  const double xlow  = std::max ( xmin , t ) ;
  const double xhigh = std::max ( xmax , t ) ;
  //
  if ( 0 == m_num ) { return m_ps2.integral ( xlow , xhigh ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpaceLeft> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpaceLeft)" ;
  //
  const auto F = s_integrator.make_function ( this) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag () , 
      &F     , 
      xlow   , xhigh ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      m_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
  // const double n      =  ( 3 * m_num - 5 ) * 0.5 ;
  // //
  // const double tlow   =  ( xlow  - t ) ;
  // const double thigh  =  ( xhigh - t ) ;
  // //
  // return ( std::pow ( thigh , n + 1 ) - 
  //          std::pow ( tlow  , n + 1 ) ) / ( n + 1 ) * std::pow ( m_scale , n ) ;
}
// ============================================================================
// set the new value for scale
// ============================================================================
bool Ostap::Math::PhaseSpaceLeft::setScale  ( const double value  )
{
  const double a = std::abs  ( value ) ;
  if ( s_equal ( a , m_scale ) ) { return false ; } //  RETURN
  m_scale = a ;
  return true ;  
}
// ============================================================================
// set the new value for threshold
// ============================================================================
bool Ostap::Math::PhaseSpaceLeft::setThreshold ( const double value )
{
  const double a = std::abs  ( value ) ;
  const double t = threshold (       ) ;
  if ( s_equal ( a , t ) ) { return false ; } // RETURN
  ///
  if ( 0 == m_num ) 
  {
    m_ps2.setM1 ( m_ps2.m1 () * ( a / t ) ) ;
    m_ps2.setM2 ( m_ps2.m2 () * ( a / t ) ) ;
    return true  ;
  }
  m_threshold = a ;
  return true ;
}
// ============================================================================
// get the tag  
// ============================================================================
std::size_t Ostap::Math::PhaseSpaceLeft::tag () const  // get the tag
{ return std::hash_combine ( m_threshold , m_num , m_scale , m_ps2.tag () ) ; }
// ============================================================================

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
  const double thigh  = m_threshold - xlow ;
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
  if ( s_equal ( x , m_threshold ) ) { return false ; } // RETURN
  m_threshold = x ;
  return true ;
}
// ============================================================================
// get the tag  
// ============================================================================
std::size_t Ostap::Math::PhaseSpaceRight::tag () const  // get the tag
{ return std::hash_combine ( m_threshold , m_N , m_L  ) ; }
// ============================================================================

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
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag ()        ,
      &F            , 
      low           , high      ,    // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      m_workspace.size () ,          // size of workspace
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
// get the tag  
// ============================================================================
std::size_t Ostap::Math::PhaseSpaceNL::tag () const  // get the tag
{ return std::hash_combine ( m_L , m_N , m_threshold1 , m_threshold2 ) ; }
// ============================================================================

// ============================================================================
// constructor from all masses
// ============================================================================
Ostap::Math::PSDalitz::PSDalitz
( const double M  , 
  const double m1 , 
  const double m2 , 
  const double m3 ) 
  : PSDalitz (  Ostap::Kinematics::Dalitz ( M , m1 , m2 , m3 ) )
{}
// ============================================================================
// constructor from Dalizt plot 
// ============================================================================
Ostap::Math::PSDalitz::PSDalitz
( const Ostap::Kinematics::Dalitz& dalitz )
  : m_dalitz     ( dalitz ) 
  , m_norm       ( -1     )
  , m_workspace  (        )
{
  m_norm = 1.0 / integral() ;
}
// ============================================================================
/*  get the value of PDF 
 *  @see Ostap::Kinematics::Dalitz::dRdm12 
 */
// ============================================================================
double Ostap::Math::PSDalitz::operator () ( const double x ) const 
{ return ( 0 < m_norm ? m_norm : 1.0 ) * m_dalitz.dRdm12 ( x ) ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::PSDalitz::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0        ; } // RETURN
  if (           low > high   ) { return - integral ( high , low  ) ; } // RETURN
  //
  const double x_min = xmin () ;
  const double x_max = xmax () ;
  //
  if ( low >= x_max || high <= x_min ) { return 0 ; }
  //
  const double xlow  = std::max ( low  , x_min ) ;
  const double xhigh = std::min ( high , x_max ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PSDalitz> s_integrator {} ;
  static char s_message[] = "Integral(PSDalitz)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( tag ()        ,
      &F            , 
      xlow          , xhigh     ,    // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      m_workspace.size () ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double  Ostap::Math::PSDalitz::integral() const
{ return 0 < m_norm ? 1.0 / m_norm : integral ( m1() + m2() , M() - m3 () ) ; }
// ============================================================================
// get the tag  
// ============================================================================
std::size_t Ostap::Math::PSDalitz::tag () const  // get the tag
{ return std::hash_combine ( M () , m1 () , m2 () , m3 () ) ; }
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
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache 
    ( tag () , 
      &F     , 
      low    , high      ,           // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION         ,          // relative precision
      m_workspace.size () ,          // size of workspace
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
// get the tag  
// ============================================================================
std::size_t Ostap::Math::PhaseSpace23L::tag () const  // get the tag
{ return std::hash_combine ( m_m1 , m_m2 , m_m3 , m_m , m_l , m_L ) ; }
// ============================================================================



// ============================================================================
/*  Get a full integrated phase space over Dalitz plot 
 *  \f$  R(s) = \int \int R(s_1,s_2) \deriv s_1 \deriv s_2 =
 *  \int _{(m_2+m_3)^2}^{ (\sqrt{s}-m_1)^2}
 *   \frac{\deriv s_2}{s_2}
 *   \lambda^{1/2}(s_2,s,m_1^2)
 *   \lambda^{1/2}(s_2,m_2^2,m_3^2)\f$ 
 */
// ============================================================================
double Ostap::Kinematics::phase_space 
( const Ostap::Kinematics::Dalitz& dalitz )  
{ return Ostap::Math::PSDalitz ( dalitz ).phasespace () ; }
// ============================================================================
//                                                                      The END 
// ============================================================================
