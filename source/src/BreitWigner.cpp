// ============================================================================
// Include files
// ============================================================================
// STD & STL  
// ============================================================================
#include <cmath>
#include <array>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PhaseSpace.h"
#include "Ostap/BreitWigner.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/MoreMath.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "local_gsl.h"
#include "Integrator1D.h"
// ============================================================================
/** @file 
 *  implementation of useful models for describing signal peaks with the natural width \
 *  - Breit-Wigner
 *  - Flatte 
 *  - LASS  (kappa) 
 *  - Bugg  (sigma-pole)
 *  - Gounaris-Sakurai
 *  - Voight &Co 
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  // Breit-Wigner & Co 
  // ==========================================================================
  /// get the complex Breit amplitude
  std::complex<double> breit_amp
  ( const double x     ,
    const double m0    ,
    const double gamma )
  {
    //
    static const std::complex<double> s_j ( 0 , 1 ) ;
    //
    const std::complex<double> v = m0 * m0 - x * x - s_j * m0 * gamma ;
    //
    // attention: normalization factors and phase space are here!
    //
    // const double d = 2 / M_PI ;
    // const double d = 2 * std::abs ( m0 * gamma  * x ) / M_PI ;
    //
    return  1.0 / v ;
  }
  // ==========================================================================
  //// calculate the current width
  double gamma_run 
  ( const double gam0    ,
    const double x       ,
    const double m1      ,
    const double m2      ,
    const double m0      ,
    const unsigned int L ,
    const Ostap::Math::FormFactor* fun  = 0 )
  {
    //
    if ( m1 + m2 >= x ) { return 0 ; }   // RETURN
    //
    const double q  = Ostap::Math::PhaseSpace2::q ( x  , m1 , m2 ) ;
    const double q0 = Ostap::Math::PhaseSpace2::q ( m0 , m1 , m2 ) ;
    //
    if ( 0 >= q || 0 >= q0 ) { return 0 ; }  // RETURN
    //
    const double r  = 0 != fun ? (*fun) ( x  , m0 , m1 , m2 ) : 1.0 ;
    const double r0 = 0 != fun ? (*fun) ( m0 , m0 , m1 , m2 ) : 1.0 ;
    //
    if ( 0 >= r0 )           { return 0 ; }  // RETURN
    //
    return gam0 * Ostap::Math::pow ( q / q0 , 2 * L + 1 ) * ( r / r0 ) ;
  }
  // ==========================================================================
  /** @var s_BUKIN
   *  useful constant (neded for pseudo-Vogt)
   *  \f$ \sqrt{ 2 \log 2 } \f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  constexpr double s_BUKIN   = std::sqrt ( 2.0 * std::log ( 2.0 ) ) ;
  // ==========================================================================
} //                                            The end of  anonymous namespace 
// ============================================================================
// Rho-functions from Jackson
// ============================================================================
/* the simplest function: constant
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_0
( double /* m  */ ,
  double /* m0 */ ,
  double /* m1 */ ,
  double /* m2 */ ) { return 1 ; }
// ============================================================================
/* the simple function for \f$ 1^- \rightarrow 0^- 0^- \f$, l = 1
 *  \f$\rho(\omega)= \omega^{-1}\f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A2
( double    m     ,
  double /* m0 */ ,
  double /* m1 */ ,
  double /* m2 */ ) { return 1./m ; }
// ============================================================================
/*  the simple function for \f$ 1^- \rightarrow 0^- 1^- \f$, l = 1
 *  \f$\rho(\omega)= \omega \f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A3
( double    m     ,
  double /* m0 */ ,
  double /* m1 */ ,
  double /* m2 */ ) { return m ; }
// ============================================================================
/*  the simple function for
 *  \f[ \frac{3}{2}^+ \rightarrow \frac{1}{2}^+ 0^- \f], l = 1
 *  $\rho(\omega)= \frac{ ( \omega + M )^2 - m^2 }{ \omega^2} \f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @param m the invariant mass
 *  @param m1 the invariant mass of the first  (spinor) particle
 *  @param m2 the invariant mass of the secodn (scalar) particle
 *  @return the value of rho-function
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A4
( double    m     ,
  double /* m0 */ ,
  double    m1    ,
  double    m2    )
{
  const double a = m + m1 ;
  //
  return ( a * a  - m2 * m2 ) / ( m * m ) ;
}
// ============================================================================
/*  the simple function for
 *  \f$ \frac{3}{2}^- \rightarrow \frac{1}{2}^+ 0^- \f$, l = 2
 *  $\rho(\omega)= \left[ ( \omega + M )^2 - m^2 \right]^{-1} \f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @param m the invariant mass
 *  @param m1 the invariant mass of the first  (spinor) particle
 *  @param m2 the invariant mass of the second (scalar) particle
 *  @return the value of rho-function
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A5
( double    m     ,
  double /* m0 */ ,
  double    m1    ,
  double    m2    )
{
  const double a = m + m1 ;
  //
  return 1 / ( a * a  - m2 * m2 ) ;
}
// ============================================================================
/*  the simple function for \f$\rho^- \rightarrow \pi^+ \pi^-\f$
 *  \f$ 1- \rightarrow 0^- 0^- \f$, l = 1
 *  $\rho(\omega)= \left[ q_0^2 + q^2 \right]^{-1}f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @param m  the invariant mass
 *  @param m0 the nominal   mass
 *  @param m1 the invariant mass of the first  particle
 *  @param m2 the invariant mass of the second particle
 *  @return the value of rho-function
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A7
( double    m   ,
  double    m0  ,
  double    m1  ,
  double    m2  )
{
  //
  const double q  = Ostap::Math::PhaseSpace2::q ( m  , m1 , m2 ) ;
  const double q0 = Ostap::Math::PhaseSpace2::q ( m0 , m1 , m2 ) ;
  //
  if ( 0 >= q && 0 >= q0 ) { return 1 ; }
  //
  return 1. / ( q * q + q0 * q0 ) ;
}
// ============================================================================
// FORMFACTORS 
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Math::FormFactor::~FormFactor (){}
// ============================================================================
// default constructor
// ============================================================================
Ostap::Math::FormFactors::Jackson::Jackson() 
  : Ostap::Math::FormFactor() 
  , m_rho ( nullptr ) 
{}
// ============================================================================
// constructor from enum
// ============================================================================
Ostap::Math::FormFactors::Jackson::Jackson
( const Ostap::Math::FormFactors::JacksonRho rho ) 
  : Ostap::Math::FormFactor() 
  , m_rho ( nullptr ) 
{
  switch ( rho )
  {
  case   Ostap::Math::FormFactors::Jackson_0  :
    m_rho = &Ostap::Math::Jackson::jackson_0  ; break ;
  case   Ostap::Math::FormFactors::Jackson_A2 :
    m_rho = &Ostap::Math::Jackson::jackson_A2 ; break ;
  case   Ostap::Math::FormFactors::Jackson_A3 :
    m_rho = &Ostap::Math::Jackson::jackson_A3 ; break ;
  case   Ostap::Math::FormFactors::Jackson_A4 :
    m_rho = &Ostap::Math::Jackson::jackson_A4 ; break ;
  case   Ostap::Math::FormFactors::Jackson_A5 :
    m_rho = &Ostap::Math::Jackson::jackson_A5 ; break ;
  case   Ostap::Math::FormFactors::Jackson_A7 :
    m_rho = &Ostap::Math::Jackson::jackson_A7 ; break ;
  default         :
    m_rho = nullptr ; 
  }
  //
}
// ============================================================================
// constructor from function itself 
// ============================================================================
Ostap::Math::FormFactors::Jackson::Jackson
( const Ostap::Math::FormFactors::rho_fun rho ) 
  : Ostap::Math::FormFactor() 
  , m_rho ( rho ) 
{ if ( !m_rho ) { m_rho = &Ostap::Math::Jackson::jackson_0 ; } }
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Math::FormFactors::Jackson::~Jackson(){}
// ============================================================================
// clone method ("virtual constructor")
// ============================================================================
Ostap::Math::FormFactors::Jackson* 
Ostap::Math::FormFactors::Jackson:: clone() const 
{ return new Ostap::Math::FormFactors::Jackson ( *this ) ; }
// ============================================================================
// the only important method 
// ============================================================================
double Ostap::Math::FormFactors::Jackson::operator() 
  ( const double m  , const double m0 ,
    const double m1 , const double m2 ) const
{ 
  return nullptr == m_rho ? 1.0 : (*m_rho)( m , m0 , m1 , m2 ) ; 
}
// ============================================================================
// Blatt-Weisskopf formfactors 
// ============================================================================
namespace
{
  // ==========================================================================
  // Coefficients for Blatt-Weisskopf formfactors 
  // ==========================================================================
  // const std::array<int,1> s_BW_0 { {                               1 } } ;
  // const std::array<int,2> s_BW_1 { {                            1, 1 } } ;
  const std::array<int,3> s_BW_2 { {                        9,  3, 1 } } ;
  const std::array<int,4> s_BW_3 { {                 225,  45,  6, 1 } } ;
  const std::array<int,5> s_BW_4 { {         11025, 1575, 135, 10, 1 } } ;
  const std::array<int,6> s_BW_5 { { 893025, 99225, 6300, 315, 15, 1 } } ;  
  // ==========================================================================
}    
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf::BlattWeisskopf
( const Ostap::Math::FormFactors::BlattWeisskopf::Case L , 
  const double                                         b )
  : Ostap::Math::FormFactor() 
  , m_L ( L ) 
  , m_b ( b )
{
  switch ( L ) 
  {
  case Zero  : break ;
  case One   : break ;
  case Two   : break ;
  case Three : break ;
  case Four  : break ;
  case Five  : break ;
  default:   
    Ostap::throwException( "Illegal Blatt-Weisskopf form factor" , "Math" ) ;
  }
}


// ============================================================================
// default constructor (needed for  serialization)
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf::BlattWeisskopf()
  : Ostap::Math::FormFactor() 
  , m_L ( Ostap::Math::FormFactors::BlattWeisskopf::Zero ) 
  , m_b ( 0.0 )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf::~BlattWeisskopf(){}
// ============================================================================
// clone method ("virtual constructor")
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf*
Ostap::Math::FormFactors::BlattWeisskopf::clone() const 
{ return new Ostap::Math::FormFactors::BlattWeisskopf(*this) ; }
// ============================================================================
// get the barrier factor 
// ============================================================================
double Ostap::Math::FormFactors::BlattWeisskopf::b 
( const double z   , 
  const double z0  ) const 
{
  if ( Zero == m_L || s_equal ( z , z0 ) ) { return 1 ; }
  //
  const long double r2 =
    //
    One   == m_L ? ( 1 + z0 ) / ( 1  + z )  :
    //
    Two   == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_2.rbegin() , s_BW_2.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_2.rbegin() , s_BW_2.rend  () , z  ).first :
    //
    Three == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_3.rbegin() , s_BW_3.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_3.rbegin() , s_BW_3.rend  () , z  ).first :
    //
    Four == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_4.rbegin() , s_BW_4.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_4.rbegin() , s_BW_4.rend  () , z  ).first :
    //
    Five == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_5.rbegin() , s_BW_5.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_5.rbegin() , s_BW_5.rend  () , z  ).first : 
    //
    1.0L ;
  //
  return std::sqrt ( r2 ) ;
}
// ============================================================================
// the only important method 
// ============================================================================
double Ostap::Math::FormFactors::BlattWeisskopf::operator() 
  ( const double m  , const double m0 ,
    const double m1 , const double m2 ) const 
{
  //
  if ( s_equal ( m , m0 ) ) { return    1   ; }
  if ( s_zero  ( m_b )    ) { return m0 / m ; }
  //
  /// get the momenta 
  const double q  = Ostap::Math::PhaseSpace2::q ( m  , m1 , m2 ) ;
  const double q0 = Ostap::Math::PhaseSpace2::q ( m0 , m1 , m2 ) ;
  ///
  const double _z  = q  * m_b ;
  const double _z0 = q0 * m_b ;
  //
  return ( m0 / m ) * b ( _z * _z , _z0 * _z0 ) ;
}
// ============================================================================

// ============================================================================
// Breit-Wigner itself 
// ============================================================================



// ============================================================================
// constructor
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner
( const double         m0   ,
  const double         gam0 ,
  const double         m1   ,
  const double         m2   ,
  const unsigned short L    )
  : m_m0         (             m0    )
  , m_gam0       ( std::abs ( gam0 ) )
  , m_m1         ( std::abs (   m1 ) )
  , m_m2         ( std::abs (   m2 ) )
  , m_L          (              L    )
  , m_formfactor ( nullptr ) 
    //
  , m_workspace  ()
    //
{}
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner
( const double                                m0   ,
  const double                                gam0 ,
  const double                                m1   ,
  const double                                m2   ,
  const unsigned short                        L    ,
  const Ostap::Math::FormFactors::JacksonRho  r    )
  : m_m0         (             m0    )
  , m_gam0       ( std::abs ( gam0 ) )
  , m_m1         ( std::abs (   m1 ) )
  , m_m2         ( std::abs (   m2 ) )
  , m_L          (              L    )
    //
  , m_formfactor ( new Ostap::Math::FormFactors::Jackson ( r ) ) 
    //
  , m_workspace  ()
    //
{}
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner
( const double                   m0   ,
  const double                   gam0 ,
  const double                   m1   ,
  const double                   m2   ,
  const unsigned short           L    ,
  const Ostap::Math::FormFactor& ff   )
  : m_m0         (             m0    )
  , m_gam0       ( std::abs ( gam0 ) )
  , m_m1         ( std::abs (   m1 ) )
  , m_m2         ( std::abs (   m2 ) )
  , m_L          (              L    )
    //
  , m_formfactor ( ff.clone() ) 
    //
  , m_workspace  ()
    //
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner 
( const Ostap::Math::BreitWigner& bw ) 
  : m_m0         ( bw.m_m0    )
  , m_gam0       ( bw.m_gam0  )
  , m_m1         ( bw.m_m1    )
  , m_m2         ( bw.m_m2    )
  , m_L          ( bw.m_L     )
    //
  , m_formfactor ( nullptr == bw.m_formfactor ? nullptr : bw.m_formfactor->clone() ) 
    //
  , m_workspace  ()
    //
{}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner 
( Ostap::Math::BreitWigner&& bw ) 
  : m_m0         ( bw.m_m0    )
  , m_gam0       ( bw.m_gam0  )
  , m_m1         ( bw.m_m1    )
  , m_m2         ( bw.m_m2    )
  , m_L          ( bw.m_L     )
    //
  , m_formfactor ( bw.m_formfactor ) 
    //
  , m_workspace  ()
    //
{
  bw.m_formfactor = nullptr ;
}
// ============================================================================
Ostap::Math::BreitWigner*
Ostap::Math::BreitWigner::clone() const
{ return new BreitWigner ( *this ) ; }  
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::BreitWigner::~BreitWigner ()
{ if ( 0 != m_formfactor ) { delete m_formfactor ; m_formfactor = nullptr ; } }
// ============================================================================
//  calculate the Breit-Wigner amplitude
// ============================================================================
std::complex<double>
Ostap::Math::BreitWigner::amplitude ( const double x ) const
{
  //
  if ( m_m1 + m_m2 >= x ) { return 0 ; }
  //
  const double g  = gamma ( x ) ;
  if ( 0 >= g ) { return 0 ; }
  //
  return std::sqrt ( m0 () * gam0 () ) * breit_amp ( x , m0() , g ) ;
}
// ============================================================================
/*  calculate the Breit-Wigner shape
 *  \f$\frac{1}{\pi}\frac{\omega\Gamma(\omega)}{ (\omega_0^2-\omega^2)^2-\omega_0^2\Gammma^2(\omega)-}\f$
 */
// ============================================================================
double Ostap::Math::BreitWigner::breit_wigner ( const double x ) const
{
  //
  if ( m_m1 + m_m2 >= x ) { return 0 ; }
  //
  const double g  = gamma ( x ) ;
  if ( 0 >= g ) { return 0 ; }
  //
  std::complex<double> a = amplitude ( x ) ;
  //
  return 2 * x * std::norm ( a )* g / gam0() / M_PI ;
  //
  // const double omega2 = m_m0 * m_m0 ;
  // const double delta = omega2        -          x * x ;
  // const double v     = delta * delta + omega2 * g * g ;
  //
  // return 2 * x * m_m0 * g / v / M_PI  ;
}
// ============================================================================
/*  calculate the Breit-Wigner shape
 *  \f$\frac{1}{\pi}\frac{\omega\Gamma(\omega)}{ (\omega_0^2-\omega^2)^2-\omega_0^2\Gammma^2(\omega)-}\f$
 */
// ============================================================================
double Ostap::Math::BreitWigner::operator() ( const double x ) const
{ return breit_wigner ( x ) ; }
// ============================================================================
// calculate the current width
// ============================================================================
double Ostap::Math::BreitWigner::gamma ( const double x ) const
{
  //
  return gamma_run ( m_gam0       ,
                     x            ,
                     m_m1         ,
                     m_m2         ,
                     m_m0         ,
                     m_L          ,
                     m_formfactor ) ;
  //
}
// ===========================================================================
// get the value of formfactor at given m 
// ============================================================================
double Ostap::Math::BreitWigner::formfactor ( const double m ) const 
{ 
  return 
    nullptr == m_formfactor ? 1. : 
    (*m_formfactor)( m , m_m0 , m_m1 , m_m2 ) ; 
}
// ============================================================================


// ============================================================================
bool Ostap::Math::BreitWigner::setM0     ( const double x )
{
  const double v       = std::abs ( x ) ;
  if ( s_equal ( v , m_m0 ) ) { return false ; } // RETURN
  m_m0   = v ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BreitWigner::setGamma0 ( const double x )
{
  const double v       = std::abs ( x ) ;
  if ( s_equal ( v , m_gam0 ) ) { return false ; } // RETURN
  m_gam0  = v ;
  return true ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::BreitWigner::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  if ( m_m1 + m_m2 >= high ) { return                              0   ; }
  if ( m_m1 + m_m2 >  low  ) { return integral  ( m_m1 + m_m2 , high ) ; }
  //
  //
  // split into reasonable sub intervals
  //
  const double x1     = m_m0 - 10 * m_gam0 ;
  const double x2     = m_m0 + 10 * m_gam0  ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  if ( low < x_low  && x_low < high )
  {
    return
      integral (   low , x_low  ) +
      integral ( x_low ,   high ) ;
  }
  if ( low <  x_high && x_high < high )
  {
    return
      integral (   low  , x_high  ) +
      integral ( x_high ,   high  ) ;
  }
  //
  // split, if interval too large
  //
  const double width = std::max ( m_gam0 , 0.0 ) ;
  if ( 0 < width &&  3 * width < high - low  )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<BreitWigner> s_integrator {} ;
  static char s_message[] = "Integral(BreitWigner)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION         ,          // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL :
      s_PRECISION       ,            // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral b
// ============================================================================
double  Ostap::Math::BreitWigner::integral () const
{
  //
  // split into reasonable sub intervals
  //
  const double x1     = std::max ( m_m1  + m_m2 , m_m0 - 10 * m_gam0 ) ;
  const double x2     = std::max ( m_m1  + m_m2 , m_m0 + 10 * m_gam0 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<BreitWigner> s_integrator {} ;
  static char s_message[] = "Integral(BreitWigner/tail)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaqiu_integrate
    ( &F , 
      x_high              ,          // low edge
      workspace ( m_workspace ) ,    // workspace
      s_PRECISION         ,          // absolute precision
      s_PRECISION_TAIL    ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result + integral ( m_m1 + m_m2 , x_high );
}




// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Rho0::Rho0
( const double m0       ,
  const double gam0     ,
  const double pi_mass  )
  : Ostap::Math::BreitWigner ( m0         ,
                               gam0       ,
                               pi_mass    ,
                               pi_mass    ,
                               1          ,
                               Ostap::Math::FormFactors::Jackson_A7 )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Rho0::~Rho0(){}

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Kstar0::Kstar0
( const double m0       ,
  const double gam0     ,
  const double k_mass   ,
  const double pi_mass  )
  : Ostap::Math::BreitWigner ( m0         ,
                               gam0       ,
                               k_mass     ,
                               pi_mass    ,
                               1          ,
                               Ostap::Math::FormFactors::Jackson_A2 )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Kstar0::~Kstar0(){}

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Phi0::Phi0
( const double m0       ,
  const double gam0     ,
  const double k_mass   )
  : Ostap::Math::BreitWigner ( m0         ,
                               gam0       ,
                               k_mass     ,
                               k_mass     ,
                               1          ,
                               Ostap::Math::FormFactors::Jackson_A2 )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Phi0::~Phi0(){}

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Rho0FromEtaPrime::Rho0FromEtaPrime
( const double m0        ,
  const double gam0      ,
  const double pi_mass   ,
  const double eta_prime )
  : Ostap::Math::Rho0 ( m0 , gam0 , pi_mass )
  , m_eta_prime ( std::abs ( eta_prime ) )
{}
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Rho0FromEtaPrime::Rho0FromEtaPrime
( const Ostap::Math::Rho0& rho       ,
  const double             eta_prime )
  : Ostap::Math::Rho0 ( rho )
  , m_eta_prime ( std::abs ( eta_prime ) )
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Math::Rho0FromEtaPrime*
Ostap::Math::Rho0FromEtaPrime::clone() const
{ return new Rho0FromEtaPrime ( *this ) ; }  
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Rho0FromEtaPrime::~Rho0FromEtaPrime(){}
// ============================================================================
// calculate the function
// ============================================================================
double Ostap::Math::Rho0FromEtaPrime::operator() ( const double x ) const
{
  //
  if ( m_eta_prime <= x ) { return 0 ; }
  //
  const double k_gamma = Ostap::Math::PhaseSpace2::q ( m_eta_prime , x , 0 ) ;
  if ( 0 >= k_gamma     ) { return 0 ; }
  //
  const double rho     = breit_wigner ( x ) ;
  if ( 0 >= rho         ) { return 0 ; }
  //
  return rho * Ostap::Math::pow ( 2 * k_gamma / m_eta_prime , 3 ) * 20 ;
  //
}
// ============================================================================
//               Flatte
// ============================================================================
/*  constructor  from three parameters
 *  @param m0    the mass
 *  @param m0g1  parameter \f$ m_0\times g_1\f$
 *  @param g2og2 parameter \f$ g2/g_1       \f$
 */
// ============================================================================
Ostap::Math::Flatte::Flatte
( const double m0    ,
  const double m0g1  ,
  const double g2og1 ,
  const double mA1   ,
  const double mA2   ,
  const double mB1   ,
  const double mB2   )
  : m_m0    ( std::fabs ( m0    ) )
  , m_m0g1  ( std::fabs ( m0g1  ) )
  , m_g2og1 ( std::fabs ( g2og1 ) )
  , m_A1    ( std::fabs ( mA1   ) )
  , m_A2    ( std::fabs ( mA2   ) )
  , m_B1    ( std::fabs ( mB1   ) )
  , m_B2    ( std::fabs ( mB2   ) )
    //
  , m_workspace ()
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Flatte::~Flatte(){}
// ============================================================================
Ostap::Math::Flatte*
Ostap::Math::Flatte::clone() const { return new Ostap::Math::Flatte ( *this ) ; }
// ============================================================================
// get the value of Flatte function
// ============================================================================
double Ostap::Math::Flatte::operator() ( const double x ) const
{ return flatte ( x ) ; }
// ============================================================================
  // get the complex Flatte amplitude
// ============================================================================
std::complex<double> Ostap::Math::Flatte::flatte_amp
( const double x     )  const
{
  //
  const std::complex<double> rho_AA =
    Ostap::Math::PhaseSpace2::q1 ( x , mA1 () , mA2 () )  ;
  const std::complex<double> rho_BB =
    Ostap::Math::PhaseSpace2::q1 ( x , mB1 () , mB2  () )  ;
  //
  static const std::complex<double> s_j ( 0 , 1 ) ;
  //
  const std::complex<double> v =
    m0() * m0 () - x * x - s_j * m0g1() * ( rho_AA + g2og1 () * rho_BB ) ;
  //
  return  1.0 / v ;
}
// ===========================================================================
// get the function for pipi-channel
// ===========================================================================
double Ostap::Math::Flatte::flatte ( const double x ) const
{
  //
  if ( thresholdA () >= x ) { return 0 ; }
  //
  // get the amplitude...
  std::complex<double> amp = flatte_amp ( x ) ;
  //
  const double ps = Ostap::Math::PhaseSpace2::phasespace ( x ,  mA1() , mA2() ) ;
  //
  return x * ps * std::norm ( amp ) * 2 / M_PI * m0g1() ;
}
// ===========================================================================
// get the function for KK-channel
// ===========================================================================
double Ostap::Math::Flatte::flatte2 ( const double x ) const
{
  //
  if ( thresholdB () >= x ) { return 0 ; }
  //
  // get the amplitude...
  std::complex<double> amp = flatte_amp ( x ) ;
  //
  const double ps = Ostap::Math::PhaseSpace2::phasespace ( x ,  mB1() , mB2() ) ;
  //
  return x * ps * std::norm ( amp ) * 2 / M_PI * m0g1 () * g2og1 () ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Flatte::integral
( const double low  ,
  const double high ) const
{
  //
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  const double a = threshold() ;
  if ( a >= high ) { return                     0 ; }
  if ( a >  low  ) { return integral ( a , high ) ; }
  //
  const double b = std::max ( thresholdA () , thresholdB () ) ;
  if ( low < b     && b    < high ) 
  { return integral ( low , b ) + integral ( b , high ) ; }
  //
  if ( low < m_m0  && m_m0 < high ) 
  { return integral ( low , m_m0 ) + integral ( m_m0 , high ) ; }
  //
  const double width =
    0 > m_m0 ? 0.0 :
    std::abs ( m_m0g1 / m_m0           ) +
    std::abs ( m_m0g1 / m_m0 * m_g2og1 ) ;
  //
  for ( unsigned int i = 0 ; ( i < 5 ) && ( 0 < width ) ; ++ i ) 
  {
    const double x1 = m_m0 + i * width ;
    if ( low < x1  && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 - i * width ;
    if ( low < x2  && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  const double x_low  = 0 < width ? m_m0 - 20 * width : low  ;
  const double x_high = 0 < width ? m_m0 + 20 * width : high ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Flatte> s_integrator {} ;
  static char s_message[] = "Integral(Flatte)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION       ,            // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL :
      s_PRECISION       ,            // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral b
// ============================================================================
double  Ostap::Math::Flatte::integral () const
{
  //
  // split into reasonable sub intervals
  //
  const double x_low  = threshold () ;
  const double x_high = m_m0
    + 15 * std::abs ( m_m0g1 / m_m0           )
    + 15 * std::abs ( m_m0g1 / m_m0 * m_g2og1 ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Flatte> s_integrator {} ;
  static char s_message[] = "Integral(Flatte/tail)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaqiu_integrate
    ( &F , 
      x_high   ,                     // high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION         ,          // absolute precision
      s_PRECISION_TAIL    ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result + integral ( x_low , x_high );
}
// ============================================================================
// set mass
// ============================================================================
bool Ostap::Math::Flatte::setM0     ( const double x )
{
  const double v = std::abs ( x ) ;
  //
  if ( s_equal ( v , m_m0 ) ) { return false ; }
  //
  m_m0 = v ;
  //
  return true ;
}
// ============================================================================
// set mass times G1
// ============================================================================
bool Ostap::Math::Flatte::setM0G1   ( const double x )
{
  const double v = std::abs ( x ) ;
  //
  if ( s_equal ( v , m_m0g1 ) ) { return false ; }
  //
  m_m0g1 = v ;
  //
  return true ;
}
// ============================================================================
// set G2 over G1
// ============================================================================
bool Ostap::Math::Flatte::setG2oG1  ( const double x )
{
  const double v = std::abs ( x ) ;
  //
  if ( s_equal ( v , m_g2og1 ) ) { return false ; }
  //
  m_g2og1 = v ;
  //
  return true ;
}
// ============================================================================

// ============================================================================
//               Flatte-2
// ============================================================================
/* constructor  from three parameters
 *  @param m0    the mass
 *  @param m0g1  parameter \f$ m_0\times g_1\f$
 *  @param g2og2 parameter \f$ g2/g_1       \f$
 */
// ============================================================================
Ostap::Math::Flatte2::Flatte2
( const double m0    ,
  const double m0g1  ,
  const double g2og1 ,
  const double mA1   ,
  const double mA2   ,
  const double mB1   ,
  const double mB2   )
  : Ostap::Math::Flatte ( m0 , m0g1 , g2og1 , mA1 , mA2 , mB1 , mB2 )
{}
// ============================================================================
// constructor  from Flatte
// ============================================================================
Ostap::Math::Flatte2::Flatte2
( const Ostap::Math::Flatte& flatte )
  : Ostap::Math::Flatte ( flatte ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Flatte2::~Flatte2(){}
// ============================================================================
Ostap::Math::Flatte2*
Ostap::Math::Flatte2::clone() const { return new Ostap::Math::Flatte2 ( *this ) ; }
// ============================================================================
// get the value of Flatte function
// ============================================================================
double Ostap::Math::Flatte2::operator() ( const double x ) const
{ return flatte2 ( x ) ; }
// ============================================================================

// ============================================================================
// LASS: Kpi S-wave 
// ============================================================================
/*  constructor from all masses and angular momenta
 *  @param m1 the mass of the first  particle
 *  @param m2 the mass of the second particle
 *  @param a  the LASS parameter
 *  @param r  the LASS parameter
 *  @param e  the LASS parameter
 */
// ============================================================================
Ostap::Math::LASS::LASS
( const double         m1 ,
  const double         m2 ,
  const double         m0 ,
  const double         g0 ,
  const double         a  ,
  const double         r  ,
  const double         e  )
  : m_m0  ( std::abs ( m0 ) )
  , m_g0  ( std::abs ( g0 ) )
  , m_a   ( std::abs ( a  ) )
  , m_r   ( std::abs ( r  ) )
  , m_e   ( std::abs ( e  ) )
// phase space
  , m_ps2 ( m1 , m2 )
//
  , m_workspace ()
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::LASS::~LASS(){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setM0 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_m0 ) ) { return false ; }
  //
  m_m0 = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setG0 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_g0 ) ) { return false ; }
  //
  m_g0 = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setA ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_a ) ) { return false ; }
  //
  m_a = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setR ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_r ) ) { return false ; }
  //
  m_r = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setE ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_e ) ) { return false ; }
  //
  m_e = v ;
  //
  return true ;
}
// ============================================================================
// get the (complex) LASS amplitude
// ============================================================================
std::complex<double>
Ostap::Math::LASS::amplitude ( const double x ) const
{
  //
  const double q  = m_ps2.q_ ( x ) ;
  if ( 0 >= q                ) { return 0 ; }  // RETURN
  //
  // get the width:
  const double gs = gamma_run ( m_g0        ,
                                x           ,
                                m_ps2.m1 () ,
                                m_ps2.m2 () ,
                                m_m0        ,
                                // K*(1430) is a scalar! 
                                0           ) * m_m0 / x  ;
  //
  // phase shift:
  const double cotB = 1.0 / ( m_a * q ) + 0.5 * m_r * q  ;
  // phase shift:
  const double cotR = ( m_m0 * m_m0 - x * x )  / m_m0 / gs ;
  //
  // const double sinB =  1.0 / std::sqrt ( 1 + cotB*cotB ) ;
  const double sinB =  1.0 / std::hypot ( 1.0 ,  cotB ) ;
  const double cosB = cotB * sinB ;
  //
  // exp( i*pi/2 )
  static const std::complex<double> i = std::complex<double>( 0 , 1 );
  //
  // exp( i*Delta_B )
  std::complex<double> deltaB ( cosB , sinB ) ;
  //
  // the amplitude
  std::complex<double> A =
    1.0 / ( cotB - i ) + m_e * deltaB * deltaB / ( cotR - i ) ;
  //
  // scale it!
  std::complex<double> T = A * ( x / q ) ;
  //
  return T ;
}
// ============================================================================
// get the phase space factor
// ============================================================================
double Ostap::Math::LASS::phaseSpace ( const double x ) const
{ return std::max ( 0.0 , m_ps2 ( x ) ) ; }
// ============================================================================
// evaluate LASS
// ============================================================================
double Ostap::Math::LASS::operator () ( const double x ) const
{
  const double result = phaseSpace  ( x ) ;
  if ( 0 >= result ) { return 0 ; }
  //
  return result * std::norm ( amplitude( x ) ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::LASS::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  if ( high <= m_ps2.lowEdge  () ) { return 0 ; }
  //
  if ( low  <  m_ps2.lowEdge  () )
  { return integral ( m_ps2.lowEdge() , high ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<LASS> s_integrator {} ;
  static char s_message[] = "Integral(LASS)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
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

// ============================================================================
// Bugg
// ============================================================================
/*  constructor from all masses and angular momenta
 *  @param M  mass of sigma (very different from the pole positon!)
 *  @param g2 width parameter g2 (4pi width)
 *  @param b1 width parameter b1  (2pi coupling)
 *  @param b2 width parameter b2  (2pi coupling)
 *  @param s1 width parameter s1  (cut-off for 4pi coupling)
 *  @param s2 width parameter s2  (cut-off for 4pi coupling)
 *  @param a  parameter a (the exponential cut-off)
 *  @param m1 the mass of the first  particle
 */
// ============================================================================
Ostap::Math::Bugg::Bugg
( const double         M  ,
  const double         g2 ,
  const double         b1 ,
  const double         b2 ,
  const double         a  ,
  const double         s1 ,
  const double         s2 ,
  const double         m1 )
//
  : m_M  ( std::abs ( M  ) )
  , m_g2 ( std::abs ( g2 ) )
  , m_b1 ( std::abs ( b1 ) )
  , m_b2 ( std::abs ( b2 ) )
  , m_s1 ( std::abs ( s1 ) )
  , m_s2 ( std::abs ( s2 ) )
  , m_a  ( std::abs ( a  ) )
// phase space
  , m_ps ( m1 , m1 )
//
  , m_workspace ()
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Bugg::~Bugg(){}
// ============================================================================
double Ostap::Math::Bugg::rho2_ratio ( const double x ) const
{
  if ( lowEdge() >= x ) { return 0 ; }
  //
  return
    Ostap::Math::PhaseSpace2::phasespace ( x    , m1() , m2 () ) /
    Ostap::Math::PhaseSpace2::phasespace ( M () , m1() , m2 () ) ;
}
// ============================================================================
std::complex<double>
Ostap::Math::Bugg::rho4_ratio ( const double x ) const
{
  //
  if ( 2 * m1() >= x ) { return 0 ; }
  //
  return rho4 ( x ) / rho4 ( M() ) ;
}
// ============================================================================
std::complex<double>
Ostap::Math::Bugg::rho4 ( const double x ) const
{
  const double s  = x * x ;
  //
  const double r2 = 1 - 16 * m1() * m1() / s ;
  //
  const double r  =
    std::sqrt ( std::abs ( r2 ) ) *
    ( 1 + std::exp ( ( s1 () - s )  / s2 () ) ) ;
  //
  return 0 <= r2 ?
    std::complex<double> ( r , 0 ) :
    std::complex<double> ( 0 , r ) ;
}
// ============================================================================
// Adler's pole
// ============================================================================
double Ostap::Math::Bugg::adler ( const double x ) const
{
  if ( lowEdge() >= x ) { return 0 ; }
  //
  const double pole = 0.5 * m1 () * m1 ()  ;
  //
  return ( x * x - pole ) / ( M2 () - pole ) ;
}
// ============================================================================
// get the running width by Bugg
// ============================================================================
std::complex<double>
Ostap::Math::Bugg::gamma ( const double x ) const
{
  //
  if ( lowEdge() >= x ) { return 0 ; }
  //
  const double s = x * x ;
  //
  const double g1 =
    b     ( x ) *
    adler ( x ) * std::exp ( -1 * ( s - M2() )  / a() ) ;
  //
  return g1 * rho2_ratio ( x ) + g2 () * rho4_ratio ( x ) ;
}
// ============================================================================
// get the amlitude  (not normalized!)
// ============================================================================
std::complex<double>
Ostap::Math::Bugg::amplitude (  const double x ) const
{
  if ( lowEdge() >= x ) { return 0 ; }
  //
  static const std::complex<double> j ( 0 , 1 ) ;
  //
  std::complex<double> d = M2() - x * x  - j * M() * gamma ( x ) ;
  //
  return 1.0 / d ;
}
// ============================================================================
// evaluate Bugg
// ============================================================================
double Ostap::Math::Bugg::pdf ( const double x ) const
{
  //
  if ( lowEdge() >= x ) { return 0 ; }
  //
  const double result = phaseSpace  ( x ) ;
  if ( 0 >= result ) { return 0 ; }
  //
  return result * std::norm ( amplitude ( x ) ) ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Bugg::setM ( const double x )
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
bool Ostap::Math::Bugg::setG2 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_g2 ) ) { return false ; }
  //
  m_g2 = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Bugg::setB1 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_b1 ) ) { return false ; }
  //
  m_b1 = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Bugg::setB2 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_b2 ) ) { return false ; }
  //
  m_b2 = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Bugg::setS1 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_s1 ) ) { return false ; }
  //
  m_s1 = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Bugg::setS2 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_s2 ) ) { return false ; }
  //
  m_s2 = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Bugg::setA ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_a ) ) { return false ; }
  //
  m_a = v ;
  //
  return true ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Bugg::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  if ( high <= lowEdge  () ) { return 0 ; }
  //
  if ( low  <  lowEdge  () )
  { return integral ( lowEdge() , high        ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Bugg> s_integrator {} ;
  static char s_message[] = "Integral(Bugg)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
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


// ============================================================================
// SWANSON CUSP 
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::Swanson::Swanson
( const double         m1     ,   // the first  real particle 
  const double         m2     ,   // the second real particle                
  const double         m1_0   ,   // the first  particle for cusp
  const double         m2_0   ,   // the second particle for cusp 
  const double         beta_0 ,   // beta_0 parameter
  const unsigned short L      )   // orbital momentum for real particles 
  : m_bw ( ( std::abs ( m1 )  + std::abs ( m2 ) ) * 2.1 , // almost arbitrary 
           ( std::abs ( m1 )  + std::abs ( m2 ) ) * 0.5 , // almost arbitrary  
           std::abs ( m1 ) , 
           std::abs ( m2 ) ,  
           L               ) 
  , m_m1         ( std::abs (   m1_0 ) )
  , m_m2         ( std::abs (   m2_0 ) )
  , m_beta0      ( std::abs ( beta_0 ) )
    //
  , m_workspace  ()
    //
{}
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::Swanson::Swanson
( const double         m1             ,   // the first  real particle 
  const double         m2             ,   // the second real particle                
  const double         m1_0           ,   // the first  particle for cusp
  const double         m2_0           ,   // the second particle for cusp 
  const double         beta_0         ,   // beta_0 parameter
  const unsigned short L              ,   // orbital momentum for real particles 
  const Ostap::Math::FormFactors::JacksonRho  r )  //  formfactor
  : m_bw ( ( std::abs ( m1 )  + std::abs ( m2 ) ) * 2.1 , // almost arbitrary 
           ( std::abs ( m1 )  + std::abs ( m2 ) ) * 0.5 , // almost arbitrary  
           std::abs ( m1 ) , 
           std::abs ( m2 ) ,  
           L  , r          )            
  , m_m1         ( std::abs (   m1_0 ) )
  , m_m2         ( std::abs (   m2_0 ) )
  , m_beta0      ( std::abs ( beta_0 ) )
    //
  , m_workspace  ()
    //
{}

// ============================================================================
// constructor
// ============================================================================
Ostap::Math::Swanson::Swanson
( const double         m1             ,   // the first  real particle 
  const double         m2             ,   // the second real particle                
  const double         m1_0           ,   // the first  particle for cusp
  const double         m2_0           ,   // the second particle for cusp 
  const double         beta_0         ,   // beta_0 parameter
  const unsigned short L              ,   // orbital momentum for real particles 
  const Ostap::Math::FormFactor&   f  )  //  formfactor
  : m_bw ( ( std::abs ( m1 )  + std::abs ( m2 ) ) * 2.1 , // almost arbitrary 
           ( std::abs ( m1 )  + std::abs ( m2 ) ) * 0.5 , // almost arbitrary  
           std::abs ( m1 ) , 
           std::abs ( m2 ) ,  
           L  , f          )            
  , m_m1         ( std::abs (   m1_0 ) )
  , m_m2         ( std::abs (   m2_0 ) )
  , m_beta0      ( std::abs ( beta_0 ) )
    //
  , m_workspace  ()
    //
{}


// ============================================================================
// constructor
// ============================================================================
Ostap::Math::Swanson::Swanson
( const Ostap::Math::BreitWigner&   bw             ,   // breit-wigner 
  const double         m1_0   ,   // the first  particle for cusp
  const double         m2_0   ,   // the second particle for cusp 
  const double         beta_0 )   // beta_0 parameter
  : m_bw     ( bw ) 
  , m_m1     ( std::abs (   m1_0 ) )
  , m_m2     ( std::abs (   m2_0 ) )
  , m_beta0  ( std::abs ( beta_0 ) )
    //
  , m_workspace  ()
    //
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::Swanson::Swanson
( const Ostap::Math::Swanson& sw ) 
  : m_bw    ( sw.m_bw    )
  , m_m1    ( sw.m_m1    )
  , m_m2    ( sw.m_m2    )
  , m_beta0 ( sw.m_beta0 )
    //
  , m_workspace  ()
    //
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Swanson::~Swanson (){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Swanson::setM1_0 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_m1 ) ) { return false ; }
  //
  m_m1 = v ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Swanson::setM2_0 ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_m2 ) ) { return false ; }
  //
  m_m2 = v ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Swanson::setBeta0( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_beta0 ) ) { return false ; }
  //
  m_beta0 = v ;
  //
  return true ;
}
// ============================================================================
//  calculate the Swanson amplitude
// ============================================================================
std::complex<double>
Ostap::Math::Swanson::amplitude ( const double x ) const
{
  //
  const double  f = - s_SQRT2PISQUAREDi*m_beta0/(1/m_m1+1/m_m2) ;
  //
  const double zf = 4 * m_m1 * m_m2 / ( m_beta0 * m_beta0 * ( m_m1 + m_m2 ) ) ;
  const double z  = zf * ( m_m1 + m_m2 - x ) ;
  //
  // above threshold, Z is negative 
  std::complex<double> iZ = 
    0 <= z ? 
    std::complex<double>(     std::sqrt (            z   ) , 0 ) :
    std::complex<double>( 0 , std::sqrt ( std::abs ( z ) )     ) ;
  //
  return f * 0.5 * s_SQRTPIHALF * ( 1.0 - s_SQRTPI * iZ * Ostap::Math::erfcx ( iZ ) ) ;
}
// ============================================================================
//  calculate the Swanson shape 
// ============================================================================
double Ostap::Math::Swanson::swanson ( const double x ) const
{
  if ( m_bw.m1() + m_bw.m2() >= x ) { return 0 ; }
  //
  const double g  = m_bw.gamma ( x ) ;
  if ( 0 >= g ) { return 0 ; }
  //
  const std::complex<double> a = amplitude ( x ) ;
  //
  return 2 * x * std::norm ( a ) * g / m_bw.gam0() / M_PI ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Swanson::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  const double x_min  = m_bw.m1() + m_bw.m2() ;
  if ( x_min >= high ) { return                        0   ; }
  if ( x_min >  low  ) { return integral  ( x_min , high ) ; }
  //
  // split into reasonable sub intervals
  //
  const double x1   = x_min +  1 * ( m_m1 + m_m2 ) ;
  const double x2   = x_min +  2 * ( m_m1 + m_m2 ) ;
  const double x5   = x_min +  5 * ( m_m1 + m_m2 ) ;
  const double x10  = x_min + 10 * ( m_m1 + m_m2 ) ;
  //
  if ( low <  x1 &&  x1 < high ) { return integral ( low ,  x1 ) + integral (  x1 , high ) ; }
  if ( low <  x2 &&  x2 < high ) { return integral ( low ,  x2 ) + integral (  x2 , high ) ; }
  if ( low <  x5 &&  x5 < high ) { return integral ( low ,  x5 ) + integral (  x5 , high ) ; }
  if ( low < x10 && x10 < high ) { return integral ( low , x10 ) + integral ( x10 , high ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Swanson> s_integrator {} ;
  static char s_message[] = "Integral(Swanson)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION         ,          // absolute precision
      ( x10  <= low  ) ? s_PRECISION_TAIL :
      s_PRECISION         ,          // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================





// ============================================================================
// Voigtian
// ============================================================================
// constructor  from all parameters
// ============================================================================
Ostap::Math::Voigt::Voigt
( const double m0    ,
  const double gamma ,
  const double sigma )
  : m_m0        ( m0 )
  , m_gamma     ( std::abs ( gamma ) )
  , m_sigma     ( std::abs ( sigma ) )
//
  , m_workspace ()
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Voigt::~Voigt(){}
// ============================================================================
// get the value of Voigt function
// ============================================================================
double Ostap::Math::Voigt::operator() ( const double x ) const
{
  //
  const double s1 = 1 / ( m_sigma * s_SQRT2   ) ;
  const double s2 = 1 / ( m_sigma * s_SQRT2PI ) ;
  //
  return Ostap::Math::faddeeva_w
    ( std::complex<double> ( x - m_m0 , m_gamma ) * s1 ).real() * s2 ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Voigt::integral
( const double low  ,
  const double high ) const
{
  //
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  const double width = std::max ( m_sigma , m_gamma ) ;
  //
  // split into reasonable sub intervals
  //
  const double x_low   = m_m0 - 4 * width ;
  const double x_high  = m_m0 + 4 * width ;
  //
  if      ( low <  x_low  && x_low  < high )
  {
    return
      integral (   low  , x_low   ) +
      integral ( x_low  ,   high  ) ;
  }
  else if ( low <  x_high && x_high < high )
  {
    return
      integral (   low  , x_high  ) +
      integral ( x_high ,   high  ) ;
  }
  //
  // split, if interval too large
  //
  if ( 0 < width && 10 * width < high - low  )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  const bool in_tail = 
    ( low  > m_m0 + 10 * width ) || ( high < m_m0 + 10 * width ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Voigt> s_integrator {} ;
  static char s_message[] = "Integral(Voigt)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Voigt::integral () const { return 1 ; }
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Voigt::setM0 ( const double x )
{
  //
  if ( s_equal ( x , m_m0 ) ) { return false ; }
  //
  m_m0 = x ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Voigt::setGamma ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_gamma ) ) { return false ; }
  //
  m_gamma = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Voigt::setSigma ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_sigma ) ) { return false ; }
  //
  m_sigma = v ;
  //
  return true ;
}
// ============================================================================
/*  full width at half maximum 
 *  @see http://en.wikipedia.org/wiki/Voigt_profile
 */
// ============================================================================
double Ostap::Math::Voigt::fwhm   () const 
{
  const double fg = 2 * m_sigma * s_BUKIN ;
  return 0.5346 * m_gamma + std::sqrt ( 0.2166 * m_gamma * m_gamma + fg * fg ) ;
}
// ============================================================================


// ============================================================================
// PseudoVoigtian
// T. Ida, M. Ando and H. Toraya
// "Extended pseudo-Voigt function for approximating the Voigt profile"
// J. Appl. Cryst. (2000). 33, 1311-1316
// doi:10.1107/S0021889800010219
// http://dx.doi.org/10.1107/S0021889800010219
// ============================================================================
// constructor  from all parameters
// ============================================================================
Ostap::Math::PseudoVoigt::PseudoVoigt
( const double m0    ,
  const double gamma ,
  const double sigma )
  : m_m0        ( m0 )
  , m_gamma     ( std::abs ( gamma ) )
  , m_sigma     ( std::abs ( sigma ) )
    //
  , m_w         ( 4 , 0 ) 
  , m_eta       ( 4 , 0 )
  , m_workspace ()
{
  update() ;  
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::PseudoVoigt::~PseudoVoigt(){}
// ============================================================================
namespace 
{
  // ==========================================================================
  /// gaussian profile 
  inline double f_gauss       ( const double dx , const double gamma ) 
  { return my_exp ( - dx * dx / ( gamma * gamma) ) / ( gamma * s_SQRTPI ) ; }
  // ==========================================================================
  /// lorenzian profile 
  inline double f_lorentzian  ( const double dx , const double gamma ) 
  { return gamma / ( ( dx*dx + gamma * gamma)  * M_PI ) ; }
  // ==========================================================================
  /// irrational profile 
  inline double f_irrational  ( const double dx , const double gamma ) 
  { return std::pow ( 1.0 +  dx*dx/(gamma*gamma) , -1.5 ) / ( 2 * gamma ) ; }
  // ==========================================================================
  /// squared sech profile 
  inline double f_sech2       ( const double dx , const double gamma ) 
  { 
    const double s = Ostap::Math::sech ( dx / gamma ) ;
    return s * s / ( 2 * gamma )  ; 
  }
  // ==========================================================================
  // parametrization data
  // ==========================================================================
  const std::array<double,7> s_Ai = {{   0.66000 ,   0.15021 ,  -1.24984 , 
                                         4.74052 ,  -9.48291 ,   8.48252 , -2.95553  }} ;
  const std::array<double,7> s_Bi = {{ -0.42179  ,  -1.25693 ,  10.30003 , 
                                       -23.45651 ,  29.14158 , -16.60453 ,  3.19974  }} ;
  const std::array<double,7> s_Ci = {{  1.19913  ,   1.43021 , -15.36331 , 
                                        47.06071 , -73.61822 ,  57.92559 , -17.80614 }} ;
  const std::array<double,7> s_Di = {{   1.10186 ,  -0.47745 ,  -0.68688 , 
                                         2.76622 ,  -4.55466 ,   4.05475 ,  -1.26571 }} ;
  const std::array<double,7> s_Fi = {{ -0.30165  ,  -1.38927 ,   9.31550 , 
                                       -24.10743 ,  34.96491 , -21.18862 ,   3.70290 }} ;
  const std::array<double,7> s_Gi = {{ 0.25437   ,  -0.14107 ,   3.23653 ,
                                       -11.09215 ,  22.10544 , -24.12407 ,   9.76947 }} ;
  const std::array<double,7> s_Hi = {{ 1.01579   ,   1.50429 ,  -9.21815 ,
                                       23.59717  , -39.71134 ,  32.83023 , -10.02142 }} ;
  // ==========================================================================
  inline double w_G ( const double rho ) 
  { return 1 - rho    *Ostap::Math::Clenshaw::monomial_sum ( s_Ai.rbegin() , 
                                                             s_Ai.rend()   , rho ).first ; }
  inline double w_L ( const double rho ) 
  { return 1 - (1-rho)*Ostap::Math::Clenshaw::monomial_sum ( s_Bi.rbegin() , 
                                                             s_Bi.rend()   , rho ).first ; }
  inline double w_I ( const double rho ) 
  { return             Ostap::Math::Clenshaw::monomial_sum ( s_Ci.rbegin() , 
                                                             s_Ci.rend()   , rho ).first ; }
  
  inline double w_P ( const double rho ) 
  { return             Ostap::Math::Clenshaw::monomial_sum ( s_Di.rbegin() ,
                                                             s_Di.rend()   , rho ).first ; }
  
  inline double eta_L ( const double rho ) 
  { return rho * ( 1 + ( 1 - rho ) * Ostap::Math::Clenshaw::monomial_sum ( s_Fi.rbegin() , 
                                                                           s_Fi.rend()   , rho ).first ) ; } 
  inline double eta_I ( const double rho ) 
  { return rho * ( 1 - rho ) * Ostap::Math::Clenshaw::monomial_sum ( s_Gi.rbegin() , 
                                                                     s_Gi.rend()   , rho ).first  ; }
  inline double eta_P ( const double rho ) 
  { return rho * ( 1 - rho ) * Ostap::Math::Clenshaw::monomial_sum ( s_Hi.rbegin() , 
                                                                     s_Hi.rend()   , rho ).first  ; }  
  // ==========================================================================
  // constants 
  // ==========================================================================
  // W_G <--> gamma_G 
  const double s_PV_cG = 1.0 / ( 2*std::sqrt ( std::log ( 2.0 ) ) ) ;
  // W_L <--> gamma_L 
  const double s_PV_cL = 0.5  ;
  // W_I <--> gamma_I 
  const double s_PV_cI = 1/(2.0*std::sqrt(std::pow(2.0,2.0/3)-1)) ;
  // W_P <--> gamma_P 
  const double s_PV_cP = 1/(2.0*std::acosh(std::sqrt(2.0))) ;
  // ==========================================================================
}
// ============================================================================
double Ostap::Math::PseudoVoigt::fwhm_gauss()  const 
{ return 2 * m_sigma * s_BUKIN ; }
// ============================================================================
void Ostap::Math::PseudoVoigt::update() 
{
  const double _rho = rho() ;
  //
  m_w  [0] =   w_G ( _rho ) * s_PV_cG ;
  m_w  [1] =   w_L ( _rho ) * s_PV_cL ;
  m_w  [2] =   w_I ( _rho ) * s_PV_cI ;
  m_w  [3] =   w_P ( _rho ) * s_PV_cP ;
  //
  m_eta[1] = eta_L ( _rho )           ;
  m_eta[2] = eta_I ( _rho )           ;
  m_eta[3] = eta_P ( _rho )           ;
  //
  m_eta[0] = 1 - m_eta[1] - m_eta[2] - m_eta[3] ;
}
// ============================================================================
// get the value of PseudoVoigt function
// ============================================================================
double Ostap::Math::PseudoVoigt::operator() ( const double x ) const
{
  //
  const double gamma_sum = fwhm_gauss() + fwhm_lorentzian() ;
  //
  const double dx =  ( x - m_m0 ) / gamma_sum ;
  //
  return 
    ( f_gauss      ( dx , m_w[0] ) * m_eta[0] + 
      f_lorentzian ( dx , m_w[1] ) * m_eta[1] +             
      f_irrational ( dx , m_w[2] ) * m_eta[2] + 
      f_sech2      ( dx , m_w[3] ) * m_eta[3]   ) / gamma_sum ;
}
// ============================================================================
// get the Gaussian component 
// ============================================================================
double Ostap::Math::PseudoVoigt::gaussian   ( const double x ) const 
{
  const double gamma_sum = fwhm_gauss() + fwhm_lorentzian() ;
  const double dx =  ( x - m_m0 ) / gamma_sum ;
  //
  return f_gauss ( dx , m_w[0] ) * m_eta[0] / gamma_sum ;
}
// ============================================================================
// get the Lorentzian component 
// ============================================================================
double Ostap::Math::PseudoVoigt::lorentzian   ( const double x ) const 
{
  const double gamma_sum = fwhm_gauss() + fwhm_lorentzian() ;
  const double dx =  ( x - m_m0 ) / gamma_sum ;
  //
  return f_lorentzian ( dx , m_w[1] ) * m_eta[1] / gamma_sum ;
}
// ============================================================================
// get the Irrational component 
// ============================================================================
double Ostap::Math::PseudoVoigt::irrational  ( const double x ) const 
{
  const double gamma_sum = fwhm_gauss() + fwhm_lorentzian() ;
  const double dx =  ( x - m_m0 ) / gamma_sum ;
  //
  return f_irrational ( dx , m_w[2] ) * m_eta[2] / gamma_sum ;
}
// ============================================================================
// get the Sech2 component 
// ============================================================================
double Ostap::Math::PseudoVoigt::sech2  ( const double x ) const 
{
  const double gamma_sum = fwhm_gauss() + fwhm_lorentzian() ;
  const double dx =  ( x - m_m0 ) / gamma_sum ;
  //
  return f_sech2 ( dx , m_w[3] ) * m_eta[3] / gamma_sum ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::PseudoVoigt::integral
( const double low  ,
  const double high ) const
{
  //
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  const double width = std::max ( m_sigma , m_gamma ) ;
  //
  // split into reasonable sub intervals
  //
  const double x_low   = m_m0 - 4 * width ;
  const double x_high  = m_m0 + 4 * width ;
  //
  if      ( low <  x_low  && x_low  < high )
  {
    return
      integral (   low  , x_low   ) +
      integral ( x_low  ,   high  ) ;
  }
  else if ( low <  x_high && x_high < high )
  {
    return
      integral (   low  , x_high  ) +
      integral ( x_high ,   high  ) ;
  }
  //
  // split, if interval too large
  //
  if ( 0 < width && 10 * width < high - low  )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  const bool in_tail = 
    ( low  > m_m0 + 10 * width ) || ( high < m_m0 + 10 * width ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PseudoVoigt> s_integrator {} ;
  static char s_message[] = "Integral(PseudoVoigt)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::PseudoVoigt::integral () const { return 1 ; }
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::PseudoVoigt::setM0 ( const double x )
{
  //
  if ( s_equal ( x , m_m0 ) ) { return false ; }
  //
  m_m0 = x ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::PseudoVoigt::setGamma ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_gamma ) ) { return false ; }
  //
  m_gamma = v ;
  //
  // recalculate data 
  update() ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::PseudoVoigt::setSigma ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_sigma ) ) { return false ; }
  //
  m_sigma = v ;
  //
  // recalculate data 
  update() ;
  //
  return true ;
}
// ============================================================================

// ============================================================================
// "2-from-3" shapes 
// ============================================================================
 

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::BW23L::BW23L
( const double         m0   ,
  const double         gam0 ,
  const double         m1   ,
  const double         m2   ,
  const double         m3   ,
  const double         m    ,
  const unsigned short L1   ,
  const unsigned short L2   )
  //
  : m_bw ( std::make_unique<BreitWigner> ( m0 , gam0 , m1  , m2 , L1 ) )
  , m_ps ( m1 , m2   , m3  , m  , L2 , L1 )
    //
  , m_workspace ()
{}
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::BW23L::BW23L
( const double                               m0   ,
  const double                               gam0 ,
  const double                               m1   ,
  const double                               m2   ,
  const double                               m3   ,
  const double                               m    ,
  const unsigned short                       L1   ,
  const unsigned short                       L2   ,
  const Ostap::Math::FormFactors::JacksonRho r    )
  //
  : m_bw ( std::make_unique<BreitWigner> ( m0 , gam0 , m1  , m2 , L1 , r ) )
  , m_ps ( m1 , m2   , m3  , m  , L2 , L1 )
//
  , m_workspace ()
{}
// ============================================================================
// constructor from BreitWigner
// ============================================================================
Ostap::Math::BW23L::BW23L
( const Ostap::Math::BreitWigner& bw ,
  const double                    m3 ,
  const double                    m  ,
  const unsigned short            L2 )
  //
  : m_bw ( bw.clone () ) 
  , m_ps ( bw.m1() , bw.m2() , m3  , m  , L2 , bw. L())
    //
  , m_workspace ()
{}
// ============================================================================
// COPY
// ============================================================================
Ostap::Math::BW23L::BW23L ( const Ostap::Math::BW23L& right ) 
  : m_bw ( right.m_bw->clone() ) 
  , m_ps ( right.m_ps )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::BW23L::~BW23L (){}
// ============================================================================
// calculate the shape
// ============================================================================
double Ostap::Math::BW23L::operator() ( const double x ) const
{
  if (  lowEdge() >= x || highEdge()  <= x ) { return 0 ; }
  //
  const double bw = std::norm ( m_bw->amplitude ( x ) )   ;
  //
  // // get the incomplete phase space factor
  // const double ps  =                   // get the incomplete phase space factor
  //   x / M_PI *
  //   // =======================================================================
  //   // the second factor is already in our BW !!!
  //   Ostap::Math::PhaseSpace2::phasespace ( x          , 
  //                                          m_bw.m1 () , 
  //                                          m_bw.m2 () , 
  //                                          m_bw.L  () ) *
  //   // =======================================================================
  //   Ostap::Math::PhaseSpace2::phasespace ( m_ps.m  () ,
  //                                          x          ,
  //                                          m_ps.m3 () ,
  //                                          m_ps.L  () ) ;
  //
  return bw * m_ps ( x ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::BW23L::integral
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
  if ( low  <  lowEdge  () )
  { return integral ( lowEdge() , high        ) ; }
  //
  if ( high >  highEdge () )
  { return integral ( low       , highEdge () ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<BW23L> s_integrator {} ;
  static char s_message[] = "Integral(BW23L)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
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
// get the integral
// ============================================================================
double  Ostap::Math::BW23L::integral () const
{ return integral ( lowEdge () , highEdge() ) ; }
// ============================================================================


// ============================================================================
// Flatte23L
// ============================================================================
/*  constructor  from all parameters
 *  @param m0    the mass
 *  @param m0g1  parameter \f$ m_0\times g_1\f$
 *  @param g2og2 parameter \f$ g2/g_1       \f$
 *  @param mK    kaon mass
 *  @param mPi   pion mass
 *  @param m3    the mass of the third particle
 *  @param m     the mass of mother particle
 *  @param L     the orbital momentum between the pair and the third particle
 */
// ============================================================================
Ostap::Math::Flatte23L::Flatte23L
( const double         m0    ,     // MeV
  const double         m0g1  ,     // MeV^2
  const double         g2og1 ,     // dimensionless
  const double         mA    ,     // MeV
  const double         mB    ,     // MeV
  const double         m3    ,     // MeV
  const double         m     ,     // MeV
  const unsigned short L     )
  //
  : m_flatte ( std::make_unique<Flatte> ( m0  , m0g1 , g2og1 , mA , mA , mB , mB ) )
  , m_ps        ( mA  , mA  , m3    , m  , L    )
    //
  , m_workspace ()
{}
// ============================================================================
/* constructor  from flatte function
 *  @param m3    the mass of the third particle
 *  @param m     the mass of mother particle
 *  @param L     the orbital momentum between the pair and the third particle
 */
// ============================================================================
Ostap::Math::Flatte23L::Flatte23L
( const Ostap::Math::Flatte& fun ,     // MeV
  const double               m3  ,     // MeV
  const double               m   ,     // MeV
  const unsigned short       L   )
//
  : m_flatte    ( fun.clone() ) 
  , m_ps        ( fun.mA1() , fun.mA2()  , m3    , m  , L    )
    //
  , m_workspace ()
{}
// ============================================================================
// copy
// ============================================================================
Ostap::Math::Flatte23L::Flatte23L
( const Ostap::Math::Flatte23L& right ) 
  : m_flatte ( right.m_flatte->clone() )
  , m_ps     ( right.m_ps )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Flatte23L::~Flatte23L (){}
// ============================================================================
// get the value of Flatte function
// ============================================================================
double Ostap::Math::Flatte23L::operator() ( const double x ) const
{
  //
  if ( lowEdge () >= x || highEdge() <= x ) { return 0 ; } // RETURN
  //
  // get the amplitude...
  std::complex<double> amp = m_flatte->flatte_amp ( x ) ;
  //
  return m_ps ( x ) * std::norm ( amp ) * 2 / M_PI * m0g1() ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Flatte23L::integral
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
  if ( low  <  lowEdge  () )
  { return integral ( lowEdge() , high        ) ; }
  //
  if ( high >  highEdge () )
  { return integral ( low       , highEdge () ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Flatte23L> s_integrator {} ;
  static char s_message[] = "Integral(Flatte23L)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION         ,          // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double  Ostap::Math::Flatte23L::integral () const
{ return integral ( lowEdge () , highEdge() ) ; }
// ============================================================================


// ============================================================================
// Gounaris & Sakurai shape
// ============================================================================
/* constructor from all masses and angular momenta
 *  @param M  mass of rho
 *  @param g0 width parameter
 *  @param m1 the mass of the first  particle (the same as the second)
 *  @param m3 the mass of the third  particle
 *  @param m  the mass of the mother particle (m>m1+m2+m3)
 *  @param L  the angular momentum between the first pair and the third
 */
// ============================================================================
Ostap::Math::Gounaris23L::Gounaris23L
( const double         M  ,  // GeV
  const double         g0 ,  // GeV
  const double         m1 ,  // MeV
  const double         m3 ,  // MeV
  const double         m  ,  // MeV
  const unsigned short L  )
//
  : m_M  ( std::abs ( M  ) )
  , m_g0 ( std::abs ( g0 ) )
//
  , m_ps ( m1 , m1 , m3 , m , L , 1 )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Gounaris23L::~Gounaris23L(){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Gounaris23L::setM ( const double x )
{
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_M ) ) { return false ; }
  m_M = v ;
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::Gounaris23L::setG0 ( const double x )
{
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_g0 ) ) { return false ; }
  m_g0 = v ;
  return true ;
}
// ============================================================================
// get h-factor
// ============================================================================
double Ostap::Math::Gounaris23L::h ( const double x ,
                                     const double k ) const
{
  if ( lowEdge() > x || highEdge() < x ) { return 0 ; }
  return 2 * k  / M_PI / x * std::log ( ( x + 2 * k ) / 2 / m1() ) ;
}
// ============================================================================
// get h-factor
// ============================================================================
double Ostap::Math::Gounaris23L::h ( const double x ) const
{
  if ( lowEdge() > x ) { return 0 ; }
  const double k = PhaseSpace2::q ( x , m1 () , m1() ) ;
  return h ( x , k ) ;
}
// ============================================================================
// get h'-factor
// ============================================================================
double Ostap::Math::Gounaris23L::h_prime ( const double x ) const
{
  if ( lowEdge() > x ) { return 0 ; }
  const double k = PhaseSpace2::q ( x , m1 () , m1() ) ;
  return h_prime ( x , k ) ;
}
// ============================================================================
// get h'-factor
// ============================================================================
double Ostap::Math::Gounaris23L::h_prime ( const double x ,
                                           const double k ) const
{
  //
  if ( lowEdge() > x ) { return 0 ; }
  //
  const double f =  ( x + 2 * k ) / ( 2  * m1 () ) ;
  //
  return k / M_PI / x / x * ( - std::log ( f ) / x  + 0.5 / m1() / f ) ;
}
// ============================================================================
// get the amlitude  (not normalized!)
// ============================================================================
std::complex<double>
Ostap::Math::Gounaris23L::amplitude (  const double x ) const
{
  //
  if ( x <= lowEdge() ) { return 0 ; }
  //
  const double k    = PhaseSpace2::q ( x    , m1 () , m1 () ) ;
  const double k0   = PhaseSpace2::q ( M () , m1 () , m1 () ) ;
  const double k03  = k0 * k0 * k0 ;
  //
  const double m0_2 = M() * M() ;
  //
  const double v1   = m0_2 - x * x ;
  //
  const double dh   = h ( x , k ) - h ( M() , k0 ) ;
  const double hp   = h_prime ( m() , k0 ) ;
  //
  const double v2 = k * k * dh + k0 * k0 * hp * ( m0_2 - x * x ) ;
  const double v3 = Ostap::Math::pow ( k/k0 , 3 ) * m0() / x ;
  //
  return
    std::sqrt ( g0 () * m0 () ) /
    std::complex<double> ( v1 + v2 * g0() * m0_2 / k03 ,
                           v3      * g0() * m0 ()      ) ;
}
// ============================================================================
// calculate the Gounaris-Sakurai shape
// ============================================================================
double Ostap::Math::Gounaris23L::operator() ( const double x ) const
{
  //
  if ( lowEdge() >= x || highEdge() <= x ) { return 0 ; }
  //
  std::complex<double> amp = amplitude ( x ) ;
  const double  ps = m_ps( x ) ;
  //
  return x * ps * std::norm ( amp ) * 2 / M_PI  ;
}

// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Gounaris23L::integral
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
  if ( low  <  lowEdge  () )
  { return integral ( lowEdge() , high        ) ; }
  //
  if ( high >  highEdge () )
  { return integral ( low       , highEdge () ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Gounaris23L> s_integrator {} ;
  static char s_message[] = "Integral(Gounaris23L)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION         ,          // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double  Ostap::Math::Gounaris23L::integral () const
{ return integral ( lowEdge () , highEdge() ) ; }
// ============================================================================


// ============================================================================
// LASS: Kpi S-wave for  X -> (K pi) Y decays..
// ============================================================================
/*  constructor from all masses and angular momenta
 *  @param m1 the mass of the first  particle
 *  @param m2 the mass of the second particle
 *  @param m3 the mass of the third  particle
 *  @param m  the mass of the mother particle (m>m1+m2+m3)
 *  @param L  the angular momentum between the first pair and the third
 *  @param a  the LASS parameter
 *  @param r  the LASS parameter
 */
// ============================================================================
Ostap::Math::LASS23L::LASS23L
( const double         m1 ,
  const double         m2 ,
  const double         m3 ,
  const double         m  ,
  const double         m0 ,
  const double         g0 ,
  const unsigned short L  ,
  const double         a  ,
  const double         r  ,
  const double         e  )
// LASS-function 
  : m_lass ( m1 , m2 , m0 , g0  , a , r , e )
// phase space
  , m_ps   ( m1 , m2 , m3 , m   , L , 0 )
//
  , m_workspace ()
{}
// ============================================================================
/*  constructor from LASS and 3-rd particle 
 *  @param lass the actual lass shape 
 *  @param m3   the mass of third particle (Y)
 *  @param m    the mass of mother particle (X)
 *  @param L    the orbital momentum between Y and (Kpi) 
 */
// ============================================================================
Ostap::Math::LASS23L::LASS23L
( const Ostap::Math::LASS& lass   , 
  const double             m3     ,
  const double             m      ,
  const unsigned short     L      ) 
// LASS-function 
  : m_lass ( lass ) 
// phase space
  , m_ps   ( lass.m1 ()  , lass.m2 () , m3 , m   , L , 0 )
//
  , m_workspace ()
{}  
// ============================================================================

// ============================================================================
// destructor
// ============================================================================
Ostap::Math::LASS23L::~LASS23L(){}
// ============================================================================
// get the (complex) LASS amplitude
// ============================================================================
std::complex<double>
Ostap::Math::LASS23L::amplitude ( const double x ) const
{ return m_lass.amplitude ( x )  ; }  
// ============================================================================
// get the phase space factor
// ============================================================================
double Ostap::Math::LASS23L::phaseSpace ( const double x ) const
{ return std::max ( 0.0 , m_ps ( x ) ) ; }
// ============================================================================
// evaluate LASS
// ============================================================================
double Ostap::Math::LASS23L::operator () ( const double x ) const
{
  const double result = phaseSpace  ( x ) ;
  if ( 0 >= result ) { return 0 ; }
  //
  return result * std::norm ( amplitude( x ) ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::LASS23L::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  if ( high <= m_ps.lowEdge  () ) { return 0 ; }
  if ( low  >= m_ps.highEdge () ) { return 0 ; }
  //
  if ( low  <  m_ps.lowEdge  () )
  { return integral ( m_ps.lowEdge() , high             ) ; }
  if ( high >  m_ps.highEdge () )
  { return integral ( low            , m_ps.highEdge () ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<LASS23L> s_integrator {} ;
  static char s_message[] = "Integral(LASS23L)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION         ,          // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double  Ostap::Math::LASS23L::integral () const
{ return integral ( m_ps.lowEdge () , m_ps.highEdge() ) ; }
// ============================================================================

// ============================================================================
// Bugg23L
// ============================================================================
/*  constructor from all masses and angular momenta
 *  @param M  mass of sigma (very different from the pole positon!)
 *  @param g2 width parameter g2 (4pi width)
 *  @param b1 width parameter b1  (2pi coupling)
 *  @param b2 width parameter b2  (2pi coupling)
 *  @param s1 width parameter s1  (cut-off for 4pi coupling)
 *  @param s2 width parameter s2  (cut-off for 4pi coupling)
 *  @param a  parameter a (the exponential cut-off)
 *  @param m1 the mass of the first  particle
 *  @param m3 the mass of the third  particle
 *  @param m  the mass of the mother particle (m>m1+m2+m3)
 *  @param L  the angular momentum between the first pair and the third
 */
// ============================================================================
Ostap::Math::Bugg23L::Bugg23L
( const double         M  ,
  const double         g2 ,
  const double         b1 ,
  const double         b2 ,
  const double         a  ,
  const double         s1 ,
  const double         s2 ,
  const double         m1 ,
  const double         m3 ,
  const double         m  ,
  const unsigned short L  )
//
  : m_bugg ( M  , g2 , b1 , b2 , a , s1 , s2 , m1 ) 
  , m_ps   ( m1 , m1 , m3 , m  , L , 0 )
//
  , m_workspace ()
{}
// ============================================================================
/*  constructor from bugg & phase space parameters 
 *  @param m3 the mass of the third  particle
 *  @param m  the mass of the mother particle (m>m1+m2+m3)
 *  @param L  the angular momentum between the first pair and the third
 */
// ============================================================================
Ostap::Math::Bugg23L::Bugg23L
( const Ostap::Math::Bugg& bugg ,
  const double             m3   ,  // MeV
  const double             m    ,  // MeV
  const unsigned short     L    ) 
//
  : m_bugg ( bugg ) 
  , m_ps   ( bugg.m1 () , bugg.m1 ()  , m3 , m  , L , 0 )
//
  , m_workspace ()
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Bugg23L::~Bugg23L(){}
// ============================================================================
// evaluate Bugg
// ============================================================================
double Ostap::Math::Bugg23L::pdf ( const double x ) const
{
  //
  if ( lowEdge() >= x || highEdge() <= x ) { return 0 ; }
  //
  const double result = phaseSpace  ( x ) ;
  if ( 0 >= result ) { return 0 ; }
  //
  return result * std::norm ( amplitude ( x ) ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::Bugg23L::integral
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
  if ( low  <  lowEdge  () )
  { return integral ( lowEdge() , high        ) ; }
  //
  if ( high >  highEdge () )
  { return integral ( low       , highEdge () ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Bugg23L> s_integrator {} ;
  static char s_message[] = "Integral(Bugg23L)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( &F , 
      low , high          ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      s_PRECISION         ,          // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double  Ostap::Math::Bugg23L::integral () const
{ return integral ( lowEdge () , highEdge() ) ; }
// ============================================================================





// ============================================================================
//                                                                      The END 
// ============================================================================

