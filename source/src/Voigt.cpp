// ============================================================================
// Include files
// =============================================================================
// STD&STL
// =============================================================================
#include <cmath>
#include <array>
// =============================================================================
// Ostap
// ============================================================================
#include "Ostap/Voigt.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/MoreMath.h"
// ============================================================================
//  local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "local_gsl.h"
#include "Integrator1D.h"
// ============================================================================
/** @file
 *  Implementation file for 
 *  - class Ostap::Math::Voigt 
 *  - class Ostap::Math::PseudoVoigt 
 */
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
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()  ,          // size of workspace
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
  const double fg = 2 * m_sigma * s_Bukin ;
  return 0.5346 * m_gamma + std::sqrt ( 0.2166 * m_gamma * m_gamma + fg * fg ) ;
}
// ============================================================================



// ============================================================================
// PseudoVoigtian
// T. Ida, M. Ando and H. Toraya
// "Extended pseudo-Voigt function for approximating the Voigt profile"
// J. Appl. Cryst. (2000). 33, 1311-1316
// doi:10.1107/S0021889800010219
// https://doi.org/10.1107/S0021889800010219
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
  const std::array<double,7> s_Ai = {{   0.66000  ,  
                                         0.15021  ,  
                                         -1.24984 , 
                                         4.74052  , 
                                         -9.48291 ,  
                                         8.48252  , 
                                         -2.95553 }} ;
  const std::array<double,7> s_Bi = {{ -0.42179   , 
                                       -1.25693   , 
                                       10.30003   , 
                                       -23.45651  , 
                                       29.14158   , 
                                       -16.60453  ,  
                                       3.19974    }} ;
  const std::array<double,7> s_Ci = {{  1.19913   , 
                                        1.43021   , 
                                        -15.36331 , 
                                        47.06071  , 
                                        -73.61822 ,  
                                        57.92559  , 
                                        -17.80614 }} ;
  const std::array<double,7> s_Di = {{   1.10186  , 
                                         -0.47745 , 
                                         -0.68688 , 
                                         2.76622  ,
                                         -4.55466 ,
                                         4.05475  , 
                                         -1.26571 }} ;
  const std::array<double,7> s_Fi = {{ -0.30165   , 
                                       -1.38927   ,
                                       9.31550    , 
                                       -24.10743  , 
                                       34.96491   , 
                                       -21.18862  ,
                                       3.70290    }} ;
  const std::array<double,7> s_Gi = {{ 0.25437    , 
                                       -0.14107   , 
                                       3.23653    ,
                                       -11.09215  , 
                                       22.10544   , 
                                       -24.12407  ,  
                                       9.76947    }} ;
  const std::array<double,7> s_Hi = {{ 1.01579    , 
                                       1.50429    , 
                                       -9.21815   ,
                                       23.59717   , 
                                       -39.71134  , 
                                       32.83023   , 
                                       -10.02142  }} ;
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
{ return 2 * m_sigma * s_Bukin ; }
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
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size () ,          // size of workspace
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
//                                                                      The END 
// ============================================================================
