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
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_log.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_psi.h"
#include "gsl/gsl_sf_zeta.h"
#include "gsl/gsl_integration.h"
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
#include "Ostap/Models.h"
// ============================================================================
//  Local 
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "local_gsl.h"
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
  /** @var s_TRUNC
   *  trunkating parameter for CrystalBall-functions
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double      s_TRUNC = 15.0 ;
  // ==========================================================================
  /** @var S_ATLAL 
   *  magic constant - integral for Atlas function 
   *  @see Ostap::Math::Atlas 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-08-21
   */   
  const double s_ATLAS = 3.052369876253939 ;
  // ==========================================================================
  // Bukin & Co
  // ==========================================================================
  /** @var s_Bukin
   *  useful constant for Bukin's function
   *  \f$ \sqrt{ 2 \log 2 } \f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  const double s_Bukin   = std::sqrt ( 2.0 * std::log ( 2.0 ) ) ;
  // ==========================================================================
  /** @var s_ln2
   *  useful constant for Bukin's function
   *  \f$ \log 2 \f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  const double s_ln2 = std::log ( 2.0 ) ;
  // ==========================================================================
  /** helper function for itegration of Bukin's function
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  double bukin_GSL ( double x , void* params )
  {
    const Ostap::Math::Bukin* bukin = (Ostap::Math::Bukin*) params ;
    return (*bukin)(x) ;
  }
  // ==========================================================================
  // Novosibirsk & Co
  // ==========================================================================
  /** @var s_Novosibirsk
   *  useful constant for evaliuation of ``Novosibirsk'' function
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  const double s_Novosibirsk = std::sqrt ( std::log ( 4.0 ) ) ;
  // ==========================================================================
  /** helper function for itegration of Novosibirsk's function
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double novosibirsk_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Novosibirsk* novosibirsk =
      (Ostap::Math::Novosibirsk*) params ;
    //
    return (*novosibirsk)(x) ;
  }
  // ==========================================================================
  /** helper function for itegration of Sigmoid function
   *  @see Ostap::Math::Sigmoid
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double sigmoid_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Sigmoid* sigmoid =
      (Ostap::Math::Sigmoid*) params ;
    //
    return (*sigmoid)(x) ;
  }
  // ==========================================================================
  /** helper function for itegration of Atlas's function
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-08-21
   */
  double atlas_GSL ( double x , void* params )
  {
    const Ostap::Math::Atlas* atlas = (Ostap::Math::Atlas*) params ;
    return (*atlas)(x) ;
  }
  // ==========================================================================
  /** evaluate the helper function \f[ f = \frac{\sinh (x) }{x} \f]
   *  it allows to calculate Novosibirsk's function in efficient and regular way
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  inline double x_sinh
  ( const double x , double precision = s_PRECISION )
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
      while ( std::fabs ( delta ) > 0.1 * precision && n < 1000 ) ;
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
  /// helper expression for erf(x)/x 
  inline double error_func_x ( const long double x ) 
  { return 0 == x || s_zero ( x ) ? s_HALFSQRTPIi : error_func( x ) / x ; }
  // ==========================================================================
  /** helper function for itegration
   *  \f[ f =  \mathrm{e}^{ \kappa x^2 + \xi x }  \f]
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  double gauss_GSL ( double x , void* params )
  {
    const double* gauss = (double*) params ;
    //
    const double kappa = gauss[0] ;
    const double xi    = gauss[1] ;
    //
    const double arg   = kappa * x * x + xi * x ;
    //
    return my_exp ( arg ) ;
  }
  // ==========================================================================
  /** get the gaussian integral numerically
   *  \f[ f = \int_a^b \mathrm{e}^ {-\alpha x^2 + \beta x } \mathrm{d}x \f]
   *  @param alpha the alpha parameter
   *  @param beta  the beta  parameter
   *  @param a     the low  integration limit
   *  @param b     the high integration limit
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double gaussian_int_num
  ( const double alpha ,
    const double beta  ,
    const double a     ,
    const double b     )
  {
    //
    if      ( s_equal ( a , b ) ) { return  0.0 ; }
    //
    // use GSL to evaluate the integral numerically 
    //
    Sentry sentry ;
    //
    gsl_function F            ;
    F.function = &gauss_GSL   ;
    double params[2]          ;
    params[0]  = -alpha       ;    // ATTENTION: minus sign here!
    params[1]  =  beta        ;
    F.params   = params       ;
    //
    //
    double result   = 1.0 ;
    double error    = 1.0 ;
    //
    Ostap::Math::WorkSpace ws ;
    const int ierror = gsl_integration_qag
      ( &F                ,            // the function
        a     ,   b       ,            // low & high edges
        s_PRECISION       ,            // absolute precision
        s_PRECISION       ,            // relative precision
        s_SIZE            ,            // size of workspace
        GSL_INTEG_GAUSS31 ,            // integration rule
        workspace ( ws )  ,            // workspace
        &result           ,            // the result
        &error            ) ;          // the error in result
    //
    if ( ierror )
    {
      gsl_error ( "Ostap::Math::gaussian_int" ,
                  __FILE__ , __LINE__ , ierror ) ;
    }
    //
    return result ;
    //
  }
  // ==========================================================================
  /** get the exponential integral 
   *  \f[ f = \int_a^b \exp { \beta x } \mathrm{d}x \f]
   *  @param beta  the beta  parameter
   *  @param a     the low  integration limit
   *  @param b     the high integration limit
   */
  double exponent_int
  ( const double beta , 
    const double a    , 
    const double b    ) 
  {
    //
    if      ( s_equal ( a , b ) ) { return 0     ; }
    else if ( s_zero  ( beta  ) ) { return b - a ; }
    //
    const double beta_b = beta * b ;
    const double beta_a = beta * a ;
    //
    double result = 0 ;
    if ( !s_zero ( b ) ) { result += b * reduced_exp ( beta_b ) ; }
    if ( !s_zero ( a ) ) { result -= a * reduced_exp ( beta_a ) ; }
    //
    return result ;
  }
  // ==========================================================================
  double gaussian_int_L ( const double alpha ,
                          const double beta  ,
                          const double b     );  
  double gaussian_int_R ( const double alpha ,
                          const double beta  ,
                          const double b     );  
  // ==========================================================================
  /** get the gaussian integral:
   *  \f[ f = \int_a^b \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
   *  @attention note the sign for alpha-term!
   *  @param alpha the alpha parameter
   *  @param beta  the beta  parameter
   *  @param a     the low  integration limit
   *  @param b     the high integration limit
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double gaussian_int
  ( const double alpha ,
    const double beta  ,
    const double a     ,
    const double b     )
  {
    // trivial case 
    if ( s_equal ( a , b ) ) { return 0 ; }
    if ( a > b             ) { return - gaussian_int ( alpha , beta , b , a ) ; }
    //
    // 1) simple exponential integral ?
    //
    if ( s_zero  ( alpha ) ) { return exponent_int ( beta , a , b ) ; }    
    //
    // just a standard error-function 
    //
    if ( s_zero  ( beta  ) && 0 < alpha )  
    {
      const double sqrt_alpha = std::sqrt ( alpha ) ;
      const double ba = b * sqrt_alpha ;
      const double aa = a * sqrt_alpha ;      
      return s_HALFSQRTPI * ( b * error_func_x ( ba ) - a * error_func_x ( aa ) ) ;
    }
    else if ( 0 < alpha ) 
    {
      const double b2a        = beta / ( 2 * alpha ) ;
      if ( a < b2a && b2a < b ) 
      {
        return 
          gaussian_int ( alpha , beta , a   , b2a ) + 
          gaussian_int ( alpha , beta , b2a , b   ) ;
      }
      const double c   = b2a * alpha * b2a    ;      
      if  ( b2a <= a && b2a <= b ) 
      {
        const double sqrt_alpha = std::sqrt ( alpha ) ;      
        const double a1  = ( a - b2a ) * sqrt_alpha ;
        const double b1  = ( b - b2a ) * sqrt_alpha ;
        return s_HALFSQRTPI / sqrt_alpha * 
          ( my_exp ( -alpha * a * a + beta * a ) * Ostap::Math::erfcx ( a1 ) - 
            my_exp ( -alpha * b * b + beta * b ) * Ostap::Math::erfcx ( b1 ) ) ;  
      }
      else if ( a <= b2a  && b<= b2a ) 
      { return gaussian_int ( alpha , beta , 2*b2a - b , 2*b2a - a ) ; }
      //
      // .. should never be here, except some testing regime
      //
      if ( c < 0.1 * GSL_LOG_DBL_MAX ) 
      { return my_exp ( c ) * gaussian_int ( alpha , 0 , a - b2a  , b - b2a ) ; } 
    }
    //
    // use GSL to evaluate the integral numerically 
    //
    return gaussian_int_num ( alpha , beta , a , b ) ;
  }
  // ==========================================================================
  /** get the gaussian integral:
   *  \f[ f = \int_a^{\inf} \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
   *  @param alpha the alpha parameter
   *  @param beta  the beta  parameter
   *  @param a     the low integration limit
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double gaussian_int_R ( const double alpha ,
                          const double beta  ,
                          const double a     )
  {
    //
    if      ( 0 > alpha ) { return s_INFINITY ; }
    else if ( s_zero ( alpha ) ) 
    { return  0 > beta ? -beta * my_exp ( beta * a ) : s_INFINITY ; }
    //
    const double sqrt_alpha = std::sqrt ( alpha ) ;      
    const double b2a        = beta / ( 2 * alpha ) ;
    if ( b2a <= a ) 
    {
      const double a1  = ( a - b2a ) * sqrt_alpha ;
      return s_HALFSQRTPI / sqrt_alpha * 
        my_exp ( -alpha * a * a + beta * a ) * Ostap::Math::erfcx ( a1 ) ;
    }
    //
    return 
      gaussian_int    ( alpha , beta , a   , b2a ) +
      gaussian_int_R  ( alpha , beta ,       b2a ) ;
  }
  // ==========================================================================
  /** get the gaussian integral:
   *  \f[ f = \int_{-\inf}^{b} \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
   *  @param alpha the alpha parameter
   *  @param beta  the beta  parameter
   *  @param b     the high integration limit
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double gaussian_int_L ( const double alpha ,
                          const double beta  ,
                          const double b     )
  { return gaussian_int_R ( alpha , -beta , -b ) ; }
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
  /** @var x_sqrt2
   *  \f$\sqrt{2}\f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_sqrt2  = s_SQRT2 ;
  // ==========================================================================
  /** helper function for integration of Gram-Charlier A function
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double gram_charlier_A_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::GramCharlierA* gca =
      (Ostap::Math::GramCharlierA*) params ;
    //
    return (*gca)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Breit-Wigner shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double breit_wigner_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::BreitWigner* bw =
      (Ostap::Math::BreitWigner*) params ;
    //
    return (*bw)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Swanson shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2016-06-11
   */
  double swanson_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Swanson* sw =
      (Ostap::Math::Swanson*) params ;
    //
    return (*sw)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Flatte shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double flatte_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Flatte* f = (Ostap::Math::Flatte*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of PhaseSpace2 shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double phase_space_2_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::PhaseSpace2* ps =
      (Ostap::Math::PhaseSpace2*) params ;
    //
    return (*ps)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of PhaseSpace shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double phase_space_3_1_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::PhaseSpace3* ps =
      (Ostap::Math::PhaseSpace3*) params ;
    //
    return ps -> ps2_aux (x) ;
  }
  // ==========================================================================
  /** helper function for itegration of PhaseSpace shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double phase_space_3_2_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::PhaseSpace3* ps =
      (Ostap::Math::PhaseSpace3*) params ;
    //
    return (*ps)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of PhaseSpaceNL shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double phase_space_NL_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::PhaseSpaceNL* ps =
      (Ostap::Math::PhaseSpaceNL*) params ;
    //
    return (*ps)(x) ;
  }
  // ==========================================================================
  /** helper function for itegration of PhaseSpacePol shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2014-05-05
   */
  double phase_space_POL_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::PhaseSpacePol* ps =
      (Ostap::Math::PhaseSpacePol*) params ;
    //
    return (*ps)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of PhaseSpace23L shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double phase_space_23L_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::PhaseSpace23L* ps23L =
      (Ostap::Math::PhaseSpace23L*) params ;
    //
    return (*ps23L)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of LASS shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2013-10--5
   */
  double LASS_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::LASS* lass = (Ostap::Math::LASS*) params ;
    //
    return (*lass)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of LASS23L shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double LASS_23L_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::LASS23L* lass = (Ostap::Math::LASS23L*) params ;
    //
    return (*lass)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Bugg23L shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2012-05-23
   */
  double Bugg_23L_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Bugg23L* bugg = (Ostap::Math::Bugg23L*) params ;
    //
    return (*bugg)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Bugg shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2012-05-23
   */
  double Bugg_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Bugg* bugg = (Ostap::Math::Bugg*) params ;
    //
    return (*bugg)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of BW23L shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2012-05-23
   */
  double BW_23L_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::BW23L* bw = (Ostap::Math::BW23L*) params ;
    //
    return (*bw)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Flatte23L shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2012-05-24
   */
  double Flatte_23L_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Flatte23L* f = (Ostap::Math::Flatte23L*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Gounaris23L shape
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2012-05-24
   */
  double Gounaris_23L_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Gounaris23L* f = (Ostap::Math::Gounaris23L*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Voigt shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  double voigt_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Voigt* f = (Ostap::Math::Voigt*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of PseudoVoigt shape
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2016-06-13
   */
  double pseudovoigt_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::PseudoVoigt* f = (Ostap::Math::PseudoVoigt*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Apolonios shape 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2013-12-1
   */
  double apolonios_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Apolonios* f = (Ostap::Math::Apolonios*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Apolonios2 shape 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2013-12-1
   */
  double apolonios2_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Apolonios2* f = (Ostap::Math::Apolonios2*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of Tsallis shape 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-07-11
   */
  double tsallis_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::Tsallis* f = (Ostap::Math::Tsallis*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
  /** helper function for integration of QGSM shape 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-07-11
   */
  double qgsm_GSL ( double x , void* params )
  {
    //
    const Ostap::Math::QGSM* f = (Ostap::Math::QGSM*) params ;
    //
    return (*f)(x) ;
  }
  // ==========================================================================
} //                                                 end of anonymous namespace
// ============================================================================



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





// ============================================================================
namespace
{
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
    const std::complex<double> v =
      m0 * m0 - x * x - s_j * m0 * gamma ;
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
} // end of the anonymous namespace
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
double Ostap::Math::BifurcatedGauss::pdf ( const double x ) const
{
  const double dx = x - m_peak ;
  //
  const double norm = s_SQRTPIHALF * ( sigmaL() + sigmaR() ) ;
  //
  return
    dx < 0 ?
    my_exp ( -0.5 * dx * dx / sigmaL () / sigmaL () ) / norm :
    my_exp ( -0.5 * dx * dx / sigmaR () / sigmaR () ) / norm ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::BifurcatedGauss::integral () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::BifurcatedGauss::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  const double norm = s_SQRTPIHALF * ( sigmaL() + sigmaR() ) ;
  //
  // left half-gaussian
  if       ( high <= m_peak )
  {
    return gaussian_int ( 0.5 / sigmaL() / sigmaL() ,
                          0                       ,
                          low  - m_peak           ,
                          high - m_peak           ) / norm ;

  }
  // right half-gaussian
  else if ( low >= m_peak )
  {
    return gaussian_int ( 0.5 / sigmaR() / sigmaR() ,
                          0                       ,
                          low  - m_peak           ,
                          high - m_peak           ) / norm ;

  }
  //
  return
    integral ( low    , m_peak ) +
    integral ( m_peak , high   ) ;
}
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
    m_gbeta2  = gsl_sf_lngamma ( 3 / beta() ) ;    
    m_gbeta2 -= gsl_sf_lngamma ( 1 / beta() ) ;
    // m_gbeta2  = gsl_sf_exp     ( m_gbeta2   ) ;
    m_gbeta2  = my_exp     ( m_gbeta2   ) ;
  }
  else 
  { 
    m_gbeta1 = 1. / gsl_sf_gamma ( 1 / beta() )            ;
    m_gbeta2 =      gsl_sf_gamma ( 3 / beta() ) * m_gbeta1 ;
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
    double result  = gsl_sf_log ( 0.5 * beta() / alpha() ) ;
    result        -= delta2 ;
    result        -= gsl_sf_lngamma ( 1 / beta() ) ;
    return my_exp ( result ) ;
  }
  //
  double result   = 0.5 * beta() / alpha() ;
  result         *= my_exp   ( -delta2 ) ;
  // result         *= gsl_sf_exp   ( -delta2 ) ;
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
  double result  =   gsl_sf_lngamma ( 5 / beta() ) ;
  result        +=   gsl_sf_lngamma ( 1 / beta() ) ;
  result        -= 2*gsl_sf_lngamma ( 3 / beta() ) ;
  //
  return gsl_sf_exp ( result ) - 3 ;
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
    -gsl_sf_log        ( 1 + delta ) / kappa() :
    -gsl_sf_log_1plusx (     delta ) / kappa() ;  
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
  const double e_   = gsl_sf_erf ( y_ * s_SQRT2i ) ;
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
  return xi() - 0.5 * alpha() * kappa() * gsl_sf_exprel ( k2 ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::variance () const 
{
  if ( s_equal ( kappa() , 0 ) ) { return alpha () * alpha() ; }
  //
  const double k2 = kappa() * kappa() ;
  //
  return alpha () * alpha() * gsl_sf_exp ( k2 ) * gsl_sf_exprel ( k2 ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::sigma () const 
{ return std::sqrt ( variance() ) ; }
// ============================================================================
double Ostap::Math::GenGaussV2::skewness () const 
{
  const double k2     = kappa () * kappa() ;
  //
  const double a1     = gsl_sf_exprel (     k2 ) ;
  const double a3     = gsl_sf_exprel ( 3 * k2 ) ;
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
  const double ek2 = gsl_sf_exp ( kappa() * kappa() ) ;
  //
  return  Ostap::Math::pow ( ek2 , 4 )  
    + 2 * Ostap::Math::pow ( ek2 , 3 ) 
    + 3 * Ostap::Math::pow ( ek2 , 2 ) - 6 ;
}
// ============================================================================
// #include "boost/math/distributions/skew_normal.hpp"
// // ============================================================================
// /*  constructor from all agruments 
//  *  @param xi     location/peak posiiton 
//  *  @param omega  "scale" parameter 
//  *  @param alpha  "shape" parameter 
//  */
// // ============================================================================
// Ostap::Math::SkewGauss::SkewGauss
// ( const double xi    ,
//   const double omega , 
//   const double alpha )  // alpha=0 correponds to gaussian  
//   : m_xi      ( xi                 ) 
//   , m_omega   ( std::abs ( omega ) ) 
//   , m_alpha   (            alpha   ) 
// {}
// // ============================================================================
// // desctructor
// // ============================================================================
// Ostap::Math::SkewGauss::~SkewGauss(){}
// // ============================================================================
// bool Ostap::Math::SkewGauss::setXi        ( const double value ) 
// {
//   if ( s_equal ( value , m_xi) ) { return false ; }
//   m_xi = value ;
//   return true ;
// }
// // ============================================================================
// bool Ostap::Math::SkewGauss::setOmega     ( const double value ) 
// {
//   const double value_ = std::abs ( value ) ;
//   if ( s_equal ( value_ , m_omega ) ) { return false ; }
//   m_omega = value_ ;
//   return true ;
// }
// // ============================================================================
// bool Ostap::Math::SkewGauss::setAlpha  ( const double value ) 
// {
//   //
//   if ( s_equal ( value , m_alpha ) ) { return false ; }
//   m_alpha = value ;
//   if ( s_equal ( 0     , m_alpha ) ) { m_alpha = 0  ; }
//   //
//   return true ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::pdf ( const double x ) const 
// {
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   return boost::math::pdf ( skew_ , x ) ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::cdf ( const double x ) const 
// {
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   return boost::math::cdf ( skew_ , x ) ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::integral ( const double low  , 
//                                           const double high ) const 
// {
//   if ( s_equal ( low , high ) ) { return 0 ; }
//   //
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   //
//   return 
//     boost::math::cdf ( skew_ , high ) - 
//     boost::math::cdf ( skew_ , low  ) ;  
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::mean () const 
// {
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   //
//   return boost::math::mean ( skew_ ) ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::variance () const 
// {
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   //
//   return boost::math::variance ( skew_ ) ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::skewness () const 
// {
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   //
//   return boost::math::skewness ( skew_ ) ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::kurtosis () const 
// {
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   //
//   return boost::math::kurtosis ( skew_ ) ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::mode  () const 
// {
//   typedef boost::math::skew_normal_distribution<>  skew ;
//   skew skew_ ( m_xi , m_omega , m_alpha ) ;
//   //
//   return boost::math::kurtosis ( skew_ ) ;
// }
// // ============================================================================
// double Ostap::Math::SkewGauss::sigma  () const 
// { return std::sqrt ( variance () ) ; }
// // ============================================================================






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
  // the left tail
  //
  if ( high <= std::min ( m_x1 , m_x2 ) )  // left tail
  {
    const double d =  m_peak - m_x1 ;
    return  0.5 * gaussian_int     ( m_rho_L * m_rho_L / ( d * d ) ,
                                     m_L     / m_sigma  ,
                                     low  - m_x1        ,
                                     high - m_x1        ) ;
  }
  //
  // the right tail:
  //
  if ( low >= std::max ( m_x1 , m_x2  ) )  // right tail
  {
    const double d = m_peak - m_x2 ;
    return 0.5 * gaussian_int    ( m_rho_R * m_rho_R / ( d * d )  ,
                                   -1 * m_R  / m_sigma ,
                                   low  - m_x2         ,
                                   high - m_x2         ) ;
  }
  //
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                ;
  F.function = &bukin_GSL ;
  F.params   = const_cast<Bukin*> ( this ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const bool in_tail = 
    ( high < m_x1 - 5 * std::abs ( m_x2 - m_x1 ) )  || 
    ( low  > m_x2 + 5 * std::abs ( m_x2 - m_x1 ) ) ; 
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    //
    gsl_error ( "Ostap::Math::Bukin::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  double result = 0 ;
  // left tail:
  {
    const double d     = m_peak - m_x1     ;
    const double alpha = m_rho_L / d / d   ;
    const double beta  = m_L     / m_sigma ;
    //
    result +=  0.5 * gaussian_int_L ( alpha , beta , 0 ) ;
  }
  // right tail
  {
    const double d     =    m_peak - m_x2     ;
    const double alpha =    m_rho_R / d / d   ;
    const double beta  =  - m_R     / m_sigma ;
    //
    result += 0.5 * gaussian_int_R ( alpha , beta  , 0 ) ;
  }
  //
  return result + integral ( m_x1 , m_x2 ) ;
  //
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
  //
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                ;
  F.function = &novosibirsk_GSL ;
  F.params   = const_cast<Novosibirsk*> ( this ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL :
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Novosibirsk::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  Sentry sentry ;
  //
  gsl_function F                ;
  F.function = &novosibirsk_GSL ;
  F.params   = const_cast<Novosibirsk*> ( this ) ;
  //
  // left tail:
  //
  double tail_l   =  0.0 ;
  double error_l  = -1.0 ;
  //
  const int ierror_l = gsl_integration_qagil
    ( &F                ,         // the function
      x_low             ,         // "high" edge
      s_PRECISION       ,         // absolute precision
      s_PRECISION_TAIL  ,         // relative precision
      s_SIZE            ,         // size of workspace
      workspace ( m_workspace ) , // workspace
      &tail_l           ,         // the result
      &error_l          ) ;        // the error in result
  //
  if ( ierror_l )
  {
    gsl_error ( "Ostap::Math::Novosibirsk::QAGIL" ,
                __FILE__ , __LINE__ , ierror_l ) ;
    tail_l = 0.0 ;
  }
  //
  //
  // right tail:
  //
  double tail_r   =  0.0 ;
  double error_r  = -1.0 ;
  //
  const int ierror_r = gsl_integration_qagiu
    ( &F                ,         // the function
      x_high            ,         // "low" edge
      s_PRECISION       ,         // absolute precision
      s_PRECISION_TAIL  ,         // relative precision
      s_SIZE            ,         // size of workspace
      workspace ( m_workspace ) , // workspace
      &tail_r           ,         // the result
      &error_r          ) ;       // the error in result
  //
  if ( ierror_r )
  {
    gsl_error ( "Ostap::Math::Novosibirsk::QAGIU" ,
                __FILE__ , __LINE__ , ierror_r ) ;
    tail_r = 0.0 ;
  }
  //
  // get the final result
  //
  m_integral = tail_l + integral ( x_low , x_high ) + tail_r ;
  //
}
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
  m_B =  0.5 * ( 1 + gsl_sf_erf ( - m_alpha * s_SQRT2i ) ) ;
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
  m_B  = 0.5 * ( 1 + gsl_sf_erf ( - m_alpha * s_SQRT2i ) ) ;
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
  { return s_SQRT2PIi * gaussian_int ( 0.5   , 0 , zlow  , zhigh ) ; }
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
  m_B  = 0.5 *  ( gsl_sf_erf (  m_alpha_R * s_SQRT2i ) - 
                  gsl_sf_erf ( -m_alpha_L * s_SQRT2i ) ) ;
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
  m_B        = 0.5 *  ( gsl_sf_erf (  m_alpha_R * s_SQRT2i ) - 
                        gsl_sf_erf ( -m_alpha_L * s_SQRT2i ) ) ;
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
  m_B        = 0.5 *  ( gsl_sf_erf (  m_alpha_R * s_SQRT2i ) - 
                        gsl_sf_erf ( -m_alpha_L * s_SQRT2i ) ) ;
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
  { return  s_SQRT2PIi * gaussian_int ( 0.5   , 0 , zlow  , zhigh ) ; }
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
    //
    // use GSL to evaluate the integral
    //
    Sentry sentry ;
    //
    gsl_function F                ;
    F.function = &apolonios_GSL ;
    F.params   = const_cast<Apolonios*> ( this ) ;
    //
    double result   = 1.0 ;
    double error    = 1.0 ;
    //
    const int ierror = gsl_integration_qag
      ( &F                ,            // the function
        low   , high      ,            // low & high edges
        s_PRECISION       ,            // absolute precision
        s_PRECISION       ,            // relative precision
        s_SIZE            ,            // size of workspace
        GSL_INTEG_GAUSS31 ,            // integration rule
        workspace ( m_workspace ) ,    // workspace
        &result           ,            // the result
        &error            ) ;          // the error in result
    //
    if ( ierror )
    {
      gsl_error ( "Ostap::Math::Apolonios::QAG" ,
                  __FILE__ , __LINE__ , ierror ) ;
    }
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
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F               ;
  F.function = &apolonios2_GSL ;
  F.params   = const_cast<Apolonios2*> ( this ) ;
  //
  const double tail = ( low >= xR || high <= xL ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //  
  const int ierror = gsl_integration_qag
    ( &F                ,                     // the function
      low   , high      ,                     // low & high edges
      tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE            ,                     // size of workspace
      GSL_INTEG_GAUSS31 ,                     // integration rule
      workspace ( m_workspace ) ,             // workspace
      &result           ,                     // the result
      &error            ) ;                   // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Apolonios2::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
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
  if      ( s_equal ( low , high ) ) { return                  0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low   ) ; } // RETURN
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
  Sentry sentry ;
  //
  gsl_function F                ;
  F.function = &gram_charlier_A_GSL ;
  F.params   = const_cast<GramCharlierA*> ( this ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL :
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::GramCharlierA::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  Sentry sentry ;
  //
  gsl_function F                  ;
  F.function = &phase_space_2_GSL ;
  const PhaseSpace2* _ps = this  ;
  F.params   = const_cast<PhaseSpace2*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,             // the function
      xlow   , xhigh    ,             // low & high edges
      s_PRECISION       ,             // absolute precision
      s_PRECISION       ,             // relative precision
      s_SIZE            ,             // size of workspace
      GSL_INTEG_GAUSS31 ,             // integration rule
      workspace ( m_workspace ) ,     // workspace
      &result           ,             // the result
      &error            ) ;           // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::PhaseSpace2::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
 *  \f$ \Phi = \frac{1}{8\pi} \frac{ \lambda^{\frac{1}{2}} \left( m^2 , m_1^2, m_2_2 \right) }{ m^2 }\f$,
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &phase_space_3_1_GSL ;
  const PhaseSpace3* _ps = this  ;
  F.params   = const_cast<PhaseSpace3*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::PhaseSpace3::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &phase_space_3_2_GSL ;
  const PhaseSpace3* _ps = this  ;
  F.params   = const_cast<PhaseSpace3*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,             // the function
      low   , high      ,             // low & high edges
      s_PRECISION       ,             // absolute precision
      s_PRECISION       ,             // relative precision
      s_SIZE            ,             // size of workspace
      GSL_INTEG_GAUSS31 ,             // integration rule
      workspace ( m_workspace2 ) ,    // workspace
      &result           ,             // the result
      &error            ) ;           // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::PhaseSpace3::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  if ( ( 3 * m_N * 0.5 - 3       * 0.5 ) < GSL_SF_GAMMA_XMAX &&
       ( 3 * m_L * 0.5 - 3       * 0.5 ) < GSL_SF_GAMMA_XMAX &&
       ( 3 * m_N * 0.5 - 3 * m_L * 0.5 ) < GSL_SF_GAMMA_XMAX )
  {
    m_norm  = gsl_sf_gamma   ( 3 * m_N * 0.5 - 3       * 0.5 ) ;
    m_norm /= gsl_sf_gamma   ( 3 * m_L * 0.5 - 3       * 0.5 ) ;
    m_norm /= gsl_sf_gamma   ( 3 * m_N * 0.5 - 3 * m_L * 0.5 ) ;
  }
  else
  {
    m_norm  = gsl_sf_lngamma ( 3 * m_N * 0.5 - 3       * 0.5 ) ;
    m_norm -= gsl_sf_lngamma ( 3 * m_L * 0.5 - 3       * 0.5 ) ;
    m_norm -= gsl_sf_lngamma ( 3 * m_N * 0.5 - 3 * m_L * 0.5 ) ;
    m_norm  = gsl_sf_exp     ( m_norm ) ;
  }
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &phase_space_NL_GSL ;
  const PhaseSpaceNL* _ps = this  ;
  F.params   = const_cast<PhaseSpaceNL*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::PhaseSpaceNL::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double  Ostap::Math::PhaseSpaceNL::integral() const
{ return integral ( m_threshold1 , m_threshold2 ) ; }
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
double Ostap::Math::PhaseSpacePol::operator () ( const double x ) const 
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
// desctructor 
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &phase_space_POL_GSL ;
  const PhaseSpacePol* _ps = this  ;
  F.params                 = const_cast<PhaseSpacePol*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      xlow   , xhigh    ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::PhaseSpacePol::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
}
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &breit_wigner_GSL ;
  const BreitWigner* _bw = this  ;
  F.params   = const_cast<BreitWigner*> ( _bw ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL :
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::BreitWigner::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  const double x1     = m_m0 - 10 * m_gam0 ;
  const double x2     = m_m0 + 10 * m_gam0  ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &breit_wigner_GSL ;
  const BreitWigner* _bw = this  ;
  F.params   = const_cast<BreitWigner*> ( _bw ) ;
  //
  //
  // right tail:
  //
  double result  =  0.0 ;
  double error   = -1.0 ;
  //
  const int ierror = gsl_integration_qagiu
    ( &F                ,         // the function
      x_high            ,         // "low" edge
      s_PRECISION       ,         // absolute precision
      s_PRECISION_TAIL  ,         // relative precision
      s_SIZE            ,         // size of workspace
      workspace ( m_workspace ) , // workspace
      &result           ,         // the result
      &error            ) ;       // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::BreitWigner::QAGIU" ,
                __FILE__ , __LINE__ , ierror ) ;
    result = 0.0 ;
  }
  //
  return result + integral ( m_m1 + m_m2 , x_high );
}
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
/* constructor  from three parameters
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &flatte_GSL ;
  const Flatte* _f = this  ;
  F.params   = const_cast<Flatte*> ( _f ) ;
  //
  double result   =  1.0 ;
  double error    = -1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      ( high   <= x_low  ) ? s_PRECISION_TAIL :
      ( x_high <=   low  ) ? s_PRECISION_TAIL :
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Flatte::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  //
  const double x_low  = threshold () ;
  const double x_high = m_m0
    + 15 * std::abs ( m_m0g1 / m_m0           )
    + 15 * std::abs ( m_m0g1 / m_m0 * m_g2og1 ) ;
  //
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &flatte_GSL ;
  const Flatte* _f = this  ;
  F.params   = const_cast<Flatte*> ( _f ) ;
  //
  // right tail:
  //
  double result  =  0.0 ;
  double error   = -1.0 ;
  //
  const int ierror = gsl_integration_qagiu
    ( &F                ,         // the function
      x_high            ,         // "low" edge
      s_PRECISION       ,         // absolute precision
      s_PRECISION_TAIL  ,         // relative precision
      s_SIZE            ,         // size of workspace
      workspace ( m_workspace ) , // workspace
      &result           ,         // the result
      &error            ) ;       // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Flatte::QAGIU" ,
                __FILE__ , __LINE__ , ierror ) ;
    result = 0.0 ;
  }
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
//               Flatte
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
Ostap::Math::Flatte2::Flatte2
( const Ostap::Math::Flatte& flatte ) 
  : Ostap::Math::Flatte ( flatte )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Flatte2::~Flatte2(){}
// ============================================================================
// get the value of Flatte function
// ============================================================================
double Ostap::Math::Flatte2::operator() ( const double x ) const
{ return flatte2 ( x ) ; }
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
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &voigt_GSL ;
  const Voigt* _f = this  ;
  F.params   = const_cast<Voigt*> ( _f ) ;
  //
  //
  double result   =  1.0 ;
  double error    = -1.0 ;
  //
  const double in_tail = 
    ( low  > m_m0 + 10 * width ) || ( high < m_m0 + 10 * width ) ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Voigt::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &pseudovoigt_GSL ;
  const PseudoVoigt* _f = this  ;
  F.params   = const_cast<PseudoVoigt*> ( _f ) ;
  //
  //
  double result   =  1.0 ;
  double error    = -1.0 ;
  //
  const double in_tail = 
    ( low  > m_m0 + 10 * width ) || ( high < m_m0 + 10 * width ) ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::PseudoVoigt::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function = &swanson_GSL ;
  const Swanson* _sw = this  ;
  F.params   = const_cast<Swanson*> ( _sw ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      ( x10  <= low  ) ? s_PRECISION_TAIL :
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Swanson::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
}
// // ============================================================================
// // get the integral b
// // ============================================================================
// double  Ostap::Math::Swanson::integral () const
// {
//   //
//   const double x_min = m_bw.m1() + m_bw.m2() ;
//   //
//   // split into reasonable sub intervals
//   //
//   const double x10  = x_min + 10 * ( m_m1 + m_m2 ) ;
//   //
//   // use GSL to evaluate the integral
//   //
//   Sentry sentry ;
//   //
//   gsl_function F                 ;
//   F.function = &swanson_GSL ;
//   const Swanson* _sw = this  ;
//   F.params   = const_cast<Swanson*> ( _sw ) ;
//   //
//   double result   = 1.0 ;
//   double error    = 1.0 ;
//   //
//   const double x_high = x10 ;
//   const int ierror = gsl_integration_qagiu
//     ( &F                ,         // the function
//       x_high            ,         // "low" edge
//       s_PRECISION       ,         // absolute precision
//       s_PRECISION_TAIL  ,         // relative precision
//       s_SIZE            ,         // size of workspace
//       workspace ( m_workspace ) , // workspace
//       &result           ,         // the result
//       &error            ) ;       // the error in result
//   //
//   if ( ierror )
//   {
//     gsl_error ( "Ostap::Math::Swanson::QAGIU" ,
//                 __FILE__ , __LINE__ , ierror ) ;
//     result = 0.0 ;
//   }
//   //
//   return result + integral ( x_min , x_high );
// }
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function         = &LASS_GSL ;
  const LASS* _ps    = this  ;
  F.params           = const_cast<LASS*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::LASS::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
}
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function               = &phase_space_23L_GSL ;
  const PhaseSpace23L* _ps = this  ;
  F.params                 = const_cast<PhaseSpace23L*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::PhaseSpace23L::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function         = &LASS_23L_GSL ;
  const LASS23L* _ps = this  ;
  F.params           = const_cast<LASS23L*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::LASS23L::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function         = &Bugg_GSL ;
  const Bugg*    _ps = this  ;
  F.params           = const_cast<Bugg*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::BUGG::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
// double  Ostap::Math::Bugg23L::integral () const
// { return integral ( lowEdge () , highEdge() ) ; }
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
  Sentry sentry ;
  //
  gsl_function F                 ;
  F.function         = &Bugg_23L_GSL ;
  const Bugg23L* _ps = this  ;
  F.params           = const_cast<Bugg23L*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::BUGG23L::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  : m_bw ( m0 , gam0 , m1  , m2 , L1      )
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
  : m_bw ( m0 , gam0 , m1  , m2 , L1 , r  )
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
  : m_bw ( bw )
  , m_ps ( bw.m1() , bw.m2() , m3  , m  , L2 , bw. L())
//
  , m_workspace ()
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
  const double bw = std::norm ( m_bw.amplitude ( x ) )   ;
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
  Sentry sentry ;
  //
  gsl_function F                   ;
  F.function         = &BW_23L_GSL ;
  const BW23L* _ps   = this  ;
  F.params           = const_cast<BW23L*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::BW23L::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
  : m_flatte    ( m0  , m0g1 , g2og1 , mA , mA , mB , mB ) 
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
  : m_flatte    ( fun )
  , m_ps        ( fun.mA1() , fun.mA2()  , m3    , m  , L    )
    //
  , m_workspace ()
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
  std::complex<double> amp = m_flatte.flatte_amp ( x ) ;
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
  Sentry sentry ;
  //
  gsl_function F                   ;
  F.function             = &Flatte_23L_GSL ;
  const Flatte23L* _ps   = this  ;
  F.params               = const_cast<Flatte23L*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::BW23L::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
bool Ostap::Math::Gounaris23L::setG0 ( const double x )
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
// get h-factor
// ============================================================================
double Ostap::Math::Gounaris23L::h ( const double x ,
                                     const double k ) const
{
  //
  if ( lowEdge() > x || highEdge() < x ) { return 0 ; }
  //
  return 2 * k  / M_PI / x * std::log ( ( x + 2 * k ) / 2 / m1() ) ;
}
// ============================================================================
// get h-factor
// ============================================================================
double Ostap::Math::Gounaris23L::h ( const double x ) const
{
  //
  if ( lowEdge() > x ) { return 0 ; }
  //
  const double k = PhaseSpace2::q ( x , m1 () , m1() ) ;
  //
  return h ( x , k ) ;
}
// ============================================================================
// get h'-factor
// ============================================================================
double Ostap::Math::Gounaris23L::h_prime ( const double x ) const
{
  //
  if ( lowEdge() > x ) { return 0 ; }
  //
  const double k = PhaseSpace2::q ( x , m1 () , m1() ) ;
  //
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
  Sentry sentry ;
  //
  gsl_function F                   ;
  F.function             = &Gounaris_23L_GSL ;
  const Gounaris23L* _ps = this  ;
  F.params               = const_cast<Gounaris23L*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION     ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Gounaris23L::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
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
namespace 
{
  //
  inline double student_cdf (  const double t , const double nu ) 
  {
    const double xt    = nu / ( t * t + nu ) ;
    const double value = 0.5 * gsl_sf_beta_inc ( 0.5 * nu , 0.5 , xt ) ;
    return t >= 0 ? 1 - value : value ;
  }
  //
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
// Student-T 
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
    m_normR  = gsl_sf_gamma ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
    m_normR /= std::sqrt    ( M_PI * v ) ;
  }
  //
  if ( s_equal ( v , m_nR ) ) { return false ; }
  //
  m_nR = v ;
  //
  m_normR  = gsl_sf_gamma ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
  m_normR /= std::sqrt    ( M_PI * v ) ;
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


// ============================================================================
// Generalaize Gamma distribtions 
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
// ============================================================================


// ============================================================================
namespace 
{
  //
  inline double shash ( const double x   , 
                        const double eps , 
                        const double dlt ) 
  {
    const double y = eps + dlt * std::asinh ( x ) ;
    return 
      (     GSL_LOG_DBL_MAX < y ) ?    s_INFINITY :
      ( -1* GSL_LOG_DBL_MAX > y ) ? -1*s_INFINITY : std::sinh ( y ) ;
  }
  //
}
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
    *    std::hypot ( 1 ,z )          / std::hypot ( 1 , y )  
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
  //
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                ;
  F.function = &atlas_GSL ;
  F.params   = const_cast<Atlas*> ( this ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const bool in_tail = ( high <= left || low >= right ) ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // absolute precision
      in_tail ? s_PRECISION_TAIL : s_PRECISION , // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    //
    gsl_error ( "Ostap::Math::Atlas::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
}
// ============================================================================
// overall integral, not exact but precise enough...
// ============================================================================
double Ostap::Math::Atlas::integral () const { return 1 ; }
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
  //
  inline double phi_ ( const double x ) 
  { return gsl_ran_gaussian_pdf  ( x , 1 ) ; }
  inline double Phi_ ( const double x ) 
  { return gsl_cdf_ugaussian_P   ( x     ) ; }
  inline double Psi_ ( const double x ) 
  { return Phi_ ( x ) - x * phi_ ( x ) - 0.5 ; } 
  //
}
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
  // use GSL to evaluate the integral
  //
  Sentry sentry ;
  //
  gsl_function F                   ;
  F.function             = &sigmoid_GSL ;
  const Sigmoid*     _ps = this  ;
  F.params               = const_cast<Sigmoid*> ( _ps ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      low   , high      ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  {
    gsl_error ( "Ostap::Math::Sigmoid::QAG" ,
                __FILE__ , __LINE__ , ierror ) ;
  }
  //
  return result ;
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
  Sentry sentry ;
  //
  gsl_function F                ;
  F.function = &tsallis_GSL ;
  F.params   = const_cast<Tsallis*> ( this ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      _low   , high     ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  { gsl_error ( "Ostap::Math::Tsallis::QAG" , __FILE__ , __LINE__ , ierror ) ; }
  //
  return result ;
  //
}





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
  Sentry sentry ;
  //
  gsl_function F                ;
  F.function = &qgsm_GSL ;
  F.params   = const_cast<QGSM*> ( this ) ;
  //
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  const int ierror = gsl_integration_qag
    ( &F                ,            // the function
      _low   , high     ,            // low & high edges
      s_PRECISION       ,            // absolute precision
      s_PRECISION       ,            // relative precision
      s_SIZE            ,            // size of workspace
      GSL_INTEG_GAUSS31 ,            // integration rule
      workspace ( m_workspace ) ,    // workspace
      &result           ,            // the result
      &error            ) ;          // the error in result
  //
  if ( ierror )
  { gsl_error ( "Ostap::Math::QGSM::QAG" , __FILE__ , __LINE__ , ierror ) ; }
  //
  return result ;
  //
}
// ============================================================================




// ============================================================================
/*  get the gaussian integral
 *  \f[ f = \int_a^b \exp { -\alpha^2 x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param low   the low  integration limit
 *  @param high  the high integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ============================================================================
double Ostap::Math::gaussian_integral
( const double alpha ,
  const double beta  ,
  const double low   ,
  const double high  ) 
{
  // note the difference in the arguments! 
  return gaussian_int ( alpha * alpha , beta , low , high ) ;
}
// ============================================================================
/*  get the gaussian integral
 *  \f[ f = \int_{a}^{_\inf} \exp { -\alpha^2 x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param low   the low  integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ============================================================================
double Ostap::Math::gaussian_integral_right
( const double alpha ,
  const double beta  ,
  const double low   ) 
{
  // note the difference in the arguments! 
  return gaussian_int_R ( alpha * alpha , beta , low ) ;
}
// ============================================================================
/*  get the gaussian integral
 *  \f[ f = \int_{-\inf}^b \exp { -\alpha^2 x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param high  the high integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ============================================================================
double Ostap::Math::gaussian_integral_left
( const double alpha ,
  const double beta  ,
  const double high  ) 
{
  // note the difference in the arguments! 
  return gaussian_int_L ( alpha * alpha , beta , high ) ;
}
// ============================================================================
// The END
// ============================================================================


