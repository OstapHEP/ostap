// ============================================================================
#ifndef LOCAL_MATH_H 
#define LOCAL_MATH_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <numbers>
#include <limits>
#include <complex>
// ============================================================================
#if defined ( __cplusplus ) && defined ( __cpp_lib_math_constants ) && ( 201907L <= __cpp_lib_math_constants )
#include <numbers>
#endif 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
/// local namespace to hide all techinical symbols
namespace 
{
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double>                s_equal  {} ; // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>                    s_zero   {} ; // zero for doubles
  /// zero for vectors 
  const Ostap::Math::Zero< std::vector<double> >     s_vzero  {} ; // zero for vectors
  /// zero for comples doubles  
  const Ostap::Math::Zero<std::complex<double> >     s_czero  {} ; // zero for complex doubles
  /// equality criteria for comples doubles
  const Ostap::Math::Equal_To<std::complex<double> > s_cequal {} ; // equality criteria for complex doubles
  // ==========================================================================
  // Limits? 
  // ==========================================================================
  static_assert ( std::numeric_limits<float> ::is_specialized           ,
                  "std::numeric_limits<float>  is not specialized"      ) ;
  static_assert ( std::numeric_limits<double>::is_specialized           ,
                  "std::numeric_limits<double> is not specialized"      ) ;
  static_assert ( std::numeric_limits<double>::has_denorm               ,
                  "std::numeric_limits<double> does not have denorm"    ) ;
  static_assert ( std::numeric_limits<double>::has_infinity             ,
                  "std::numeric_limits<double> does not have infinity"  ) ;
  static_assert ( std::numeric_limits<double>::has_quiet_NaN            ,
                  "std::numeric_limits<double> does not have quiet NaN" ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned short> ::is_specialized        ,
                  "std::numeric_limits<unsigned short>  is not specialized" ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned int>   ::is_specialized        ,
                  "std::numeric_limits<unsigned int>  is not specialized"   ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned long>  ::is_specialized        ,
                  "std::numeric_limits<unsigned long> is not specialized"   ) ;
  // ==========================================================================
  /** @var s_QUIETNAN
   *  quite NaN
   */
  constexpr double s_QUIETNAN = std::numeric_limits<double>::quiet_NaN () ;
  // ==========================================================================
  /** @var s_EPSILON
   *  epsilon 
   */
  constexpr double s_EPSILON = std::numeric_limits<double>::epsilon () ;
  static_assert ( 0 < s_EPSILON , "s_EPSILON is not positive!"   ) ;
  // ==========================================================================
  /** @var s_POSINF
   *  True positive infinity 
   */
  constexpr double  s_POSINF =  std::numeric_limits<double>::infinity () ;
  static_assert ( 0 < s_POSINF , "+infty/s_POSINF is not positive!"   ) ;
  // ==========================================================================
  /** @var s_NEGINF
   *  True negative infinity 
   */
  constexpr double s_NEGINF = -std::numeric_limits<double>::infinity () ;
  static_assert ( 0 > s_NEGINF , "-infty/s_NEGINF is not negative!"   ) ;  
  // ==========================================================================
  /** @var s_INFINITY
   *  representation of the almost maximal double 
   */
  constexpr double s_INFINITY  =  0.95 * std::numeric_limits<double>::max ()  ;
  constexpr double s_POSHUGE   =  0.95 * std::numeric_limits<double>::max ()  ;
  constexpr double s_NEGHUGE   = -0.95 * std::numeric_limits<double>::max ()  ;
  // ==========================================================================
  /** @var s_SMALL
   *  representation of positive "small"
   */
  constexpr double s_SMALL  = 10 * std::numeric_limits<double>::denorm_min () ;
  constexpr double s_SMALL2 =  2 * std::numeric_limits<double>::min        () ;
  // ==========================================================================
  static_assert ( 0       < s_SMALL  , "`s_SMALL' is not positive" ) ;
  static_assert ( s_SMALL < s_SMALL2 , "`S_SMALL' is too large"    ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned int>::is_specialized , 
                  "std::numeric_limits<usigned int> is not specialized" ) ;
  constexpr unsigned long s_UL_max = std::numeric_limits<unsigned int>::max () - 1 ;
  static_assert ( 1 < s_UL_max ,  "s_UL_max is not large enough!" ) ;
  // ==========================================================================
  /// small value 
  const Ostap::Math::Small<long double> s_small
  ( 2.0L * std::numeric_limits<double>::epsilon() ) ;
  // ==========================================================================
  // value that is (as float) different from zero for mULPS
  const double s_NONZERO = Ostap::Math::next_float ( 0.0f , Ostap::Math::mULPS_float + 1 ) ; 
  // ==========================================================================
  /** @var s_INFINITY_LOG_POS
  *  representation of positive INFINITY_LOG 
  */
  const double  s_INFINITY_LOG_POS = std::log ( s_INFINITY ) ;
  // ==========================================================================
  /** @var s_INFINITY_LOG_NEG
  *  representation of negative INFINITY_LOG
  */
  const double  s_INFINITY_LOG_NEG = std::log ( 2 * s_SMALL2 ) ;
  // ==========================================================================
  /** @var S_EXP_OVERFLOW 
   *  For IEEE-compatible type double, overflow is guaranteed if 709.8 < num
   *  and underflow is guaranteed if num < -708.4.
   */
  const double s_EXP_OVERFLOW = 709.8 ;
  // ===========================================================================
  /** @var S_EXP_UNDERFLOW 
   *  For IEEE-compatible type double, overflow is guaranteed if 709.8 < num
   *  and underflow is guaranteed if num < -708.4.
   */
  const double s_EXP_UNDERFLOW = -708.4 ;
  // ==========================================================================
 /** @var S_GAUSS_UNDERFLOW 
   *  For IEEE-compatible type double, overflow is guaranteed if 709.8 < num
   *  and underflow is guaranteed if num < -708.4.
   */
  const double s_GAUSS_UNDERFLOW = std::sqrt ( 2 * std::abs ( s_EXP_UNDERFLOW ) ) ;
  // ==========================================================================
  /** @var S_INFINITY_ERFC_UNDERFLOW 
   *  For the IEEE-compatible type double, underflow is guaranteed if num > 26.55.
   */
  const double s_ERFC_UNDERFLOW = 26.55 ;
  // ==========================================================================
  /// imaginary unit 
  const std::complex<double> s_j { 0.0 , 1.0 } ;
  // ==========================================================================  
  // Constants
  // ==========================================================================
#if defined ( __cplusplus ) && defined ( __cpp_lib_math_constants ) && ( 201907L <= __cpp_lib_math_constants )
  // ==========================================================================
  /// @var s_pi   constant pi 
  constexpr long double s_pi         { std::numbers::pi_v<long double>      } ;
  /// @var s_1_pi    \f$ \frac{1}{\pi} \f$  
  constexpr long double s_1_pi       { std::numbers::inv_pi_v<long double>  } ;  
  /// @var s_E   constant e  
  constexpr long double s_e          { std::numbers::e_v<long double>      } ;
  /// @var s_Mascheroni Euler-Mascheroni constant \f$ \gamma_E \f$
  constexpr long double s_Mascheroni { std::numbers::egamma_v<long double> } ;  
  /// @var s_GammaE     Euler-Mascheroni constant \f$ \gamma_E \f$
  constexpr long double s_GammaE     { std::numbers::egamma_v<long double> } ;    
  /// @var s_ln10 \f$\log 10\f$ 
  constexpr long double s_ln10       { std::numbers::ln10_v<long double>   } ;  
  /// @var s_ln2 \f$\log 2 \f$ 
  constexpr long double s_ln2        { std::numbers::ln2_v<long double>    } ;
  /// @var s_sqrt2  \f$ \sqrt{2} \f$
  constexpr long double s_sqrt2      { std::numbers::sqrt2_v<long double>  } ;
  /// @var s_sqrt3  \f$ \sqrt{3} \f$
  constexpr long double s_sqrt3      { std::numbers::sqrt3_v<long double>  } ;
  // ==========================================================================    
#else // ======================================================================
  // ==========================================================================
  /// @var s_pi  constant pi 
  const long double s_pi         = 3.141592653589793238462643383279502884L /* pi */ ; 
  /// @var s_1_pi    \f$ \frac{1}{\pi} \f$  
  const long double s_1_pi       = 1.0L / s_pi ;
  /// @var s_e   constant e  
  const long double s_3          = 2.718281828459045235360287471352L       /* e */  ; 
  /// @var s_Mascheroni Euler-Mascheroni constant \f$ \gamma_E \f$
  const long double s_Mascheroni = 0.57721566490153286060651209008240243104215933593992L ;
  /// @var s_GammaE     Euler-Mascheroni constant \f$ \gamma_E \f$
  const long double s_GammaE     = s_Mascheroni       ;
  /// @var s_ln10 \f$\log 10\f$ 
  const long double s_ln10       = std::log ( 10.0L ) ; 
  /// @var s_ln10 \f$\log 2 \f$ 
  const long double s_ln2        = std::log (  2.0L ) ; 
  /// @var s_sqrt2  \f$ \sqrt{2} \f$
  const long double s_sqrt2      = std::sqrt ( 2.0L ) ;
  /// @var s_sqrt3  \f$ \sqrt{3} \f$
  const long double s_sqrt3      = std::sqrt ( 3.0L ) ; 
  // ==========================================================================
#endif // =====================================================================
  // ==========================================================================
  /// @var s_2pi   \f$  2\pi\f$
  const long double s_2pi        = s_pi * 2 ;
  /// @var s_pi2   \f$  \pi^2 \f$
  const long double s_pi2        = s_pi  * s_pi ;
  /// @var s_pi3   \f$  \pi^3 \f$
  const long double s_pi3        = s_pi2 * s_pi ;
  /// @var s_pi4   \f$  \pi^4 \f$
  const long double s_pi4        = s_pi2 * s_pi2 ;
  /// @var s_pi_2  \f$ \frac{\pi}{2}\f$
  const long double s_pi_2       = s_pi / 2 ;  
  /// @var s_3pi_2  \f$ \frac{3\pi}{2}\f$
  const long double s_3pi_2      = s_pi * 1.5L ;
  /// @var s_pi_3  \f$ \frac{\pi}{3}\f$
  const long double s_pi_3       = s_pi / 3 ;
  /// @var s_pi_4  \f$ \frac{\pi}{4}\f$
  const long double s_pi_4       = s_pi / 4 ;
  /// @var s_pi_5  \f$ \frac{\pi}{5}\f$
  const long double s_pi_5       = s_pi / 5 ;
  /// @var s_2_pi  \f$ \frac{2}{\pi} \f$
  const long double s_2_pi       = s_1_pi * 2 ;  
  /// @var s_4_pi  \f$ \frac{4}{\pi} \f$
  const long double s_4_pi       = s_1_pi * 4 ;
  /// @var s_8_pi  \f$ \frac{8}{\pi} \f$
  const long double s_8_pi       = s_1_pi * 8 ;
  /// @var s_1_2pi \f$ \frac{1}{2\pi} \f$
  const long double s_1_2pi      = s_1_pi / 2  ;
  /// @var s_1_4pi \f$ \frac{1}{4\pi} \f$
  const long double s_1_4pi      = s_1_pi / 4  ;
  /// @var s_1_8pi \f$ \frac{1}{8\pi} \f$
  const long double s_1_8pi      = s_1_pi / 8  ;  
  /// @var s_1_sqrt2  \f$  \frac{1}{\sqrt{2}}\f$ 
  const long double s_1_sqrt2    = 1.0L / s_sqrt2 ;
  /// @var s_sqrt_pi   \f$ \sqrt { \pi } \f$
  const long double s_sqrt_pi    = std::sqrt ( s_pi ) ; 
  /// @var s_sqrt_2pi  \f$ \sqrt { 2\pi } \f$
  const long double s_sqrt_2pi   = std::sqrt ( 2.0L * s_pi ) ; 
  /// @var s+sqrt_ip_2 \f$ \sqrt { \frac { \pi }{ 2 } \f$
  const long double s_sqrt_pi_2  = std::sqrt ( s_pi_2  ) ; 
  /// @var s_sqrt_1_pi  \f$ \frac{1}{ \sqrt { \pi }\f$
  const long double s_sqrt_1_pi  = std::sqrt ( s_1_pi  ) ; 
  /// @var s_sqrt_1_2pi \f$ \frac{1}{ \sqrt { 2\pi }\f$
  const long double s_sqrt_1_2pi = 1.0L / s_sqrt_2pi ;
  /// @var s_sqrt_2_pi \f$ \sqrt{\frac{2}{\pi} }\f$
  const long double s_sqrt_2_pi  = s_sqrt2 / s_sqrt_pi ;  
  /// @var s_sqrt_1_8pi \f$ \frac{1}{ \sqrt { 8\pi }\f$
  const long double s_sqrt_1_8pi = 0.5L / s_sqrt_2pi;
  /// @var s_1_pi3   \f$  \frac{1}{\pi^2}\f$
  const long double s_1_pi2      = 1.0L / s_pi2 ;
  /// @var s_log_2pi \f$ \log 2\pi \f$
  const long double s_log_2pi    = std::log ( s_2pi ) ;
  // ==========================================================================

  // ==========================================================================
  /// precomputed value of ln(2) squared 
  const long double s_ln2_sq = s_ln2 * s_ln2 ;  
  /// precomputed value of 1/ln(10) 
  const long double s_1_ln10 = 1.0L / s_ln10 ;
  /// precomputed value of 1/ln(2) 
  const long double s_1_ln2  = 1.0L / s_ln2  ;
  // ==========================================================================

  // ==========================================================================
  // some old-fashioned names 
  // ==========================================================================  
  ///  @var s_PIHALF  \f$ \frac{\pi}{2} \f$ 
  const long double s_PIHALF     = s_pi_2 ; 
  // ==========================================================================
  /** @var s_HALFSQRTPI
   *  helper constant \f$ \frac{\sqrt{\pi}}{2}\f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  const long double s_HALFSQRTPI = 0.5L * s_sqrt_pi ;
  // ==========================================================================
  /** @var s_HALFSQRTPIi
   *  helper constant \f$ \frac{2}{\sqrt{\pi}}\f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  const long double s_HALFSQRTPIi = 1.0L / s_HALFSQRTPI  ;
  // ==========================================================================
  /** @var s_HALFSQRTPI_log
  *  helper constant \f$ \log \frac{\sqrt{\pi}}{2}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  const long double s_HALFSQRTPI_log  = std::log ( 0.5L * s_sqrt_pi ) ;
  // ==========================================================================
  /** @var s_SQRT2PISQUARED 
   *  helper constant  \f$   \sqrt{2}\pi^2\f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2016-06-11
   */
  const long double s_SQRT2PISQUARED  = s_sqrt2 * s_pi * s_pi ;
  // ==========================================================================
  /** @var s_SQRT2PISQUAREDi
   *  helper constant  \f$   \frac{1}{\sqrt{2}\pi^2}\f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2016-06-11
   */
  const long double s_SQRT2PISQUAREDi = 1.0 / s_SQRT2PISQUARED ;
  // ==========================================================================
  /** @var s_SQRT3overPI 
   *  helper constant \f$ \frac{\sqrt{3}}{\pi} \f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2016-06-14
   */
  const long double s_SQRT3overPI = s_sqrt3 / s_pi ;
  // ==========================================================================
  // Bukin & Co
  // ==========================================================================
  /** @var s_Bukin
   *  useful constant for Bukin's function
   *  \f$ \sqrt{ 2 \log 2 } \f$
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  const long double s_Bukin   = std::sqrt ( 2.0L * std::log ( 2.0L ) ) ;
  // ==========================================================================
  // Novosibirsk & Co
  // ==========================================================================
  /** @var s_Novosibirsk
   *  useful constant for evaluation of `Novosibirsk' function
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-04-19
   */
  const long double s_Novosibirsk = std::sqrt ( std::log ( 4.0L ) ) ;
  // ==========================================================================
  /** @var s_WMODE 
   *  width of the window between mean and mode:
   * |mean -mode|< sqrt(3) * sigma
   */
  const double s_WMODE            = s_sqrt3 * 1.05 ;
  // ===========================================================================
  const double s_INFINITY_LOG     = s_INFINITY_LOG_POS ;
  // ==========================================================================
  /** the protected exponent
   *  @author Vanya BELYAEV
   */
  inline long double my_exp ( const long double arg )
  {
    return 
      arg > s_INFINITY_LOG_POS ? s_INFINITY :
      arg < s_INFINITY_LOG_NEG ? s_SMALL2   : std::exp ( arg ) ;
  }
  // ==========================================================================
  /** the protected logarithm
   *  @author Vanya BELYAEV
   */
  inline long double my_log ( const long double arg )
  {
    return 
      arg <= 0          ?  -s_INFINITY_LOG :
      arg >  s_INFINITY ?   s_INFINITY_LOG :  std::log ( arg ) ;
  }
  // ==========================================================================
  /** the simple wrapper for the standard error function
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  inline long double error_func ( const long double x )
  { return 500 < x * x ? ( 0 < x ? 1 : -1 ) : std::erf ( x ) ; }
  // ==========================================================================  
  /// helper expression for erf(x)/x 
  inline double error_func_x ( const long double x ) 
  { return 0 == x || s_zero ( x ) ? s_HALFSQRTPIi : error_func( x ) / x ; }
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LOCAL_MATH_H
// ============================================================================
