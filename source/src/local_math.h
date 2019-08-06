// ============================================================================
#ifndef LOCAL_MATH_H 
#define LOCAL_MATH_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
/// local namespace to hide all tehcnical symbols 
#if  __cplusplus == 201703L
  #ifndef CONSTEXP
    #define CONSTEXP  constexpr 
  #endif
#else 
  #ifndef CONSTEXP
    #define CONSTEXP  const
  #endif
#endif

namespace 
{
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double>            s_equal{} ; // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>                s_zero {} ; // zero for doubles
  /// zero for vectors 
  const Ostap::Math::Zero< std::vector<double> > s_vzero{} ; // zero for vectors
  // ==========================================================================
  // Limits? 
  // ==========================================================================
  static_assert ( std::numeric_limits<float> ::is_specialized        ,
                  "std::numeric_limits<float>  is not specialized"   ) ;
  static_assert ( std::numeric_limits<double>::is_specialized        ,
                  "std::numeric_limits<double> is not specialized"   ) ;
  static_assert ( std::numeric_limits<double>::has_denorm            ,
                  "std::numeric_limits<double> doed not have denorm" ) ;
  // ==========================================================================
  /** @var s_INFINITY
   *  representation of positive INFINITY
   */
  constexpr double s_INFINITY  = 0.8 * std::numeric_limits<double>::max ()  ;
  // ==========================================================================
  /** @var s_SMALL
   *  representation of positive "small"
   */
  constexpr double s_SMALL  = 10 * std::numeric_limits<double>::denorm_min () ;
  constexpr double s_SMALL2 =  2 * std::numeric_limits<double>::min        () ;
  // ==========================================================================
  static_assert ( 0       < s_SMALL  , "``s_SMALL'' is not positive" ) ;
  static_assert ( s_SMALL < s_SMALL2 , "``S_SMALL'' is too large"    ) ;
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
  /** @var s_INFINITY_LOG_POS
  *  representation of positive INFINITY_LOG 
  */
  CONSTEXP double  s_INFINITY_LOG_POS = std::log ( s_INFINITY ) ;
  /** @var s_INFINITY_LOG
  *  representation of positive INFINITY_LOG 
  */
  /** @var s_INFINITY_LOG_NEG
  *  representation of negative INFINITY_LOG
  */
  CONSTEXP double  s_INFINITY_LOG_NEG = std::log ( 2 * s_SMALL2 ) ;

  // ==========================================================================
  // Constants
  // ==========================================================================
  /** @var s_LN10
  *  \f$\ln(10)\f$ 
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-05-23
  */
  CONSTEXP double       s_LN10 = std::log ( 10 ) ;
  // ==========================================================================
  /** @var s_SQRTPIHALF
  *  helper constant \f$ \sqrt{\frac{\pi}{2}}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double  s_SQRTPIHALF = std::sqrt( M_PI_2 ) ;
  // ==========================================================================
  /** @var s_SQRTPI
  *  helper constant \f$ \sqrt{\pi}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
  *  @date 2016-06-11
  */
  CONSTEXP double  s_SQRTPI = std::sqrt( M_PI ) ;
  // ==========================================================================
  /** @var s_SQRT2PI
  *  helper constant \f$ \sqrt{2\pi}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double  s_SQRT2PI    =       std::sqrt ( 2 * M_PI ) ;
  // ===========================================================================
  /** @var s_SQRT2PIi
  *  helper constant \f$ \frac{1}{\sqrt{2\pi}}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double  s_SQRT2PIi    =      1./s_SQRT2PI ;
  // ===========================================================================
  /** @var s_SQRT2 
  *  helper constant \f$\sqrt{2}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2013-08-25
  */
  CONSTEXP double  s_SQRT2      =       std::sqrt ( 2.0 )      ;
  // ===========================================================================
  /** @var s_SQRT2i 
  *  helper constant \f$\frac{1}{\sqrt{2}}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2013-08-25
  */
  CONSTEXP double  s_SQRT2i      =       1/std::sqrt ( 2.0 )    ;
  // ===========================================================================
  /** @var s_HALFSQRTPI
  *  helper constant \f$ \frac{\sqrt{\pi}}{2}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double  s_HALFSQRTPI = 0.5 * std::sqrt(     M_PI ) ;
  // ==========================================================================
  /** @var s_HALFSQRTPIi
  *  helper constant \f$ \frac{2}{\sqrt{\pi}}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double  s_HALFSQRTPIi = 1/s_HALFSQRTPI  ;
  // ==========================================================================
  /** @var s_SQRT3 
  *  helper constant \f$ \sqrt{3} \f$ 
  *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
  *  @date 2015-08-21
  */
  CONSTEXP double  s_SQRT3 = std::sqrt ( 3.0 ) ;
  // ==========================================================================
  /** @var s_HALFSQRTPI_log
  *  helper constant \f$ \log \frac{\sqrt{\pi}}{2}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double s_HALFSQRTPI_log  = std::log ( 0.5 * std::sqrt(     M_PI )  ) ;
  // ==========================================================================
  /** @var s_SQRT2PISQUARED 
  *  helper constant  \f$   \sqrt{2}\pi^2\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2016-06-11
  */
  CONSTEXP double s_SQRT2PISQUARED  = std::sqrt(2.0)*M_PI*M_PI ;
  // ==========================================================================
  /** @var s_SQRT2PISQUAREDi
  *  helper constant  \f$   \frac{1}{\sqrt{2}\pi^2}\f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2016-06-11
  */
  CONSTEXP double s_SQRT2PISQUAREDi = 1.0/(std::sqrt(2.0)*M_PI*M_PI) ;
  // ==========================================================================
  /** @var s_SQRT3overPI 
  *  helper constant \f$ \frac{\sqrt{3}}{\pi} \f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2016-06-14
  */
  CONSTEXP double s_SQRT3overPI = std::sqrt(3.0)/M_PI ;
  // ==========================================================================
  // Bukin & Co
  // ==========================================================================
  /** @var s_Bukin
  *  useful constant for Bukin's function
  *  \f$ \sqrt{ 2 \log 2 } \f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double s_Bukin   = std::sqrt ( 2.0 * std::log ( 2.0 ) ) ;
  // ==========================================================================
  /** @var s_ln2
  *  useful constant for Bukin's function
  *  \f$ \log 2 \f$
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  CONSTEXP double s_ln2 = std::log ( 2.0 ) ;
 

  // ==========================================================================
  // Novosibirsk & Co
  // ==========================================================================
  /** @var s_Novosibirsk
  *  useful constant for evaliuation of ``Novosibirsk'' function
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-04-19
  */
  const double s_Novosibirsk = std::sqrt ( std::log ( 4.0 ) ) ;

  const     double s_INFINITY_LOG     = s_INFINITY_LOG_POS ;

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
  /** evaluate the helper function  \f[ f = \frac{\log(1+x)}{x} \f]
   *  it allows to calculate Bukin' function in efficient and regular way
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  inline long double x_log ( const long double x  ) // , double precision = s_PRECISION )
  {
    //
    if      ( s_equal   ( x , 0 )       ) { return 1 ; } // RETURN
    else if ( x <= -1                   ) { return 0 ; } // RETURN ?
    else if ( s_equal   ( x , -1 )      ) { return 0 ; } // RETURN ?
    else if ( x < 0.9                   ) { return std::log1p ( x ) / x ; }
    //
    // the generic evaluation
    //
    return std::log ( 1.0 + x ) / x ;
  }
  // ==========================================================================
  /** get "reduced" exponent
   *  \f[ f(x) = \frac{ e^{x}-1 }{x} \f]
   */
  inline long double reduced_exp ( const long double x ) 
  { return ( 0 == x || s_zero ( x ) ) ? 1.0 : std::expm1 ( x ) / x ; }
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LOCAL_MATH_H
// ============================================================================
