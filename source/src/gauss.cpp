// ============================================================================
// Include files 
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_exp.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
// ============================================================================
// local
// ============================================================================
#include "local_math.h"
#include "local_gsl.h"
#include "gauss.h"
// ============================================================================
namespace 
{
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
} //                                             The end of anonymous namespace 
// ============================================================================
/*  get the gaussian integral numerically
 *  \f[ f = \int_a^b \mathrm{e}^ {-\alpha x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param a     the low  integration limit
 *  @param b     the high integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ============================================================================
double Ostap::Math::details::gaussian_int_num
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
  { gsl_error ( "Ostap::Math::gaussian_int" , __FILE__ , __LINE__ , ierror ) ; }
  //
  return result ;
}
// ==========================================================================
/*  get the gaussian integral:
 *  \f[ f = \int_a^b \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
 *  @attention note the sign for alpha-term!
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param a     the low  integration limit
 *  @param b     the high integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ==========================================================================
double Ostap::Math::details::gaussian_int
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
/*  get the gaussian integral:
 *  \f[ f = \int_a^{\inf} \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
 *  @param alpha the alpha parameter
 *  @param beta  the beta  parameter
 *  @param a     the low integration limit
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-05-23
 */
// ==========================================================================
double Ostap::Math::details::gaussian_int_R 
( const double alpha ,
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
/*  get the exponential integral 
 *  \f[ f = \int_a^b \exp { \beta x } \mathrm{d}x \f]
 *  @param beta  the beta  parameter
 *  @param a     the low  integration limit
 *  @param b     the high integration limit
 */
// ==========================================================================
double Ostap::Math::details::exponent_int
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
// ============================================================================
//                                                                      The END  
// ============================================================================
