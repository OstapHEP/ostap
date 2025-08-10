

//$Id:%
// ============================================================================
// STD& STL
// ============================================================================
#include <cmath>
#include <limits>
#include <array>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_min.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Binomial.h"
#include "Ostap/MoreMath.h"
#include "Ostap/LinAlg.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_math.h"
#include "local_gsl.h"
#include "GSL_sentry.h"
// ============================================================================
/** @file
 *  implementation of functions from file Ostap/Binomial.h
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
namespace
{
  // ==========================================================================
  typedef std::array<double,3> PARAMS ;
  // ==========================================================================
  //
  double beta_minus_alpha
  ( const double x , void* p )
  {
    //
    const PARAMS* params = static_cast<PARAMS*> ( p ) ;
    //
    const double alpha = ( *params ) [ 0 ] ;
    const double beta  = ( *params ) [ 1 ] ;
    const double P     = ( *params ) [ 2 ] ;
    // 
    const double CDF = Ostap::Math::beta_cdf      ( x       , alpha , beta ) ;
    const double Q   = Ostap::Math::beta_quantile ( CDF + P , alpha , beta ) ;
    //
    return Q - x ;
  } ;    
  // ==========================================================================  
}
// ============================================================================
/*  Bayes' theorem-based interval
 *  @see M.Paterno, "Calculating efficiencies and their uncertainties", 
 *                   FERMILAB-TM-2286-CD
 *  @see https://inspirehep.net/literature/669498
 *  @see DOI: 10.2172/15017262
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::bayes_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  if ( 0 == accepted && 0 == rejected              ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  if ( 1 <= conflevel || s_equal ( conflevel , 1 ) ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  const double a = accepted ;
  const double e = a / ( a + rejected ) ;
  if ( 0 >= conflevel || s_zero ( conflevel ) ) { return std::make_pair ( e , e ) ; }
  //
  if ( !accepted ) { return std::make_pair ( 0.0 , 1 - std::pow ( 1 - conflevel , 1.0 / ( rejected + 1 ) )       ) ; }
  if ( !rejected ) { return std::make_pair (           std::pow ( 1 - conflevel , 1.0 / ( accepted + 1 ) ) , 1.0 ) ; }  
  ///
  const double par_alpha = accepted + 1 ;
  const double par_beta  = rejected + 1 ;
  //
  // The interval is:
  //   0<= alpha <= e <= beta <= 1, such that  CDF ( beta ) - CDF ( alpha ) == cnpflevel 
  //
  /// CDF for Beta-distribution
  auto CDF = [par_alpha,par_beta] ( const double x ) -> double
  { return Ostap::Math::beta_cdf      ( x , par_alpha , par_beta ) ; } ;
  /// Quantile for Beta-distribution  
  auto Q   = [par_alpha,par_beta] ( const double p ) -> double
  { return Ostap::Math::beta_quantile ( p , par_alpha , par_beta   ) ; } ;       
  ///


  // use GSL: 
  // Ostap::Math::GSL::GSL_Error_Handler sentry ;
  Ostap::Utils::GslIgnore  sentry { true } ;

  //
  // The [alpha,beta] interval is:
  //   0 <= alpha <= e <= beta <= 1,
  //   such that  CDF ( beta ) - CDF ( alpha ) == cnpflevel
  // due to this constrain for each (fixed) alpha beta is a functiopn of alpha
  // and we need to find such alphs taht minimizes the length of [alpha,beta] interval     
  //
  
  const double e_CDF     = CDF ( e ) ;
  const double alpha_min = e_CDF <= conflevel ? 0.0 : Q ( e_CDF - conflevel ) ;
  const double alpha_max = std::min ( e , Q ( 1 - conflevel  ) ) ; 
  ///
  // Jump to GSL-world 
  //
  PARAMS params { par_alpha  , par_beta , conflevel } ;
  // ==========================================================================  
  // 
  const gsl_min_fminimizer_type    *T ;
  gsl_min_fminimizer               *s ;
  //
  gsl_function F ;
  F.function   = &beta_minus_alpha    ; 
  F.params     = &params              ;
  //
  T = gsl_min_fminimizer_brent        ;
  s = gsl_min_fminimizer_alloc ( T )  ;
  //
  gsl_min_fminimizer_set ( s, &F,  0.5 * ( alpha_min + alpha_max ) , alpha_min , alpha_max );
  //
  const double min_e     = 0.001 * std::min ( e         , 1 - e         ) ;
  const double min_c     = 0.001 * std::min ( conflevel , 1 - conflevel ) ;
  static const double DX = std::clamp ( std::min ( min_e , min_c ) , 1.e-8 , 1.e-7 ) ; 
  //
  double low   = alpha_min ;
  double high  = alpha_max ;
  double alpha = 0.5 * ( alpha_min + alpha_max ) ;
  for ( unsigned short i = 0 ; i < 50 ; ++i )
    {
      int status = gsl_min_fminimizer_iterate   ( s ) ;
      //
      if ( status ) { gsl_error ( "Error from gsl_min_fminimizer_iterate" , __FILE__ , __LINE__ , status ) ; }
      //
      alpha      = gsl_min_fminimizer_x_minimum ( s ) ;
      low        = gsl_min_fminimizer_x_lower   ( s ) ;
      high       = gsl_min_fminimizer_x_upper   ( s ) ;
      //
      if ( GSL_SUCCESS == gsl_min_test_interval ( low , high , DX , DX ) ) { break ; } 
      else if ( std::abs ( high - low ) <= DX                            ) { break ; }          
      //
    }
  //
  gsl_min_fminimizer_free( s ) ;
  // upper edge of the interval 
  const double beta = Q ( CDF ( alpha ) + conflevel ) ;
  // result 
  return std::make_pair ( alpha , beta ) ;  
}
// ============================================================================
//                                                                      The END 
// ============================================================================


