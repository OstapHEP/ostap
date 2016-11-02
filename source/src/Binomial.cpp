//$Id:%
// ============================================================================
// STD& STL
// ============================================================================
#include <cmath>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_gamma.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Binomial.h"
#include "Ostap/MoreMath.h"
// ============================================================================
/** @file
 *  implementation of funcntions from file Ostap/Binomial.h
 * 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *
 *  Version           $Revision$
 *  Last modification $Date$
 *                 by $Author$
 */
// ============================================================================
namespace
{
  // ===========================================================================
  /** get quantile function for standard normal distribution
   *  @see http://en.wikipedia.org./wiki/Probit
   */
  inline double  _probit_ ( const double alpha  )
  { return gsl_cdf_ugaussian_Pinv ( alpha ) ; }
  // ==========================================================================
  inline double ibeta_inv  ( const double a ,
                             const double b ,
                             const double p )
  { return gsl_cdf_beta_Pinv      ( p , a , b ) ; } 
  // ==========================================================================
}
// ============================================================================
/* get quantile function for standard normal distribution
 *  @see http://en.wikipedia.org./wiki/Probitq
 *  @param alpha argument    \f$  0<\alpha<1 \f$  
 *  @return quantile value 
 */
// ============================================================================
double Ostap::Math::probit ( const double alpha  )
{
  return
    alpha <= 0 ? -1.0/0.0  :
    alpha >= 1 ? +1.0/0.0  : _probit_ ( alpha  ) ;  
}
// ============================================================================
/* normal approximation interval for binomial proportion/efficiency 
 *  ( "Wald test")
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::wald_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  //
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )     ; }
  //
  const double alpha = 1 - conflevel ;
  const double  z    = _probit_ ( 1 - 0.5  * alpha ) ;
  //
  // adjust cases of p=0 and p=1
  //
  const double p1 =
    0 == accepted ?   ( a + 1 ) / ( a + rejected + 1 ) :
    0 == rejected ?   ( a     ) / ( a + rejected + 1 ) : p ;
  //
  const double  dx = std::sqrt ( p1 * ( 1 - p1 ) / ( accepted + rejected ) ) ;
  //
  const double low  = std::max ( 0.0 , p - z * dx ) ;
  const double high = std::min ( 1.0 , p + z * dx ) ;
  //
  return std::make_pair ( low , high ) ;
  //
}
// ============================================================================
/*  Wilson score interval for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::wilson_score_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )     ; }
  //
  const double alpha = 1 - conflevel ;
  const double  z  = _probit_ ( 1 - 0.5  * alpha ) ;
  //
  const long double n = accepted + rejected ;
  const double dx  = std::sqrt ( p * ( 1 - p ) / n + z * z / ( 4 * n * n ) ) ;
  const double f1  = p + z * z / ( 2 * n ) ;
  const double f2  = 1 / ( 1 + z * z / n ) ;
  //
  const double low  = f2 * ( f1 - z * dx ) ;
  const double high = f2 * ( f1 + z * dx ) ;
  //
  return std::make_pair ( low  , high ) ;  
}
// ============================================================================
/*  Wilson score interval with continuity correction for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::wilson_score_continuity_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )     ; }
  //
  const double alpha = 1 - conflevel ;
  const double  z  = _probit_ ( 1 - 0.5  * alpha ) ;
  //
  const double      np =  accepted ;
  const long double n  = accepted + rejected ;
  //
  const double dx  = z * std::sqrt ( z*z - 1 / n + 4 * np * ( 1 - p ) + 4*p-2  ) + 1 ;
  //
  const double f1  =  2 * np + z * z ;
  const double f2  =  1 / ( 2 * ( n + z*z ) ) ;
  //
  const double low  = std::max ( 0.0 , f2 * ( f1 - dx ) ) ;
  const double high = std::min ( 1.0 , f2 * ( f1 + dx ) ) ;
  //
  return std::make_pair ( low  , high ) ;  
}
// ============================================================================
/*  ArcSin interval with continuity correction for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::arcsin_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )     ; }
  //
  const double alpha = 1 - conflevel ;
  const double  z  = _probit_ ( 1 - 0.5  * alpha ) ;
  //
  // const double p1  =
  //  0 == accepted ? ( a + 1 ) / ( a + rejected     ) :
  //  0 == rejected ? ( a     ) / ( a + rejected + 1 ) : p ;
  const double p1  =
    ( a + 3.0 / 8 ) / ( a + rejected + 3.0 / 4 ) ;
  // 0 == accepted ? ( a + 1 ) / ( a + rejected     ) :
  // 0 == rejected ? ( a     ) / ( a + rejected + 1 ) : p ;
  //
  const double asp   = std::asin ( std::sqrt ( p1 ) ) ;
  //
  const double n  = accepted + rejected ;
  //
  const double dx  = z / ( 0.5 * std::sqrt( n ) ) ;
  //
  const double l1  =
    0 == accepted ? 0.0 : std::sin ( asp - dx ) ;
  const double h1  =
    0 == rejected ? 1.0 : std::sin ( asp + dx ) ;
  //
  return std::make_pair ( l1 * l1  , h1 * h1  ) ;  
}
// ============================================================================
/*  Agresti-Coull interval with continuity correction for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::agresti_coull_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )     ; }
  //
  const double alpha = 1 - conflevel ;
  const double  z  = _probit_ ( 1 - 0.5  * alpha ) ;
  //
  const double n1 = accepted + rejected + z * z ;
  const double p1 = ( accepted + 0.5 * z * z ) / n1 ;
  //
  const double dx   = z * std::sqrt ( p1 * ( 1 - p1 ) / n1 ) ;
  //
  const double low  = std::max ( 0.0 , p1 - dx ) ;
  const double high = std::min ( 1.0 , p1 + dx ) ;
  //
  return std::make_pair ( low , high ) ;  
}
// ============================================================================
/*  Jeffreys interval for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::jeffreys_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  //
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )     ; }
  //
  const double alpha = 1 - conflevel ;
  //
  const double low   =
    0 == accepted ? 0.0 : ibeta_inv (  accepted + 0.5 , rejected + 0.5 ,     0.5 * alpha ) ;
  const double high  =
    0 == rejected ? 1.0 : ibeta_inv (  accepted + 0.5 , rejected + 0.5 , 1 - 0.5 * alpha ) ;
  //
  return std::make_pair ( low  , high ) ;
}
// ============================================================================
/* Clopper-Pearson interval for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::clopper_pearson_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  //
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )     ; }
  //
  const double alpha = 1 - conflevel ;
  //
  const double low   =
    0 == accepted ? 0.0 : ibeta_inv (  accepted     , rejected + 1 ,     0.5 * alpha ) ;
  const double high  =
    0 == rejected ? 1.0 : ibeta_inv (  accepted + 1 , rejected     , 1 - 0.5 * alpha ) ;
  //
  return std::make_pair ( low  , high ) ;
}


// ============================================================================
// The END 
// ============================================================================

