// ===========================================================================
// Include files
// ===========================================================================
// StD&STL
// ===========================================================================
#include <array>
#include <limits>
// ===========================================================================
// Ostap
// ===========================================================================
#include "Ostap/MoreMath.h"
#include "Ostap/CornishFisher.h"
// ===========================================================================
// local
// ===========================================================================
#include "local_math.h"
// ===========================================================================
/** @file 
 *  Implementation of Cornish-Fisher asymprotic expansion
 *  @see Ostap/CornishFisher.h
 *  @see Ostap::Math::cornish_fisher 
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 *  @date 2025-07-24
 */
// ===========================================================================
/* An asymptotic Cornish-Fisher expansion  
 * of the quantile function in terms of cumulants  
 * 
 * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
 * 
 * @param p  probablity, \f$ 0 < p < 1 \f$
 * @param mu mean value (1st cumulant)
 * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
 * @return approximate quantile for probabilty p   
 */
// ===========================================================================
double Ostap::Math::cornish_fisher 
( const double p        , 
  const double mu       , 
  const double sigma    ) 
  {
    if ( s_zero  ( p )     ) { return -std::numeric_limits<double>::max()       ; }
    if ( s_equal ( p , 1 ) ) { return  std::numeric_limits<double>::max()       ; }
    if ( p < 0 || 1 < p    ) { return  std::numeric_limits<double>::quiet_NaN() ; }
    //
    const double x = Ostap::Math::probit ( p ) ;
    //
    const double wp = x ;
    //
    return mu + std::abs ( sigma ) * wp ;
  }
// ===========================================================================
/* An asymptotic Cornish-Fisher expansion  
 * of the quantile function in terms of cumulants  
 * 
 * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
 * 
 * @param p  probablity, \f$ 0 < p < 1 \f$
 * @param mu mean value (1st cumulant)
 * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
 * @param skewness  skewness  , 3rd cumulant 
 * @return approximate quantile for probabilty p   
 */
// ===========================================================================
double Ostap::Math::cornish_fisher 
( const double p        , 
  const double mu       , 
  const double sigma    , 
  const double skewness ) 
  {
    if ( s_zero  ( p )     ) { return -std::numeric_limits<double>::max()       ; }
    if ( s_equal ( p , 1 ) ) { return  std::numeric_limits<double>::max()       ; }
    if ( p < 0 || 1 < p    ) { return  std::numeric_limits<double>::quiet_NaN() ; }
   //
    const double x = Ostap::Math::probit ( p ) ;
    //
    const std::array<double,3> s_He { 1.0 , x , x * x - 1 } ;
    //
    const double h1 = s_He[2] / 6 ;
    //
    const double g1 = skewness ;
    //
    double  wp  = x ; 
    wp         += g1 * h1 ;
    //
    return mu + std::abs ( sigma ) * wp ;
  }
// ===========================================================================
/* An asymptotic Cornish-Fisher expansion  
 * of the quantile function in terms of cumulants  
 * 
 * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
 * 
 * @param p  probablity, \f$ 0 < p < 1 \f$
 * @param mu mean value (1st cumulant)
 * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
 * @param skewness  skewness  , 3rd cumulant 
 * @param kurtosis   excess kurtosis, 4thcumulant  
 * @return approximate quantile for probabilty p   
 */
// ===========================================================================
double Ostap::Math::cornish_fisher 
( const double p        , 
  const double mu       , 
  const double sigma    , 
  const double skewness , 
  const double kurtosis ) 
  {
    if ( s_zero  ( p )     ) { return -std::numeric_limits<double>::max()       ; }
    if ( s_equal ( p , 1 ) ) { return  std::numeric_limits<double>::max()       ; }
    if ( p < 0 || 1 < p    ) { return  std::numeric_limits<double>::quiet_NaN() ; }
    // 
    const double x = Ostap::Math::probit ( p ) ;
    //
    std::array<double,4> s_He ;
    s_He [ 0 ] = 1.0 ;
    s_He [ 1 ] = x   ;
    for ( unsigned short n = 1 ; n < 3 ; ++n ) { s_He [ n + 1 ] = x  * s_He [ n ] - n * s_He [ n ] ;  }
    //
    const double h1  =         s_He [ 2 ]               /  6 ;    
    const double h2  =         s_He [ 3 ]               / 24 ;
    const double h11 = - ( 2 * s_He [ 3 ] + s_He [1 ] ) / 36 ; 
    //
    const double g1 = skewness ;
    const double g2 = kurtosis ; 
    //
    double wp  = x ; 
    wp        += g1 * h1 ; 
    wp        += g2 * h2 + g1 * g1 * h11 ;
    //
    return mu + std::abs ( sigma ) * wp ;   
    //
  }
// ===========================================================================
/*  An asymptotic Cornish-Fisher expansion  
 * of the quantile function in terms of cumulants  
 * 
 * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
 * 
 * @param p  probablity, \f$ 0 < p < 1 \f$
 * @param mu mean value (1st cumulant)
 * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
 * @param skewness  skewness  , 3rd cumulant 
 * @param kurtosis  excess kurtosis, 4th cumulant
 * @param kappa5    5th cumulant   
 * @return approximate quantile for probabilty p   
 */
// ===========================================================================
double Ostap::Math::cornish_fisher 
( const double p        , 
  const double mu       , 
  const double sigma    , 
  const double skewness , 
  const double kurtosis , 
  const double kappa5   ) 
  {
    if ( s_zero  ( p )     ) { return -std::numeric_limits<double>::max()       ; }
    if ( s_equal ( p , 1 ) ) { return  std::numeric_limits<double>::max()       ; }
    if ( p < 0 || 1 < p    ) { return  std::numeric_limits<double>::quiet_NaN() ; }
    //
    const double x = Ostap::Math::probit ( p ) ;

    std::array<double,5> s_He ;
    s_He [ 0 ] = 1.0 ;
    s_He [ 1 ] = x   ;
    for ( unsigned short n = 1 ; n < 4 ; ++n ) { s_He [ n + 1 ] = x  * s_He [ n ] - n * s_He [ n ] ;  }
    //
    const double h1   =          s_He [ 2 ]                     /   6 ;    
    const double h2   =          s_He [ 3 ]                     /  24 ;
    const double h11  = - (  2 * s_He [ 3 ] +      s_He [ 1 ] ) /  36 ; 
    const double h3   =          s_He [ 4 ]                     / 120 ; 
    const double h12  = - (      s_He [ 4 ] +      s_He [ 2 ] ) /  24 ; 
    const double h111 =   ( 12 * s_He [ 4 ] + 19 * s_He [ 2 ] ) / 324 ; 
    //
    const double g1 = skewness ;
    const double g2 = kurtosis ;
    const double g3 = kappa5 / sigma ;
    //
    double wp  = x ; 
    wp        +=  g1 * h1 ; 
    wp        +=  g2 * h2 + g1 * g1 * h11 ;
    wp        +=  g3 * h3 + g1 * g2 * h12 + g1 * g1 * g1 * h111 ;  
    //
    return mu + std::abs ( sigma ) * wp ;   
    //
  }
// ===========================================================================
//                                                                     The END
// =========================================================================== 