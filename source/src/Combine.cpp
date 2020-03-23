// $Id$ 
// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/GenericVectorTypes.h"
#include "Ostap/Combine.h"
// ============================================================================
/** @file
 *  Implementation file for utilitied fromm LHCbMath/Comiiner.h file  
 *  @see Ostap::Math::Combine
 *  Helper utility to combine   correlated measurements 
 *  @see P.Avery "Combining measurements with correlated errors", CBX 95 55
 *  @see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
 *  @see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
 *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-09-28
 */
// ============================================================================
/*  combine two measurements <code>x</code> and <code>y</code>
 *  with covarinace matrix <code>cov</code>
 *  @param x (INPUT) the first  measurement 
 *  @param y (INPUT) the second measurement 
 *  @param cov2 (INPUT) covariance matrix 
 *  @return combined result
 *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-09-28
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::combine  
( const double               x   , 
  const double               y   , 
  const Ostap::SymMatrix2x2& cov ) 
{
  Ostap::Vector2 data ( x,y) ;
  Ostap::Math::Combine<2> combiner  ( data , cov ) ;
  return combiner.result() ;
}
// ============================================================================
/* combine two measurements <code>x1</code> and <code>x2</code>
 *  using correlation coefficient <code>rho</code>:  \f$-1\le\rho\le1\f$
 *  @param x1  (INPUT) the first  measurement 
 *  @param x2  (INPUT) the second measurement 
 *  @param rho (INPUT) correlation coefficient 
 *  @return combined result
 *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-09-28
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::combine  
( const Ostap::Math::ValueWithError& x1  ,
  const Ostap::Math::ValueWithError& x2  , 
  const double                       rho ) 
{
  //
  Ostap::Vector2      data ( x1.value() , x2.value() ) ;
  Ostap::SymMatrix2x2 cov  ;
  cov ( 0 , 0 ) =     x1.cov2() ;
  cov ( 0 , 1 ) = rho*std::sqrt( x1.cov2() * x2.cov2() ) ;
  cov ( 1 , 1 ) =     x2.cov2() ;
  //
  Ostap::Math::Combine<2> combiner  ( data , cov ) ;
  return combiner.result() ;
}
// =========================================================================
/*  combine two measurements <code>x1</code> and <code>x2</code>
 *  using theie "statistical" uncertainties (assumed to be uncorrelated) 
 *  and a covariance matrix of "systematic" uncertainnties
 *  @param x1   (INPUT) the first  measurement 
 *  @param x2   (INPUT) the second measurement 
 *  @param syst (INPUT) covariance matrix of systematic uncertainnties  
 *  @return combined result
 *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-09-28
 */
// =========================================================================
Ostap::Math::ValueWithError 
Ostap::Math::combine  
( const Ostap::Math::ValueWithError& x1   ,
  const Ostap::Math::ValueWithError& x2   , 
  const Ostap::SymMatrix2x2&         syst ) 
{
  std::array<double,2> data { x1.value() , x2.value() } ;
  Ostap::SymMatrix2x2 cov  ( syst ) ;
  cov ( 0 , 0 ) += x1.cov2() ;
  cov ( 1 , 1 ) += x2.cov2() ;
  //
  Ostap::Math::Combine<2> combiner  ( data , cov ) ;
  return combiner.result() ;
}
// =========================================================================
/* combine three measurements <code>x1</code>, <code>x2</code> and <code>x3</code>
 *  using their "statistical" uncertainties (assumed to be uncorrelated) 
 *  and a covariance matrix of "systematic" uncertainties
 *  @param x1   (INPUT) the first  measurement 
 *  @param x2   (INPUT) the second measurement 
 *  @param x3   (INPUT) the third  measurement 
 *  @param syst (INPUT) covariance matrix of systematic uncertainties  
 *  @return combined result
 *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-09-28
 */
// =========================================================================
Ostap::Math::ValueWithError 
Ostap::Math::combine  
( const Ostap::Math::ValueWithError& x1   ,
  const Ostap::Math::ValueWithError& x2   , 
  const Ostap::Math::ValueWithError& x3   , 
  const Ostap::SymMatrix3x3&         syst ) 
{
  std::array<double,3> data { x1.value() , x2.value() , x3.value () }  ;
  Ostap::SymMatrix3x3 cov  ( syst ) ;
  cov ( 0 , 0 ) += x1.cov2() ;
  cov ( 1 , 1 ) += x2.cov2() ;
  cov ( 2 , 2 ) += x3.cov2() ;
  //
  Ostap::Math::Combine<3> combiner  ( data , cov ) ;
  return combiner.result() ;
}
// ============================================================================
/*  combine four measurements:
 *  - <code>x1</code>, 
 *  - <code>x2</code>,
 *  - <code>x3</code> and 
 *  - <code>x4</code>
 *  using their "statistical" uncertainties (assumed to be uncorrelated) 
 *  and a covariance matrix of "systematic" uncertainties
 *  @param x1   (INPUT) the first  measurement 
 *  @param x2   (INPUT) the second measurement 
 *  @param x3   (INPUT) the third  measurement 
 *  @param x4   (INPUT) the fourth measurement 
 *  @param syst (INPUT) covariance matrix of systematic uncertainties  
 *  @return combined result
 *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-09-28
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::combine  
( const Ostap::Math::ValueWithError& x1   ,
  const Ostap::Math::ValueWithError& x2   , 
  const Ostap::Math::ValueWithError& x3   , 
  const Ostap::Math::ValueWithError& x4   , 
  const Ostap::SymMatrix4x4&         syst ) 
{
  std::array<double,4> data { x1.value() , x2.value() , x3.value () , x4.value() } ;
  Ostap::SymMatrix4x4 cov  ( syst ) ;
  cov ( 0 , 0 ) += x1.cov2() ;
  cov ( 1 , 1 ) += x2.cov2() ;
  cov ( 2 , 2 ) += x3.cov2() ;
  cov ( 3 , 3 ) += x4.cov2() ;
  //
  Ostap::Math::Combine<4> combiner  ( data , cov ) ;
  return combiner.result() ;

}
// ============================================================================
//                                                                      The END 
// ============================================================================
