// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <complex>
#include <cmath>
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/SMatrix.h"
// ============================================================================
// Ostap
// =============================================================================
#include "Ostap/GenericMatrixTypes.h"
#include "Ostap/MoreMath.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/MatrixMath.h"
// ============================================================================
/* Moebius transformation
 * \f[ f(x) = \frac{ax+b}{cx+d}\f]
 *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
 */
// ============================================================================
std::complex<double>
Ostap::Math::moebius
( const Ostap::Matrix2x2 &     m , 
  const std::complex<double>&  x )
{
  static const std::complex<double> s_1 { 1 , 0 } ; 
  return moebius ( x                 ,
		   s_1 * m ( 0 , 0 ) ,
		   s_1 * m ( 0 , 1 ) ,
		   s_1 * m ( 1 , 0 ) ,
		   s_1 * m ( 1 , 1 ) ) ;
}
// ========================================================================
/* Moebius transformation
 * \f[ f(x) = \frac{ax+b}{cx+d}\f]
 *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
 */
// ========================================================================
Ostap::Math::ValueWithError 
Ostap::Math::moebius
( const Ostap::Matrix2x2 &            m , 
  const Ostap::Math::ValueWithError&  x )
{ return moebius ( x ,
		   m ( 0 , 0 ) ,
		   m ( 0 , 1 ) ,
		   m ( 1 , 0 ) ,
		   m ( 1 , 1 ) ) ; }
// ========================================================================
/*  Moebius transformation
 * \f[ f(x) = \frac{ax+b}{cx+d}\f]
 *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
 */
// ========================================================================
double
Ostap::Math::moebius
( const Ostap::Matrix2x2 &            m , 
  const double                        x ) 
{ return moebius ( x ,
		   m ( 0 , 0 ) ,
		   m ( 0 , 1 ) ,
		   m ( 1 , 0 ) ,
		   m ( 1 , 1 ) ) ; }

// =============================================================================
//                                                                      The END 
// =============================================================================
