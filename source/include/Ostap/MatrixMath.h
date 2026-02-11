// ============================================================================
#ifndef OSTAP_MATRIXMATH_H 
#define OSTAP_MATRIXMATH_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <complex>
#include <cmath>
// ============================================================================
// Ostap
// =============================================================================
#include "Ostap/SymmetricMatrixTypes.h"
#include "Ostap/GenericMatrixTypes.h"
#include "Ostap/ValueWithError.h"
// ============================================================================
/** @file Ostap/MatrixMath.h
 *  collection of various helper math functions  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** Moebius transformation
     * \f[ f(x) = \frac{ax+b}{cx+d}\f]
     *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
     */
    std::complex<double>
    moebius
    ( const Ostap::Matrix2x2&      m , 
      const std::complex<double>&  x ) ;
    // ========================================================================
    /** Moebius transformation
     * \f[ f(x) = \frac{ax+b}{cx+d}\f]
     *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
     */
    Ostap::Math::ValueWithError 
    moebius
    ( const Ostap::Matrix2x2 &            m , 
      const Ostap::Math::ValueWithError&  x ) ;
    // ========================================================================
    /** Moebius transformation
     * \f[ f(x) = \frac{ax+b}{cx+d}\f]
     *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
     */
    double
    moebius
    ( const Ostap::Matrix2x2 &            m , 
      const double                        x ) ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MATRIXMATH_H
// ============================================================================
