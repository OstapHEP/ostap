// ============================================================================
#ifndef OSTAP_GENERICMATRIXTYPES_H 
#define OSTAP_GENERICMATRIXTYPES_H 1
// ============================================================================
// Include files
// ============================================================================
#include "Math/SMatrix.h"
// ============================================================================
/** @file Ostap/GenericMatrixTypes.h
 *  useful typedefs 
 */
// ============================================================================
namespace Ostap
{
  // ============================================================================
  typedef ROOT::Math::SMatrix<double, 1, 1> Matrix1x1;   ///< Generic 1x1 matrix (double)
  typedef ROOT::Math::SMatrix<double, 2, 2> Matrix2x2;   ///< Generic 2x2 matrix (double)
  typedef ROOT::Math::SMatrix<double, 3, 3> Matrix3x3;   ///< Generic 3x3 matrix (double)
  typedef ROOT::Math::SMatrix<double, 4, 4> Matrix4x4;   ///< Generic 4x4 matrix (double)
  typedef ROOT::Math::SMatrix<double, 5, 5> Matrix5x5;   ///< Generic 5x5 matrix (double)
  typedef ROOT::Math::SMatrix<double, 6, 6> Matrix6x6;   ///< Generic 6x6 matrix (double)
  typedef ROOT::Math::SMatrix<double, 7, 7> Matrix7x7;   ///< Generic 7x7 matrix (double)
  typedef ROOT::Math::SMatrix<double, 8, 8> Matrix8x8;   ///< Generic 8x8 matrix (double)
  typedef ROOT::Math::SMatrix<double, 9, 9> Matrix9x9;   ///< Generic 9x9 matrix (double)
  // ============================================================================
  typedef ROOT::Math::SMatrix<double, 1, 3> Matrix1x3;   ///< Generic 1x3 matrix (double)
  typedef ROOT::Math::SMatrix<double, 1, 5> Matrix1x5;   ///< Generic 1x5 matrix (double)
  typedef ROOT::Math::SMatrix<double, 1, 6> Matrix1x6;   ///< Generic 1x6 matrix (double)
  typedef ROOT::Math::SMatrix<double, 4, 3> Matrix4x3;   ///< Generic 4x3 matrix (double)
  typedef ROOT::Math::SMatrix<double, 3, 4> Matrix3x4;   ///< Generic 3x4 matrix (double)
  typedef ROOT::Math::SMatrix<double, 3, 5> Matrix3x5;   ///< Generic 3x5 matrix (double)
  typedef ROOT::Math::SMatrix<double, 3, 6> Matrix3x6;   ///< Generic 3x6 matrix (double)
  typedef ROOT::Math::SMatrix<double, 2, 3> Matrix2x3;   ///< Generic 2x3 matrix (double)
  typedef ROOT::Math::SMatrix<double, 3, 2> Matrix3x2;   ///< Generic 3x2 matrix (double)
// ============================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_GENERICMATRIXTYPES_H
// ============================================================================
