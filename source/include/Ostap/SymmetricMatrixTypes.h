// ============================================================================
#ifndef OSTAP_SYMMETRICMATRIXTYPES_H 
#define OSTAP_SYMMETRICMATRIXTYPES_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/SMatrix.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  typedef ROOT::Math::SMatrix<double, 1, 1, 
                              ROOT::Math::MatRepSym<double,1> > SymMatrix1x1; ///< Symmetrix 1x1 matrix (double)
  typedef ROOT::Math::SMatrix<double, 2, 2,
                              ROOT::Math::MatRepSym<double,2> > SymMatrix2x2; ///< Symmetrix 2x2 matrix (double)
  typedef ROOT::Math::SMatrix<double, 3, 3,
                              ROOT::Math::MatRepSym<double,3> > SymMatrix3x3; ///< Symmetrix 3x3 matrix (double)
  typedef ROOT::Math::SMatrix<double, 4, 4,
                              ROOT::Math::MatRepSym<double,4> > SymMatrix4x4; ///< Symmetrix 4x4 matrix (double)
  typedef ROOT::Math::SMatrix<double, 5, 5,
                              ROOT::Math::MatRepSym<double,5> > SymMatrix5x5; ///< Symmetrix 5x5 matrix (double)
  typedef ROOT::Math::SMatrix<double, 6, 6,
                              ROOT::Math::MatRepSym<double,6> > SymMatrix6x6; ///< Symmetrix 6x6 matrix (double)
  typedef ROOT::Math::SMatrix<double, 7, 7,
                              ROOT::Math::MatRepSym<double,7> > SymMatrix7x7; ///< Symmetrix 7x7 matrix (double)
  typedef ROOT::Math::SMatrix<double, 8, 8,
                              ROOT::Math::MatRepSym<double,8> > SymMatrix8x8; ///< Symmetrix 8x8 matrix (double)
  typedef ROOT::Math::SMatrix<double, 9, 9,
                              ROOT::Math::MatRepSym<double,9> > SymMatrix9x9; ///< Symmetrix 9x9 matrix (double)
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_SYMMETRICMATRIXTYPES_H
// ============================================================================
