// ============================================================================
#ifndef OSTAP_GSL_LINALGUTILS_H 
#define OSTAP_GSL_LINALGUTILS_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TVectorT.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/MatrixAsBuffer.h" 
#include "Ostap/LinAlg.h" 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Math
  {
    // =======================================================================
    namespace GSL
    {
      // =====================================================================
      // SMatrix -> GSL matrix
      // =====================================================================
      
      // =====================================================================
      /// create GSL matrix from the SMatrix     
      template <class T, unsigned int D1, unsigned int D2,class R>      
      inline Matrix matrix ( const ROOT::Math::SMatrix<T,D1,D2,R>& m )
      { return Matrix { D1 , D2  , Ostap::Utils::buffer ( m ) } ; }      
      // ====================================================================
      
      // =====================================================================
      /// convert T-matrix into GSL matrix 
      Matrix matrix ( const TMatrixT<float>&     m ) ;    
      /// convert T-matrix into GSL matrix 
      Matrix matrix ( const TMatrixT<double>&    m ) ;      
      // =====================================================================
      /** convert TMatrixTSym into GSL matrix
       *  @attention It reads only the upper triangle part and fill both upper & lower parts
       */
      Matrix matrix ( const TMatrixTSym<float>&  m ) ;
      // =====================================================================
      /** convert TMatrixTSym into GSL matrix
       *  @attention It reads only the upper triangle part and fill both upper & lower parts
       */
      Matrix matrix ( const TMatrixTSym<double>& m ) ;
      // =====================================================================
      
      // =====================================================================    
      /// S-vector -> GSL Vector 
      template <class T, unsigned int D>
      inline Vector vector ( const ROOT::Math::SVector<T,D>& v )
      { return Vector ( v.begin () , v.end () ) ; }
      // =====================================================================    
      /// T-vector -> GSL Vector 
      template <class T>
      inline Vector vector ( const TVectorT<T>& v )
      { return Vector ( v.begin () , v.end () ) ; }
      // =====================================================================    
      
      // =====================================================================          
    } //                                 The end of namesapce Ostap::Math::GSL
    // =======================================================================
  } //                                        The end of namesapce Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_GSL_LINALGUTILS_H
// ============================================================================
//                                                                      The END 
// ============================================================================
