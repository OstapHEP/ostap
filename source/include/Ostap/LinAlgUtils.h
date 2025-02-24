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
#include "Ostap/LinAlg.h" 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace GSL
  {
    // ========================================================================
    // SMatrix -> GSL matrix 
    // =======================================================================
    /// create GSL matrix from the S-     
    template <class T, unsigned int D1, unsigned int D2,class R>      
      inline Matrix matrix ( const ROOT::Math::SMatrix<T,D1,D2,R>& m )
    {
      Matrix result { D1 , D2 } ;
      for ( unsigned int i = 0 ; i < D1 ; ++i )
        {
          for ( unsigned j = 0 ; j < D2 ; ++j )
            {
              const double mij = m ( i , j ) ;
              result.set ( i , j , mij ) ;
            }
        }
      return result ;
    }
    // =======================================================================
    /// Symmetric S-matrix -> GSL matrix 
    template <class T, unsigned int D>      
      inline Matrix matrix
      ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&  m )
    {
      Matrix result { D } ;
      for ( unsigned int i = 0 ; i < D ; ++i )
        {
          result.set ( i , i , m ( i , i ) ) ;
          for ( unsigned j = i + 1  ; j < D ; ++j )
            {
              const double mij = m ( i , j ) ;
              result.set ( i , j , mij ) ;
              result.set ( j , i , mij ) ;
            }
        }
      return result ;
    }
    // ========================================================================
    /// template creator of vector from sequence of data 
    template <class ITERATOR>
      inline Vector vector
      ( ITERATOR begin ,
        ITERATOR end   )
        {
          const unsigned int N = std::distance ( begin , end ) ;
          Vector result { N } ;
          for ( unsigned int index = 0 ; begin != end ; ++begin, ++index )
            { result.set ( index , *begin ) ; }
          return result ;
        }
    // =========================================================================
    /// S-vector -> GSL Vector 
    template <class T, unsigned int D>
      inline Vector vector ( const ROOT::Math::SVector<T,D>& v )
    { return vector ( v.begin() , v.end() ) ; }
    // ========================================================================    
    /// convert T-matrix into GSL matrix 
    Matrix matrix ( const TMatrixT<float>&          m ) ;    
    /// convert T-matrix into GSL matrix 
    Matrix matrix ( const TMatrixT<double>&         m ) ;
    /// convert T-matrix into GSL matrix 
    Matrix matrix ( const TMatrixT<long double>&    m ) ;    
    /// convert T-matrix into GSL matrix 
    Matrix matrix ( const TMatrixTSym<float>&       m ) ;
    /// convert T-matrix into GSL matrix 
    Matrix matrix ( const TMatrixTSym<double>&      m ) ;
    /// convert T-matrix into GSL matrix 
    Matrix matrix ( const TMatrixTSym<long double>& m ) ;    
    // convert T-vector into GSL vector 
    Vector vector ( const TVectorT<float>&          v ) ;
    // convert T-vector into GSL vector 
    Vector vector ( const TVectorT<double>&         v ) ;
    // convert T-vector into GSL vector 
    Vector vector ( const TVectorT<long double>&    v ) ;    
    // ========================================================================
  } //                                          The end of namesapce Ostap::GSL
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_GSL_LINALGUTILS_H
// ============================================================================
//                                                                      The END 
// ============================================================================
