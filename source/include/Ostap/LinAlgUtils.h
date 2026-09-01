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
      // ====================================================================

      // ====================================================================
      // Symmetric S-matrix -> GSL matrix
      // ====================================================================
      
      // =====================================================================
      /// create GSL matrix from symmetric SMatrix     
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
      // =====================================================================
      
      // =====================================================================
      /// convert T-matrix into GSL matrix 
      Matrix matrix ( const TMatrixT<float>&          m ) ;    
      /// convert T-matrix into GSL matrix 
      Matrix matrix ( const TMatrixT<double>&         m ) ;   
      /// convert T-matrix into GSL matrix 
      Matrix matrix ( const TMatrixTSym<float>&       m ) ;
      /// convert T-matrix into GSL matrix 
      Matrix matrix ( const TMatrixTSym<double>&      m ) ;
      // =====================================================================
      
      // =====================================================================
      /// template creator of GSL-vector from sequence of data 
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename category   = typename std::iterator_traits<ITERATOR>::iterator_category,
                typename std::enable_if< std::is_convertible<value_type, double>::value &&
                                         std::is_base_of<std::forward_iterator_tag, category>::value,int>::type = 0>      
      inline Vector vector
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        const auto dist     = std::distance ( begin , end ) ;
        const std::size_t N = 0 <= dist ? static_cast<std::size_t>( dist ) : 0 ;
        Vector result { N } ;
        for ( unsigned int index = 0 ; begin != end ; ++begin, ++index )
        { result.set ( index , *begin ) ; }
        return result ;
      }
      // =====================================================================    
      /// S-vector -> GSL Vector 
      template <class T, unsigned int D>
      inline Vector vector ( const ROOT::Math::SVector<T,D>& v )
      { return vector ( v.begin () , v.end () ) ; }
      // =====================================================================    
      /// T-vector -> GSL Vector 
      template <class T>
      inline Vector vector ( const TVectorT<T>& v )
      { return vector ( v.begin () , v.end () ) ; }
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
