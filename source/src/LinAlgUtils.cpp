// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/LinAlg.h"
#include "Ostap/LinAlgUtils.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_linalg.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/LinAlg.h"
#include "Ostap/Buffer.h"
// ============================================================================
// local
// ============================================================================
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for helper  GSL classes 
 *  @date 2020-09-22 
 *  @author Vanya BELYAEV IvanBelyaev@iter.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================
  // convert TMatrixT into GSL matrix 
  // ============================================================================
  template <typename FLOAT>
  inline Ostap::Math::GSL::Matrix
  _matrix_ ( const TMatrixT<FLOAT>& m )
  {
    //
    Ostap::Assert ( m.IsValid() && 1 <= m.GetNrows() && 1 <= m.GetNcols () , 
                    "T-matrix is not valid!" , 
                    "Ostap::Math::GSL::matrix" ,
                    INVALID_TMATRIX , __FILE__ , __LINE__ ) ;
    
    const std::size_t NR = m.GetNrows () ;
    const std::size_t NC = m.GetNcols () ;
    //
    return Ostap::Math::GSL::Matrix { NR , NC , Ostap::Utils::Buffer<FLOAT> ( m.GetMatrixArray() , NR * NC ) } ;
  }
  // ==========================================================================

  // ==========================================================================
  /** convert TMatrixTSym into GSL matrix
   *  @attention It reads only the upper triangle part and fill both upper&lower parts
   */
  // ============================================================================
  template <typename FLOAT>
  inline Ostap::Math::GSL::Matrix
  _matrix_ ( const TMatrixTSym<FLOAT>& m )
  {
    Ostap::Assert ( m.IsValid() && 1 <= m.GetNrows() && m.GetNrows() == m.GetNcols () , 
                    "T-matrix is not valid!" , 
                    "Ostap::Math::GSL::matrix" ,
                    INVALID_TMATRIX , __FILE__ , __LINE__ ) ;
    
    const std::size_t N = m.GetNrows () ;
    //
    Ostap::Math::GSL::Matrix result { N } ;
    for ( std::size_t i = 0 ; i < N ; ++i )
    {
      result.set ( i , i , m ( i , i ) ) ;      
      for ( std::size_t j = i + 1 ; j < N ; ++j )
      {
        const double mij = m ( i , j ) ;
        result.set ( i , j , mij ) ;
        result.set ( j , i , mij ) ;
      }
    }
    //
    return result ;
  }
  // ==========================================================================  
}
// ============================================================================
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::matrix
( const TMatrixT<float>&  m ) { return _matrix_ ( m ) ; }
// ============================================================================
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::matrix
( const TMatrixT<double>&  m ) { return _matrix_ ( m ) ; }
// ============================================================================
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::matrix
( const TMatrixTSym<float>& m ) { return _matrix_ ( m ) ; }
// ============================================================================
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::matrix
( const TMatrixTSym<double>& m ) { return _matrix_ ( m ) ; }

// ============================================================================
//                                                                      The END 
// ============================================================================
