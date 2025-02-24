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
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::GSL::Matrix
Ostap::GSL::matrix
( const TMatrixT<float>&          m )
{
  Ostap::Assert ( m.IsValid() && m.GetNrows() && m.GetNcols () , 
                  "T-matrix is not valid!" , 
                  "Ostap::GSL::matrix"     ,
                  INVALID_TMATRIX          , __FILE__ , __LINE__ ) ;
  
  const unsigned int nr = m.GetNrows () ;
  const unsigned int nc = m.GetNcols () ;
  //
  Matrix result { nr , nc }  ;
  for ( unsigned long i = 0 ; i < nr ; ++i )
    { for ( unsigned long j = 0 ; j < nc ; ++j )
        { result.set ( i , j , m ( i , j ) ) ; } }
  //
  return result ;
}
// ============================================================================
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::GSL::Matrix
Ostap::GSL::matrix
( const TMatrixT<double>&          m )
{
  Ostap::Assert ( m.IsValid() && m.GetNrows() && m.GetNcols () , 
                  "T-matrix is not valid!" , 
                  "Ostap::GSL::matrix"     ,
                  INVALID_TMATRIX          , __FILE__ , __LINE__ ) ;
  
  const unsigned int nr = m.GetNrows () ;
  const unsigned int nc = m.GetNcols () ;
  //
  Matrix result { nr , nc }  ;
  for ( unsigned long i = 0 ; i < nr ; ++i )
    { for ( unsigned long j = 0 ; j < nc ; ++j )
        { result.set ( i , j , m ( i , j ) ) ; } }
  //
  return result ;
}
// ============================================================================
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::GSL::Matrix
Ostap::GSL::matrix
( const TMatrixTSym<float>&          m )
{
  Ostap::Assert ( m.IsValid() && m.GetNrows() && m.GetNrows() == m.GetNcols () , 
                  "T-matrix is not valid!" , 
                  "Ostap::GSL::matrix"     ,
                  INVALID_TMATRIX          , __FILE__ , __LINE__ ) ;
  
  const unsigned int N  = m.GetNrows () ;
  //
  Matrix result { N }  ;
  for ( unsigned long i = 0 ; i < N  ; ++i )
    {
      result.set ( i , i , m ( i , i ) ) ;      
      for ( unsigned long j = i + 1  ; j < N ; ++j )
        {
          const double mij = m ( i , j ) ;
          result.set ( i , j , mij ) ;
          result.set ( j , i , mij ) ;
        }
    }
  //
  return result ;
}
// ============================================================================
// convert T-matrix into GSL matrix 
// ============================================================================
Ostap::GSL::Matrix
Ostap::GSL::matrix
( const TMatrixTSym<double>&         m )
{
  Ostap::Assert ( m.IsValid() && m.GetNrows() && m.GetNrows() == m.GetNcols () , 
                  "T-matrix is not valid!" , 
                  "Ostap::GSL::matrix"     ,
                  INVALID_TMATRIX          , __FILE__ , __LINE__ ) ;
  
  const unsigned int N  = m.GetNrows () ;
  //
  Matrix result { N }  ;
  for ( unsigned long i = 0 ; i < N  ; ++i )
    {
      result.set ( i , i , m ( i , i ) ) ;      
      for ( unsigned long j = i + 1  ; j < N ; ++j )
        {
          const double mij = m ( i , j ) ;
          result.set ( i , j , mij ) ;
          result.set ( j , i , mij ) ;
        }
    }
  //
  return result ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
