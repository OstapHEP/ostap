// ============================================================================
// Include files 
// ============================================================================
// STD/STL
// ============================================================================
#include <cstddef>
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/LinAlg.h"
#include "Ostap/StatusCode.h"
#include "Ostap/Buffer.h"
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_version.h"
// ============================================================================
// local
// ============================================================================
#include "GSL_sentry.h"
#include "format.h"
#include "status_codes.h"
#include "local_math.h"
// ============================================================================
/** @file 
 *  Implementation file for helper  GSL classes 
 *  @date 2020-09-22 
 *  @author Vanya BELYAEV IvanBelyaev@iter.ru
 */
// ============================================================================
// GSL version major 
// ============================================================================
std::size_t Ostap::Math::GSL::GSL_version_major () 
{ return GSL_MAJOR_VERSION ; }
// ============================================================================
// GSL version minor
// ============================================================================
std::size_t Ostap::Math::GSL::GSL_version_minor () 
{ return GSL_MINOR_VERSION ; }
// ============================================================================
// GSL versionmajor  x 1000 + GAL version minor  
// ============================================================================
std::size_t Ostap::Math::GSL::GSL_version_int   () 
{ return 1000 * GSL_MAJOR_VERSION + GSL_MINOR_VERSION ; }
// ============================================================================
// GSL version as string
// ============================================================================
std::string Ostap::Math::GSL::GSL_version () { return gsl_version ; }
// ============================================================================

// ===========================================================================
// Matrix class 
// ===========================================================================

template <class TYPE>
void Ostap::Math::GSL::Matrix::fill_impl
( const TYPE*       data ,
  const std::size_t size )
{
  //
  const std::size_t rows = nRows () ;
  const std::size_t cols = nCols () ;
  //
  Ostap::Assert ( nullptr != data                       , 
                  "Matrix: null buffer pointer"         , 
                  "Ostap::Math::GSL::Matrix::fill_impl" , 
                  INVALID_BUFFER                        , __FILE__ , __LINE__ ) ;
  
  // =========================================================================
  // (1) symmetric representation (lower triangular matrix)
  // =========================================================================

  if ( rows == cols && rows * ( rows + 1 ) / 2 == size )
  {
    std::size_t idx = 0;
    for ( std::size_t i = 0; i < rows ; ++i )
    {
      for ( std::size_t j = 0; j <= i ; ++j )
      {
        const double val = static_cast<double> ( data [ idx++ ] ) ; // Elements go one by one: (0,0), (1,0), (1,1), (2,0)...
        gsl_matrix_set ( m_matrix , i , j , val ) ;
        if ( i != j ) { gsl_matrix_set ( m_matrix , j , i , val ) ; } 
      }
    }
    return ;
  }
  
  // =========================================================================
  // (2) Generic representation
  // =========================================================================
  
  Ostap::Assert ( rows * cols == size                   , 
                  "Matrix: invalid buffer size"         , 
                  "Ostap::Math::GSL::Matrix::fill_impl" , 
                  INVALID_BUFFER                        , __FILE__ , __LINE__  ) ;
  //
  for ( std::size_t i = 0 ; i < rows ; ++i )
  { for ( std::size_t j = 0 ; j < cols ; ++j )
    { gsl_matrix_set ( m_matrix , i , j , static_cast<double>( data[i * cols + j] ) ) ; } }
  //
}
// ============================================================================
template void Ostap::Math::GSL::Matrix::fill_impl<float>       ( const       float* , const std::size_t ) ;
template void Ostap::Math::GSL::Matrix::fill_impl<double>      ( const      double* , const std::size_t ) ;
template void Ostap::Math::GSL::Matrix::fill_impl<long double> ( const long double* , const std::size_t ) ;
// ============================================================================

// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const std::size_t N1 , 
  const std::size_t N2 )
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{
  Ostap::Assert ( m_matrix                         ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::Matrix" ,
                  MATRIX_ALLOCATION_FAILURE  , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const std::size_t N1    , 
  const std::size_t N2    , 
  const double       value )
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{
  Ostap::Assert ( m_matrix                         ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::Matrix"       ,
                  MATRIX_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( std::isfinite ( value )          ,
                  "Cannot use !std::isfinite"      ,
                  "Ostap::Math::GSL::Matrix"       ,
                  INVALID_SCALE                    , __FILE__ , __LINE__ ) ;
  gsl_matrix_set_all ( m_matrix , value ) ;
}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const std::size_t N1    , 
  const std::size_t N2    , 
  const Ostap::Math::GSL::Matrix::Zero /* zero */ ) 
  : m_matrix ( gsl_matrix_calloc ( N1 , N2 ) )  // NB: calloc is here!!!
{
  Ostap::Assert ( m_matrix                         ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::Matrix"       ,
                  MATRIX_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const std::size_t N1    , 
  const std::size_t N2    , 
  const Ostap::Math::GSL::Matrix::Id /* zero */ ) 
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{
  Ostap::Assert ( m_matrix                         ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::Matrix"       ,
                  MATRIX_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;  
  gsl_matrix_set_identity ( m_matrix ) ;
}
// ============================================================================
// allocate square GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const std::size_t N  )
  : Matrix ( N , N )
{}
// ============================================================================
// allocate square GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const std::size_t                    N     , 
  const Ostap::Math::GSL::Matrix::Zero zero  )
  : Matrix ( N , N , zero )
{}
// ============================================================================
// allocate square GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const std::size_t                 N  , 
  const Ostap::Math::GSL::Matrix::Id  id ) 
  : Matrix ( N , N , id )
{}
// ============================================================================
// allocate square permutation GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix
( const Ostap::Math::GSL::Permutation& p ) 
  : m_matrix ( gsl_matrix_calloc ( p.size () , p.size () ) ) 
{
  //
  Ostap::Assert ( m_matrix                         ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::Matrix"       ,
                  MATRIX_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;  
  //
  Ostap::Assert ( p.valid()                        ,
                  "(GSL)Permutation is invalid!"   , 
                  "Ostap::Math::GSL::Matrix"       ,
                  INVALID_PERMUTATION              , __FILE__ , __LINE__ ) ;
  //
  for ( std::size_t j = 0 ; j < nRows() ; ++j )
  {
    const std::size_t k = p.get ( j ) ; 
    set ( j , k , 1 ) ;
  }
}
// ==========================================================================
Ostap::Math::GSL::Matrix::Matrix
( const Ostap::Math::GSL::Vector & v ) 
  : m_matrix ( gsl_matrix_calloc ( v.size () , v.size () ) ) 
{
  //
  Ostap::Assert ( m_matrix                         ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::Matrix"       ,
                  MATRIX_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;  
  //
  const std::size_t N = v.size() ;
  for ( std::size_t i = 0 ; i < N ; ++i )
    { set ( i , i , v.get ( i ) ) ; }
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix  
( const Ostap::Math::GSL::Matrix&  right ) 
  : m_matrix ( gsl_matrix_alloc ( right.m_matrix->size1 , 
                                  right.m_matrix->size2 ) )  
{
  //
  Ostap::Assert ( m_matrix                         ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::Matrix"       ,
                  MATRIX_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
  //
  gsl_matrix_memcpy ( m_matrix , right.m_matrix ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::GSL::Matrix::Matrix  
(       Ostap::Math::GSL::Matrix&&  right ) 
  : m_matrix ( right.m_matrix )  
{
  right.m_matrix = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-matrix 
// ============================================================================
Ostap::Math::GSL::Matrix::~Matrix () 
{
  if ( nullptr != m_matrix ) 
  {
    gsl_matrix_free ( m_matrix ) ; 
    m_matrix = nullptr ; 
  }
}
// ============================================================================
// copy assignement! 
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::operator=
( const Ostap::Math::GSL::Matrix&  right ) 
{
  if ( &right == this ) { return *this ; }
  resize ( right.nRows() , right.nCols() ) ;
  gsl_matrix_memcpy ( m_matrix , right.m_matrix ) ;
  return *this ;
}
// ============================================================================
// move assignement! 
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::operator=
( Ostap::Math::GSL::Matrix&& right ) 
{
  if ( &right == this ) { return *this ; }
  std::swap ( m_matrix , right.m_matrix ) ;
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::resize
( const std::size_t n1 ,
  const std::size_t n2 )
{
  if ( n1 != m_matrix->size1 || n2 != m_matrix->size2 )
  {
    gsl_matrix_free ( m_matrix ) ; 
    m_matrix = gsl_matrix_alloc ( n1 , n2 ) ;
    //
    Ostap::Assert ( m_matrix                           ,
                    "(GSL)Matrix allocation failure"   ,
                    "Ostap::Math::GSL::Matrix::resize" ,
                    MATRIX_ALLOCATION_FAILURE          , __FILE__ , __LINE__ ) ;
  }
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::resize
( const std::size_t n1    ,
  const std::size_t n2    ,
  const double      value )
{
  if ( !value ) { return resize ( n1 , n2 , Zero() ) ; }
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannot use !std::isfinite"         ,
                  "Ostap::Math::GSL::Matrix::resize"  ,
                  INVALID_SCALE                       , __FILE__ , __LINE__ ) ;
  resize ( n1 , n2 ) ;
  gsl_matrix_set_all ( m_matrix , value ) ;
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::resize
( const std::size_t n1    ,
  const std::size_t n2    ,
  const Ostap::Math::GSL::Matrix::Zero /* zero */ ) 
{
  if ( n1 != m_matrix->size1 || n2 != m_matrix->size2 )
    {
      gsl_matrix_free ( m_matrix ) ; 
      m_matrix = gsl_matrix_calloc ( n1 , n2 ) ;      // CALLOC!
      Ostap::Assert ( m_matrix                           ,
                      "(GSL)Matrix allocation failure"   ,
                      "Ostap::Math::GSL::Matrix::resize" ,
                      MATRIX_ALLOCATION_FAILURE          , __FILE__ , __LINE__ ) ;      
    }
  else
    { gsl_matrix_set_all ( m_matrix , 0 ) ; }
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::resize
( const std::size_t n1    ,
  const std::size_t n2    ,
  const Ostap::Math::GSL::Matrix::Id /* zero */ ) 
{
  resize ( n1 , n2 ) ;
  gsl_matrix_set_identity ( m_matrix ) ;
  return *this ;
}
// ============================================================================
// get matrix row by value)
// ============================================================================
Ostap::Math::GSL::Vector
Ostap::Math::GSL::Matrix::row
( const std::size_t n ) const  
{
  Ostap::Assert ( n < nRows ()                    ,
                  "(GSL)Invalid row index"        ,
                  "Ostap::Math::GSL::Matrix::row" ,
                  INVALID_INDEX                   , __FILE__ , __LINE__ ) ;
  //
  Vector r { nCols() } ; 
  gsl_matrix_get_row ( r.vector() , m_matrix , n ) ;
  return r ;
}
// ============================================================================
// get matrix column by value)
// ============================================================================
Ostap::Math::GSL::Vector
Ostap::Math::GSL::Matrix::column
( const std::size_t n ) const  
{
  Ostap::Assert ( n < nCols ()                       ,
                  "(GSL)Invalid column index"        ,
                  "Ostap::Math::GSL::Matrix::column" ,
                  INVALID_INDEX                      , __FILE__ , __LINE__ ) ;
  //
  Vector c { nRows() } ; 
  gsl_matrix_get_col ( c.vector () , m_matrix , n ) ;
  return c ;
}
// ============================================================================
// Numerical equality of two GSL matrices 
// ============================================================================
bool Ostap::Math::GSL::Matrix::equal
( const Ostap::Math::GSL::Matrix& r ) const 
{
  if ( &r == this       ) { return true  ; }
  //
  const std::size_t NR { nRows () } ;
  if ( NR != r.nRows () ) { return false ; }
  //
  const std::size_t NC { nCols () } ;
  if ( NC != r.nCols () ) { return false ; }
  //
  for ( std::size_t i = 0 ; i < NR ; ++i )
  { for ( std::size_t j = 0 ; j < NC ; ++j )
    {
      const double e1 =   get ( i , j ) ;
      if ( !std::isfinite ( e1 ) ) { return false ; }
      const double e2 = r.get ( i , j ) ;
      if ( !std::isfinite ( e2 ) ) { return false ; }        
      if ( !s_equal ( e1 , e2 )  ) { return false ; }
    }
  }
  //
  return true ;
}
// ============================================================================
// swap two matrices
// ============================================================================
void Ostap::Math::GSL::Matrix::swap
( Ostap::Math::GSL::Matrix& right )
{ std::swap ( m_matrix , right.m_matrix ) ; }
// ============================================================================
// swap two rows in matrix 
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::swap_rows
( const std::size_t i1 ,
  const std::size_t i2 )
{
  if ( i1 == i2 ) { return *this ; }
  Ostap::Assert ( i1 < nRows () && i2 <= nRows ()       , 
                  "Invalid row index!"                  ,
                  "Ostap::Math::GSL::Martix::swap_rows" ,
                  INVALID_ROWINDEX                      , __FILE__ , __LINE__ ) ;
  gsl_matrix_swap_rows    ( m_matrix , i1 , i2 ) ;
  return *this ; 
}
// ============================================================================
// swap two columns in matrix 
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::swap_cols
( const std::size_t i1 ,
  const std::size_t i2 )
{
  if ( i1 == i2 ) { return *this ; }
  Ostap::Assert ( i1 < nCols () && i2 <= nCols ()       , 
                  "Invalid column index!"               ,
                  "Ostap::Math::GSL::Martix::swap_cols" ,
                  INVALID_COLINDEX                      , __FILE__ , __LINE__ ) ;
  gsl_matrix_swap_columns ( m_matrix , i1 , i2 ) ;
  return *this ; 
}
// ============================================================================
// permute the rows     of the ematrix according to permutation 
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::permute_rows
( const Ostap::Math::GSL::Permutation& p )
{
  Ostap::Assert ( nRows () == p.size ()                    ,
                  "Inconsistent Permutation structure"     , 
                  "Ostap::Math::GSL::Matrix::permute_rows" ,
                  INVALID_PERMUTATION                      , __FILE__ , __LINE__ ) ;
  //
  for ( std::size_t i = 0 ; i < nRows() ; ++i )
  {
    const std::size_t k = p ( i ) ;
    if ( i != k ) { gsl_matrix_swap_rows ( m_matrix , i , k ) ; }      
  }
  return *this ;
}
// ============================================================================
// permute the columns of the ematrix according to permutation 
// ============================================================================
Ostap::Math::GSL::Matrix&
Ostap::Math::GSL::Matrix::permute_cols
( const Ostap::Math::GSL::Permutation& p )
{
  Ostap::Assert ( nCols () == p.size ()                    ,
                  "Inconsistent Permutation structure"     , 
                  "Ostap::Math::GSL::Matrix::permute_cols" ,
                  INVALID_PERMUTATION                      , __FILE__ , __LINE__ ) ;
  //
  for ( std::size_t i = 0 ; i < nCols() ; ++i )
  {
    const std::size_t k = p ( i ) ;
    if ( i != k ) { gsl_matrix_swap_columns ( m_matrix , i , k ) ; }      
  }
  return *this ;
}
// ============================================================================
// Are all elements numerically equal to zero?      
// ============================================================================
bool Ostap::Math::GSL::Matrix::iszero   () const
{
  for ( std::size_t i = 0 ; i < nRows() ; ++i )
  { for ( std::size_t j = 0 ; j < nCols() ; ++j )
    { if ( !s_zero ( get ( i ,j ) ) ) { return false ; } } }
  return true ;
}
// ============================================================================
// Are all elements finite ? 
// ============================================================================
bool Ostap::Math::GSL::Matrix::isfinite () const
{
  for ( std::size_t i = 0 ; i < nRows() ; ++i )
  { for ( std::size_t j = 0 ; j < nCols() ; ++j )
    { if ( !std::isfinite ( get ( i ,j ) ) ) { return false ; } } }
  return true ;  
}
// ============================================================================
// scale matrix 
// ============================================================================
Ostap::Math::GSL::Matrix&                             
Ostap::Math::GSL::Matrix::imul
( const double value )
{
  if ( 1 == value ) { return *this ; }
  Ostap::Assert ( std::isfinite ( value )           ,
                  "Cannot add !std::isfinite"       ,
                  "Ostap::Math::GSL::Matrix::iadd"  ,
                  INVALID_SCALE                     , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( std::isfinite ( value )           ,
                  "Cannto scale but !std::isfinite" ,
                  "Ostap::Math::GSL::Matrix::imul"  ,
                  INVALID_SCALE                     , __FILE__ , __LINE__ ) ;
  //
  gsl_matrix_scale ( m_matrix , value ) ;
  return *this;
}
// ============================================================================
// add&subtract matrix 
// ============================================================================
Ostap::Math::GSL::Matrix&                             
Ostap::Math::GSL::Matrix::iadd
( const Ostap::Math::GSL::Matrix& right )
{
  Ostap::Assert ( this->nRows () == right.nRows() &&
                  this->nCols () == right.nCols()                ,
                  "Cannot add matrix of incompatible structure"  ,
                  "Ostap::Math::GSL::Martix::iadd"               ,
                  INVALID_GMATRIX                                , __FILE__ , __LINE__ ) ;
  //
  gsl_matrix_add ( m_matrix , right.m_matrix ) ;
  return *this;
}
// ============================================================================
// add&subtract matrix 
// ============================================================================
Ostap::Math::GSL::Matrix&                             
Ostap::Math::GSL::Matrix::isub
( const Ostap::Math::GSL::Matrix& right )
{
  Ostap::Assert ( this->nRows () == right.nRows() &&
                  this->nCols () == right.nCols()               ,
                  "Cannot sub matrix of incompatible structure" ,
                  "Ostap::Math::GSL::Martix::iasub"             ,
                  INVALID_GMATRIX                               , __FILE__ , __LINE__ ) ;
  //
  gsl_matrix_sub ( m_matrix , right.m_matrix ) ;
  return *this;
}
// ============================================================================
// add identity matrix 
// ============================================================================
Ostap::Math::GSL::Matrix&                             
Ostap::Math::GSL::Matrix::iadd
( const double value )
{
  if ( !value ) { return *this ; }
  Ostap::Assert ( std::isfinite ( value )          ,
                  "Cannot add  !std::isfinite"     ,
                  "Ostap::Math::GSL::Matrix::iadd" ,
                  INVALID_SCALE                    , __FILE__ , __LINE__ ) ;
  //
  const std::size_t N = std::min ( nRows () , nCols() );
  for ( std::size_t i = 0 ; i < N ; ++i )
  { set ( i , i , get ( i , i ) + value ) ; }
  return *this;
}
// ============================================================================
// multiply matrices  using CBLAS dgemm function 
// ============================================================================
Ostap::Math::GSL::Matrix                            
Ostap::Math::GSL::Matrix::multiply
( const Ostap::Math::GSL::Matrix& right ) const
{
  Ostap::Assert ( this->nCols () == right.nRows()                      ,
                  "Cannot multiply matrices of incompatible structure" ,
                  "Ostap::Math::GSL::Martix::multiply"                 ,
                  INVALID_GMATRIX                                      , __FILE__ , __LINE__ ) ;
  
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  Matrix result { nRows() , right.nCols() } ;
  // 
  const int status = gsl_blas_dgemm ( CblasNoTrans    ,
                                      CblasNoTrans    ,
                                      1.0             ,
                                      this ->matrix() ,
                                      right .matrix() ,
                                      0.0             ,
                                      result.matrix() ) ;
  //
  Ostap::Assert ( !status ,
                  "Error from gsl_blas_dgemm function" ,
                  "Ostap::Math::GSL::Matrix::multiply" , 
                  ERROR_GSL + status                   , __FILE__ , __LINE__ ) ; 
  //
  return result ;
}
// ============================================================================
// multiply matrices  using CBLAS dgemm function 
// ============================================================================
Ostap::Math::GSL::Matrix&                           
Ostap::Math::GSL::Matrix::imul
( const Ostap::Math::GSL::Matrix& right )
{
  Ostap::Assert ( this->nCols () == right.nRows () ,
                  "Cannot multiply matrices of incompatible structure" ,
                  "Ostap::Math::GSL::Martix::multiply"  ,
                  INVALID_GMATRIX , __FILE__ , __LINE__ ) ;
  
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  Matrix result { nRows() , right.nCols() } ;
  //
  const int status = gsl_blas_dgemm ( CblasNoTrans    ,
                                      CblasNoTrans    ,
                                      1.0             ,
                                      this ->matrix() ,
                                      right .matrix() ,
                                      0.0             ,
                                      result.matrix() ) ;
  //
  Ostap::Assert ( !status ,
                  "Error from gsl_blas_dgem function"  ,
                  "Ostap::Math::GSL::Matrix::imul"     , 
                  ERROR_GSL + status                   , __FILE__ , __LINE__ ) ; 
  //
  this->swap ( result ) ;
  //
  return *this ;
}
// ============================================================================
// multiply matrix and vector using CBLAS dgemv function 
// ============================================================================
Ostap::Math::GSL::Vector
Ostap::Math::GSL::Matrix::multiply
( const Ostap::Math::GSL::Vector& right ) const
{
  Ostap::Assert ( this->nCols () == right.size()                            ,
                  "Cannot multiply matrix&vector of incompatible structure" ,
                  "Ostap::Math::GSL::Martix::multiply"                      ,
                  INVALID_GMATRIX                                           , __FILE__ , __LINE__ ) ;
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //  
  Vector result { nRows() } ;
  //
  const int status = gsl_blas_dgemv ( CblasNoTrans     ,
                                      1.0              ,
                                      this ->matrix () ,
                                      right .vector () ,
                                      0.0              ,
                                      result.vector () ) ;
  Ostap::Assert ( !status ,
                  "Error from gsl_blas_dgemv function" ,
                  "Ostap::Math::GSL::Matrix::multiply" , 
                  ERROR_GSL + status                   , __FILE__ , __LINE__ ) ; 
  //
  return result ;
}
// ============================================================================
// multiply matrix abd vector using CBLAS dgemv function 
// ============================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::Matrix::multiply
( const Ostap::Math::GSL::Permutation& right ) const
{
  Ostap::Assert ( nCols() == right.size()                      ,
                  "Mismatch for permutation/matrix structure!" ,
                  "Ostap::GLS::Matrix::multiply"               , 
                  INVALID_PERMUTATION  , __FILE__ , __LINE__   ) ;
  return (*this) * Matrix ( right ) .T() ;
}
// ===========================================================================
// transpose the matrix
// ===========================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::Matrix::T () const
{
  // prepare result 
  Matrix result { nCols () , nRows () } ;
  gsl_matrix_transpose_memcpy ( result.matrix () , matrix () ) ;
  return result; 
}
// ============================================================================



// ============================================================================
/*  Matrix multiplication:
 *  \f$ C = a^{(T_a)} \times b^{(T_b)}\f$
 *  @param a (INPUT) the first input matrix
 *  @param Ta (INPUT) transpose the first matrix?
 *  @param b (INPUT) the second input matrix
 *  @param Tb (INPUT) transpose the second matrix?
 *  @param c (OUTPUT) the output matrix
 *  @return status code      
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::MM
( const Ostap::Math::GSL::Matrix& a  ,
  const bool                      Ta ,
  const Ostap::Math::GSL::Matrix& b  ,
  const bool                      Tb , 
  Ostap::Math::GSL::Matrix&       c  ) 
{
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //

  // Determine transpose operations based on flags
  const CBLAS_TRANSPOSE_t transA = Ta ? CblasTrans : CblasNoTrans;
  const CBLAS_TRANSPOSE_t transB = Tb ? CblasTrans : CblasNoTrans;

  //
  const std::size_t rows_A = Ta ? a.nCols() : a.nRows() ;
  const std::size_t cols_A = Ta ? a.nRows() : a.nCols() ;
  const std::size_t rows_B = Tb ? b.nCols() : b.nRows() ;
  const std::size_t cols_B = Tb ? b.nRows() : b.nCols() ;
  //
  if ( cols_A != rows_B )
  {
    const int status = GSL_EBADLEN ; 
    gsl_error ( "Ostap::Math::GSL::MM: Incompatible matrix dimensions for multiplication" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ; 
  }

  // resize output if the size is wrong
  if ( ( &a != &c ) &&  ( &b != &c ) && ( c.nRows() != rows_A || c.nCols() != cols_B ) )
  { c = Matrix ( rows_A , cols_B ) ; }
  
  // aliasing ? 
  if ( &a == &c || &b == &c )
  {
    // use temporary matrix to avoid aliasing issues
    Matrix tmp { rows_A , cols_B } ;
    Ostap::StatusCode sc = MM ( a , Ta , b , Tb , tmp ) ;
    if ( sc.isSuccess() ) { c.swap ( tmp ) ; } 
    return sc ; 
  }

  // OPTIMIZATION: Symmetric product (A^T * A or A * A^T) when A == B and T1 != T2
  // Uses dsyrk which only computes the upper triangle, doubling the performance.
  if ( &a == &b && Ta != Tb ) 
  {
    // T1 = 1, T2 = 0 -> A^T * A (trans = CblasTrans)
    // T1 = 0, T2 = 1 -> A * A^T (trans = CblasNoTrans)
    const CBLAS_TRANSPOSE_t trans = Ta ? CblasTrans : CblasNoTrans;

    // 1. Compute upper triangle of C: C = 1.0 * op(A) * op(A)^T + 0.0 * C
    const int status = gsl_blas_dsyrk ( CblasUpper , trans , 1.0 , a.matrix () , 0.0 , c.matrix () ) ;
    if ( status ) 
    { 
      gsl_error ( "Ostap::Math::GSL::MM:   Error from gsl_blas_dsyrk function" , __FILE__ , __LINE__ , status ) ;
      return ERROR_GSL + status ;
    }

    // 2. Mirror upper triangle to lower triangle to complete the symmetric matrix
    const std::size_t n = c.nRows  () ;
    gsl_matrix*       m = c.matrix ();
    for ( std::size_t i = 0; i < n; i++) 
    {
      for ( std::size_t j = i + 1 ; j < n; j++ ) 
      {
        const double val = gsl_matrix_get ( m , i , j ) ;
        gsl_matrix_set ( m , j , i , val ) ;
      }
    }
    //
    return Ostap::StatusCode::SUCCESS;
  }
  // END OF OPTIMIZATION  

  //
  // general case: C = op(A) * o
  // Standard BLAS Level 3 general matrix multiplication: C = 1.0 * op(A) * op(B) + 0.0 * C
  const int status = gsl_blas_dgemm ( transA , transB , 1.0 , a.matrix () , b.matrix () , 0.0 , c.matrix () ) ;
  if ( status ) 
  { 
    gsl_error ( "Ostap::Math::GSL::MM:   Error from gsl_blas_dgemm function" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  return Ostap::StatusCode::SUCCESS;
}

// ============================================================================
/*  Matrix multiplication of 3 matrices:
 *  \f$ D = a^{(T_a)} \times b^{(T_b)} \times c^{(T_c)}\f$
 *  @param a  (INPUT) the first input matrix
 *  @param Ta (INPUT) transpose the first matrix?
 *  @param b  (INPUT) the second input matrix
 *  @param Tb (INPUT) transpose the second matrix?
 *  @param c  (INPUT) the third input matrix
 *  @param Tc (INPUT) transpose the third matrix?
 *  @param d  (OUTPUT) the output matrix
 *  @return status code      
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::MMM
( const Ostap::Math::GSL::Matrix& a  ,
  const bool                      Ta ,
  const Ostap::Math::GSL::Matrix& b  ,
  const bool                      Tb , 
  const Ostap::Math::GSL::Matrix& c  ,
  const bool                      Tc , 
  Ostap::Math::GSL::Matrix&       d  ) 
{
  // use GSL error handler sentry
  Ostap::Math::GSL::GSL_Error_Handler sentry ;

  // Calculate dimensions for transposed/non-transposed representations
  const std::size_t m      = Ta ? a.nCols() : a.nRows() ; // Rows of op(A)
  const std::size_t k      = Ta ? a.nRows() : a.nCols() ; // Cols of op(A)
  const std::size_t rows_B = Tb ? b.nCols() : b.nRows() ;
  const std::size_t p      = Tb ? b.nRows() : b.nCols() ; // Cols of op(B)
  const std::size_t rows_C = Tc ? c.nCols() : c.nRows() ;
  const std::size_t n      = Tc ? c.nRows() : c.nCols() ; // Cols of op(C)

  // Validate internal dimensions consistency
  if ( k != rows_B || p != rows_C )
  {
    const int status = GSL_EBADLEN ; 
    gsl_error ( "Ostap::Math::GSL::MMM: Incompatible matrix dimensions for 3-matrix multiplication" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ; 
  }

  // Handle Aliasing (&a == &d || &b == &d || &c == &d):
  // Compute into a temporary matrix first to prevent corrupting inputs
  if ( &a == &d || &b == &d || &c == &d )
  {
    Matrix tmp { m , n } ;
    Ostap::StatusCode sc = MMM ( a , Ta , b , Tb , c , Tc , tmp ) ;
    if ( sc.isSuccess() ) { d.swap ( tmp ) ; } 
    return sc ; 
  }

  // Ensure matrix d has the correct size (m x n)
  if ( d.nRows() != m || d.nCols() != n ) { d = Matrix ( m , n ) ;  }

  // Cost analysis to choose optimal multiplication chain order:
  // Cost 1: (A * B) * C = 2 * m * p * (k + n)
  // Cost 2: A * (B * C) = 2 * k * n * (p + m)
  const std::size_t cost_AB_C = m * p * ( k + n ) ;
  const std::size_t cost_A_BC = k * n * ( p + m ) ;

  if ( cost_AB_C <= cost_A_BC )
  {
    // Strategy 1: Compute AB_tmp = op(A) * op(B), then d = AB_tmp * op(C)
    Matrix AB_tmp { m , p } ;
    Ostap::StatusCode sc = MM ( a , Ta , b , Tb , AB_tmp ) ;
    if ( sc.isFailure () ) { return sc ; }
    //
    return MM ( AB_tmp , false , c , Tc , d ) ;
  }
  else
  {
    // Strategy 2: Compute BC_tmp = op(B) * op(C), then d = op(A) * BC_tmp
    Matrix BC_tmp { k , n } ;
    Ostap::StatusCode sc = MM ( b , Tb , c , Tc , BC_tmp ) ;
    if ( sc.isFailure () ) { return sc ; }
    //
    return MM ( a , Ta , BC_tmp , false , d ) ;
  }
}

// ============================================================================
/*  Matrix multiplication of 4 matrices:
 *  \f$ E = a^{(T_a)} \times b^{(T_b)} \times c^{(T_c)} \times d^{(T_d)}\f$
 *  @param a  (INPUT) the first input matrix
 *  @param Ta (INPUT) transpose the first matrix?
 *  @param b  (INPUT) the second input matrix
 *  @param Tb (INPUT) transpose the second matrix?
 *  @param c  (INPUT) the third input matrix
 *  @param Tc (INPUT) transpose the third matrix?
 *  @param d  (INPUT) the fourth input matrix
 *  @param Td (INPUT) transpose the fourth matrix?
 *  @param e  (OUTPUT) the output matrix
 *  @return status code      
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::MMMM
( const Ostap::Math::GSL::Matrix& a  ,
  const bool                      Ta ,
  const Ostap::Math::GSL::Matrix& b  ,
  const bool                      Tb , 
  const Ostap::Math::GSL::Matrix& c  ,
  const bool                      Tc , 
  const Ostap::Math::GSL::Matrix& d  ,
  const bool                      Td , 
  Ostap::Math::GSL::Matrix&       e  ) 
{
  // use GSL error handler sentry
  Ostap::Math::GSL::GSL_Error_Handler sentry ;

  // Dimensions of the 4 matrices (rows x cols) after accounting for transposes
  // op(A): m x k
  // op(B): k x p
  // op(C): p x q
  // op(D): q x n
  const std::size_t m = Ta ? a.nCols() : a.nRows() ; // Rows of op(A)
  const std::size_t k = Ta ? a.nRows() : a.nCols() ; // Cols of op(A)
  const std::size_t rows_B = Tb ? b.nCols() : b.nRows() ;
  const std::size_t p      = Tb ? b.nRows() : b.nCols() ; // Cols of op(B)
  const std::size_t rows_C = Tc ? c.nCols() : c.nRows() ;
  const std::size_t q      = Tc ? c.nRows() : c.nCols() ; // Cols of op(C)
  const std::size_t rows_D = Td ? d.nCols() : d.nRows() ;
  const std::size_t n      = Td ? d.nRows() : d.nCols() ; // Cols of op(D)

  // Validate internal matrix dimension compatibility
  if ( k != rows_B || p != rows_C || q != rows_D )
  {
    const int status = GSL_EBADLEN ; 
    gsl_error ( "Ostap::Math::GSL::MMMM: Incompatible matrix dimensions for 4-matrix multiplication" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ; 
  }

  // Handle Aliasing (&a == &e || &b == &e || &c == &e || &d == &e):
  // Compute into a temporary matrix first to prevent overwriting inputs
  if ( &a == &e || &b == &e || &c == &e || &d == &e )
  {
    Matrix tmp { m , n } ;
    Ostap::StatusCode sc = MMMM ( a , Ta , b , Tb , c , Tc , d , Td , tmp ) ;
    if ( sc.isSuccess() ) { e.swap ( tmp ) ; } 
    return sc ; 
  }

  // Ensure matrix e has correct dimensions (m x n)
  if ( e.nRows() != m || e.nCols() != n ) { e = Matrix ( m , n ) ;  }

  // Cost calculation for all 5 possible evaluation trees (ignoring multiplier 2 for simplicity):
  // 1. ((AB)C)D : AB(m*k*p) + (AB)C(m*p*q) + ((AB)C)D(m*q*n)
  const std::size_t cost1 = m*k*p + m*p*q + m*q*n ;
  
  // 2. (A(BC))D : BC(k*p*q) + A(BC)(m*k*q) + (A(BC))D(m*q*n)
  const std::size_t cost2 = k*p*q + m*k*q + m*q*n ;
  
  // 3. (AB)(CD) : AB(m*k*p) + CD(p*q*n) + (AB)(CD)(m*p*n)
  const std::size_t cost3 = m*k*p + p*q*n + m*p*n ;
  
  // 4. A((BC)D) : BC(k*p*q) + (BC)D(k*q*n) + A((BC)D)(m*k*n)
  const std::size_t cost4 = k*p*q + k*q*n + m*k*n ;
  
  // 5. A(B(CD)) : CD(p*q*n) + B(CD)(k*p*n) + A(B(CD))(m*k*n)
  const std::size_t cost5 = p*q*n + k*p*n + m*k*n ;

  // Find minimum cost among the 5 strategies
  std::size_t min_cost = cost1 ;
  int strategy = 1 ;

  if ( cost2 < min_cost ) { min_cost = cost2 ; strategy = 2 ; }
  if ( cost3 < min_cost ) { min_cost = cost3 ; strategy = 3 ; }
  if ( cost4 < min_cost ) { min_cost = cost4 ; strategy = 4 ; }
  if ( cost5 < min_cost ) { min_cost = cost5 ; strategy = 5 ; }

  // Execute optimal strategy utilizing our robust MMM or MM subroutines
  switch ( strategy )
  {
    case 1: 
    {
      // ((A * B) * C) * D  ==>  MMM(A, B, C) * D
      Matrix ABC_tmp { m , q } ;
      Ostap::StatusCode sc = MMM ( a , Ta , b , Tb , c , Tc , ABC_tmp ) ;
      if ( sc.isFailure () ) { return sc ; }
      return MM ( ABC_tmp , false , d , Td , e ) ;
    }
    case 2: 
    {
      // (A * (B * C)) * D  ==>  MMM(A, B, C) * D
      Matrix ABC_tmp { m , q } ;
      Ostap::StatusCode sc = MMM ( a , Ta , b , Tb , c , Tc , ABC_tmp ) ;
      if ( sc.isFailure () ) { return sc ; }
      return MM ( ABC_tmp , false , d , Td , e ) ;
    }
    case 3: 
    {
      // (A * B) * (C * D)  ==>  MM(A, B) * MM(C, D)
      Matrix AB_tmp { m , p } ;
      Matrix CD_tmp { p , n } ;
      Ostap::StatusCode sc1 = MM ( a , Ta , b , Tb , AB_tmp ) ;
      if ( sc1.isFailure () ) { return sc1 ; }
      Ostap::StatusCode sc2 = MM ( c , Tc , d , Td , CD_tmp ) ;
      if ( sc2.isFailure () ) { return sc2 ; }
      return MM ( AB_tmp , false , CD_tmp , false , e ) ;
    }
    case 4: 
    {
      // A * ((B * C) * D)  ==>  A * MMM(B, C, D)
      Matrix BCD_tmp { k , n } ;
      Ostap::StatusCode sc = MMM ( b , Tb , c , Tc , d , Td , BCD_tmp ) ;
      if ( sc.isFailure () ) { return sc ; }
      return MM ( a , Ta , BCD_tmp , false , e ) ;
    }
    case 5: 
    {
      // A * (B * (C * D))  ==>  A * MMM(B, C, D)
      Matrix BCD_tmp { k , n } ;
      Ostap::StatusCode sc = MMM ( b , Tb , c , Tc , d , Td , BCD_tmp ) ;
      if ( sc.isFailure () ) { return sc ; }
      return MM ( a , Ta , BCD_tmp , false , e ) ;
    }
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}

// ============================================================================
/*  Matrix multiplication with a diagonal matrix in the middle:
 *  \f$ C = a^{(T_a)} \times \text{diag}(d) \times b^{(T_b)}\f$
 *
 *  @param a  (INPUT) the first input matrix
 *  @param Ta (INPUT) transpose the first matrix?
 *  @param d  (INPUT) diagonal elements of the middle matrix
 *  @param b  (INPUT) the second input matrix
 *  @param Tb (INPUT) transpose the second matrix?
 *  @param c  (OUTPUT) the output matrix
 *  @return status code      
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::MDM
( const Ostap::Math::GSL::Matrix& a  ,
  const bool                      Ta ,
  const Ostap::Math::GSL::Vector& d  ,
  const Ostap::Math::GSL::Matrix& b  ,
  const bool                      Tb ,
  Ostap::Math::GSL::Matrix&       c  )
{
  // use GSL error handler sentry
  Ostap::Math::GSL::GSL_Error_Handler sentry ;

  // Dimensions of op(A): rows_A x cols_A
  const std::size_t rows_A = Ta ? a.nCols() : a.nRows() ;
  const std::size_t cols_A = Ta ? a.nRows() : a.nCols() ;

  // Dimensions of op(B): rows_B x cols_B
  const std::size_t rows_B = Tb ? b.nCols() : b.nRows() ;
  const std::size_t cols_B = Tb ? b.nRows() : b.nCols() ;

  // Size of the diagonal vector d
  const std::size_t dim_d = d.size() ;

  // Check dimension compatibility: cols(op(A)) == dim(d) && dim(d) == rows(op(B))
  if ( cols_A != dim_d || rows_B != dim_d )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::MDM: Incompatible matrix/vector dimensions" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }

  // Handle Aliasing (&a == &c || &b == &c):
  // Compute into a temporary matrix first to prevent overwriting inputs
  if ( &a == &c || &b == &c )
  {
    Matrix tmp { rows_A , cols_B } ;
    Ostap::StatusCode sc = MDM ( a , Ta , d , b , Tb , tmp ) ;
    if ( sc.isSuccess() ) { c.swap ( tmp ) ; }
    return sc ;
  }

  // Ensure output matrix c has correct dimensions (rows_A x cols_B)
  if ( c.nRows() != rows_A || c.nCols() != cols_B )
  { c = Matrix ( rows_A , cols_B ) ; }

  // OPTIMIZATION for Symmetric case: A * diag(d) * A^T or A^T * diag(d) * A
  // Uses dsyrk via column scaling with sqrt(d) if all d_k >= 0
  if ( &a == &b && Ta != Tb && gsl_vector_isnonneg ( d.vector() ) )
  {
    // 1. Prepare scaled matrix: A_scaled = op(A) * diag(sqrt(d))
    Matrix A_scaled { rows_A , dim_d } ;
    for ( std::size_t k = 0 ; k < dim_d ; ++k )
    {
      const double sqrt_d_k = std::sqrt ( d ( k ) ) ;
      if ( Ta )
      { for ( std::size_t i = 0 ; i < rows_A ; ++i )
        { A_scaled.set ( i , k , a ( k , i ) * sqrt_d_k ) ; } }
      else
        { for ( std::size_t i = 0 ; i < rows_A ; ++i )
          { A_scaled.set ( i , k , a ( i , k ) * sqrt_d_k ) ; } }
    }

    // 2. Compute upper triangle via dsyrk: C = 1.0 * A_scaled * A_scaled^T
    const int status = gsl_blas_dsyrk ( CblasUpper , CblasNoTrans , 1.0 , A_scaled.matrix() , 0.0 , c.matrix() ) ;
    if ( status )
    {
      gsl_error ( "Ostap::Math::GSL::MDM: Error from gsl_blas_dsyrk function" , __FILE__ , __LINE__ , status ) ;
      return ERROR_GSL + status ;
    }

    // 3. Mirror upper triangle to lower triangle
    gsl_matrix* m_c = c.matrix () ;
    for ( std::size_t i = 0 ; i < rows_A ; ++i )
    {
      for ( std::size_t j = i + 1 ; j < rows_A ; ++j )
      {
        const double val = gsl_matrix_get ( m_c , i , j ) ;
        gsl_matrix_set ( m_c , j , i , val ) ;
      }
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  }

  // General Case:
  // Step 1: Scale columns of op(A) by vector d: A_scaled = op(A) * diag(d)
  Matrix A_scaled { rows_A , dim_d } ;
  gsl_matrix* as = A_scaled.matrix() ;
  for ( std::size_t k = 0 ; k < dim_d ; ++k )
  {
    const double d_k = d ( k ) ;
    if ( Ta )
    { for ( std::size_t i = 0 ; i < rows_A ; ++i )
      { gsl_matrix_set ( as , i , k , a ( k , i ) * d_k ) ; } }
    else
    { for ( std::size_t i = 0 ; i < rows_A ; ++i )
      { gsl_matrix_set ( as , i , k , a ( i , k ) * d_k ) ; } }
  }

  // Step 2: Multiply A_scaled by op(B) using optimized general MM
  return MM ( A_scaled , false , b , Tb , c ) ;
}

// ===========================================================================
/* Element-wise multiplication of two diagonal matrices: d = a * b
 *  @param a [in]  first diagonal vector
 *  @param b [in]  second diagonal vector
 *  @param d [out] resulting diagonal vector
 *  @return status code (Ostap::StatusCode::SUCCESS on success)
 */
// ===========================================================================
Ostap::StatusCode 
Ostap::Math::GSL::DD 
( const Ostap::Math::GSL::Vector& a , 
  const Ostap::Math::GSL::Vector& b , 
  Ostap::Math::GSL::Vector&       d )
{
  if ( a.size () != b.size() )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::DD: vector size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  // 1. Create a copy of vector 'a' into a temporary object 'tmp'
  //    to guarantee safety against parameter aliasing (d == a or d == b)
  Vector tmp { a } ;
  //
  // 2. Delegate element-wise multiplication to GSL C-API: tmp_i = tmp_i * b_i
  const int status = gsl_vector_mul ( tmp.vector () , b.vector () ) ;
  if ( status ) 
  { 
    gsl_error ( "Ostap::Math::GSL::DD: gsl_vector_mul execution failed" , __FILE__ , __LINE__ , status ) ; 
    return ERROR_GSL + status ;
  }
  //     
  // 3. Move temporary result to output parameter d (O(1) pointer swap)
  d.swap ( tmp ) ; 
  //
  return Ostap::StatusCode::SUCCESS ;
}


// ===========================================================================
// Symmetric matrix scaling shortcut: r = d * m * d
// ===========================================================================
/* Symmetric scaling: r = d * m * d
 *  @param d [in]  diagonal vector
 *  @param m [in]  input matrix
 *  @param r [out] resulting matrix
 *  @return status code
 */
// ===========================================================================
Ostap::StatusCode 
Ostap::Math::GSL::DMD 
( const Ostap::Math::GSL::Vector& d , 
  const Ostap::Math::GSL::Matrix& m , 
  Ostap::Math::GSL::Matrix&       r )
{ return DMD ( d , m , d , r ) ; }


// ===========================================================================
// General asymmetric matrix scaling: r = a * m * b
// ===========================================================================
/* Asymmetric scaling: r = a * m * b
 *  @param a [in]  left diagonal vector (scales rows)
 *  @param m [in]  input matrix
 *  @param b [in]  right diagonal vector (scales columns) 
 *  @param r [out] resulting matrix
 *  @return status code
 */
// ===========================================================================
Ostap::StatusCode 
Ostap::Math::GSL::DMD 
( const Ostap::Math::GSL::Vector& a , 
  const Ostap::Math::GSL::Matrix& m , 
  const Ostap::Math::GSL::Vector& b , 
  Ostap::Math::GSL::Matrix&       r )
{
  const std::size_t rows = m.nRows () ;
  const std::size_t cols = m.nCols () ;
  //
  if ( a.size() != rows )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::DMD: vector/matrix  size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
   if ( b.size() != cols )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::DMD: matrix/vector size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  // 1. Calculate into a temporary matrix to prevent aliasing (e.g. r == m)
  Matrix tmp { rows , cols } ;
  gsl_matrix* t = tmp.matrix() ;
  for ( std::size_t i = 0 ; i < rows ; ++i ) 
  {
    const double scale_i = a [ i ] ;
    for ( std::size_t j = 0 ; j < cols ; ++j ) 
    {
      const double scale = scale_i * b [ j ] ;
      gsl_matrix_set ( t , i , j , scale * m ( i , j ) ) ;
    }
  }
  //
  // 2. Safely move the temporary result to output matrix r
  r.swap ( tmp ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}

// ===========================================================================
/** Apply permutation transformation to a square matrix: R = P * M * P^T
 *  
 *  @param p [in]  Permutation matrix P
 *  @param m [in]  Square matrix M (N x N)
 *  @param r [out] Resulting permuted matrix (N x N)
 *  @return Status code (Ostap::StatusCode::SUCCESS on success)
 *
 *  @note Safe against argument aliasing (e.g., PMPt(p, m, m))
 */
// ===========================================================================
Ostap::StatusCode Ostap::Math::GSL::PMP
( const Ostap::Math::GSL::Permutation& P ,
  const Ostap::Math::GSL::Matrix&      m ,
  Ostap::Math::GSL::Matrix&            r ) 
{ return PMP ( P , m , P , r ) ; } 

// ===========================================================================
/*  General asymmetric permutation transformation: R = Pl * M * Pr^T
 *  
 *  @param pl [in]  Left permutation matrix Pl (size matching M.k1())
 *  @param m  [in]  Input matrix M (M x N)
 *  @param pr [in]  Right permutation matrix Pr (size matching M.k2())
 *  @param r  [out] Resulting permuted matrix (M x N)
 *  @return Status code (Ostap::StatusCode::SUCCESS on success)
 *
 *  @note Safe against argument aliasing (e.g., PLMPrT(pl, m, pr, m))
 */
// ========================================================================
Ostap::StatusCode 
Ostap::Math::GSL::PMP
( const Ostap::Math::GSL::Permutation& PL ,
  const Ostap::Math::GSL::Matrix&      m  ,
  const Ostap::Math::GSL::Permutation& PR ,
  Ostap::Math::GSL::Matrix&            r  ) 
{
  //
  const std::size_t rows = m.nRows () ;
  const std::size_t cols = m.nCols () ;
  //
  if ( PL.size() != rows )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::PMP: permutation/matrix  size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
   if ( PR.size() != cols )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::PMP: matrix/permutation size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  // 1. Temporary buffer for aliasing safety (e.g., PLMPrT(pl, m, pr, m))
  Matrix tmp { rows , cols } ;

  // 2. Direct indexing: R(i, j) = M( Pl[i], Pr[j] )
  const gsl_permutation* pl = PL .permutation () ;
  const gsl_permutation* pr = PR .permutation () ;
  gsl_matrix*            t  = tmp.matrix      () ;  
  const gsl_matrix*      a  = m  .matrix      () ;
  for ( std::size_t i = 0 ; i < rows ; ++i )
  {
    const std::size_t pli = pl->data[i] ;
    for ( std::size_t j = 0 ; j < cols ; ++j )
    {
      const std::size_t prj = pr->data[j] ;
      const double      val = gsl_matrix_get ( a , pli , prj ) ;
      gsl_matrix_set ( t , i , j , val ) ;
    }
  }

  // 3. Move result to output matrix
  r.swap ( tmp ) ;
  // 
  return Ostap::StatusCode::SUCCESS ;
}

// ===========================================================================
/** General asymmetric combination: 
 *  \f$ R = P_l \cdot D_1 \cdot M \cdot D_2 \cdot P_r^T \f$
 *  
 *  @param pl [in]  Left permutation matrix Pl
 *  @param d1 [in]  Left diagonal vector D1
 *  @param m  [in]  Input matrix A (M x N)
 *  @param d2 [in]  Right diagonal vector D2
 *  @param pr [in]  Right permutation matrix Pr
 *  @param r  [out] Resulting matrix R (M x N)
 *  @return Status code
 */
// =========================================================================== 
Ostap::StatusCode Ostap::Math::GSL::PDM 
( const Ostap::Math::GSL::Permutation& pl ,
  const Ostap::Math::GSL::Vector&      d1 ,
  const Ostap::Math::GSL::Matrix&      m  ,
  const Ostap::Math::GSL::Vector&      d2 ,
  const Ostap::Math::GSL::Permutation& pr ,
  Ostap::Math::GSL::Matrix&            r  ) 
{
  //
  const std::size_t rows = m.nRows () ;
  const std::size_t cols = m.nCols () ;
  //
  if ( pl.size() != rows )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::PDM: permutation/matrix  size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  if ( pr.size() != cols )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::PDM: matrix/permutation size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  if ( pl.size() != d1.size()  )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::PDM: permutation/vector size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  if ( pr.size() != d2.size() )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "Ostap::Math::GSL::PDM: vector/permutation size mismatch" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //

  // 1. Temporary matrix to prevent argument aliasing (e.g. r == a)
  Matrix tmp { rows , cols } ;

  // 2. Direct single-pass index & scaling mapping
  const gsl_permutation* lp = pl .permutation () ;
  const gsl_permutation* rp = pr .permutation () ;
  const gsl_matrix*      a  = m  .matrix      () ;
  gsl_matrix*            t  = tmp.matrix      () ;
  const gsl_vector*      s1 = d1 .vector      () ;
  const gsl_vector*      s2 = d2 .vector      () ;      
  for ( std::size_t i = 0 ; i < rows ; ++i )
  {
    const std::size_t pli     = lp -> data [ i ] ;
    const double      scale_i = gsl_vector_get ( s1 , pli ) ;
    for ( std::size_t j = 0 ; j < cols ; ++j )
    {
      const std::size_t prj     = rp -> data [ j ] ;
      const double      scale_j = gsl_vector_get ( s2  , prj ) ;
      const double      val     = scale_i * gsl_matrix_get ( a , pli , prj ) * scale_j ;
      gsl_matrix_set ( t , i , j , val );
    }
  }

  // 3. Move temporary buffer to output matrix
  r.swap ( tmp ) ; 

  return Ostap::StatusCode::SUCCESS ;
}

// ===========================================================================
/** Symmetric combination of permutation and diagonal scaling: 
 *  \f$ R = P \cdot D \cdot M \cdot D \cdot P^T \f$
 *  
 *  @param p [in]  Permutation matrix P
 *  @param d [in]  Diagonal vector D
 *  @param m [in]  Input square matrix A (N x N)
 *  @param r [out] Resulting matrix R (N x N)
 *  @return Status code (Ostap::StatusCode::SUCCESS on success)
 *
 *  @note Safe against argument aliasing (e.g., PDADPt(p, d, a, a))
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::PDM 
( const Ostap::Math::GSL::Permutation& p  ,
  const Ostap::Math::GSL::Vector&      d  ,
  const Ostap::Math::GSL::Matrix&      m  ,
  Ostap::Math::GSL::Matrix&            r  ) 
{ return PDM  ( p , d , m , d , p , r ) ; } 
 
// ============================================================================
/*  Compute Moore-Penrose Pseudoinverse using SVD & MDM:
 *  \f$ A^+ = V \times \text{diag}(\sigma^+) \times U^T \f$
 *
 *  @param a     (INPUT)  Input matrix A (m x n)
 *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
 *  @param tol   (INPUT)  Tolerance for zeroing small singular values (< 0 for default)
 *  @return status code
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::PINV
( const Ostap::Math::GSL::Matrix& a      ,
  Ostap::Math::GSL::Matrix&       a_pinv ,
  double                          tol    )
{
  Ostap::Math::GSL::GSL_Error_Handler sentry ;

  const std::size_t m = a.nRows() ;
  const std::size_t n = a.nCols() ;

  // Handle Aliasing (&a == &a_pinv)
  if ( &a == &a_pinv )
  {
    Matrix tmp { n , m } ;
    Ostap::StatusCode sc = PINV ( a , tmp , tol ) ;
    if ( sc.isSuccess() ) { a_pinv.swap ( tmp ) ; }
    return sc ;
  }

  // Ensure target dimensions (n x m)
  if ( a_pinv.nRows() != n || a_pinv.nCols() != m )
  { a_pinv = Matrix ( n , m ) ; }

  // SVD requires working copy of A (m x n), matrix V (n x n), vector S (n)
  // GSL computes A = U * S * V^T, where A_copy holds U on output
  Matrix A_copy { a } ; // Holds U after SVD
  Matrix V      { n , n } ;
  Vector S      { n } ;
  Vector work   { n } ;

  // Compute SVD
  int status = gsl_linalg_SV_decomp ( A_copy.matrix() , V.matrix() , S.vector() , work.vector() ) ;
  if ( status )
  {
    gsl_error ( "Ostap::Math::GSL::PINV: SVD decomposition failed" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }

  // Determine cutoff threshold for singular values
  const double max_s = S ( 0 ) ;
  if ( tol < 0.0 )
  {
    const double eps = 2.2204460492503131e-16 ; // machine epsilon for double
    tol = std::max ( m , n ) * max_s * eps ;
  }

  // Build the inverted singular vector d = sigma^+
  Vector d { n } ;
  gsl_vector* vd = d.vector() ;
  for ( std::size_t i = 0 ; i < n ; ++i )
  {
    const double s_i = S ( i ) ;
    gsl_vector_set ( vd , i , ( s_i > tol ) ? ( 1.0 / s_i ) : 0.0 ) ;
  }

  // Compute A^+ = V * diag(d) * U^T  using our optimized MDM function!
  // V (n x n), d (n), A_copy holds U (m x n), transposed Tu = true -> U^T (n x m)
  return MDM ( V , false , d , A_copy , true , a_pinv ) ;
}

// ============================================================================
// Vector 
// ============================================================================

// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::Math::GSL::Vector::Vector
( const std::size_t N  ) 
  : m_vector ( gsl_vector_alloc ( N ) )
{
  //
  Ostap::Assert ( m_vector                         ,
                  "(GSL)Vector allocation failure" ,
                  "Ostap::Math::GSL::Vector" ,
                  VECTOR_ALLOCATION_FAILURE  , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::Math::GSL::Vector::Vector
( const std::size_t  N     ,   
  const double       value )
  : m_vector ( gsl_vector_alloc ( N ) )
{
  //
  Ostap::Assert ( m_vector                         ,
                  "(GSL)Vector allocation failure" ,
                  "Ostap::Math::GSL::Vector"       ,
                  VECTOR_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( std::isfinite ( value )          ,
                  "Cannot use !std::isfinite"      ,
                  "Ostap::Math::GSL::Vector"       ,
                  INVALID_SCALE                    , __FILE__ , __LINE__ ) ;
  gsl_vector_set_all ( m_vector , value ) ;
}
// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::Math::GSL::Vector::Vector
( const std::size_t                N       ,   
  const Ostap::Math::GSL::Vector::Zero /* zero */ ) 
  : m_vector ( gsl_vector_calloc ( N ) )  //    NB! calloc here! 
{
  //
  Ostap::Assert ( m_vector                         ,
                  "(GSL)Vector allocation failure" ,
                  "Ostap::Math::GSL::Vector"       ,
                  VECTOR_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::GSL::Vector::Vector  
( const Ostap::Math::GSL::Vector&  right ) 
  : m_vector ( gsl_vector_alloc ( right.m_vector->size ) )  
{
  //
  Ostap::Assert ( m_vector                         ,
                  "(GSL)Vector allocation failure" ,
                  "Ostap::Math::GSL::Vector"       ,
                  VECTOR_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
  //
  gsl_vector_memcpy ( m_vector , right.m_vector ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::GSL::Vector::Vector  
(       Ostap::Math::GSL::Vector&&  right ) 
  : m_vector ( right.m_vector )  
{
  right.m_vector = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-vector
// ============================================================================
Ostap::Math::GSL::Vector::~Vector () 
{
  if ( nullptr != m_vector ) 
  {
    gsl_vector_free ( m_vector ) ; 
    m_vector = nullptr ; 
  }
}
// ============================================================================
// copy assignement! 
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::operator=
( const Ostap::Math::GSL::Vector&  right ) 
{
  if ( &right == this ) { return *this ; }
  resize ( right.m_vector->size ) ;
  gsl_vector_memcpy ( m_vector , right.m_vector ) ;
  return *this ;
}
// ============================================================================
// move assignement! 
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::operator=
( Ostap::Math::GSL::Vector&& right ) 
{
  if ( &right == this ) { return *this ; }
  std::swap ( m_vector , right.m_vector ) ;
  return *this ;
}
// ============================================================================
// resize the vector
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::resize
( const std::size_t n )
{
  if ( n != m_vector->size )
  {
    gsl_vector_free ( m_vector ) ; 
    m_vector = gsl_vector_alloc ( n ) ;
    //
    Ostap::Assert ( m_vector                           ,
                    "(GSL)Vector allocation failure"   ,
                    "Ostap::Math::GSL::Vector::resize" ,
                    VECTOR_ALLOCATION_FAILURE          , __FILE__ , __LINE__ ) ;
    //
  }
  return *this ;
}
// ============================================================================
// resize the vector
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::resize
( const std::size_t n     ,
  const double       value ) 
{
  if ( !value ) { return resize ( n , Zero() ) ; }
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannot use !std::isfinite"         ,
                  "Ostap::Math::GSL::Vector::resize"  ,
                  INVALID_SCALE                       , __FILE__ , __LINE__ ) ;
  resize ( n ) ;
  gsl_vector_set_all ( m_vector , value ) ;
  return *this ;
}
// ============================================================================
// resize the vector
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::resize
( const std::size_t n     ,
  const Ostap::Math::GSL::Vector::Zero /* zero */ )  
{
  if ( n != m_vector->size )
    {
      gsl_vector_free ( m_vector ) ; 
      m_vector = gsl_vector_calloc ( n ) ; // CALLOC HERE
      //
      Ostap::Assert ( m_vector                           ,
                      "(GSL)Vector allocation failure"   ,
                      "Ostap::Math::GSL::Vector::resize" ,
                      VECTOR_ALLOCATION_FAILURE          , __FILE__ , __LINE__ ) ;
      //
    }
  else { gsl_vector_set_all ( m_vector , 0 ) ; }
  //
  return *this ;
}
// ============================================================================
// Numerical equality of two vectors 
// ============================================================================
bool Ostap::Math::GSL::Vector::equal
( const Ostap::Math::GSL::Vector& r ) const 
{
  if ( &r       == this       ) { return true  ; }
  //
  const std::size_t N { size ()} ;
  if ( N != r.size () ) { return false ; }
  //
  for ( std::size_t i = 0 ; i < N ; ++i )
  {
    const double e1 =   get ( i ) ;
    if ( !std::isfinite ( e1 ) ) { return false ; }
    const double e2 = r.get ( i ) ;
    if ( !std::isfinite ( e2 ) ) { return false ; }        
    if ( !s_equal ( e1 , e2 )  ) { return false ; }
  }
  //
  return true ;
}
// ============================================================================
// swap two vectors 
// ============================================================================
void Ostap::Math::GSL::Vector::swap
( Ostap::Math::GSL::Vector& right )
{ std::swap ( m_vector , right.m_vector ) ; }
// ============================================================================
// Are all elements numerically equal to zero?      
// ============================================================================
bool Ostap::Math::GSL::Vector::iszero   () const
{
  for ( std::size_t i = 0 ; i < size ()  ; ++i )
    { if ( !s_zero ( get ( i ) ) ) { return false ; } }
  return true ;
}
// ============================================================================
// Are all elements finite ? 
// ============================================================================
bool Ostap::Math::GSL::Vector::isfinite () const
{
  for ( std::size_t i = 0 ; i < size ()  ; ++i )
    { if ( !std::isfinite ( get ( i ) ) ) { return false ; } }
  return true ;
}
// ============================================================================
// dot product of two vect
// ============================================================================
double
Ostap::Math::GSL::Vector::dot
( const Ostap::Math::GSL::Vector& value ) const 
{
  Ostap::Assert ( this->size() == value.size () ,
                  "Cannot dot vectors of incompatible structure" ,
                  "Ostap::Math::GSL::Vector::dot"  ,
                  INVALID_GVECTOR, __FILE__ , __LINE__ ) ;
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //  
  double result = 0 ;
  //
  const int status = gsl_blas_ddot ( m_vector , value.m_vector , &result ) ;
  Ostap::Assert ( !status ,
                  "Error from gsl_blas_ddot function"  ,
                  "Ostap::Math::GSL::Matrix::multiply" , 
                  ERROR_GSL + status                   , __FILE__ , __LINE__ ) ; 
  //
  return result ;
}
// ============================================================================
// cross/tensor product of two vect
// ============================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::Vector::cross
( const Ostap::Math::GSL::Vector& value ) const 
{
  Matrix result { size() , value.size() } ;
  for ( std::size_t i = 0 ; i < size()  ; ++i )
  { for ( std::size_t j = 0 ; j < value.size() ; ++j )
    { result.set ( i , j , get ( i ) * value ( j ) ) ; } }
  //
  return result ;
}
// ============================================================================
// add a vector 
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::iadd
( const Ostap::Math::GSL::Vector& value )
{
  Ostap::Assert ( size()  == value.size ()                       ,
                  "Cannot add vectors of incompatible structure" ,
                  "Ostap::Math::GSL::Vector::iadd"               ,
                  INVALID_GVECTOR                                , __FILE__  , __LINE__ ) ;
  gsl_vector_add ( m_vector , value.m_vector ) ;
  return *this ;
}
// ============================================================================
// add a cboisrant 
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::iadd
( const double value )
{
  Ostap::Assert ( std::isfinite ( value )          ,
                  "Cannot add !std::isfinite"      ,
                  "Ostap::Math::GSL::Vector::iadd" ,
                  INVALID_SCALE                    , __FILE__ , __LINE__ ) ;
  gsl_vector_add_constant  ( m_vector , value ) ;
  return *this ;
}
// ============================================================================
// subtract vector 
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::isub
( const Ostap::Math::GSL::Vector& value )
{
  Ostap::Assert ( this->size() == value.size ()                       ,
                  "Cannot subtract vectors of incompatible structure" ,
                  "Ostap::Math::GSL::Vector::isub"                    ,
                  INVALID_GVECTOR                                     , __FILE__ , __LINE__ ) ;
  gsl_vector_sub ( m_vector , value.m_vector ) ;
  return *this ;
}
// ============================================================================
// scale vector  
// ============================================================================
Ostap::Math::GSL::Vector&                             
Ostap::Math::GSL::Vector::imul
( const double value )
{
  if ( 1 == value ) { return *this ; }
  Ostap::Assert ( std::isfinite ( value )          ,
                  "Cannot scale by !std::isfinite" ,
                  "Ostap::Math::GSL::Vector::imul" ,
                  INVALID_SCALE                    , __FILE__ , __LINE__ ) ;
  //
  gsl_vector_scale ( m_vector , value ) ;
  return *this;
}
// ============================================================================
// multiply by matrix 
// ============================================================================
Ostap::Math::GSL::Vector&
Ostap::Math::GSL::Vector::imul 
( const Ostap::Math::GSL::Matrix& value )
{
  Ostap::Assert ( size() == value.nRows()                                   ,
                  "Cannot multiply vector&matrix of incompatible structure" ,
                  "Ostap::Math::GSL::Vector::imul"                          ,
                  INVALID_GMATRIX                                           , __FILE__ , __LINE__ ) ;
  
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //  
  Vector result { value.nCols() } ;
  //
  const int status = gsl_blas_dgemv ( CblasTrans       , // ATTENTIO!!! 
                                      1.0              ,
                                      value .matrix () ,
                                      m_vector         ,
                                      0.0              ,
                                      result.vector () ) ;
  //
  Ostap::Assert ( !status ,
                  "Error from gsl_blas_dgemv function" ,
                  "Ostap::Math::GSL::Matrix::imul"     , 
                  ERROR_GSL + status                   , __FILE__ , __LINE__ ) ; 
  //
  this->swap ( result ) ;
  //
  return *this ;
}

// ============================================================================
// multiply by matrix 
// ============================================================================
Ostap::Math::GSL::Vector
Ostap::Math::GSL::Vector::multiply  
( const Ostap::Math::GSL::Matrix& value ) const
{
  Ostap::Assert ( this->size() == value.nRows()                             ,
                  "Cannot multiply vector&matrix of incompatible structure" ,
                  "Ostap::Math::GSL::Vector::isub"                          ,
                  INVALID_GMATRIX                                           , __FILE__ , __LINE__ ) ;
  
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //  
  Vector result { value.nCols() } ;
  //
  const int status = gsl_blas_dgemv ( CblasTrans       , // ATTENTIO!!! 
                                      1.0              ,
                                      value .matrix () ,
                                      m_vector         ,
                                      0.0              ,
                                      result.vector () ) ;
  //
  Ostap::Assert ( !status ,
                  "Error from gsl_blas_dgemv function" ,
                  "Ostap::Math::GSL::Matrix::multiply" , 
                  ERROR_GSL + status                   , __FILE__ , __LINE__ ) ; 
  //
  return result ;
}

// ============================================================================
// constructor: allocate the permutation 
// ============================================================================
Ostap::Math::GSL::Permutation::Permutation
( const std::size_t N ) 
  : m_permutation ( gsl_permutation_calloc ( N ) ) 
{
  Ostap::Assert ( m_permutation                         ,
                  "(GSL)Permutation allocation failure" ,
                  "Ostap::Math::GSL::Permutation"       ,
                  PERMUTATION_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::GSL::Permutation::Permutation 
( const Ostap::Math::GSL::Permutation&  right ) 
  : m_permutation  ( gsl_permutation_alloc ( right.m_permutation->size ) )
{
  //
  Ostap::Assert ( m_permutation                         ,
                  "(GSL)Permutation allocation failure" ,
                  "Ostap::Math::GSL::Permutation"       ,
                  PERMUTATION_ALLOCATION_FAILURE        , __FILE__ , __LINE__ ) ;
  //
  gsl_permutation_memcpy ( m_permutation , right.m_permutation ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::GSL::Permutation::Permutation 
( Ostap::Math::GSL::Permutation&& right ) 
  : m_permutation ( right.m_permutation )  
{
  right.m_permutation = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-permutation
// ============================================================================
Ostap::Math::GSL::Permutation::~Permutation () 
{
  if ( nullptr != m_permutation ) 
  {
    gsl_permutation_free ( m_permutation ) ; 
    m_permutation = nullptr ; 
  }
}
// ============================================================================
// copy assignement! 
// ============================================================================
Ostap::Math::GSL::Permutation&
Ostap::Math::GSL::Permutation::operator=
( const Ostap::Math::GSL::Permutation&  right ) 
{
  if ( &right == this ) { return *this ; }
  if ( m_permutation->size != right.m_permutation->size )
  {
    gsl_permutation_free ( m_permutation ) ;
    m_permutation = gsl_permutation_alloc ( right.m_permutation->size ) ;
    Ostap::Assert ( m_permutation                              ,
                    "(GSL)Permutation allocation failure"      ,
                    "Ostap::Math::GSL::Permutation::operator=" ,
                    PERMUTATION_ALLOCATION_FAILURE             , __FILE__ , __LINE__ ) ;
    //      
  }
  gsl_permutation_memcpy ( m_permutation , right.m_permutation ) ;
  return *this ;
}
// ============================================================================
// move assignement! 
// ============================================================================
Ostap::Math::GSL::Permutation&
Ostap::Math::GSL::Permutation::operator=
( Ostap::Math::GSL::Permutation&& right ) 
{
  if ( &right == this ) { return *this ; }
  std::swap ( m_permutation , right.m_permutation ) ;
  return *this ;
}
// ============================================================================
// resize the permutation
// ============================================================================
Ostap::Math::GSL::Permutation&
Ostap::Math::GSL::Permutation::resize
( const std::size_t n )
{
  // no action ? 
  if ( size() == n ) { return *this ; }
  //  
  gsl_permutation_free ( m_permutation ) ;
  m_permutation = gsl_permutation_alloc ( n ) ;
  Ostap::Assert ( m_permutation                           ,
                  "(GSL)Permutation allocation failure"   ,
                  "Ostap::Math::GSL::Permutation::resize" ,
                  PERMUTATION_ALLOCATION_FAILURE          , __FILE__ , __LINE__ ) ;
  //      
  return *this ;
}
// ============================================================================
// swap two permutation 
// ============================================================================
void Ostap::Math::GSL::Permutation::swap
( Ostap::Math::GSL::Permutation& right )
{ std::swap ( m_permutation , right.m_permutation ) ; }
// ============================================================================
// valid permutation ?
// ============================================================================
bool
Ostap::Math::GSL::Permutation::valid () const
{ return GSL_SUCCESS == gsl_permutation_valid ( m_permutation ) ; }

// ============================================================================
// apply permutation to the matrix 
// ============================================================================
Ostap::Math::GSL::Matrix
Ostap::Math::GSL::Permutation::apply 
( const Ostap::Math::GSL::Matrix& value ) const 
{
  Ostap::Assert ( size() == value.nRows()                      ,
                  "Mismatch for permutation/matrix structure!" ,
                  "Ostap::Math::GLS::Permutation::apply"       , 
                  INVALID_PERMUTATION                          , __FILE__ , __LINE__   ) ;
  //
  // Matrix result { value } ;
  // result.permute_rows ( *this  ) ;
  //
  Matrix result { *this } ;
  return result * value ;
}
// ============================================================================
// print matrix to the stream
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const Ostap::Math::GSL::Matrix& m ,
  std::ostream&              s )
{ return toStream ( *m.matrix() , s ) ; }
// ============================================================================
// print vector to the stream 
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const Ostap::Math::GSL::Vector& v ,
  std::ostream&              s )
{ return toStream ( *v.vector () , s ) ; }
// ============================================================================
// print permutation to the stream 
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const Ostap::Math::GSL::Permutation& p ,
  std::ostream&                   s )
{ return toStream ( *p.permutation () , s ) ; }
// ============================================================================
/* print GSL-vector to the stream 
 *  @param v the vector 
 *  @param s the stream 
 *  @return the stream 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2012-05-28
 */
// ============================================================================
std::ostream& 
Ostap::Utils::toStream 
( const gsl_vector&  v , 
  std::ostream&      s ) 
{
  s << "[ " ;
  for ( std::size_t i = 0 ; i < v.size ; ++i ) 
  {
    if ( 0 != i ) { s << " , " ; } 
    s << gsl_vector_get ( &v , i ) ;
  }
  s << "]" ;
  return s ;
}
// ============================================================================
/*  print GSL-matrix to the stream 
 *  @param m the matrix 
 *  @param s the stream 
 *  @return the stream 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2012-05-28
 */
// ============================================================================
std::ostream& 
Ostap::Utils::toStream 
( const gsl_matrix&  m , 
  std::ostream&      s ) 
{
  for ( std::size_t i = 0 ; i < m.size1 ; ++i ) 
  {
    s << " | " ;
    for ( std::size_t j = 0 ; j < m.size2 ; ++j ) 
    {
      if ( 0 != j ) { s << ", " ; }
      s << Ostap::format ( "|%11.5g|" , gsl_matrix_get ( &m , i , j ) ) ;
    }
    s << " | "<< std::endl ;
  }
  return s ;  
}
// ===========================================================================
/*  print GSL-permutation to the stream 
 *  @param p the permutation
 *  @param s the stream 
 *  @return the stream 
 */    
// ===========================================================================
std::ostream& Ostap::Utils::toStream 
( const gsl_permutation& p , 
  std::ostream&          s )
{
  s << "( " ;
  for ( std::size_t i = 0 ; i < p.size ; ++i ) 
  {
    if ( 0 != i ) { s << " , " ; } 
    s << gsl_permutation_get ( &p , i ) ;
  }
  s << ")" ;
  return s ;
}
// ============================================================================

// ============================================================================
// few more utilties 
// ============================================================================
// get the max element
// ============================================================================
double Ostap::Math::max_element ( const Ostap::Math::GSL::Vector& v )
{ return gsl_vector_max ( v.vector () ) ; }  
// ============================================================================
// get the min element
// ============================================================================
double Ostap::Math::min_element ( const Ostap::Math::GSL::Vector& v )
{ return gsl_vector_max ( v.vector () ) ; }  
// ============================================================================
// get the max element
// ============================================================================
double Ostap::Math::max_element ( const Ostap::Math::GSL::Matrix& m )
{ return gsl_matrix_max ( m.matrix () ) ; }    
// ============================================================================
// get the min element
// ============================================================================
double Ostap::Math::min_element ( const Ostap::Math::GSL::Matrix& m ) 
{ return gsl_matrix_min ( m.matrix () ) ; }    

// ============================================================================
// get the element with maximal absolute value 
// ============================================================================
double Ostap::Math::maxabs_element
( const Ostap::Math::GSL::Matrix& m )
{
  double result = -1 ;
  for ( std::size_t i = 0 ; i < m.nRows () ; ++i )
  { for ( std::size_t j = 0 ; j < m.nCols () ; ++j )
    { result = std::max ( result , std::abs ( m ( i , j ) ) ) ; } }
  return result ;
}
// ============================================================================
// get the element with maximal absolute value 
// ============================================================================
double Ostap::Math::maxabs_element
( const Ostap::Math::GSL::Vector& v )
{
  double result = -1 ;
  for ( std::size_t i = 0 ; i < v.size() ; ++i )
  { result = std::max ( result , std::abs ( v ( i ) ) ) ; }
  return result ;
}
// ============================================================================

// ============================================================================
// Is this matrix symmetric ?
// ============================================================================
bool Ostap::Math::symmetric ( const Ostap::Math::GSL::Matrix& m )
{
  const std::size_t N = m.nRows () ;
  const std::size_t M = m.nCols () ;
  //
  if ( N != M ) { return false ; }
  //
  for ( std::size_t i = 0 ; i < N ; ++i )
  {
    const double cii = m.get ( i , i ) ;
    if ( !std::isfinite ( cii ) ) { return false ; } 
    for ( std::size_t j = i + 1 ; j < M ; ++j )
    {
      const double cij = m.get ( i , j ) ;
      if ( !std::isfinite ( cij ) ) { return false ; } 
      const double cji = m.get ( j , i ) ;
      if ( !std::isfinite ( cji ) ) { return false ; } 
      if ( !s_equal ( cij , cji ) ) { return false ; }
    }
  }
  return true ;
}
// ============================================================================
/*  Can this matrix be symmetric & positive-definite ?
 *  - Finite 
 *  - Diagonal elements are finite and positive
 *  - Symmetric
 *  - have CholeskyDecomposition
 */
// ============================================================================
bool Ostap::Math::symmetric_positive_definite
( const Ostap::Math::GSL::Matrix& m )
{
  //
  const std::size_t N = m.nRows () ;
  const std::size_t M = m.nCols () ;
  //
  // (1) Square 
  if ( N != M ) { return false ; } 
  //
  for ( std::size_t i = 0 ; i < N ; ++i )
  {
    const double cii = m.get ( i , i ) ;
    if ( !std::isfinite ( cii )        ) { return false ; } // diagonal is finite 
    if ( cii <= 0  || s_zero ( cii )   ) { return false ; } // diagonal is positive 
    
    for ( std::size_t j = 0 ; j < i ; ++j )
    {
      // already checked in "i"-loop 
      const double cjj = m.get ( j , j ) ;
      //      
      const double cij = m.get ( i , j ) ;
      if ( !std::isfinite ( cij )      ) { return false ; } // element is finite
      //
      const double cji = m.get ( j , i ) ;
      if ( !std::isfinite ( cji )      ) { return false ; } // element is finite
      //
      if ( !s_equal ( cij , cji )      ) { return false ; } // matrix is symmetric
      //
    }
    // 
  }

  // 
  // The final shot: try Cholesky decomposition 
  //
  Ostap::Math::GSL::Matrix      A { m } ;
  Ostap::Math::GSL::Permutation P { N } ;
  Ostap::Math::GSL::Vector      S { N } ;
  //
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Ignore sentry { true } ;
  const int status = gsl_linalg_pcholesky_decomp2 ( A.matrix      () ,
                                                    P.permutation () ,
                                                    S.vector      () ) ;  
  return GSL_SUCCESS == status ;
}


// ============================================================================
/*  Can this matrix be a covariance matrix?
 *  - Square
 *  - Finite 
 *  - Symmetric 
 *  - Diagonal elements are finite and positive
 *  - Off-diagonal elements are finite and not too large 
 *  - have CholeskyDecomposition
 */
// ============================================================================
bool Ostap::Math::covariance_matrix 
( const Ostap::Math::GSL::Matrix& m )
{
  //
  const std::size_t N = m.nRows () ;
  const std::size_t M = m.nCols () ;
  //
  // (1) Square 
  if ( N != M ) { return false ; } 
  //
  for ( std::size_t i = 0 ; i < N ; ++i )
  {
    const double cii = m.get ( i , i ) ;
    if ( !std::isfinite ( cii )        ) { return false ; } // diagonal is finite 
    if ( cii <= 0  || s_zero ( cii )   ) { return false ; } // diagonal is positive 
    
    for ( std::size_t j = 0 ; j < i ; ++j )
    {
      // already checked in "i"-loop 
      const double cjj = m.get ( j , j ) ;
      //      
      const double cij = m.get ( i , j ) ;
      if ( !std::isfinite ( cij )      ) { return false ; } // element is finite
      //
      const double cji = m.get ( j , i ) ;
      if ( !std::isfinite ( cji )      ) { return false ; } // element is finite
      //
      if ( !s_equal ( cij , cji )      ) { return false ; } // matrix is symmetric
      //
      const double cc = cii * cjj ;
      //
      if  ( cc < cij * cij             ) { return false ; }  // off-diagonal element is too large 
    }
    // 
  }

  // 
  // The final shot: try Cholesky decomposition 
  //
  Ostap::Math::GSL::Matrix      A { m } ;
  Ostap::Math::GSL::Permutation P { N } ;
  Ostap::Math::GSL::Vector      S { N } ;
  //
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Ignore sentry { true } ;
  const int status = gsl_linalg_pcholesky_decomp2 ( A.matrix      () ,
                                                    P.permutation () ,
                                                    S.vector      () ) ;  
  return GSL_SUCCESS == status ;
}



// ============================================================================
// Actual Linear Algebra starts here 
// ============================================================================

// ============================================================================
/* "in-place" LU decomposition  
 *  \f$ PA = LU \f$, where 
 *   - A is \f$ M \times N \f$ matrix 
 *   - P is \f$ M \times M \f$ permutation matrix  
 *   - L is \f$ M \times \min (M,N)\f$  lower trianhular matris 
 *   - U is \f$ \min (M,N) \times N \f$ upper trianhular matris 
 * 
 * For square matrices:
 *   - L is a lower unit triangular matrix
 *   - U is upper triangular
 * 
 * For \f$ M>N \f$: 
 *   - L is a unit lower trapezoidal matrix of size \f$ M\timex N \f$ 
 * 
 * For \f$ M < N \f$: 
 *  - U is upper trapezoidal of size \f$ M \times M \f$  
 *
 *  For square matrices this decomposition can be used to convert the linear 
 *  system \f$ Ax=b\f$  into a pair of triangular systems, 
 *  \f$ Lu=Pb\f$ and  \f$ Ux=y\f$, which can be solved by forward and 
 *   back-substitution. 
 *  Note that the LU decomposition is also valid for singular matrices.
 *
 *  @parameter A  (update) input/update MxN  marix 
 *  @return    M-permutations 
 * 
 *  The matrix at the end contans two matrices: 
 * 
 *  On output the diagonal and upper triangular (or trapezoidal) part of the 
 *  input matrix A contain the matrix U. 
 *  The lower triangular (or trapezoidal) part of the input matrix (excluding 
 *  the diagonal) contains L. The diagonal elements of U are unity, and are not stored.
 *
 *  The permutation matrix P is encoded in the permutation p on output.
 *  The j-th column of the matrix P is given by the j-th column of the 
 *   identity matrix, where \f$ k = p_j \f%  
 *  the j-th element of the permutation vector. 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::PLU 
( Ostap::Math::GSL::Matrix&      A ,
  Ostap::Math::GSL::Permutation& P )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //    
  P.resize ( A.nRows() ) ;
  int signum = 0 ;
  //    
  int status = gsl_linalg_LU_decomp ( A.matrix() , P.permutation() , &signum ) ;
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_LU_decomp" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================
/*  perfom  LU decomposition  
 *  @param  A  (INOUT)         input matrix 
 *  @param  LU (UPDATE/OUTPUT) output LU matrix 
 *  @return M-permutations 
 *  @see gsl_linalg_LU_decomp 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::PLU
( const Ostap::Math::GSL::Matrix& A  ,
  Ostap::Math::GSL::Permutation&  P  , 
  Ostap::Math::GSL::Matrix&       LU )
{
  LU = A ;
  return PLU ( LU , P ) ; 
}
// ============================================================================
/*  perfom LU decomposition  
 *  @param  A   (INOUT)         input matrix 
 *  @param  L   (UPDATE/OUTPUT) lower triangular matrix 
 *  @param  U   (UPDATE/OUTPUT) upper triangular matrix 
 *  @return M-permutation 
 *  @see gsl_linalg_LU_decomp 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::PLU
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Permutation&  P ,   
  Ostap::Math::GSL::Matrix&       L ,
  Ostap::Math::GSL::Matrix&       U )
{
  //
  Ostap::Math::GSL::Matrix  LU { A } ;  
  Ostap::StatusCode sc = PLU ( LU , P ) ;
  if ( sc.isFailure() ) { return sc ; }
  //
  const std::size_t M { A.nRows () } ;
  const std::size_t N { A.nCols () } ;
  const std::size_t K { std::min ( M , N ) } ; 
  //
  L.resize ( M , K , Ostap::Math::GSL::Matrix::Id   () ) ; 
  U.resize ( K , N , Ostap::Math::GSL::Matrix::Zero () ) ; 
  //
  // =========================================================================
  // Fill L-matrix 
  // =========================================================================
  for ( std::size_t i = 0 ; i < K ; ++i )
  { for ( std::size_t j = 0 ; j < i  ; ++j )
    { L.set ( i , j , LU ( i , j ) ) ; } } 
  for ( std::size_t i = K ; i < M ; ++i )
  { for ( std::size_t j = 0 ; j < K ; ++j )
    { L.set ( i , j , LU ( i , j ) ) ; } }  
  // =========================================================================
  // Fill U-matrix 
  // =========================================================================
  for ( std::size_t i = 0 ; i < K ; ++i )
    { for ( std::size_t j = i ; j < N ; ++j ) { U.set ( i , j , LU ( i , j ) ) ; } }
  //
  return sc ;
}
// ============================================================================

// ============================================================================
// QR decomposition with column pivoting 
// ============================================================================

// ============================================================================
namespace 
{
  // ==========================================================================
  /** make QR Decomposion of matrix A : \f$ AP = QR\f$ where 
   *  - A is input                 MxN matrix  
   *  - P is permuutation matrix   NxN 
   *  - Q is orthogonal matrix     MxM 
   *  - R is right triaular matrix MxN 
   *  - r is the reciprocal condition number of R
   *  
   *  @param A  (input) the matrix to decopose 
   *  @param P  (output/update) permutation matrix P
   *  @param Q  (output/update) orthogonal matrix Q 
   *  @param R  (output/update) rigth triangular matrix R 
   *  @param r  (output/update) reciprocal condition number of R
   *  @return status code 
   */
  Ostap::StatusCode _PQR_
  ( const Ostap::Math::GSL::Matrix& A ,
    Ostap::Math::GSL::Permutation&  P , 
    Ostap::Math::GSL::Matrix&       Q ,
    Ostap::Math::GSL::Matrix&       R , 
    double*                         r = nullptr )
  {
    // use GSL: 
    Ostap::Math::GSL::GSL_Error_Handler sentry ;
    //    
    const std::size_t M = A.nRows() ;
    const std::size_t N = A.nCols() ;
    const std::size_t K = std::min ( M , N ) ; 
    //
    P.resize ( N     ) ;
    //
    Q.resize ( M , M ) ;
    R.resize ( M , N , Ostap::Math::GSL::Matrix::Zero() ) ;
    //
    Ostap::Math::GSL::Vector      tau  { K } ;
    Ostap::Math::GSL::Vector      norm { N } ;
    //
    int signum = 0 ;
    int status = gsl_linalg_QRPT_decomp2
      ( A.matrix      () ,
        Q.matrix      () ,
        R.matrix      () ,
        tau.vector    () ,
        P.permutation () , &signum , norm.vector() ) ;
    //
    if ( status )
    {
      gsl_error ( "PQR: Error from gsl_linalg_QRPT_decomp2" , __FILE__ , __LINE__ , status ) ;
      return ERROR_GSL + status ;
    }
    //
    if ( r )
    {
      Ostap::Math::GSL::Vector ws { 3 * N } ;
      status = gsl_linalg_QRPT_rcond ( R.matrix () , r , ws.vector() ) ;
      if ( status )
      {
        gsl_error ( "PQR: Error from gsl_linalg_QRPT_rcond" , __FILE__ , __LINE__ , status ) ;
        return ERROR_GSL + status ;
      }
    }
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
}
// ============================================================================

// ============================================================================
/*  mape QR Decomposion of matrix A : \f$ AP = QR\f$ where 
 *  - A is input                 MxN matrix  
 *  - P is permuutation matrix   NxN 
 *  - Q is orthogonal matrix     MxM 
 *  - R is right triaular matrix MxN 
 *  
 *  @param A  (input) the matrix to decopose 
 *  @param P  (output/update) permutation matrix P
 *  @param Q  (output/update) orthogonal matrix Q 
 *  @param R  (output/update) rigth triangular matrix R 
 *  @return status code  
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::PQR
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Permutation&  P , 
  Ostap::Math::GSL::Matrix&       Q ,
  Ostap::Math::GSL::Matrix&       R )
{ return _PQR_ ( A , P , Q , R , nullptr ) ; }
// ============================================================================
/*  make QR Decomposion of matrix A : \f$ AP = QR\f$ where 
 *  - A is input                 MxN matrix  
 *  - P is permuutation matrix   NxN 
 *  - Q is orthogonal matrix     MxM 
 *  - R is right triaular matrix MxN 
 *  - r is the reciprocal condition number of R
 *  
 *  @param A  (input) the matrix to decopose 
 *  @param P  (output/update) permutation matrix P
 *  @param Q  (output/update) orthogonal matrix Q 
 *  @param R  (output/update) rigth triangular matrix R 
 *  @param r  (output/update) reciprocal condition number of R
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::PQR
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Permutation&  P , 
  Ostap::Math::GSL::Matrix&       Q ,
  Ostap::Math::GSL::Matrix&       R , 
  double&                         r )
{ return _PQR_ ( A , P , Q , R , &r ) ; }
// ============================================================================


// ============================================================================
// LQ decomposition
// ============================================================================
/* LQ decomposition of matrix A: \f$ A = LQ\f$, where 
 *  - L is lower trapezoidal MxN 
 *  - Q is orthogonal NxN 
 */ 
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::LQ
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Matrix&       L ,
  Ostap::Math::GSL::Matrix&       Q )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //    
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  const std::size_t K = std::min ( M , N ) ;
  //
  L.resize ( M , N , Ostap::Math::GSL::Matrix::Zero() ) ;
  Q.resize ( N , N ) ;
  //
  Vector tau { K } ;
  Matrix R   { A } ;
  //
  int status = gsl_linalg_LQ_decomp ( R.matrix() , tau.vector() ) ;
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_LQ_decomp function" , __FILE__ , __LINE__ , status ) ; 
    return ERROR_GSL + status ;
  }
  //
  status = gsl_linalg_LQ_unpack ( R.matrix   () ,
                                  tau.vector () ,
                                  Q.matrix   () , 
                                  L.matrix   () ) ;
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_LQ_unpack function" , __FILE__ , __LINE__ , status ) ; 
    return ERROR_GSL + status ;
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================

// ============================================================================
// QL decomposition
// ============================================================================
/*  QL decomposition of matrix A: \f$ A = QL\f$, where 
 *  - Q is orthogonal MxM
 *  - L is lower trapezoidal MxN 
 */ 
// ============================================================================

#define OSTAP_GSL_VERSION(a,b) (1000*(a)+(b))
#define OSTAP_GSL_CODE_VERSION OSTAP_GSL_VERSION(GSL_MAJOR_VERSION, GSL_MINOR_VERSION)

Ostap::StatusCode Ostap::Math::GSL::QL
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Matrix&       Q ,
  Ostap::Math::GSL::Matrix&       L )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
#if OSTAP_GSL_CODE_VERSION < OSTAP_GSL_VERSION ( 2 , 7 )
  // 
  int status = GSL_ERROR ;
  gsl_error ( "For A=QL decomposition GSL version>2.6 is needed", __FILE__ , __LINE__ , status ) ;
  return GSL_ERROR + GSL_VERSION_IS_TOO_OLD ;
  //
#else
  //
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  //
  Q.resize ( M , M ) ;
  L.resize ( M , N , Ostap::Math::GSL::Matrix::Zero() ) ;
  //
  Vector tau { N } ;
  Matrix R   { A } ;
  //
  int status = gsl_linalg_QL_decomp ( R.matrix() , tau.vector() ) ;
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_QL_decomp"  , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  status = gsl_linalg_QL_unpack ( R.matrix   () ,
                                  tau.vector () ,
                                  Q.matrix   () , 
                                  L.matrix   () ) ;
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_QL_unpack"  , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  return Ostap::StatusCode::SUCCESS ;
  //
#endif 
}
// ============================================================================

// ============================================================================
// COD decomposition
// ============================================================================
/*  COD - Complete Orthogonal Decomposion
 *  \f$ AP = Q R Z^T \f$ 
 *  - A input MxN matrix 
 *  - P is permutation matrix 
 *  - Q is MxM orthogonal matrix 
 *  - R is 2x2 block matrix with top-left blobck being right triangular matrix and
 *    other blocks are zeroes 
 *  - Z is NxN orthogonal matrix 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Math::GSL::COD
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Permutation&  P , 
  Ostap::Math::GSL::Matrix&       Q ,
  Ostap::Math::GSL::Matrix&       R ,
  Ostap::Math::GSL::Matrix&       Z )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  const std::size_t K = std::min ( M , N ) ;
  //
  P.resize ( N ) ; 
  //
  Q.resize ( M , M ) ;
  R.resize ( M , N , Ostap::Math::GSL::Matrix::Zero() ) ;
  Z.resize ( N , N ) ;
  //
  Matrix D     { A } ;
  Vector tau_Q { K } ; 
  Vector tau_Z { K } ;
  Vector work  { N } ;
  //
  std::size_t  rank ;
  int status = gsl_linalg_COD_decomp
    ( D.matrix      () ,
      tau_Q.vector  () , 
      tau_Z.vector  () ,
      P.permutation () ,
      &rank            ,
      work.vector()  ) ;
  //
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_COD_decomp" ,__FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ; 
  }
  //
  status = gsl_linalg_COD_unpack
    ( D.matrix() ,
      tau_Q.vector () ,
      tau_Z.vector () ,
      rank            ,
      Q.matrix     () ,
      R.matrix     () ,
      Z.matrix     () ) ;
  //
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_COD_unpack" ,__FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ; 
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
// SVD decomposition
// ============================================================================
/* SVD : Singular Value Decomposition  \f$ A = U S V^T\f$   
 *  - A input MxN matrix 
 *  - K = min ( M , N ) : 
 *  - U MxK orthogonal matrix 
 *  - S KxK Diagonal matrix of singular values 
 *  - V NxK orthogonal matrix 
 *  @param A (input)  input matrix A 
 *  @param U (update) orthogonal matrix U 
 *  @param V (update) orthogonal matrix V 
 *  @param golub (input) use Golub or Jacobi algorithm? 
 *  @return vector of singular values 
 * -  Jacobi algorithm is more prrcise  and Golub algorithm is more CPU efficient 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::SVD
( const Ostap::Math::GSL::Matrix& A     ,
  Ostap::Math::GSL::Vector&       S     ,
  Ostap::Math::GSL::Matrix&       U     ,
  Ostap::Math::GSL::Matrix&       V     ,
  const bool                      golub )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  //
  if ( M < N ) { return SVD ( A.T() , S , V , U ) ; }
  //
  // replace U with A 
  U = A ;
  V.resize ( N , N ) ;
  //
  S.resize ( N ) ;
  //
  /// Use one sided Jacobi orthogonalization 
  if ( !golub )
  {
    int status = gsl_linalg_SV_decomp_jacobi 
      ( U.matrix () ,
        V.matrix () ,
        S.vector () ) ;
    if ( status )
    {
      gsl_error ( "Error from gsl_linalg_SV_decomp_jacobi" , __FILE__ , __LINE__ , status ) ;
      return ERROR_GSL + status ;
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  }
  //
  // workspace 
  Vector work { N } ;
  //
  /// standard Golub' algorithms 
  if ( M < 4 * N )
  {
    int status = gsl_linalg_SV_decomp
      ( U.matrix    () ,
        V.matrix    () ,
        S.vector    () ,
        work.vector () ) ;
    if ( status )
    {
      gsl_error ( "Error from gsl_linalg_SV_decomp" , __FILE__ , __LINE__ , status ) ;
      return ERROR_GSL + status ;
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  }
  
  /// additional workspace 
  Matrix X { N , N } ;
  
  /// modified Golub algorithm for M>>N
  int status = gsl_linalg_SV_decomp_mod 
    ( U.matrix    () ,
      X.matrix    () , 
      V.matrix    () ,
      S.vector    () ,
      work.vector () ) ;
  if ( status )
  {
    gsl_error ( "Error from gsl_linalg_SV_decomp_mod" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}


// ===========================================================================
/** LLT : Cholesky decomposition of positive definite matrix \f$ A = L L^T\f$, 
 *  Only lower triangular part of the matrix A is used.
 *  @param A (input)  input MxM matrix
 *  @param L (update) lower triangular matrix
 *  @return status code
 */  
// ===========================================================================
Ostap::StatusCode Ostap::Math::GSL::LLT
( const Ostap::Math::GSL::Matrix& A , 
  Ostap::Math::GSL::Matrix&       L ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  // 
  if ( N != M ) 
  {
    gsl_error ( "LLT: matrix is not square" , __FILE__ , __LINE__ , GSL_EBADLEN ) ;
    return MATRIX_IS_NOT_SQUARE ;
  }
  //
  if ( L.nRows() != M || L.nCols() != N ) { L.resize ( M , N , 0 ) ; }
  //
  Matrix aux { A } ;
  //
  const int status = gsl_linalg_cholesky_decomp1 ( aux.matrix() ) ;
  if ( status ) 
  {
    gsl_error ( "LLT: Error from gsl_linalg_cholesky_decomp1" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }  
  /// copy results to L 
  gsl_matrix* l = L.matrix() ;
  for ( std::size_t i = 0 ; i < M ; ++i ) 
  {
    for ( std::size_t j = 0 ; j <= i ; ++j ) 
    { gsl_matrix_set ( l , i , j , gsl_matrix_get ( aux.matrix() , i , j ) ) ; }
    for ( std::size_t j = i + 1 ; j < N ; ++j ) 
    { gsl_matrix_set ( l , i , j , 0.0 ) ; }  
  }
  //
   return Ostap::StatusCode::SUCCESS ; 
}

// ===========================================================================
/* LDLT : Cholesky decomposition of positive definite matrix 
 * \f$ PSASP^T = L D L^T\f$, 
 *  Only lower triangular part of the matrix A is used.
 *  @param A (input)  input MxM matrix
 *  @param S (output/update) scale vector/diagonal matrix 
 *  @param P (output/update) permutation 
 *  @param L (utput/update)  lower triangular matrix
 *  @param D (output/update) vector/diagonal matrix   
 *  @return status code
 */
// ===========================================================================  
Ostap::StatusCode Ostap::Math::GSL::LDLT
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Vector&       S ,
  Ostap::Math::GSL::Permutation&  P , 
  Ostap::Math::GSL::Matrix&       L , 
  Ostap::Math::GSL::Vector&       D ) 
{

  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  const std::size_t M = A.nRows () ;
  const std::size_t N = A.nCols () ;
  //
  if ( N != M ) 
  {
    gsl_error ( "LDLT: matrix is not square" , __FILE__ , __LINE__ , GSL_EBADLEN ) ;
    return MATRIX_IS_NOT_SQUARE ;
  }
  //
  if ( L.nRows() != M || L.nCols() != N ) { L.resize ( M , N , 0 ) ; }
  //
  Matrix      a { A } ;
  Permutation p { N } ;
  if ( S.size() != N ) { S.resize ( N ) ; } 
  if ( D.size() != N ) { D.resize ( N ) ; } 
  const int status = gsl_linalg_pcholesky_decomp2 ( a.matrix      () , 
                                                    p.permutation () ,  
                                                    S.vector      () ) ;
  if ( status ) 
  {
    gsl_error ( "LLT: Error from gsl_linalg_pcholesky_decomp2" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }  

  /// copy results to L and D  
  gsl_matrix* l = L.matrix () ;
  gsl_vector* d = D.vector () ;
  gsl_matrix* m = a.matrix () ; 
  for ( std::size_t i = 0 ; i < M ; ++i ) 
  {
    // lower triangle 
    for ( std::size_t j = 0 ; j < i ; ++j ) 
    { gsl_matrix_set ( l , i , j , gsl_matrix_get ( m , i , j ) ) ; }
    // unit diagonal 
    const double aii = gsl_matrix_get( m , i , i ) ;
    gsl_vector_set ( d , i , aii   ) ;
    gsl_matrix_set ( l , i , i , 1 ) ;
    // above diagonal
    for ( std::size_t j = i + 1 ; j < N ; ++j ) 
    { gsl_matrix_set ( l , i , j , 0.0 ) ; }  
  }
  //
  P.swap ( p ) ;
  // 
  return Ostap::StatusCode::SUCCESS ;
}


// ============================================================================
/*  Polar decompositon of the square matrix A: \f$ A = UP \f$
 *  - U ius orthogonal 
 *  - P is positiev semi-definitive 
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::POLAR
( const Ostap::Math::GSL::Matrix& A ,
  Ostap::Math::GSL::Matrix      & U ,
  Ostap::Math::GSL::Matrix      & P )
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  // Polar decomposition exists only for square matrices!
  if ( A.nRows() != A.nCols() ) { return MATRIX_IS_NOT_SQUARE ; }
  //
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  const std::size_t K = std::min ( M , N ) ;
  //
  Matrix auxu { M , K } ;
  Matrix auxv { N , K } ;
  Vector S    { K }     ; // vector of singular values 
  //
  Ostap::StatusCode sc = SVD ( A , S , auxu , auxv ) ;
  //
  U = auxu                * auxv . T () ;
  ///
  /// TO DO :  NEED TO MAKE IT EFFICIENT!! MAtrix(S) is a diagonal!!!
  P = auxv * Matrix ( S ) * auxv . T () ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ===============================================================================
namespace Ostap
{
  // =============================================================================
  namespace Math
  {
    // ===========================================================================
    namespace GSL
    {
      // =========================================================================
      class SchurWorkspace 
      {
      public :
        explicit SchurWorkspace 
        ( const std::size_t M )
          : m_ws ( gsl_eigen_nonsymm_alloc  ( M ) )
        {}
        ~SchurWorkspace ()  
        { if ( m_ws ) { gsl_eigen_nonsymm_free( m_ws ) ; m_ws = nullptr ;} }
      public:
        gsl_eigen_nonsymm_workspace* workspace () const { return m_ws ; }
      private:
        gsl_eigen_nonsymm_workspace* m_ws { nullptr } ;
      };
      // =======================================================================
      class ComplexVector
      {
      public:
        ComplexVector 
        ( const std::size_t N )
          : m_v ( gsl_vector_complex_alloc ( N ) )
        {}
        ~ComplexVector () 
        {
          if ( m_v ) { gsl_vector_complex_free ( m_v ) ; m_v = nullptr ; }
        } 
      public:
        gsl_vector_complex* vector() { return m_v ; }
      private:
        gsl_vector_complex* m_v { nullptr } ;
      };
      // =========================================================================
    }
    // ===========================================================================
  }
  // =============================================================================
}
// ===============================================================================
/* Schur decomposition of square matrix \f$ A = Z T Z^T\f$, where 
 *  - A is inpur MxM (square) matrix
 *  - T is Schur form of matix  
 *  - Z is orthogonam matrix 
 */
// ==============================================================================
Ostap::StatusCode Ostap::Math::GSL::SCHUR 
( const Ostap::Math::GSL::Matrix&  A ,  
  Ostap::Math::GSL::Matrix&        Z , 
  Ostap::Math::GSL::Matrix&        T ) 
{
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  // Schur decomposition exists only for square matrices!
  if ( A.nRows() != A.nCols() ) { return MATRIX_IS_NOT_SQUARE ; }
  //
  const std::size_t N { A.nRows() } ; 
  //
  T = A ;
  Z.resize ( N , N ) ;
  Ostap::Math::GSL::SchurWorkspace ws   { N } ;
  Ostap::Math::GSL::ComplexVector  eval { N } ;
  //
  gsl_eigen_nonsymm_params ( 1 , 0 , ws.workspace () ) ;
  int status = gsl_eigen_nonsymm_Z ( T.matrix     () , 
                                     eval.vector  () , 
                                     Z.matrix     () , 
                                     ws.workspace () ) ;

  if ( status )
  {
    gsl_error ( "Error from gsl_eigen_nonsymm_Z" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;
  }
  // need to clean the lower left part of T
  // gsl_vector_complex_fprintf  ( stderr , eval.vector() , "%+.4g") ; 
  //
  for  ( std::size_t j = 0 ; j < N ; ++j   )
    { for ( std::size_t i = j + 2 ; i < N ; ++i ) { T.set ( i , j , 0.0 ) ; } }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
