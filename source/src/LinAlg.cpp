// ============================================================================
// Include files 
// ============================================================================
// STD/STL
// ============================================================================
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/LinAlg.h"
#include "Ostap/StatusCode.h"
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
std::size_t Ostap::GSL::GSL_version_major () 
{ return GSL_MAJOR_VERSION ; }
// ============================================================================
// GSL version minor
// ============================================================================
std::size_t Ostap::GSL::GSL_version_minor () 
{ return GSL_MINOR_VERSION ; }
// ============================================================================
// GSL versionmajor  x 1000 + GAL version minor  
// ============================================================================
std::size_t Ostap::GSL::GSL_version_int   () 
{ return 1000 * GSL_MAJOR_VERSION + GSL_MINOR_VERSION ; }
// ============================================================================
// GSL version as string
// ============================================================================
std::string Ostap::GSL::GSL_version () { return gsl_version ; }
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const std::size_t N1 , 
  const std::size_t N2 )
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const std::size_t N1    , 
  const std::size_t N2    , 
  const double       value )
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannot use !std::isfinite"         ,
                  "Ostap::GSL::Matrix"                ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  gsl_matrix_set_all ( m_matrix , value ) ;
}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const std::size_t N1    , 
  const std::size_t N2    , 
  const Ostap::GSL::Matrix::Zero /* zero */ ) 
  : m_matrix ( gsl_matrix_calloc ( N1 , N2 ) )  // NB: calloc is here!!!
{}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const std::size_t N1    , 
  const std::size_t N2    , 
  const Ostap::GSL::Matrix::Id /* zero */ ) 
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{
  gsl_matrix_set_identity ( m_matrix ) ;
}
// ============================================================================
// allocate square GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const std::size_t N  )
  : Matrix ( N , N )
{}
// ============================================================================
// allocate square GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const std::size_t             N     , 
  const Ostap::GSL::Matrix::Zero zero  )
  : Matrix ( N , N , zero )
{}
// ============================================================================
// allocate square GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const std::size_t                 N  , 
  const Ostap::GSL::Matrix::Id  id ) 
  : Matrix ( N , N , id )
{}
// ============================================================================
// allocate square permutation GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::Matrix
( const Ostap::GSL::Permutation& p ) 
  : m_matrix ( gsl_matrix_calloc ( p.size () , p.size () ) ) 
{
  //
  Ostap::Assert ( p.valid()                                 ,
                  "(GSL)Permutation is invalid!"            , 
                  "Ostap::GSL::Matrix"                      ,
                  INVALID_PERMUTATION , __FILE__ , __LINE__ ) ;
  //
  for ( std::size_t j = 0 ; j < nRows() ; ++j )
    {
      const std::size_t k = p.get ( j ) ; 
      set ( j , k , 1 ) ;
    }
}
// ==========================================================================
Ostap::GSL::Matrix::Matrix
( const Ostap::GSL::Vector & v ) 
  : m_matrix ( gsl_matrix_calloc ( v.size () , v.size () ) ) 
{
  const std::size_t N = v.size() ;
  for ( std::size_t i = 0 ; i < N ; ++i )
    { set ( i , i , v.get ( i ) ) ; }
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::GSL::Matrix::Matrix  
( const Ostap::GSL::Matrix&  right ) 
  : m_matrix ( gsl_matrix_alloc ( right.m_matrix->size1 , 
                                  right.m_matrix->size2 ) )  
{
  gsl_matrix_memcpy ( m_matrix , right.m_matrix ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::GSL::Matrix::Matrix  
(       Ostap::GSL::Matrix&&  right ) 
  : m_matrix ( right.m_matrix )  
{
  right.m_matrix = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-matrix 
// ============================================================================
Ostap::GSL::Matrix::~Matrix () 
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
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::operator=
( const Ostap::GSL::Matrix&  right ) 
{
  if ( &right == this ) { return *this ; }
  resize ( right.nRows() , right.nCols() ) ;
  gsl_matrix_memcpy ( m_matrix , right.m_matrix ) ;
  return *this ;
}
// ============================================================================
// move assignement! 
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::operator=
( Ostap::GSL::Matrix&& right ) 
{
  if ( &right == this ) { return *this ; }
  std::swap ( m_matrix , right.m_matrix ) ;
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::resize
( const std::size_t n1 ,
  const std::size_t n2 )
{
  if ( n1 != m_matrix->size1 || n2 != m_matrix->size2 )
    {
      gsl_matrix_free ( m_matrix ) ; 
      m_matrix = gsl_matrix_alloc ( n1 , n2 ) ;
    }
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::resize
( const std::size_t n1    ,
  const std::size_t n2    ,
  const double      value )
{
  if ( !value ) { return resize ( n1 , n2 , Zero() ) ; }
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannot use !std::isfinite"         ,
                  "Ostap::GSL::Matrix"                ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  resize ( n1 , n2 ) ;
  gsl_matrix_set_all ( m_matrix , value ) ;
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::resize
( const std::size_t n1    ,
  const std::size_t n2    ,
  const Ostap::GSL::Matrix::Zero /* zero */ ) 
{
  if ( n1 != m_matrix->size1 || n2 != m_matrix->size2 )
    {
      gsl_matrix_free ( m_matrix ) ; 
      m_matrix = gsl_matrix_calloc ( n1 , n2 ) ;      // CALLOC! 
    }
  else
    { gsl_matrix_set_all ( m_matrix , 0 ) ; }
  return *this ;
}
// ============================================================================
// resize/reset matrx
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::resize
( const std::size_t n1    ,
  const std::size_t n2    ,
  const Ostap::GSL::Matrix::Id /* zero */ ) 
{
  resize ( n1 , n2 ) ;
  gsl_matrix_set_identity ( m_matrix ) ;
  return *this ;
}
// ============================================================================
// swap two matrices
// ============================================================================
void Ostap::GSL::Matrix::swap
( Ostap::GSL::Matrix& right )
{ std::swap ( m_matrix , right.m_matrix ) ; }
// ============================================================================
// swap two rows in matrix 
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::swap_rows
( const std::size_t i1 ,
  const std::size_t i2 )
{
  if ( i1 == i2 ) { return *this ; }
  Ostap::Assert ( i1 < nRows () && i2 <= nRows () , 
                  "Invalid row index!"            ,
                  "Ostap::GSL::Martix::swap_rows" ,
                  INVALID_ROWINDEX , __FILE__ , __LINE__ ) ;
  gsl_matrix_swap_rows    ( m_matrix , i1 , i2 ) ;
  return *this ; 
}
// ============================================================================
// swap two columns in matrix 
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::swap_cols
( const std::size_t i1 ,
  const std::size_t i2 )
{
  if ( i1 == i2 ) { return *this ; }
  Ostap::Assert ( i1 < nCols () && i2 <= nCols ()  , 
                  "Invalid coumn index!"          ,
                  "Ostap::GSL::Martix::swap_cols" ,
                  INVALID_COLINDEX , __FILE__ , __LINE__ ) ;
  gsl_matrix_swap_columns ( m_matrix , i1 , i2 ) ;
  return *this ; 
}
// ============================================================================
// permute the rows     of the ematrix according to permutation 
// ============================================================================
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::permute_rows
( const Ostap::GSL::Permutation& p )
{
  Ostap::Assert ( nRows () == p.size ()                     ,
                  "Inconsistent Permutation structure"      , 
                  "Ostap::GSL::Matrix::permute_rows"        ,
                  INVALID_PERMUTATION , __FILE__ , __LINE__ ) ;
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
Ostap::GSL::Matrix&
Ostap::GSL::Matrix::permute_cols
( const Ostap::GSL::Permutation& p )
{
  Ostap::Assert ( nCols () == p.size ()                     ,
                  "Inconsistent Permutation structure"      , 
                  "Ostap::GSL::Matrix::permute_cols"        ,
                  INVALID_PERMUTATION , __FILE__ , __LINE__ ) ;
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
bool Ostap::GSL::Matrix::iszero   () const
{
  for ( std::size_t i = 0 ; i < nRows() ; ++i )
    { for ( std::size_t j = 0 ; j < nCols() ; ++j )
        { if ( !s_zero ( get ( i ,j ) ) ) { return false ; } } }
  return true ;
}
// ============================================================================
// Are all elements finite ? 
// ============================================================================
bool Ostap::GSL::Matrix::isfinite () const
{
  for ( std::size_t i = 0 ; i < nRows() ; ++i )
    { for ( std::size_t j = 0 ; j < nCols() ; ++j )
        { if ( !std::isfinite ( get ( i ,j ) ) ) { return false ; } } }
  return true ;  
}
// ============================================================================
// scale matrix 
// ============================================================================
Ostap::GSL::Matrix&                             
Ostap::GSL::Matrix::imul
( const double value )
{
  if ( 1 == value ) { return *this ; }
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannto add !std::isfinite"         ,
                  "Ostap::GSL::Matrix::iadd"          ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannto scale but !std::isfinite"   ,
                  "Ostap::GSL::Matrix::imul"          ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  //
  gsl_matrix_scale ( m_matrix , value ) ;
  return *this;
}
// ============================================================================
// add&subtract matrix 
// ============================================================================
Ostap::GSL::Matrix&                             
Ostap::GSL::Matrix::iadd
( const Ostap::GSL::Matrix& right )
{
  Ostap::Assert ( this->nRows () == right.nRows() &&
                  this->nCols () == right.nCols() ,
                  "Cannot add matrix of incompatible structure" ,
                  "Ostap::GSL::Martix::iadd" ,
                  INVALID_GMATRIX , __FILE__ , __LINE__ ) ;
  //
  gsl_matrix_add ( m_matrix , right.m_matrix ) ;
  return *this;
}
// ============================================================================
// add&subtract matrix 
// ============================================================================
Ostap::GSL::Matrix&                             
Ostap::GSL::Matrix::isub
( const Ostap::GSL::Matrix& right )
{
  Ostap::Assert ( this->nRows () == right.nRows() &&
                  this->nCols () == right.nCols() ,
                  "Cannot sub matrix of incompatible structure" ,
                  "Ostap::GSL::Martix::iasub" ,
                  INVALID_GMATRIX , __FILE__ , __LINE__ ) ;
  //
  gsl_matrix_sub ( m_matrix , right.m_matrix ) ;
  return *this;
}
// ============================================================================
// add identity matrix 
// ============================================================================
Ostap::GSL::Matrix&                             
Ostap::GSL::Matrix::iadd
( const double value )
{
  if ( !value ) { return *this ; }
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannto add  !std::isfinite"        ,
                  "Ostap::GSL::Matrix::iadd"          ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  //
  const std::size_t N = std::min ( nRows () , nCols() );
  for ( std::size_t i = 0 ; i < N ; ++i )
    { set ( i , i , get ( i , i ) + value ) ; }
  return *this;
}
// ============================================================================
// multiply matrices  using CBLAS dgemm function 
// ============================================================================
Ostap::GSL::Matrix                            
Ostap::GSL::Matrix::multiply
( const Ostap::GSL::Matrix& right ) const
{
  Ostap::Assert ( this->nCols () == right.nRows() ,
                  "Cannot multiply matrices of incompatible structure" ,
                  "Ostap::GSL::Martix::multiply"  ,
                  INVALID_GMATRIX , __FILE__ , __LINE__ ) ;
  
  Matrix result { nRows() , right.nCols() } ;
  // 
  gsl_blas_dgemm ( CblasNoTrans    ,
                   CblasNoTrans    ,
                   1.0             ,
                   this ->matrix() ,
                   right .matrix() ,
                   0.0             ,
                   result.matrix() );
  //
  return result ;
}
// ============================================================================
// multiply matrices  using CBLAS dgemm function 
// ============================================================================
Ostap::GSL::Matrix&                           
Ostap::GSL::Matrix::imul
( const Ostap::GSL::Matrix& right )
{
  Ostap::Assert ( this->nCols () == right.nRows () ,
                  "Cannot multiply matrices of incompatible structure" ,
                  "Ostap::GSL::Martix::multiply"  ,
                  INVALID_GMATRIX , __FILE__ , __LINE__ ) ;
  
  Matrix result { nRows() , right.nCols() } ;
  // 
  gsl_blas_dgemm ( CblasNoTrans    ,
                   CblasNoTrans    ,
                   1.0             ,
                   this ->matrix() ,
                   right .matrix() ,
                   0.0             ,
                   result.matrix() );
  //
  this->swap ( result ) ;
  //
  return *this ;
}
// ============================================================================
// multiply matrix abd vector using CBLAS dgemv function 
// ============================================================================
Ostap::GSL::Vector
Ostap::GSL::Matrix::multiply
( const Ostap::GSL::Vector& right ) const
{
  Ostap::Assert ( this->nCols () == right.size() ,
                  "Cannot multiply matrix&vector of incompatible structure" ,
                  "Ostap::GSL::Martix::multiply"  ,
                  INVALID_GMATRIX , __FILE__ , __LINE__ ) ;
  //
  Vector result { nRows() } ;
  //
  gsl_blas_dgemv ( CblasNoTrans     ,
                   1.0              ,
                   this ->matrix () ,
                   right .vector () ,
                   0.0              ,
                   result.vector () ) ;
  //
  return result ;
}
// ============================================================================
// multiply matrix abd vector using CBLAS dgemv function 
// ============================================================================
Ostap::GSL::Matrix
Ostap::GSL::Matrix::multiply
( const Ostap::GSL::Permutation& right ) const
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
Ostap::GSL::Matrix
Ostap::GSL::Matrix::T () const
{
  // preare result 
  Matrix result { nCols() , nRows() } ;
  gsl_matrix_transpose_memcpy ( result.matrix() , matrix () ) ;
  return result; 
}
// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::GSL::Vector::Vector
( const std::size_t N  ) 
  : m_vector ( gsl_vector_alloc ( N ) )
{}
// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::GSL::Vector::Vector
( const std::size_t  N     ,   
  const double        value )
  : m_vector ( gsl_vector_alloc ( N ) )
{
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannot use !std::isfinite"         ,
                  "Ostap::GSL::Vector"                ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  gsl_vector_set_all ( m_vector , value ) ;
}
// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::GSL::Vector::Vector
( const std::size_t                N       ,   
  const Ostap::GSL::Vector::Zero /* zero */ ) 
  : m_vector ( gsl_vector_calloc ( N ) )  //    NB! calloc here! 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::GSL::Vector::Vector  
( const Ostap::GSL::Vector&  right ) 
  : m_vector ( gsl_vector_alloc ( right.m_vector->size ) )  
{
  gsl_vector_memcpy ( m_vector , right.m_vector ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::GSL::Vector::Vector  
(       Ostap::GSL::Vector&&  right ) 
  : m_vector ( right.m_vector )  
{
  right.m_vector = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-vector
// ============================================================================
Ostap::GSL::Vector::~Vector () 
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
Ostap::GSL::Vector&
Ostap::GSL::Vector::operator=
( const Ostap::GSL::Vector&  right ) 
{
  if ( &right == this ) { return *this ; }
  resize ( right.m_vector->size ) ;
  gsl_vector_memcpy ( m_vector , right.m_vector ) ;
  return *this ;
}
// ============================================================================
// move assignement! 
// ============================================================================
Ostap::GSL::Vector&
Ostap::GSL::Vector::operator=
( Ostap::GSL::Vector&& right ) 
{
  if ( &right == this ) { return *this ; }
  std::swap ( m_vector , right.m_vector ) ;
  return *this ;
}
// ============================================================================
// resize the vector
// ============================================================================
Ostap::GSL::Vector&
Ostap::GSL::Vector::resize
( const std::size_t n )
{
  if ( n != m_vector->size )
    {
      gsl_vector_free ( m_vector ) ; 
      m_vector = gsl_vector_alloc ( n ) ;
    }
  return *this ;
}
// ============================================================================
// resize the vector
// ============================================================================
Ostap::GSL::Vector&
Ostap::GSL::Vector::resize
( const std::size_t n     ,
  const double       value ) 
{
  if ( !value ) { return resize ( n , Zero() ) ; }
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannot use !std::isfinite"         ,
                  "Ostap::GSL::Vector"                ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  resize ( n ) ;
  gsl_vector_set_all ( m_vector , value ) ;
  return *this ;
}
// ============================================================================
// resize the vector
// ============================================================================
Ostap::GSL::Vector&
Ostap::GSL::Vector::resize
( const std::size_t n     ,
  const Ostap::GSL::Vector::Zero /* zero */ )  
{
  if ( n != m_vector->size )
    {
      gsl_vector_free ( m_vector ) ; 
      m_vector = gsl_vector_calloc ( n ) ; // CALLOC HERE 
    }
  else { gsl_vector_set_all ( m_vector , 0 ) ; }
  //
  return *this ;
}
// ============================================================================
// swap two vectors 
// ============================================================================
void Ostap::GSL::Vector::swap
( Ostap::GSL::Vector& right )
{ std::swap ( m_vector , right.m_vector ) ; }
// ============================================================================
// Are all elements numerically equal to zero?      
// ============================================================================
bool Ostap::GSL::Vector::iszero   () const
{
  for ( std::size_t i = 0 ; i < size ()  ; ++i )
    { if ( !s_zero ( get ( i ) ) ) { return false ; } }
  return true ;
}
// ============================================================================
// Are all elements finite ? 
// ============================================================================
bool Ostap::GSL::Vector::isfinite () const
{
  for ( std::size_t i = 0 ; i < size ()  ; ++i )
    { if ( !std::isfinite ( get ( i ) ) ) { return false ; } }
  return true ;
}
// ============================================================================
// dot product of two vect
// ============================================================================
double
Ostap::GSL::Vector::dot
( const Ostap::GSL::Vector& value ) const 
{
  Ostap::Assert ( this->size() == value.size () ,
                  "Cannot dot vectors of incompatible structure" ,
                  "Ostap::GSL::Vector::dot"  ,
                  INVALID_GVECTOR, __FILE__ , __LINE__ ) ;
  //
  double result = 0 ;
  gsl_blas_ddot ( m_vector , value.m_vector , &result ) ;
  return result ;
}
// ============================================================================
// cross product of two vect
// ============================================================================
Ostap::GSL::Matrix
Ostap::GSL::Vector::cross
( const Ostap::GSL::Vector& value ) const 
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
Ostap::GSL::Vector&
Ostap::GSL::Vector::iadd
( const Ostap::GSL::Vector& value )
{
  Ostap::Assert ( size()  == value.size ()   ,
                  "Cannot add vectors of incompatible structure" ,
                  "Ostap::GSL::Vector::iadd" ,
                  INVALID_GVECTOR, __FILE__  , __LINE__ ) ;
  gsl_vector_add ( m_vector , value.m_vector ) ;
  return *this ;
}
// ============================================================================
// add a cboisrant 
// ============================================================================
Ostap::GSL::Vector&
Ostap::GSL::Vector::iadd
( const double value )
{
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannot use !std::isfinite"         ,
                  "Ostap::GSL::Vector"                ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  gsl_vector_add_constant  ( m_vector , value ) ;
  return *this ;
}
// ============================================================================
// subtract vector 
// ============================================================================
Ostap::GSL::Vector&
Ostap::GSL::Vector::isub
( const Ostap::GSL::Vector& value )
{
  Ostap::Assert ( this->size()  == value.size () ,
                  "Cannot subtract vectors of incompatible structure" ,
                  "Ostap::GSL::Vector::isub"  ,
                  INVALID_GVECTOR, __FILE__ , __LINE__ ) ;
  gsl_vector_sub ( m_vector , value.m_vector ) ;
  return *this ;
}


// ============================================================================
// scale vector  
// ============================================================================
Ostap::GSL::Vector&                             
Ostap::GSL::Vector::imul
( const double value )
{
  if ( 1 == value ) { return *this ; }
  Ostap::Assert ( std::isfinite ( value )             ,
                  "Cannto scale but !std::isfinite"   ,
                  "Ostap::GSL::Vector::imul"          ,
                  INVALID_SCALE , __FILE__ , __LINE__ ) ;
  //
  gsl_vector_scale ( m_vector , value ) ;
  return *this;
}
// ============================================================================
// multiply by matrix 
// ============================================================================
Ostap::GSL::Vector&
Ostap::GSL::Vector::imul 
( const Ostap::GSL::Matrix& value )
{
  Ostap::Assert ( size()  == value.nRows() ,
                  "Cannot multiply vector&matrix of incompatible structure" ,
                  "Ostap::GSL::Vector::isub"  ,
                  INVALID_GMATRIX, __FILE__ , __LINE__ ) ;
  
  Vector result { value.nCols() } ;
  //
  gsl_blas_dgemv ( CblasTrans       , // ATTENTIO!!! 
                   1.0              ,
                   value .matrix () ,
                   m_vector         ,
                   0.0              ,
                   result.vector () ) ;
  //
  this->swap ( result ) ;
  //
  return *this ;
}

// ============================================================================
// multiply by matrix 
// ============================================================================
Ostap::GSL::Vector
Ostap::GSL::Vector::multiply  
( const Ostap::GSL::Matrix& value ) const
{
  Ostap::Assert ( this->size()  == value.nRows() ,
                  "Cannot multiply vector&matrix of incompatible structure" ,
                  "Ostap::GSL::Vector::isub"  ,
                  INVALID_GMATRIX, __FILE__ , __LINE__ ) ;
  
  Vector result { value.nCols() } ;
  //
  gsl_blas_dgemv ( CblasTrans       , // ATTENTIO!!! 
                   1.0              ,
                   value .matrix () ,
                   m_vector         ,
                   0.0              ,
                   result.vector () ) ;
  //
  return result ;
}


// ============================================================================
// constructor: allocate the permutation 
// ============================================================================
Ostap::GSL::Permutation::Permutation
( const std::size_t N ) 
  : m_permutation ( gsl_permutation_calloc ( N ) ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::GSL::Permutation::Permutation 
( const Ostap::GSL::Permutation&  right ) 
  : m_permutation  ( gsl_permutation_alloc ( right.m_permutation->size ) )
{
  gsl_permutation_memcpy ( m_permutation , right.m_permutation ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::GSL::Permutation::Permutation 
( Ostap::GSL::Permutation&& right ) 
  : m_permutation ( right.m_permutation )  
{
  right.m_permutation = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-permutation
// ============================================================================
Ostap::GSL::Permutation::~Permutation () 
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
Ostap::GSL::Permutation&
Ostap::GSL::Permutation::operator=
( const Ostap::GSL::Permutation&  right ) 
{
  if ( &right == this ) { return *this ; }
  if ( m_permutation->size != right.m_permutation->size )
    {
      gsl_permutation_free ( m_permutation ) ;
      m_permutation = gsl_permutation_alloc ( right.m_permutation->size ) ;
    }
  gsl_permutation_memcpy ( m_permutation , right.m_permutation ) ;
  return *this ;
}
// ============================================================================
// move assignement! 
// ============================================================================
Ostap::GSL::Permutation&
Ostap::GSL::Permutation::operator=
( Ostap::GSL::Permutation&& right ) 
{
  if ( &right == this ) { return *this ; }
  std::swap ( m_permutation , right.m_permutation ) ;
  return *this ;
}
// ============================================================================
// swap two permutation 
// ============================================================================
void Ostap::GSL::Permutation::swap
( Ostap::GSL::Permutation& right )
{ std::swap ( m_permutation , right.m_permutation ) ; }
// ============================================================================
// valid permutation ?
// ============================================================================
bool
Ostap::GSL::Permutation::valid () const
{ return GSL_SUCCESS == gsl_permutation_valid ( m_permutation ) ; }

// ============================================================================
// apply permutation to the matrix 
// ============================================================================
Ostap::GSL::Matrix
Ostap::GSL::Permutation::apply 
( const Ostap::GSL::Matrix& value ) const 
{
  Ostap::Assert ( size() == value.nRows()                      ,
                  "Mismatch for permutation/matrix structure!" ,
                  "Ostap::GLS::Permutaton::apply"              , 
                  INVALID_PERMUTATION  , __FILE__ , __LINE__   ) ;
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
( const Ostap::GSL::Matrix& m ,
  std::ostream&              s )
{ return toStream ( *m.matrix() , s ) ; }
// ============================================================================
// print vector to the stream 
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const Ostap::GSL::Vector& v ,
  std::ostream&              s )
{ return toStream ( *v.vector () , s ) ; }
// ============================================================================
// print permutation to the stream 
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const Ostap::GSL::Permutation& p ,
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

// ============================================================================
// get the element with maxina absolute value 
// ============================================================================
double Ostap::Math::maxabs_element
( const Ostap::GSL::Matrix& m )
{
  double result = -1 ;
  for ( std::size_t i = 0 ; i < m.nRows () ; ++i )
    { for ( std::size_t j = 0 ; j < m.nCols () ; ++j )
        { result = std::max ( result , std::abs ( m ( i , j ) ) ) ; } }
  return result ;
}
// ============================================================================
// get the element with maxina absolute value 
// ============================================================================
double Ostap::Math::maxabs_element
( const Ostap::GSL::Vector& v )
{
  double result = -1 ;
  for ( std::size_t i = 0 ; i < v.size() ; ++i )
    { result = std::max ( result , std::abs ( v ( i ) ) ) ; }
  return result ;
}
// ============================================================================
// get the element with maxina absolute value 
// ============================================================================
double Ostap::Math::maxabs_element
( const Ostap::GSL::Permutation& v ) { return v.size() ; }
// ============================================================================
// Actual Linar Algebra start here 
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
Ostap::GSL::Permutation
Ostap::GSL::PLU ( Ostap::GSL::Matrix& A )
{
  Permutation P { A.nRows() } ;
  int signum = 0 ;
  int status = gsl_linalg_LU_decomp ( A.matrix() , P.permutation() , &signum ) ;
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_LU_decomp"        ,
		  "Ostap::GSL::PLU"                        ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;  
  //
  return P ;
}
// ============================================================================
/*  perfom LU decomposition  
 *  @param  A   (INOUT)         input matrix 
 *  @param  LU  (UPDATE/OUTPUT) output LU matrix 
 *  @return M-permutations 
 *  @see gsl_linalg_LU_decomp 
 */
// ============================================================================
Ostap::GSL::Permutation
Ostap::GSL::PLU
( const Ostap::GSL::Matrix& A  ,
  Ostap::GSL::Matrix&       LU )
{
  LU = A ;
  return PLU ( LU ) ;
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
Ostap::GSL::Permutation
Ostap::GSL::PLU
( const Ostap::GSL::Matrix& A ,
  Ostap::GSL::Matrix&       L ,
  Ostap::GSL::Matrix&       U )
{
  //
  Ostap::GSL::Matrix      LU { A } ;
  Ostap::GSL::Permutation P  { PLU ( LU ) } ;
  //
  const std::size_t M { A.nRows () } ;
  const std::size_t N { A.nCols () } ;
  const std::size_t K { std::min ( M , N ) } ; 
  //
  L.resize ( M , K , Ostap::GSL::Matrix::Id   () ) ; 
  U.resize ( K , N , Ostap::GSL::Matrix::Zero () ) ; 
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
  return P ;
}
// ============================================================================


// ============================================================================
// QR decomposition with column pivoting 
// ============================================================================
/*  mape QR Decomposion of matrix A : \f$ AP = QR\f$ where 
 *  - A is input                 MxN matrix  
 *  - P is permuutation matrix   NxN 
 *  - Q is orthogonal matrix     MxM 
 *  - R is right triaular matrix MxN 
 *  
 *  @param A  (input) the matrix to decopose 
 *  @param Q  (outpt/update) orthogonal matrix Q 
 *  @param R  (outpt/update) rigth triangular matrix R 
 *  @return permutation P 
 */
// ============================================================================
Ostap::GSL::Permutation
Ostap::GSL::PQR
( const Ostap::GSL::Matrix& A ,
  Ostap::GSL::Matrix&       Q ,
  Ostap::GSL::Matrix&       R )
{
  const std::size_t M = A.nRows() ;
  const std::size_t N = A.nCols() ;
  const std::size_t K = std::min ( M , N ) ; 
  //
  //
  Q.resize ( M , M ) ;
  R.resize ( M , N , Ostap::GSL::Matrix::Zero() ) ;
  //
  Permutation P    { N } ;
  Vector      tau  { K } ;
  Vector      norm { N } ;
  //
  int signum = 0 ;
  int status = gsl_linalg_QRPT_decomp2
    ( A.matrix()   , Q.matrix() , R.matrix() ,
      tau.vector() , P.permutation() , &signum , norm.vector() ) ;
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_QRPT_decomp2"     ,
		  "Ostap::GSL::PQR"                        ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  
  //
  return P ;
}
// ============================================================================

// ============================================================================
// LQ decomposition
// ============================================================================
/* LQ decomposition of matrix A: \f$ A = LQ\f$, where 
 *  - L is lower trapezoidal MxN 
 *  - Q is orthogonal NxN 
 */ 
// ============================================================================
void Ostap::GSL::LQ
( const Ostap::GSL::Matrix& A ,
  Ostap::GSL::Matrix&       L ,
  Ostap::GSL::Matrix&       Q )
{
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  const std::size_t K = std::min ( M , N ) ;
  //
  L.resize ( M , N , Ostap::GSL::Matrix::Zero() ) ;
  Q.resize ( N , N ) ;
  //
  Vector tau { K } ;
  Matrix R   { A } ;
  //
  int status = gsl_linalg_LQ_decomp ( R.matrix() , tau.vector() ) ;
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_LQ_decomp"        ,
		  "Ostap::GSL::LQ"                         ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  //
  status = gsl_linalg_LQ_unpack ( R.matrix   () ,
				  tau.vector () ,
				  Q.matrix   () , 
				  L.matrix   () ) ; 
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_LQ_unpack"        ,
		  "Ostap::GSL::LQ"                         ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  //
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

void Ostap::GSL::QL
( const Ostap::GSL::Matrix& A ,
  Ostap::GSL::Matrix&       Q ,
  Ostap::GSL::Matrix&       L )
{
#if OSTAP_GSL_CODE_VERSION < OSTAP_GSL_VERSION ( 2 , 7 )
  Ostap::Assert ( false  , 
    "For A=QL decomposition need GSL version>2.6", 
    "Ostap::GSL::QL"      , 
    ERROR_GSL    , __FILE__ , __LINE__ ) ;
#else 
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  //
  Q.resize ( M , M ) ;
  L.resize ( M , N , Ostap::GSL::Matrix::Zero() ) ;
  //
  Vector tau { N } ;
  Matrix R   { A } ;
  //
  int status = gsl_linalg_QL_decomp ( R.matrix() , tau.vector() ) ;
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_QL_decomp"        ,
		  "Ostap::GSL::QL"                         ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  //
  status = gsl_linalg_QL_unpack ( R.matrix   () ,
				  tau.vector () ,
				  Q.matrix   () , 
				  L.matrix   () ) ; 
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_QL_unpack"        ,
		  "Ostap::GSL::QL"                         ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
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
Ostap::GSL::Permutation
Ostap::GSL::COD
( const Ostap::GSL::Matrix& A ,
  Ostap::GSL::Matrix& Q ,
  Ostap::GSL::Matrix& R ,
  Ostap::GSL::Matrix& Z )
{
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  const std::size_t K = std::min ( M , N ) ;
  //
  Permutation P { N } ;
  //
  Q.resize ( M , M ) ;
  R.resize ( M , N , Ostap::GSL::Matrix::Zero() ) ;
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
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_COD_decomp"       ,
		  "Ostap::GSL::COD"                        ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  //
  status = gsl_linalg_COD_unpack
    ( D.matrix() ,
      tau_Q.vector () ,
      tau_Z.vector () ,
      rank            ,
      Q.matrix     () ,
      R.matrix     () ,
      Z.matrix     () ) ;
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_COD_unpack"       ,
		  "Ostap::GSL::COD"                        ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  //
  return P ;
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
Ostap::GSL::Vector
Ostap::GSL::SVD
( const Ostap::GSL::Matrix& A     ,
  Ostap::GSL::Matrix&       U     ,
  Ostap::GSL::Matrix&       V     ,
  const bool                golub )
{
  //
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  //
  if ( M < N ) { return SVD ( A.T() , V , U ) ; }
  //
  // replace U with A 
  U = A ;
  V.resize ( N , N ) ;
  //
  Vector S    { N } ;
  //
  /// Use one sided JAcobi orthogonalization 
  if ( !golub )
    {
      int status = gsl_linalg_SV_decomp_jacobi 
	( U.matrix () ,
	  V.matrix () ,
	  S.vector () ) ; 
      Ostap::Assert ( GSL_SUCCESS == status                    ,
		      "Error from gsl_linalg_SV_decomp_jacobi" ,
		      "Ostap::GSL::SVD"                        ,
		      ERROR_GSL + status , __FILE__ , __LINE__ ) ;
      
      return S ;  // RETURN 
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
      Ostap::Assert ( GSL_SUCCESS == status                    ,
		      "Error from gsl_linalg_SV_decomp"        ,
		      "Ostap::GSL::SVD"                        ,
		      ERROR_GSL + status , __FILE__ , __LINE__ ) ;
      return S ;  
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
  Ostap::Assert ( GSL_SUCCESS == status                    ,
		  "Error from gsl_linalg_SV_decomp_mod"    ,
		  "Ostap::GSL::SVD"                        ,
		  ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  return S ; 
}
// ============================================================================
/*  Polar decompositon of the square matrix A: \f$ A = UP \f$
 *  - U ius orthogonal 
 *  - P is positiev semi-definitive 
 */
// ============================================================================
void Ostap::GSL::POLAR
( const Ostap::GSL::Matrix& A ,
  Ostap::GSL::Matrix      & U ,
  Ostap::GSL::Matrix      & P )
{
  //
  const std::size_t M = A.nRows  () ;
  const std::size_t N = A.nCols  () ;
  const std::size_t K = std::min ( M , N ) ;
  //
  Ostap::Assert ( M == N ,
		  "Polar decomposition exists only for square matrices!" ,
		  "Ostap::GSL::POLAR" ,
		  INVALID_GMATRIX     , __FILE__ , __LINE__ ) ;
  //
  Matrix auxu { M  , K } ;
  Matrix auxv { N  , K } ;
  //
  Matrix S    { SVD ( A , auxu , auxv ) }  ;
  //
  U = auxu     * auxv.T() ;
  P = auxv * S * auxv.T() ;
}


namespace Ostap
{
  namespace GSL
    {
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
      // ========================================================================
    }
}
// ===============================================================================
/* Schur decomposition of square matrix \f$ A = Z T Z^T\f$, where 
  *  - A is inpur MxM (square) matrix
  *  - T is Schur form of matix  
  *  - Z is orthogonam matrix 
  */
 // ==============================================================================
 #include <iostream>

void Ostap::GSL::SCHUR 
( const Ostap::GSL::Matrix&  A ,  
  Ostap::GSL::Matrix&        Z , 
  Ostap::GSL::Matrix&        T ) 
  {
    Ostap::Assert( A.nRows() == A.nCols() , 
                  "Schur decomposiiton is only for square matrices!" ,
                  "Ostap::GSL::SCHUR" , 
                  INVALID_GMATRIX , __FILE__  , __LINE__ ) ;

    const std::size_t N { A.nRows() } ; 

    T = A ;
    Z.resize ( N , N ) ;
    Ostap::GSL::SchurWorkspace ws   { N } ;
    Ostap::GSL::ComplexVector  eval { N } ;
    //
    gsl_eigen_nonsymm_params ( 1 , 0 , ws.workspace () ) ;
    int status = gsl_eigen_nonsymm_Z ( T.matrix     () , 
                                       eval.vector  () , 
                                       Z.matrix     () , 
                                       ws.workspace () ) ; 
    Ostap::Assert ( GSL_SUCCESS == status                     , 
                    "Error from gsl_eigwn_nonsymm_Z"          , 
                    "Ostap::GSL::SCHUR"                       , 
                     ERROR_GSL + status , __FILE__ , __LINE__ ) ;
  // need to clean the lower left part of T
  // ,, 
  gsl_vector_complex_fprintf  ( stderr , eval.vector() , "%.4g") ; 
  //
  bool prev_cmplx = false ; 
  for  ( std::size_t i = 0 ; i < N ; ++i   )
  {
    std:size_t k = N - i - 1 ; 
    const double re = GSL_REAL ( gsl_vector_complex_get ( eval.vector() , k ) ) ;
    const double im = GSL_IMAG ( gsl_vector_complex_get ( eval.vector() , k ) ) ; 
    //
    std::cout << " i " 
              << " re: " << re 
              << " im: " << im 
              << " T:  " << T  ( i , i ) << std::endl ;

    if ( !im  || prev_cmplx ) 
    {
      for ( std::size_t j = i + 1 ; j < N ; ++ j ) { T.set ( j , i , 0.0 ) ; }
      prev_cmplx = false ; 
    }
    else
    {
      for ( std::size_t j = i + 2 ; j < N ; ++ j ) { T.set ( j , i , 0.0 ) ; }
      prev_cmplx = true  ; 
    }
  }
}



// ===x=========================================================================
//                                                                      The END 
// ============================================================================
