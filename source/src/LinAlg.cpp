// ============================================================================
// Include files 
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
  int sign = 0 ;
  gsl_linalg_LU_decomp ( A.matrix() , P.permutation() , &sign ) ;
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
//                                                                      The END 
// ============================================================================
