// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/GSL_utils.h"
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_linalg.h"
// ============================================================================
// local
// ============================================================================
#include "GSL_helpers.h"
// ============================================================================
/** @file 
 *  Implementation file for helper  GSL classes 
 *  @date 2020-09-22 
 *  @author Vanya BELYAEV IvanBelyaev@iter.ru
 */
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL_Matrix::GSL_Matrix
( const unsigned short N1 , 
  const unsigned short N2 )
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL_Matrix::GSL_Matrix
( const unsigned short N1    , 
  const unsigned short N2    , 
  const double         value )
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{
  gsl_matrix_set_all ( m_matrix , value ) ;
}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL_Matrix::GSL_Matrix
( const unsigned short N1    , 
  const unsigned short N2    , 
  const Ostap::GSL_Matrix::Zero /* zero */ ) 
  : m_matrix ( gsl_matrix_calloc ( N1 , N2 ) )  // NB: calloc is here!!!
{}
// ============================================================================
// allocate GSL-matrix 
// ============================================================================
Ostap::GSL_Matrix::GSL_Matrix
( const unsigned short N1    , 
  const unsigned short N2    , 
  const Ostap::GSL_Matrix::Identity /* zero */ ) 
  : m_matrix ( gsl_matrix_alloc ( N1 , N2 ) )
{
  gsl_matrix_set_identity ( m_matrix ) ;
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::GSL_Matrix::GSL_Matrix  
( const Ostap::GSL_Matrix&  right ) 
  : m_matrix ( gsl_matrix_alloc ( right.m_matrix->size1 , 
                                  right.m_matrix->size2 ) )  
{
  gsl_matrix_memcpy ( m_matrix , right.m_matrix ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::GSL_Matrix::GSL_Matrix  
(       Ostap::GSL_Matrix&&  right ) 
  : m_matrix ( right.m_matrix )  
{
  right.m_matrix = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-matrix 
// ============================================================================
Ostap::GSL_Matrix::~GSL_Matrix () 
{
  if ( nullptr != m_matrix ) 
  {
    gsl_matrix_free ( m_matrix ) ; 
    m_matrix = nullptr ; 
  }
}

  



// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::GSL_Vector::GSL_Vector
( const unsigned short N  ) 
  : m_vector ( gsl_vector_alloc ( N ) )
{}
// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::GSL_Vector::GSL_Vector
( const unsigned short N     ,   
  const double         value )
  : m_vector ( gsl_vector_alloc ( N ) )
{
  gsl_vector_set_all ( m_vector , value ) ;
}
// ============================================================================
// allocate GSL-Vector 
// ============================================================================
Ostap::GSL_Vector::GSL_Vector
( const unsigned short N     ,   
  const Ostap::GSL_Vector::Zero /* zero */ ) 
  : m_vector ( gsl_vector_calloc ( N ) )  //    NB! calloc here! 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::GSL_Vector::GSL_Vector  
( const Ostap::GSL_Vector&  right ) 
  : m_vector ( gsl_vector_alloc ( right.m_vector->size ) )  
{
  gsl_vector_memcpy ( m_vector , right.m_vector ) ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::GSL_Vector::GSL_Vector  
(       Ostap::GSL_Vector&&  right ) 
  : m_vector ( right.m_vector )  
{
  right.m_vector = nullptr ;
}
// ============================================================================
///  destructor: free  GSL-vector
// ============================================================================
Ostap::GSL_Vector::~GSL_Vector () 
{
  if ( nullptr != m_vector ) 
  {
    gsl_vector_free ( m_vector ) ; 
    m_vector = nullptr ; 
  }
}
  
// ============================================================================
// constructor: allocate the permutation 
// ============================================================================
Ostap::GSL_Permutation::GSL_Permutation
 ( const unsigned short N ) 
   : m_permutation ( gsl_permutation_alloc ( N ) ) 
{}
// ============================================================================
///  destructor: free  GSL-permutation
// ============================================================================
Ostap::GSL_Permutation::~GSL_Permutation () 
{
  if ( nullptr != m_permutation ) 
  {
    gsl_permutation_free ( m_permutation ) ; 
    m_permutation = nullptr ; 
  }
}
// ============================================================================
/// print operator 
std::ostream& operator<<( std::ostream&            s ,
			  const Ostap::GSL_Matrix& m )
{
  if ( m.matrix() == nullptr ) { s << "GSL_Matrix{nullptr}" ; return s ; }
  return Ostap::Utils::toStream ( *m.matrix() , s ) ;
}
// ============================================================================
/// print operator 
std::ostream& operator<<( std::ostream&            s ,
			  const Ostap::GSL_Vector& v ) 
{
  if ( v.vector() == nullptr ) { s << "GSL_Vector{nullptr}" ; return s ; }
  return Ostap::Utils::toStream ( *v.vector() , s ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
