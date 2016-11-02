// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <vector>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_sort_vector.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/Exception.h"
#include "Ostap/EigenSystem.h"
#include "Ostap/EigenSystem.icpp"
// ============================================================================

// ============================================================================
/** @file 
 *  Implementation fiel for class Gaudi::Math::GSL::EigenSystem
 *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
 *  @date 2005-05-24
 */
// ============================================================================

// ============================================================================
// Standard constructor, initializes variables
// ============================================================================
Ostap::Math::GSL::EigenSystem::EigenSystem() 
  : m_dim1   ( 0 )
  , m_dim2   ( 0 )
  , m_work1  ( 0 ) 
  , m_work2  ( 0 ) 
  , m_matrix ( 0 ) 
  , m_evec   ( 0 ) 
  , m_vector ( 0 )  
{} 
// ============================================================================
// Destructor
// ============================================================================
Ostap::Math::GSL::EigenSystem::~EigenSystem() 
{
  if ( 0 != m_matrix ) { gsl_matrix_free      ( m_matrix ) ; }
  if ( 0 != m_evec   ) { gsl_matrix_free      ( m_evec   ) ; }
  if ( 0 != m_vector ) { gsl_vector_free      ( m_vector ) ; }
  if ( 0 != m_work1  ) { gsl_eigen_symm_free  ( m_work1  ) ; }
  if ( 0 != m_work2  ) { gsl_eigen_symmv_free ( m_work2  ) ; }
} 
// ============================================================================
// check/adjust the  internal structures 
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::EigenSystem::_check
( const unsigned int D ) const 
{
  // check existing GSL matrix, (re)allocate if needed  
  if ( 0 != m_matrix && ( D != m_matrix->size1 || D != m_matrix->size2 ) ) 
  { gsl_matrix_free ( m_matrix ) ; m_matrix = 0 ; }
  if ( 0 == m_matrix   ) { m_matrix = gsl_matrix_alloc ( D , D ) ; }
  if ( 0 == m_matrix   ) { return Ostap::StatusCode ( MatrixAllocationFailure    ) ; }
  
  // check existing GSL matrix, (re)allocate if needed 
  if ( 0 != m_evec && ( D != m_evec->size1 || D != m_evec->size2 ) ) 
  { gsl_matrix_free ( m_evec ) ; m_evec = 0 ; }
  if ( 0 == m_evec     ) { m_evec = gsl_matrix_alloc   ( D , D ) ; }
  if ( 0 == m_evec     ) { return Ostap::StatusCode ( MatrixAllocationFailure    ) ; }
  
  // check the GSL vector, (re)allocate if needed 
  if ( 0 != m_vector && D != m_vector->size )
  { gsl_vector_free ( m_vector ) ; m_vector = 0 ; }
  if ( 0 == m_vector   ) { m_vector = gsl_vector_alloc ( D ) ; }
  if ( 0 == m_vector   ) { return Ostap::StatusCode ( VectorAllocationFailure    ) ; }
 
  return Ostap::StatusCode::SUCCESS ;
} 
// ============================================================================
// find the eigenvalues   (& sort them if needed ) 
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::EigenSystem::_fun1 
( const bool sorted ) const
{
  // check the working space, (re)allocate if needed  
  if ( 0 != m_work1 && m_dim1 != m_vector->size ) 
  { gsl_eigen_symm_free  ( m_work1 ) ; m_work1 = 0 ; }
  if ( 0 == m_work1 ) 
  { m_dim1 = m_vector->size ; m_work1 = gsl_eigen_symm_alloc ( m_dim1 ) ; ; }
  if ( 0 == m_work1 ) { return Ostap::StatusCode ( WorkspaceAllocationFailure ) ; }
  
  const int result = gsl_eigen_symm  ( m_matrix , m_vector , m_work1 ) ;
  if ( result ) { return Ostap::StatusCode ( ErrorFromGSL + result ) ; }
  if ( sorted ) { gsl_sort_vector   ( m_vector )  ; }
  return Ostap::StatusCode::SUCCESS ;
} 
// ============================================================================
// find the eigenvalues&eigenvectors (& sort them if needed ) 
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::EigenSystem::_fun2 
( const bool sorted ) const
{
  // check the working space, (re)allocate if needed  
  if ( 0 != m_work2 &&  m_dim2 != m_vector->size ) 
  { gsl_eigen_symmv_free ( m_work2 ) ; m_work2 = 0 ; }  
  if ( 0 == m_work2 ) 
  { m_dim2 = m_vector->size ; m_work2 = gsl_eigen_symmv_alloc ( m_dim2 ) ; ; }
  if ( 0 == m_work2 ) { return Ostap::StatusCode ( WorkspaceAllocationFailure ) ; }
  
  const int result = gsl_eigen_symmv ( m_matrix , m_vector , m_evec , m_work2 ) ;
  if ( result ) { return Ostap::StatusCode ( ErrorFromGSL + result ) ; }
  if ( sorted ) 
  { gsl_eigen_symmv_sort ( m_vector , m_evec , GSL_EIGEN_SORT_VAL_ASC )  ; }

  return Ostap::StatusCode::SUCCESS ;
} 
// ============================================================================
// thrown the exception 
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::EigenSystem::Exception
( const Ostap::StatusCode& sc ) const 
{
  Ostap::throwException
    ("Error from Ostap::Math::GSL::EigenSystem" , "*GSL*" , sc ) ;
  return sc ;
}
// ============================================================================

// ============================================================================
// The END 
// ============================================================================


