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
#include "Ostap/EigenSystem.h"
#include "Ostap/EigenSystem.icpp"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
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
  //
  Ostap::Assert ( nullptr != m_matrix ,
                  "(GSL)Matrix allocation failure" ,
                  "Ostap::Math::GSL::EigenSystem"  ,
                  MatrixAllocationFailure , __FILE__ , __LINE__ ) ;
  //
  // check existing GSL matrix, (re)allocate if needed 
  if ( 0 != m_evec && ( D != m_evec->size1 || D != m_evec->size2 ) ) 
    { gsl_matrix_free ( m_evec ) ; m_evec = 0 ; }
  if ( 0 == m_evec     ) { m_evec = gsl_matrix_alloc   ( D , D ) ; }
  //
  Ostap::Assert ( nullptr != m_evec ,
                  "(GSL)Vector allocation failure" ,
                  "Ostap::Math::GSL::EigenSystem"  ,
                  VectorAllocationFailure , __FILE__ , __LINE__ ) ;
  //
  // check the GSL vector, (re)allocate if needed 
  if ( 0 != m_vector && D != m_vector->size )
    { gsl_vector_free ( m_vector ) ; m_vector = 0 ; }
  if ( 0 == m_vector   ) { m_vector = gsl_vector_alloc ( D ) ; }
  //
  Ostap::Assert ( nullptr != m_vector  ,
                  "(GSL)Vector allocation failure" ,
                  "Ostap::Math::GSL::EigenSystem"  ,
                  VectorAllocationFailure , __FILE__ , __LINE__ ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
} 
// ============================================================================
// find the eigenvalues   (& sort them if needed ) 
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::EigenSystem::_fun1 
( const bool sorted ) const
{
  // check the working space, (re)allocate if needed  
  if ( 0 != m_work1 && m_dim1 != m_vector->size ) 
    { gsl_eigen_symm_free  ( m_work1 ) ; m_work1 = 0 ; }
  if ( 0 == m_work1 ) 
    { m_dim1 = m_vector->size ; m_work1 = gsl_eigen_symm_alloc ( m_dim1 ) ; ; }
  // 
  Ostap::Assert ( nullptr != m_work1                  ,
                  "(GSL)Workspace allocation failure" ,
                  "Ostap::Math::GSL::EigenSystem"     ,
                  WorkspaceAllocationFailure , __FILE__ , __LINE__ ) ;
  //
  const int result = gsl_eigen_symm  ( m_matrix , m_vector , m_work1 ) ;
  Ostap::Assert ( !result                          ,
                  "Error from gsl_eigen_symm/+200" ,
                  "Ostap::Math::GSL::EigenSystem"  ,
                  ErrorFromGSL + 1 , __FILE__ , __LINE__ ) ;
  //
  if ( sorted ) { gsl_sort_vector   ( m_vector )  ; }
  //
  return Ostap::StatusCode::SUCCESS ;
} 
// ============================================================================
// find the eigenvalues&eigenvectors (& sort them if needed ) 
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::EigenSystem::_fun2 
( const bool sorted    ,
  const bool ascending ) const
{
  // check the working space, (re)allocate if needed  
  if ( 0 != m_work2 &&  m_dim2 != m_vector->size ) 
  { gsl_eigen_symmv_free ( m_work2 ) ; m_work2 = 0 ; }  
  if ( 0 == m_work2 ) 
  { m_dim2 = m_vector->size ; m_work2 = gsl_eigen_symmv_alloc ( m_dim2 ) ; ; }
  //
  Ostap::Assert ( nullptr != m_work2                  ,
                  "(GSL)Workspace allocation failure" ,
                  "Ostap::Math::GSL::EigenSystem"     ,
                  WorkspaceAllocationFailure , __FILE__ , __LINE__ ) ;
  //
  const int result = gsl_eigen_symmv ( m_matrix , m_vector , m_evec , m_work2 ) ;
  Ostap::Assert ( !result                           ,
                  "Error from gsl_eigen_symmv/+200" ,
                  "Ostap::Math::GSL::EigenSystem"   ,
                  ErrorFromGSL + 1 , __FILE__ , __LINE__ ) ;
  //
  if ( sorted ) { gsl_eigen_symmv_sort
      ( m_vector  ,
        m_evec    ,
        ascending ? GSL_EIGEN_SORT_VAL_ASC : GSL_EIGEN_SORT_VAL_DESC ) ; }
  //
  return Ostap::StatusCode::SUCCESS ;
} 
// ============================================================================
// thrown the exception 
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::EigenSystem::Exception
( const Ostap::StatusCode& sc   ,
  const char*              file ,
  const std::size_t        line ) const 
{
  Ostap::Assert ( sc.isSuccess ()                 ,
                  "Eigensystem error"             , 
                  "Ostap::Math::GSL::EigenSystem" , sc , file , line ) ;
  return sc ;
}
// ============================================================================

// ============================================================================
// The END 
// ============================================================================


