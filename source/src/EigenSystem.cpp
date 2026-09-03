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
#include "Ostap/LinAlg.h"
#include "Ostap/EigenSystem.h"
// ============================================================================
// Local
// ============================================================================
#include "GSL_sentry.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::GSL::EigenSystem
 *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
 *  @date 2005-05-24
 */
// ============================================================================
Ostap::Math::GSL::EigenSystem::EigenSystem
( const unsigned short N ) 
  : m_ws_symm  { nullptr } 
  , m_ws_symmv { nullptr }
{
  if ( 0 < N )
  {
    m_ws_symm  = gsl_eigen_symm_alloc  ( N ) ;
    Ostap::Assert ( m_ws_symm ,
                    "(GSL)EigenSYMM  workspace allocation failure" ,
                    "Ostap::Math::GSL::EigenSystem        "       ,
                    WORKSPACE_ALLOCATION_FAILURE                  , __FILE__ , __LINE__ ) ;
    //
    m_ws_symmv = gsl_eigen_symmv_alloc ( N ) ;
    Ostap::Assert ( m_ws_symmv ,
                    "(GSL)EigenSYMMV workspace allocation failure" ,
                    "Ostap::Math::GSL::EigenSystem        "       ,
                    WORKSPACE_ALLOCATION_FAILURE                  , __FILE__ , __LINE__ ) ;  
  }
}
// ===========================================================================
// copy constructor
// ===========================================================================
Ostap::Math::GSL::EigenSystem::EigenSystem
( const Ostap::Math::GSL::EigenSystem& /* right */ )
  : m_ws_symm  { nullptr } 
  , m_ws_symmv { nullptr }
{}
// ===========================================================================
// move constructor
// ===========================================================================
Ostap::Math::GSL::EigenSystem::EigenSystem
( Ostap::Math::GSL::EigenSystem&& right  )
  : m_ws_symm  { nullptr } 
  , m_ws_symmv { nullptr }
{
  this->swap ( right ) ;
}
// ============================================================================
// Destructor
// ============================================================================
Ostap::Math::GSL::EigenSystem::~EigenSystem() { release () ; } 
// ============================================================================
// assignement operator
// ============================================================================
Ostap::Math::GSL::EigenSystem&
Ostap::Math::GSL::EigenSystem::operator=
( const Ostap::Math::GSL::EigenSystem& /* right */ )
{ return *this ; }
// ============================================================================
// assignement operator
// ============================================================================
Ostap::Math::GSL::EigenSystem&
Ostap::Math::GSL::EigenSystem::operator=
( Ostap::Math::GSL::EigenSystem&& right  )
{
  if ( &right == this ) { return *this ; }
  this->swap ( right ) ;
  return *this ;
}
// ============================================================================
// assignement operator
// ============================================================================
void Ostap::Math::GSL::EigenSystem::swap 
( Ostap::Math::GSL::EigenSystem& right  )
{
  if ( &right != this ) { return ; } 
  std::swap ( m_ws_symm  , right.m_ws_symm  ) ;
  std::swap ( m_ws_symmv , right.m_ws_symmv ) ;       
}
// ============================================================================
// release allocated workspaces 
// ============================================================================
void Ostap::Math::GSL::EigenSystem::release ()
{
  release1 () ;
  release2 () ;
}  
// ============================================================================
// release allocated SYMM  workspace 
// ============================================================================
void Ostap::Math::GSL::EigenSystem::release1 () const 
{ if ( m_ws_symm  ) { gsl_eigen_symm_free  ( static_cast<gsl_eigen_symm_workspace*>  ( m_ws_symm  ) ) ; m_ws_symm  = nullptr ; } }  
// ============================================================================
// release allocated SYMMV workspace
// ============================================================================
void Ostap::Math::GSL::EigenSystem::release2 () const 
{ if ( m_ws_symmv ) { gsl_eigen_symmv_free ( static_cast<gsl_eigen_symmv_workspace*> ( m_ws_symmv ) ) ; m_ws_symmv = nullptr ; } }  
// ============================================================================
/*  Get the eigenvalues of symmetrical matrix
 *  @param mtrx   (INPUT) input matrix
 *  @param values (UPDATE) output vector of eigenvalues
 *  @param ascending (INPUT)  sorting order 
 *  @param sorted (INPUT) get the eigenvaleus sorted ?
 *  @return Status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::GSL::EigenSystem::eigenValues
( const Ostap::Math::GSL::Matrix& mtrx      ,
  Ostap::Math::GSL::Vector&       values    ,
  const bool                      sorted    , 
  const bool                      ascending ) const 
{
  // square matrix? 
  if ( mtrx.nRows() != mtrx.nCols() )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "eigenValues: input matrix is not square!" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;    
  }
  const std::size_t N = mtrx.nRows() ;
  //
  // workspace exists but too small
  if (  m_ws_symm && static_cast<gsl_eigen_symm_workspace*> ( m_ws_symm ) -> size < N ) { release1 () ; } 
  //
  if ( !m_ws_symm )
  {
    m_ws_symm  = gsl_eigen_symm_alloc ( N ) ;    
    Ostap::Assert ( m_ws_symm ,
                    "(GSL)EigenSYMM  workspace allocation failure" ,
                    "Ostap::Math::GSL::EigenSystem::eigenValues"   ,
                    WORKSPACE_ALLOCATION_FAILURE                   , __FILE__ , __LINE__ ) ;
  }
  //
  /// working matrix
  Ostap::Math::GSL::Matrix A { mtrx } ; 
  /// resize output vector 
  if ( values.size() != N ) { values.resize ( N ) ;}
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  const int status = gsl_eigen_symm ( A     .matrix () ,
                                      values.vector () ,
                                      static_cast<gsl_eigen_symm_workspace*> ( m_ws_symm ) ) ;
  if ( status )
  {
    gsl_error ( "Ostap::Math::GSL::EigenSystem::eigenValues: error from gsl_eigen_symm" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ; 
  }
  //
  if ( sorted )
  {
    gsl_sort_vector   ( values.vector () )  ;
    if ( !ascending ) { gsl_vector_reverse  ( values.vector () ) ; }
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ====================================================================
/*  Get the eigenvalues & eigenvectors of symmetrical matrix
 *  @param mtrx      (INPUT)  input matrix
 *  @param values    (UPDATE) output vector of eigenvalues
 *  @param vectors   (UPDATE) output matrix where each column is eigen-vector
 *  @param sorted    (INPUT)  get the eigenvalues sorted ?
 *  @param ascending (INPUT)  sorting order 
 *  @return Status code 
 */
// ====================================================================
Ostap::StatusCode
Ostap::Math::GSL::EigenSystem::eigenVectors
( const Ostap::Math::GSL::Matrix& mtrx      ,
  Ostap::Math::GSL::Vector&       values    ,
  Ostap::Math::GSL::Matrix&       vectors   ,
  const bool                      sorted    ,
  const bool                      ascending ) const 
{
  // square matrix? 
  if ( mtrx.nRows() != mtrx.nCols() )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "eigenVectors: input matrix is not square!" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;    
  }
  //
  const std::size_t N = mtrx.nRows () ;
  //
  // workspace exists but too small
  if (  m_ws_symmv && static_cast<gsl_eigen_symmv_workspace*> ( m_ws_symmv ) -> size < N ) { release2 () ; } 
  //
  if ( !m_ws_symmv )
  {
    m_ws_symmv  = gsl_eigen_symmv_alloc ( N ) ;    
    Ostap::Assert ( m_ws_symm ,
                    "(GSL)EigenSYMMS workspace allocation failure" ,
                    "Ostap::Math::GSL::EigenSystem::eigenVectors"  ,
                    WORKSPACE_ALLOCATION_FAILURE                   , __FILE__ , __LINE__ ) ;
  }
  /// working matrix 
  Ostap::Math::GSL::Matrix A { mtrx } ; 
  if ( values.size   () != N                          ) { values  .resize ( N     ) ; }
  if ( vectors.nRows () != N || vectors.nCols () != N ) { vectors .resize ( N , N ) ; }
  //
  // use GSL: 
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  const int status = gsl_eigen_symmv ( A       .matrix () ,
                                       values  .vector () ,
                                       vectors .matrix () ,
                                       static_cast<gsl_eigen_symmv_workspace*> ( m_ws_symmv ) ) ;
  if ( status )
  {
    gsl_error ( "Ostap::Math::GSL::EigenSystem::eigenVectors: error from gsl_eigen_symmv" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ; 
  }
  //
  if ( sorted ) { gsl_eigen_symmv_sort ( values . vector () ,
                                         vectors .matrix () , 
                                         ascending          ? GSL_EIGEN_SORT_VAL_ASC : GSL_EIGEN_SORT_VAL_DESC ) ; }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ========================================================================================
/*  Get the eigenvalues & eigenvectors of symmetrical matrix
 *  @param mtrx      (INPUT)  input matrix
 *  @param values    (UPDATE) output vector of eigenvalues
 *  @param vectors   (UPDATE) output matrix where each column is eigen-vector
 *  @param sorted    (INPUT)  get the eigenvaleus sorted ?
 *  @param ascending (INPUT)  sorting order 
 *  @return Status code 
 */
// =====================================================================
Ostap::StatusCode
Ostap::Math::GSL::EigenSystem::eigenVectors
( const Ostap::Math::GSL::Matrix&        mtrx      ,
  Ostap::Math::GSL::Vector&              values    ,
  std::vector<Ostap::Math::GSL::Vector>& vectors   ,
  const bool                             sorted    ,
  const bool                             ascending ) const 
{
  // square matrix? 
  if ( mtrx.nRows() != mtrx.nCols() )
  {
    const int status = GSL_EBADLEN ;
    gsl_error ( "eigenVectors: input matrix is not square!" , __FILE__ , __LINE__ , status ) ;
    return ERROR_GSL + status ;    
  }
  //
  const std::size_t N = mtrx.nRows () ;
  // eigenvectors 
  Ostap::Math::GSL::Matrix vctrs { N } ;
  // delegate actual calculation to another method
  Ostap::StatusCode sc = eigenVectors ( mtrx , values , vctrs , sorted , ascending ) ;
  if ( sc.isFailure() ) { vectors.clear() ; return sc ; }
  ///
  vectors.clear   () ;
  vectors.reserve ( N ) ;
  for ( std::size_t k = 0 ; k < N ; ++k )
  { vectors.emplace_back ( vctrs.column ( k ) ) ; }
  //
  return Ostap::StatusCode::SUCCESS ;
}
 
// ============================================================================
//                                                                      The END 
// ============================================================================


