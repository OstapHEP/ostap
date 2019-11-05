// ============================================================================
// Include files 
// ============================================================================
// STD  & STL
// ============================================================================
#include <utility>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_errno.h"
#include "gsl/gsl_integration.h"
// ============================================================================
// Ostap
#include "Ostap/Workspace.h"
// ============================================================================
// Local
// ============================================================================
#include "local_gsl.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Workspace
 */
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace ( const std::size_t size ) 
  : m_workspace ( nullptr )
  , m_size      ( size )
{}
// ============================================================================
// (fictive) copy constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace
( const Ostap::Math::WorkSpace& /* right */ )
  : m_workspace ( nullptr )
  , m_size      ( 0       )
{}
// ============================================================================
// (move constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace
(       Ostap::Math::WorkSpace&& right  )
  : m_workspace ( nullptr )
  , m_size      ( 0       )
{ 
  std::swap ( m_workspace , right.m_workspace ) ; 
  std::swap ( m_size      , right.m_size      ) ; 
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::WorkSpace::~WorkSpace ()
{
  if ( nullptr != m_workspace )
  {
    gsl_integration_workspace * _ws = (gsl_integration_workspace*) m_workspace ;
    m_workspace = nullptr ;
    gsl_integration_workspace_free ( _ws );
  }
}
// ============================================================================
// get the integration workspace
// ============================================================================
void* Ostap::Math::WorkSpace::workspace () const
{
  if ( nullptr == m_workspace ) 
  { m_workspace = (char*) gsl_integration_workspace_alloc ( 0 == m_size ? s_SIZE : m_size ) ; }
  return m_workspace ;
}
// ============================================================================
// fictive assignement operator
// ============================================================================
Ostap::Math::WorkSpace&
Ostap::Math::WorkSpace::operator=
  ( const Ostap::Math::WorkSpace& /* right */ ) { return *this ; }
// ============================================================================
// move assignement operator
// ============================================================================
Ostap::Math::WorkSpace&
Ostap::Math::WorkSpace::operator=
(       Ostap::Math::WorkSpace&& right ) 
{ 
  std::swap  ( m_workspace ,  right.m_workspace ) ;
  std::swap  ( m_size      ,  right.m_size      ) ;
  return *this ;
}
// ============================================================================
// swap
// ============================================================================
void Ostap::Math::WorkSpace::swap ( Ostap::Math::WorkSpace& right ) 
{ 
  std::swap  ( m_workspace ,  right.m_workspace ) ;
  std::swap  ( m_size      ,  right.m_size      ) ;
}
// ============================================================================
// resize the workspace 
// ============================================================================
std::size_t Ostap::Math::WorkSpace::resize ( const std::size_t newsize ) 
{
  if      ( m_size  == newsize     ) { return m_size ; }
  else if ( nullptr != m_workspace ) 
  {
    //
    // free existing workspace:
    //
    gsl_integration_workspace * _ws = (gsl_integration_workspace*) m_workspace ;
    m_workspace = nullptr ;
    gsl_integration_workspace_free ( _ws );
  }
  //
  m_size = newsize ;
  return m_size ;
}
// ============================================================================


// ============================================================================
//  The END 
// ============================================================================
