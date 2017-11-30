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
/** @file implementation file for class Ostap::Math::Workspace
 */
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace () : m_workspace ( nullptr ){}
// ============================================================================
// (fictive) copy constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace
( const Ostap::Math::WorkSpace& /* right */ )
  : m_workspace ( nullptr )
{}
// ============================================================================
// (move constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace
(       Ostap::Math::WorkSpace&& right  )
  : m_workspace ( nullptr )
{ 
  std::swap ( m_workspace , right.m_workspace ) ; 
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
  { m_workspace = (char*) gsl_integration_workspace_alloc ( s_SIZE ) ; }
  //
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
{ std::swap  ( m_workspace ,  right.m_workspace ) ;}
// ============================================================================
// swap
// ============================================================================
void Ostap::Math::WorkSpace::swap ( Ostap::Math::WorkSpace& right ) 
{ std::swap  ( m_workspace ,  right.m_workspace ) ; }
// ============================================================================


// ============================================================================
//  The END 
// ============================================================================
