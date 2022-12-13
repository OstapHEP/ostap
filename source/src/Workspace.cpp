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
Ostap::Math::WorkSpace::WorkSpace 
( const std::size_t    size         , 
  const unsigned short size_cquad   , 
  const unsigned short size_romberg ) 
  : m_workspace         ( nullptr      )
  , m_workspace_cquad   ( nullptr      )
  , m_workspace_romberg ( nullptr      )
  , m_size              ( size         )
  , m_size_cquad        ( size_cquad   )
  , m_size_romberg      ( size_romberg )
{}
// ============================================================================
// (almost fictive) copy constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace
( const Ostap::Math::WorkSpace& right          )
  : m_workspace         ( nullptr              )
  , m_workspace_cquad   ( nullptr              )
  , m_workspace_romberg ( nullptr              )
  , m_size              ( right.m_size         )
  , m_size_cquad        ( right.m_size_cquad   )
  , m_size_romberg      ( right.m_size_romberg )
{}
// ============================================================================
// move constructor
// ============================================================================
Ostap::Math::WorkSpace::WorkSpace
(       Ostap::Math::WorkSpace&& right )
  : m_workspace         ( nullptr )
  , m_workspace_cquad   ( nullptr )
  , m_workspace_romberg ( nullptr )
  , m_size              ( 0       )
  , m_size_cquad        ( 0       )
  , m_size_romberg      ( 0       )
{ 
  std::swap ( m_workspace         , right.m_workspace         ) ; 
  std::swap ( m_workspace_cquad   , right.m_workspace_cquad   ) ; 
  std::swap ( m_workspace_romberg , right.m_workspace_romberg ) ; 
  std::swap ( m_size              , right.m_size              ) ; 
  std::swap ( m_size_cquad        , right.m_size_cquad        ) ; 
  std::swap ( m_size_romberg      , right.m_size_romberg      ) ; 
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::WorkSpace::~WorkSpace ()
{
  //
  if ( nullptr != m_workspace )
  {
    gsl_integration_workspace * _ws = 
      (gsl_integration_workspace*) m_workspace ;
    m_workspace = nullptr ;
    gsl_integration_workspace_free ( _ws );
  }
  //
  if ( nullptr != m_workspace_cquad )
  {
    gsl_integration_cquad_workspace * _ws = 
      (gsl_integration_cquad_workspace*) m_workspace_cquad ;
    m_workspace_cquad = nullptr ;
    gsl_integration_cquad_workspace_free ( _ws );
  }
  //
  if ( nullptr != m_workspace_romberg )
  {
    gsl_integration_romberg_workspace * _ws = 
      (gsl_integration_romberg_workspace*) m_workspace_romberg ;
    m_workspace_romberg = nullptr ;
    gsl_integration_romberg_free ( _ws );
  }
  //
}
// ============================================================================
// get the main integration workspace
// ============================================================================
void* Ostap::Math::WorkSpace::workspace () const
{
  if ( nullptr == m_workspace ) 
  { m_workspace = gsl_integration_workspace_alloc ( 0 == m_size ? s_SIZE : m_size ) ; }
  return m_workspace ;
}
// ============================================================================
// get the integration workspace for CQUAD integrtaor 
// ============================================================================
void* Ostap::Math::WorkSpace::workspace_cquad () const
{
  if ( nullptr == m_workspace_cquad ) 
  { m_workspace_cquad = gsl_integration_cquad_workspace_alloc
      ( 0 == m_size_cquad ? s_SIZE_CQUAD : m_size_cquad ) ; }
  return m_workspace_cquad ;
}
// ============================================================================
// get the integration workspace for Romberg integrtaor 
// ============================================================================
void* Ostap::Math::WorkSpace::workspace_romberg () const
{
  if ( nullptr == m_workspace_romberg ) 
  { m_workspace_romberg = gsl_integration_romberg_alloc
      ( 0 == m_size_romberg ? s_SIZE_ROMBERG : m_size_romberg ) ; }
  return m_workspace_romberg ;
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
(       Ostap::Math::WorkSpace&& /* right */ ) { return *this ; }
// ============================================================================
// swap
// ============================================================================
void Ostap::Math::WorkSpace::swap ( Ostap::Math::WorkSpace& right ) 
{ 
  std::swap  ( m_workspace         ,  right.m_workspace         ) ;
  std::swap  ( m_workspace_cquad   ,  right.m_workspace_cquad   ) ;
  std::swap  ( m_workspace_romberg ,  right.m_workspace_romberg ) ;
  std::swap  ( m_size              ,  right.m_size              ) ;
  std::swap  ( m_size_cquad        ,  right.m_size_cquad        ) ;
  std::swap  ( m_size_romberg      ,  right.m_size_romberg      ) ;
}
// ============================================================================
// resize the main integration workspace 
// ============================================================================
std::size_t Ostap::Math::WorkSpace::resize ( const std::size_t newsize ) 
{
  if      ( m_size  == newsize     ) { return m_size ; }
  else if ( nullptr != m_workspace ) 
  {
    // free existing workspace:
    gsl_integration_workspace * _ws = (gsl_integration_workspace*) m_workspace ;
    m_workspace = nullptr ;
    gsl_integration_workspace_free ( _ws );
  }
  //
  m_size = newsize ;
  return m_size ;
}
// ============================================================================
// resize the integration workspace for CQUAD integrator 
// ============================================================================
std::size_t Ostap::Math::WorkSpace::resize_cquad 
( const std::size_t newsize ) 
{
  if      ( m_size_cquad  == newsize     ) { return m_size_cquad ; }
  else if ( nullptr != m_workspace_cquad ) 
  {
    // free existing workspace:
    gsl_integration_cquad_workspace * _ws = 
      (gsl_integration_cquad_workspace*) m_workspace_cquad  ;
    m_workspace_cquad = nullptr ;
    gsl_integration_cquad_workspace_free ( _ws );
  }
  //
  m_size_cquad = newsize ;
  return m_size_cquad ;
}
// ============================================================================
// resize the integration workspace for Romberg integrator 
// ============================================================================
std::size_t Ostap::Math::WorkSpace::resize_romberg 
( const std::size_t newsize ) 
{
  if      ( m_size_romberg  == newsize     ) { return m_size_romberg ; }
  else if ( nullptr != m_workspace_romberg ) 
  {
    // free existing workspace:
    gsl_integration_romberg_workspace * _ws = 
      (gsl_integration_romberg_workspace*) m_workspace_romberg  ;
    m_workspace_romberg = nullptr ;
    gsl_integration_romberg_free ( _ws );
  }
  //
  m_size_romberg = newsize ;
  return m_size_romberg ;
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
