// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Parameters.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Parameters
 *  @see Ostap::Math::Parameters
 *  @author Vanya BELYAEV
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<double>::is_specialized           , 
                  "mumeric_limits are not specialized for doubles"      ) ;
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double> s_equal {} ;       // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>     s_zero  {} ;       // zero for doubles
  /// zero fo vectors 
  const Ostap::Math::Zero< std::vector<double> > s_vzero {} ; // zero for vectors
  // ==========================================================================
}
// ============================================================================
// PARAMETERS
// ============================================================================
// constructor from number of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
( const std::size_t np ) 
  : m_pars ( np , 0.0 ) 
{}
// ============================================================================
// constructor from  the list of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
( const std::vector<double>&  pars   ) 
  : m_pars ( pars ) 
{}
// ============================================================================
// constructor from the list of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
(       std::vector<double>&& pars   ) 
  : m_pars ( std::forward<std::vector<double> > ( pars ) ) 
{}
// ============================================================================
// all zero ?
// ============================================================================
bool Ostap::Math::Parameters::zero  () const { return s_vzero ( m_pars ) ; }
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Parameters::_setPar 
( const std::size_t k     , 
  const double      value ) 
{
  if ( m_pars.size() <= k               ) { return false ; }
  if ( s_equal ( m_pars [ k ] , value ) ) { return false ; }
  m_pars [ k ] = value ;
  return true ;
}
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::Parameters::swap (  Ostap::Math::Parameters& right ) 
{ std::swap ( m_pars ,  right.m_pars ); }


// ============================================================================
//                                                                      The END 
// ============================================================================
