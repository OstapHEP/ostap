// ============================================================================
// Include files 
// =============================================================================
// Incldue files 
// =============================================================================
// STD& STL 
// =============================================================================
#include <cmath>
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Laplace.h"
// =============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Laplace
 *  @date 2024-10-27 
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 */
// =============================================================================
double Ostap::Math::Laplace::operator() ( const double x ) const
{
  auto the_fun = [x,this] ( const double z ) -> double
  { return this->m_func ( z ) * std::exp ( - x * z ) ; };
  //
  const std::size_t ntag =
    ( 0 == m_tag ) ? m_tag : Ostap::Utils::hash_combiner ( m_tag , x ) ; 
  return m_integrator.integrate_to_infinity 
    ( std::cref ( the_fun ) ,
      0.0                   , // xmin    
      ntag                  , // tag 
      m_aprecision          ,
      m_rprecision          ) ;
}
// =============================================================================
//                                                                       The END 
// =============================================================================


