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
#include "Ostap/Hilbert.h"
// =============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Hilbert
 *  @date 2024-10-27 
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 */
// =============================================================================
/* constructor from the function
 * @param func   the function
 * @param tag     unique tag/label for cache 
 * @param rescale rescale function for better numerical precision 
 * @param aprecision absolute precision 
 * @param rprecision relative precision 
 * @param size    size of integration workspace  
 */
// =============================================================================
Ostap::Math::Hilbert::Hilbert 
( Ostap::Math::Hilbert::function1 func       ,
  const std::size_t               tag        ,
  const unsigned short            rescale    ,
  const double                    aprecision ,
  const double                    rprecision ,
  const double                    width      ,           
  const std::size_t               size       )
  : m_func       ( func       )
  , m_tag        ( tag        ) 
  , m_rescale    ( rescale    )
  , m_aprecision ( aprecision ) 
  , m_rprecision ( rprecision )
  , m_width      ( width      )
  , m_integrator ( size       )
{} 
// =============================================================================
double Ostap::Math::Hilbert::operator() ( const double x ) const
{ return m_integrator.cauchy_pv_infinity 
    ( std::cref ( m_func ) ,
      x                    ,
      m_tag                ,
      m_rescale            ,
      m_aprecision         ,
      m_rprecision         ,
      m_width              ) / M_PI ; }
// =============================================================================
//                                                                       The END 
// =============================================================================


