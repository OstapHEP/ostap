// ============================================================================
// Include files 
// =============================================================================
// Incldue files 
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/KramersKronig.h"
// =============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::KramersKronig
 *  @date 2020-09-01 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// =============================================================================
/* the only important method
 *   \f[ \chi_{\omega } = 
 *    s \frac{s^n}{\pi} {\mathcal{P}} \int\limit^{+\infty}{\omega_0}
 *     \frac{ \rho (\omega^\prime} 
 *     { \omega^{\prime n} \left( \omega^{\prime} - \omega \right) }  
 *     d \omega^{\prime} \f]
 * @param x value of \f$ \omega \f$
 * @see Ostap::Math::Integrator
 * @see Ostap::Math::Integrator::kramers_kronig
 */
// =============================================================================
double Ostap::Math::KramersKronig::operator() ( const double x ) const
{ return m_scale * m_integrator.kramers_kronig 
    ( std::cref ( m_rho ) ,
      x                   ,
      m_omega0            ,
      m_n                 ,
      m_tag               ,
      m_rescale           , 
      m_aprecision        , 
      m_rprecision        , 
      m_width             ) ; }
// =============================================================================


// =============================================================================
//                                                                       The END 
// =============================================================================


