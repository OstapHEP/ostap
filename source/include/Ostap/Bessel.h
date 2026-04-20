// ============================================================================
#ifndef OSTAP_BESSEL_H 
#define OSTAP_BESSEL_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <cstdint>
#include <vector>
#include <complex>
#include <cmath>
#include <numeric>
#include <algorithm>
// ============================================================================
// Ostap
// =============================================================================
#include "Ostap/Math.h"
// ============================================================================
/** @file Ostap/Bessel.h
 *  Bessel's functions 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    // Bessel functions 
    // ========================================================================

    // ========================================================================
    /** regular Bessel function of the first kind
     *  \f$ J_n(x)\f$ 
     * @see https://en.wikipedia.org/wiki/Bessel_function
     * @see gsl_sf_bessel_J0
     * @see gsl_sf_bessel_J1
     * @see gsl_sf_bessel_Jn
     */
    double bessel_Jn
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** Irregular Bessel function of the first kind
     *  \f$ Y_n(x)\f$ 
     * @see https://en.wikipedia.org/wiki/Bessel_function
     * @see gsl_sf_bessel_Y0
     * @see gsl_sf_bessel_Y1
     * @see gsl_sf_bessel_Yn
     */
    double bessel_Yn
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** regular Bessel function of the first kind
     *  \f$ J_{\nu}(x)\f$ 
     * @see https://en.wikipedia.org/wiki/Bessel_function
     * @see gsl_sf_bessel_Jnu
     */
    double bessel_Jnu
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
    /** Irregular Bessel function of the first kind
     *  \f$ Y_{\nu}(x)\f$ 
     * @see https://en.wikipedia.org/wiki/Bessel_function
     * @see gsl_sf_bessel_Y0
     * @see gsl_sf_bessel_Y1
     * @see gsl_sf_bessel_Yn
     */
    double bessel_Ynu
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
    /** modified Bessel function of the fist kind  
     *  \f$ I_n(x) \f$ for \f$ x>0 \f$
     * @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1 
     * @see gsl_sf_bessel_I0
     * @see gsl_sf_bessel_I1
     * @see gsl_sf_bessel_In
     */
    double bessel_In
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** modified Bessel function of the second kind  
     *  \f$ K_n(x) \f$ for \f$ x>0 \f$
     *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
     *  @see gsl_sf_bessel_K0_e 
     *  @see gsl_sf_bessel_K1_e 
     *  @see gsl_sf_bessel_Kn_e 
     */
    double bessel_Kn    
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** modified Bessel function of the second kind  
     *  \f$ I_{\nu}(x) \f$ for \f$ x>0, \nu>0 \f$
     *  @see https://en.wikipedia.org/wiki/Bessel_function
     *  @see gsl_sf_bessel_Inu_e 
     */
    double bessel_Inu   
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
    /** modified Bessel function of the second kind  
     *  \f$ K_{\nu}(x) \f$ for \f$ x>0, \nu>0 \f$
     *  @see https://en.wikipedia.org/wiki/Bessel_function
     *  @see gsl_sf_bessel_Knu_e 
     */
    double bessel_Knu   
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
    /** scaled modified Bessel function of the second kind 
     *  \f$ \mathrm{e}^x K_n(x) \f$ for \f$ x>0 \f$
     *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
     *  @see gsl_sf_bessel_K0_scaled_e 
     *  @see gsl_sf_bessel_K1_scaled_e 
     *  @see gsl_sf_bessel_Kn_scaled_e 
     */
    double bessel_Kn_scaled  
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** scaled modified Bessel function of the second kind 
     *  \f$ \mathrm{e}^x K_{\nu}(x) \f$ for \f$ x>0, \nu>0 \f$
     *  @see https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I%CE%B1,_K%CE%B1
     *  @see gsl_sf_bessel_Knu_scaled_e 
     */
    double bessel_Knu_scaled
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
       
    // ========================================================================
    // derivatives for Bessel functions 
    // ========================================================================
    
    // ========================================================================
    /** derivative for the  regular Bessel function of the first kind
     *  - \f$ J_0^{\prime}(x) =  - J_1(s)\f$
     *  - \f$ J_{n}^{\prime}(x) = \frac{1}{2}\left(J_{n-1}(x) - J_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_Jn
     */
    double der_bessel_Jn
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** derivative for the  regular Bessel function of the first kind
     *  - \f$ J_0^{\prime}(x) =  - J_1(s)\f$
     *  - \f$ J_{n}^{\prime}(x) = \frac{1}{2}\left(J_{n-1}(x) - J_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_Jnu
     */
    double der_bessel_Jnu
    ( const double nu , 
      const double x  ) ;
    // =========================================================================
    /** derivative for the  irregular Bessel function of the first kind
     *  - \f$ Y_0^{\prime}(x) =  - Y_1(s)\f$
     *  - \f$ Y_{n}^{\prime}(x) = \frac{1}{2}\left(Y_{n-1}(x) - Y_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_Yn
     */
    double der_bessel_Yn
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** derivative for the irregular Bessel function of the first kind
     *  - \f$ Y_0^{\prime}(x) =  - Y_1(s)\f$
     *  - \f$ Y_{n}^{\prime}(x) = \frac{1}{2}\left(Y_{n-1}(x) - Y_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_Ynu
     */
    double der_bessel_Ynu
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
    /** derivative for the  modified Bessel function
     *  - \f$ I_0^{\prime}(x) =  I_1(s)\f$
     *  - \f$ I_{n}^{\prime}(x) = \frac{1}{2}\left(I_{n-1}(x) + I_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_In
     */
    double der_bessel_In
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** derivative for the modified Bessel function
     *  - \f$ I_0^{\prime}(x) =  I_1(s)\f$
     *  - \f$ I_{n}^{\prime}(x) = \frac{1}{2}\left(I_{n-1}(x) + I_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_Inu
     */
    double der_bessel_Inu
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
    /** derivative for the  modified Bessel function
     *  - \f$ K_0^{\prime}(x) = -K_1(s)\f$
     *  - \f$ K_{n}^{\prime}(x) = -\frac{1}{2}\left(K_{n-1}(x) + K_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_Kn
     */
    double der_bessel_Kn
    ( const int    n , 
      const double x ) ;
    // ========================================================================
    /** derivative for the modified Bessel function
     *  - \f$ K_0^{\prime}(x) = -K_1(s)\f$
     *  - \f$ K_{n}^{\prime}(x) = -\frac{1}{2}\left(K_{n-1}(x) + K_{n+1}(x)\right)\f$ 
     *  @see Ostap::Math::bessel_Knu
     */
    double der_bessel_Knu
    ( const double nu , 
      const double x  ) ;
    // ========================================================================
    
    
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_BESSEL_H
// ============================================================================
