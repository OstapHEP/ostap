// ============================================================================
#ifndef OSTAP_KRAMERSKRONIG_H 
#define OSTAP_KRAMERSKRONIG_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Integrator.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace  Math
  {
    // ========================================================================
    /** @class KramersKronig Ostap/KramersKronig.h
     *  Simple clss to implement Kramers-Kronig relations 
     *  @see  https://en.wikipedia.org/wiki/Kramers%E2%80%93Kronig_relations
     *
     *   \f[ \chi_{\omega } = 
     *    \frac{s^n}{\pi} {\mathcal{P}} \int\limit^{+\infty}{\omega_0}
     *     \frac{ \rho (\omega^\prime} 
     *     { \omega^{\prime n} \left( \omega^{\prime} - \omega \right) }  
     *     d \omega^{\prime} \f]
     *  - Note  the sign! 
     *  @see Ostap::Math::Integrtor
     *  @author Vanya Belyaev@itep.ru
     *  @date   2020-08-01
     */
    class KramersKronig
    {
    public:
      // ======================================================================
      /**  templated contructor from the function, low integration edge, 
       *   number of subtractions and the scale factor
       *   @param rho the function
       *   @param omega0 low intergation edge  
       *   @param n      number of subtractions  
       *   @param scale  scale factor (e.g. sign)
       *   @param size   size of integrtaion workspace  
       */
      template <class FUNCTION>
      KramersKronig ( FUNCTION             rho       ,
                      const double         omega0    ,
                      const unsigned short n     = 0 ,
                      const double         scale = 1 , 
                      const std::size_t    size  = 0 )
        : m_rho        ( rho    )
        , m_omega0     ( omega0 )
        , m_n          ( n      )
        , m_scale      ( scale  )          
        , m_integrator ( size   )
      {}
      // ======================================================================
    public:
      // ======================================================================
      /**  create the object from the function, low integration edge, number of 
       *   subtractions and scale  factor 
       *   @param rho the function
       *   @param omega0 low intergation edge  
       *   @param n      number of subtractions  
       *   @param scale  scale factor (e.g. sign)
       *   @param size   size of integrtaion workspace  
       */
      template <class FUNCTION>
      inline static KramersKronig
      create  ( FUNCTION             rho        ,
                const double         omega0     ,
                const unsigned short n      = 0 ,
                const double         scale  = 1 , 
                const std::size_t    size   = 0 )
      { return KramersKronig ( rho , omega0 , n , scale , size ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** the only important method
       *   \f[ \chi_{\omega } = 
       *    s \frac{s^n}{\pi} {\mathcal{P}} \int\limit^{+\infty}{\omega_0}
       *     \frac{ \rho (\omega^\prime} 
       *     { \omega^{\prime n} \left( \omega^{\prime} - \omega \right) }  
       *     d \omega^{\prime} \f]
       * @param x value of \f$ \omega \f$
       * @see Ostap::Math::Integrator
       * @see Ostap::Math::Integrator::kramers_kronig
       */
      double operator() ( const double x ) const
      { return m_scale * m_integrator.kramers_kronig ( std::cref ( m_rho ) , x , m_omega0 , m_n ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the function 
      std::function<double(double)>  m_rho    ; // the function
      /// the low integrtaion limit
      double                         m_omega0 ; // low integration limit 
      /// number of subtractions
      unsigned short                 m_n      ; // number of subtractions
      /// scale factor (e.g. sign) 
      unsigned short                 m_scale  ; // scale factor (e.g. sign) 
      /// Integrator
      Ostap::Math::Integrator        m_integrator ; // integrator 
      // ======================================================================
    };
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END  
// ============================================================================
#endif // OSTAP_KRAMERSKRONIG_H
// ============================================================================
