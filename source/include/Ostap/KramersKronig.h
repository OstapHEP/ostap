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
     *  Simple class to implement Kramers-Kronig relations 
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
       *   @param tag    unique tag/label for cacheing 
       *   @param rescale rescale function for better numerical precison 
       *   @param size   size of integration workspace  
       */
      template <class FUNCTION>
      KramersKronig 
      ( FUNCTION             rho         ,
        const double         omega0      ,
        const unsigned short n       = 0 ,
        const double         scale   = 1 ,                      
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 , 
        const std::size_t    size    = 0 )
        : m_rho        ( rho     )
        , m_omega0     ( omega0  )
        , m_n          ( n       )
        , m_scale      ( scale   )          
        , m_tag        ( tag     ) 
        , m_rescale    ( rescale ) 
        , m_integrator ( size    )
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
       *   @param tag    unique tag/label for cacheing 
       *   @param size   size of integration workspace  
       */
      template <class FUNCTION>
      inline static KramersKronig
      create
      ( FUNCTION             rho         ,
        const double         omega0      ,
        const unsigned short n       = 0 ,
        const double         scale   = 1 , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ,
        const std::size_t    size    = 0 )
      { return KramersKronig ( rho , omega0 , n , scale , tag , rescale , size ) ; }
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
      double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ====================================================================== 
      /// get the value of \f$ \varrho \f$ function  
      double         rho     ( const double x ) const { return m_rho ( x ) ; }
      /// number of subtractions 
      unsigned short n       () const { return m_n      ; }
      /// scale  factor 
      double         scale   () const { return m_scale  ; }
      /// low integrtaion edge 
      double         lowEdge () const { return m_omega0 ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the function 
      std::function<double(double)>  m_rho        ; // the function
      /// the low integration eddge 
      double                         m_omega0     ; // low integration edge 
      /// number of subtractions
      unsigned short                 m_n          ; // number of subtractions
      /// scale factor (e.g. sign) 
      double                         m_scale      ; // scale factor (e.g. sign) 
      /// unique tag/label 
      std::size_t                    m_tag        ; // unique tag/label 
      /// rescale fnuction for better numerical precision 
      unsigned short                 m_rescale    ; // #rescale points 
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
