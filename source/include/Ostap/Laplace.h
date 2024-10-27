// ============================================================================
#ifndef OSTAP_LAPLACE_H 
#define OSTAP_LAPLACE_H 1
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
    /** @class Laplace Ostap/Laplace.h
     *  Simple class to implement Hilbert transform
     *  https://en.wikipedia.org/wiki/Hilbert_transform
     */
    class Laplace
    {
      // ======================================================================
      /**  templated constructor from the function
       *   @param func   the function
       *   @param tag     unique tag/label for cache 
       *   @param size    size of integration workspace  
       */
      template <class FUNCTION>
        Laplace
        ( FUNCTION             func           ,
          const std::size_t    tag        = 0 ,
          const double         aprecision = 0 ,
          const double         rprecision = 0 ,
          const std::size_t    size       = 0 )
        : m_func       ( func       )
        , m_tag        ( tag        ) 
        , m_aprecision ( aprecision ) 
        , m_rprecision ( rprecision )
        , m_integrator ( size       )
        {} ;
      // ======================================================================
    public:
      // ======================================================================
      /**  templated creator from the function
       *   @param func   the function
       *   @param tag     unique tag/label for cache 
       *   @param size    size of integration workspace  
       */
      template <class FUNCTION>
        inline static Laplace
        create
        ( FUNCTION             func           ,
          const std::size_t    tag        = 0 ,
          const double         aprecision = 0 ,
          const double         rprecision = 0 ,
          const std::size_t    size       = 0 )
          { return Laplace ( func , tag , aprecision , rprecision , size ) ; }        
      // ======================================================================
    public:
      // ======================================================================
      /// Get the value of Laplace transform
      double operator() ( const double x ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of the unction  
      inline double func ( const double x ) const { return m_func ( x ) ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// the function 
      std::function<double(double)>  m_func       ; // the function
      /// unique tag/label 
      std::size_t                    m_tag        ; // unique tag/label 
      /// absolte precison
      double                         m_aprecision ; // absoluet preciison
      /// relative precison
      double                         m_rprecision ; // relative precision
      /// Integrator
      Ostap::Math::Integrator        m_integrator ; // integrator 
      // ======================================================================
    }; //                                 The end of class Ostap::Math::Laplace
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namepsace Ostap
// ============================================================================
//                                                                      The END  
// ============================================================================
#endif // OSTAP_LAPLACE_H
// ============================================================================
