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
     *  Simple class to implement Laplace transform
     *  https://en.wikipedia.org/wiki/Laplace_transform
     */
    class Laplace
    {
      // ======================================================================
    public:
      // ======================================================================
      /// the actual function type 
      typedef std::function<double(double)>                         function1 ;
      // ======================================================================
    public:
      // ======================================================================
      /**  templated constructor from the function
       *   @param func   the function
       *   @param tag     unique tag/label for cache 
       *   @param aprecision absolute precision 
       *   @param rprecision relative precision 
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
      /**  constructor from the function
       *   @param func   the function
       *   @param tag     unique tag/label for cache 
       *   @param aprecision absolute precision 
       *   @param rprecision relative precision 
       *   @param size    size of integration workspace  
       */
      Laplace
      ( function1            func           ,
        const std::size_t    tag        = 0 ,
        const double         aprecision = 0 ,
        const double         rprecision = 0 ,
        const std::size_t    size       = 0 ) ;
      // ======================================================================      
    public:
      // ======================================================================
      /**  templated creator from the function
       *   @param func   the function
       *   @param tag    unique tag/label for cache 
       *   @param aprecision absolute precision 
       *   @param rprecision relative precision 
       *   @param size   size of integration workspace  
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
      /// get the value of the original function  
      inline double func ( const double x ) const { return m_func ( x ) ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// the function 
      function1                      m_func       ; // the function
      /// unique tag/label 
      std::size_t                    m_tag        ; // unique tag/label 
      /// absolute precision
      double                         m_aprecision ; // absolute preciison
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
