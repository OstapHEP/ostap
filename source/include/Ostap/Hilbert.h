// ============================================================================
#ifndef OSTAP_HILBERT_H 
#define OSTAP_HILBERT_H 1
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
    /** @class Hilbert Ostap/Hilbert.h
     *  Simple class to implement Hilbert transform
     *  https://en.wikipedia.org/wiki/Hilbert_transform
     */
    class Hilbert
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
       *   @param rescale rescale function for better numerical precision 
       *   @param aprecision absolute precision 
       *   @param rprecision relative precision 
       *   @param size    size of integration workspace  
       */
      template <class FUNCTION>
      Hilbert 
      ( FUNCTION             func           ,
	const std::size_t    tag        = 0 ,
	const unsigned short rescale    = 0 ,
	const double         aprecision = 0 ,
	const double         rprecision = 0 ,
	const double         width      = 0 ,           
	const std::size_t    size       = 0 )
        : m_func       ( func       )
        , m_tag        ( tag        ) 
        , m_rescale    ( rescale    )
        , m_aprecision ( aprecision ) 
        , m_rprecision ( rprecision )
        , m_width      ( width      )
        , m_integrator ( size       )
      {} ;
      // ======================================================================
      /** constructor from the function
       *   @param func   the function
       *   @param tag     unique tag/label for cache 
       *   @param rescale rescale function for better numerical precision 
       *   @param aprecision absolute precision 
       *   @param rprecision relative precision 
       *   @param size    size of integration workspace  
       */
      Hilbert 
      ( function1            func           ,
	const std::size_t    tag        = 0 ,
	const unsigned short rescale    = 0 ,
	const double         aprecision = 0 ,
	const double         rprecision = 0 ,
	const double         width      = 0 ,           
	const std::size_t    size       = 0 ) ;
      // ======================================================================      
    public:
      // ======================================================================
      /**  templated creator from the function
       *   @param func   the function
       *   @param tag     unique tag/label for cache 
       *   @param rescale rescale function for better numerical precision 
       *   @param aprecision absolute precision 
       *   @param rprecision relative precision 
       *   @param size    size of integration workspace  
       */
      template <class FUNCTION>
      inline static Hilbert
      create
      ( FUNCTION             func           ,
	const std::size_t    tag        = 0 ,
	const unsigned short rescale    = 0 ,
	const double         aprecision = 0 ,
	const double         rprecision = 0 ,
	const double         width      = 0 ,                     
	const std::size_t    size       = 0 )
      { return Hilbert ( func , tag , rescale , aprecision , rprecision , width , size ) ; }        
      // ======================================================================
    public:
      // ======================================================================
      /// Get the value of Hilbert transform
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
      function1                      m_func       ; // the function
      /// unique tag/label 
      std::size_t                    m_tag        ; // unique tag/label 
      /// rescale function for better numerical precision 
      unsigned short                 m_rescale    ; // #rescale points
      /// absolute precison
      double                         m_aprecision ; // absolute preciison
      /// relative precison
      double                         m_rprecision ; // relative precision
      /// width 
      double                         m_width      ; // width       
      /// Integrator
      Ostap::Math::Integrator        m_integrator ; // integrator 
      // ======================================================================
    }; //                                 The end of class Ostap::Math::Hilbert 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namepsace Ostap
// ============================================================================
//                                                                      The END  
// ============================================================================
#endif // OSTAP_HILBERT_H
// ============================================================================
