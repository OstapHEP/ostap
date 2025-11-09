// ============================================================================
#ifndef OSTAP_INTEGRATOR_H 
#define OSTAP_INTEGRATOR_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <functional>
#include <utility>
#include <string>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class Integrator Ostap/Integrator.h 
     *  simple numerical integrator for 1D&,2D&3D-cases 
     * 
     *  It constains method for 
     *   - regular numeric integration using GAG from GSL 
     *   - numeric integration for infinite and semi-infinite intervals
     *   - Cauchy principal value integrals, inclusing semi-infinite intervals  
     *   - Kramers-Kronig dispersive integrals (including subtractions)  
     *   - integration of function with singular values 
     *  
     *  Methods above are based on GSL adaptive numerical integration routines 
     *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
     *  @see https://www.gnu.org/software/gsl/doc/html/integration.html#qag-adaptive-integration
     *  @see https://www.gnu.org/software/gsl/doc/html/integration.html#qagi-adaptive-integration-on-infinite-intervals 
     *  @see https://www.gnu.org/software/gsl/doc/html/integration.html#qawc-adaptive-integration-for-cauchy-principal-values
     *  @see https://www.gnu.org/software/gsl/doc/html/integration.html#qaws-adaptive-integration-for-singular-functions
     *
     *  In addition there is explict integrator using 
     *  using double adaptive CQUAD algorithm 
     *  @see https://www.gnu.org/software/gsl/doc/html/integration.html#cquad-doubly-adaptive-integration
     *  and numercial integration using Romberg algorithm 
     *  @https://www.gnu.org/software/gsl/doc/html/integration.html#romberg-integration
     *
     *  Also one has integrators for 2D and 3D functions:
     *  - integration for 2D functions using Genz-Malik's cubature method 
     *  - partial integration for 2D functions 1D-integration methods 
     *  - integration for 3D functions using Genz-Malik's cubature method 
     *  - partial integration for 3D functions using 1D and 2D integration methods 
     * 
     *
     *  The interface  contains following parts:  
     *  - template methods for integration of 1D-functions                           
     *    (Workspace is provided by the Integrator itself)
     *  - template methods for integration of 2D-functions                           
     *    (no workspace is needed) 
     *  - template methods for 1D (partial) integration of 2D-functions             
     *    (Workspace is provided by the Integrator itself)   
     *  - template methods for integration of 3D-functions                           
     *    (no workspace is needed)
     *  - template methods for 2D (partial) integration of 3D-functions              
     *    (no workspace is needed)             
     *  - template methods for 1D (partial) integration of 3D-functions              
     *    (Workspace is provided by the Integrator itself)   
     *  - template static methods for integratiro of 1D-functions                    
     *    (Workspace needs to be provdeid by user)
     *  - static   methods for integration of <code>std::function<double(double)></code>
     *    (Workspace needs to be provided by user)
     *  - static   methods for 2D-integratio of <code>std::function<double(double,double)></code>     
     *    (no workspace is needed) 
     *  - static   methods for 1D (partial) integration of <code>std::function<double(double,double)></code>  
     *    (Workspace needs to be provided by user)
     *  - static   methods for integration of <code>std::function<double(double,double,double)></code> 
     *    (no workspace is needed)
     *  - static   methods for 2D (partial) integration of <code>std::function<double(double,double,double)></code>  
     *    (no workspace is needed)
     *  - static   methods for 1D (partial) integration of <code>std::function<double(double,double,double)></code>  
     *    (Workspace needs to be provided by user)
     */
    class Integrator 
    {
    public:
      // ======================================================================
      typedef std::function<double(double)>               function1 ;
      typedef std::function<double(double,double)>        function2 ;
      typedef std::function<double(double,double,double)> function3 ;
      typedef std::pair<double,double>                    result    ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with integration workspace size
      Integrator
      ( const std::size_t    size         = 0 , 
        const unsigned short size_cquad   = 0 ,
        const unsigned short size_romberg = 0 ) ;
      // ======================================================================
      /// constructor with integration workspace
      Integrator
      ( const Ostap::Math::WorkSpace& ws        ) ;
      // ======================================================================
      /// constructor with the fictive name 
      Integrator ( const std::string& /* name */ ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param  xmin lower  integration edge 
       *  @param  xmax upper  integration edge
       *  @param  tag  unique label/tag  
       *  @param  rescale rescale function for better numerical precision
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision 
       *  @param  rule      interation rule
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double integrate
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const 
      { return integrate_
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 < rule       ? rule       : m_qag_rule          ).first ; }
      // ======================================================================
      /** calculate the integral with uncertainty
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower  integration edge 
       *  @param xmax upper  integration edge
       *  @param  tag  unique label/tag  
       *  @param  rescale rescale function for better numerical precision
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision 
       *  @param  rule      interation rule
       *  @return the value of the integral & uncertainty 
       */
      template <class FUNCTION1>
      inline double integrate_err 
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const 
      { return integrate_
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 < rule       ? rule       : m_qag_rule          ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param tag unique tag  (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double integrate_infinity
      ( FUNCTION1         f1             , 
        const std::size_t tag        = 0 , 
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_infinity_ 
          ( std::cref    ( f1 ) , 
            m_workspace , tag   , 
            0 < aprecision ? aprecision : m_abs_precision_qagi , 
            0 < rprecision ? rprecision : m_rel_precision_qagi ).first  ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param tag unique tag  (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result integrate_infinity_err 
      ( FUNCTION1         f1             , 
        const std::size_t tag        = 0 , 
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_infinity_ 
          ( std::cref    ( f1 ) , 
            m_workspace , tag   , 
            0 < aprecision ? aprecision : m_abs_precision_qagi , 
            0 < rprecision ? rprecision : m_rel_precision_qagi ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param xmin lower integration edge
       *  @param tag unique tag (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double integrate_to_infinity
      ( FUNCTION1         f1             , 
        const double      xmin           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_to_infinity_
          ( std::cref ( f1 ) , xmin , 
            m_workspace      , tag  ,
            0 < aprecision ? aprecision : m_abs_precision_qagiu , 
            0 < rprecision ? rprecision : m_rel_precision_qagiu ).first ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param xmin lower integration edge
       *  @param tag unique tag (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result integrate_to_infinity_err
      ( FUNCTION1         f1             , 
        const double      xmin           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_to_infinity_
          ( std::cref ( f1 ) , xmin , 
            m_workspace      , tag  ,
            0 < aprecision ? aprecision : m_abs_precision_qagiu , 
            0 < rprecision ? rprecision : m_rel_precision_qagiu ) ; }
      // ========================================================================
    public: 
      // ========================================================================
      /** calculate the integral 
       *  \f[ r = \mathcal{P}\int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @param tag     unique label/tag 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double integrate_from_infinity
      ( FUNCTION1         f1             , 
        const double      xmax           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const
      { return integrate_from_infinity_ 
          ( std::cref ( f1 ) , xmax , 
            m_workspace , tag       , 
            0 < aprecision ? aprecision : m_abs_precision_qagil , 
            0 < rprecision ? rprecision : m_rel_precision_qagil ).first ; }
      // ========================================================================
      /** calculate the integral 
       *  \f[ r = \mathcal{P}\int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @param tag     unique label/tag 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result integrate_from_infinity_err
      ( FUNCTION1         f1             , 
        const double      xmax           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const
      { return integrate_from_infinity_ 
          ( std::cref ( f1 ) , xmax , 
            m_workspace , tag       , 
            0 < aprecision ? aprecision : m_abs_precision_qagil , 
            0 < rprecision ? rprecision : m_rel_precision_qagil ) ; }
      // ======================================================================
    public : 
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1      the function 
       *  @param c       the parameter 
       *  @param xmin    lower integration edge 
       *  @param xmax    upper integration edge 
       *  @param tag     unique label/tag       
       *  @param rescale rescale function for better numerical treatmenrt
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision               
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double cauchy_pv 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           , 
        const double         xmax           , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const  
      { return cauchy_pv_ 
          ( std::cref ( f1 ) , c , xmin , xmax , 
            m_workspace , tag , rescale  , 
            0 < aprecision ? aprecision : m_abs_precision_cpv , 
            0 < rprecision ? rprecision : m_rel_precision_cpv ).first  ; }

      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1      the function 
       *  @param c       the parameter 
       *  @param xmin    lower integration edge 
       *  @param xmax    upper integration edge 
       *  @param tag     unique label/tag       
       *  @param rescale rescale function for better numerical treatmenrt
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision               
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result cauchy_pv_err 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           , 
        const double         xmax           , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const  
      { return cauchy_pv_ 
          ( std::cref ( f1 ) , c , xmin , xmax , 
            m_workspace , tag , rescale  , 
            0 < aprecision ? aprecision : m_abs_precision_cpv , 
            0 < rprecision ? rprecision : m_rel_precision_cpv ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision               
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double cauchy_pv_to_infinity  
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           ,  
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const
      { return cauchy_pv_to_infinity_
          ( std::cref ( f1 ) , c , xmin , 
            m_workspace , tag , rescale ,
            0 < aprecision ? aprecision : m_abs_precision_cpvi , 
            0 < rprecision ? rprecision : m_rel_precision_cpvi , width ).first ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision               
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result cauchy_pv_to_infinity_err   
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           ,  
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const
      { return cauchy_pv_to_infinity_
          ( std::cref ( f1 ) , c , xmin , 
            m_workspace , tag , rescale ,
            0 < aprecision ? aprecision : m_abs_precision_cpvi , 
            0 < rprecision ? rprecision : m_rel_precision_cpvi , width ) ; }
      // ======================================================================
    public : 
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *
       *  @param  f1 the function 
       *  @param  c  the parameter 
       *  @param  xmax upper integration edge 
       *  @param  tag     unique label/tag 
       *  @param  rescale rescale function for better numerical treatmenrt 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision               
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double cauchy_pv_from_infinity  
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmax           , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const
      { return cauchy_pv_from_infinity_ 
          ( std::cref ( f1 ) , c , xmax , 
            m_workspace , tag , rescale , 
            0 < aprecision ? aprecision : m_abs_precision_cpvi , 
            0 < rprecision ? rprecision : m_rel_precision_cpvi , width ).first ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *
       *  @param  f1 the function 
       *  @param  c  the parameter 
       *  @param  xmax upper integration edge 
       *  @param  tag     unique label/tag 
       *  @param  rescale rescale function for better numerical treatmenrt 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision               
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result cauchy_pv_from_infinity_err 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmax           , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const
      { return cauchy_pv_from_infinity_ 
          ( std::cref ( f1 ) , c , xmax , 
            m_workspace , tag , rescale , 
            0 < aprecision ? aprecision : m_abs_precision_cpvi , 
            0 < rprecision ? rprecision : m_rel_precision_cpvi , width ) ; }
      // ======================================================================
    public: 
      // ======================================================================      
      /** get Cauchy principal value integral for the infinite range 
       *  \f[ g(c) = \mathcal{P} \int_{-\infty}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *
       *  Integral is calculated as sum of three components
       *   - \f$ \int_{-\infty}^{a}\frac{f(x)}{x-c}dx \f$
       *   - \f$ \mathcal{P} \int_{a}^{b}\frac{f(x)}{x-c}dx \f$
       *   - \f$ int_{b}^{+\infty}\frac{f(x)}{x-c}dx \f$
       *
       *  @param  f1      the function 
       *  @param  c       the parameter 
       *  @param  ws      integration workspace 
       *  @param  tag     unique label/tag 
       *  @param  rescale rescale function for better numerical precision  
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CPVI is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CPVI is used) 
       *  @param  width   width parameter
       *  @return value of the integral and the estimate of the uncertainty
       *
       *  @see Ostap::Math::Integrator::cauchy_pv
       *  @see Ostap::Math::Integrator::integrate_to_infinity 
       *  @see Ostap::Math::Integrator::integrate_from_infinity 
       *  @see Ostap::Math::Integrator::cauchy_pv_b
       *  @see Ostap::Math::Integrator::cauchy_pv_b
       */
      template <class FUNCTION1>
      inline double cauchy_pv_infinity 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const 
      { return cauchy_pv_infinity_ 
          ( std::cref ( f1 ) , c , 
            m_workspace , tag , rescale , 
            0 < aprecision ? aprecision : m_abs_precision_cpvi , 
            0 < rprecision ? rprecision : m_rel_precision_cpvi , width ).first ; }

      // ======================================================================      
      /** get Cauchy principal value integral for the infinite range 
       *  \f[ g(c) = \mathcal{P} \int_{-\infty}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *
       *  Integral is calculated as sum of three components
       *   - \f$ \int_{-\infty}^{a}\frac{f(x)}{x-c}dx \f$
       *   - \f$ \mathcal{P} \int_{a}^{b}\frac{f(x)}{x-c}dx \f$
       *   - \f$ int_{b}^{+\infty}\frac{f(x)}{x-c}dx \f$
       *
       *  @param  f1      the function 
       *  @param  c       the parameter 
       *  @param  ws      integration workspace 
       *  @param  tag     unique label/tag 
       *  @param  rescale rescale function for better numerical precision  
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CPVI is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CPVI is used) 
       *  @param  width   width parameter
       *  @return value of the integral and the estimate of the uncertainty
       *
       *  @see Ostap::Math::Integrator::cauchy_pv
       *  @see Ostap::Math::Integrator::integrate_to_infinity 
       *  @see Ostap::Math::Integrator::integrate_from_infinity 
       *  @see Ostap::Math::Integrator::cauchy_pv_b
       *  @see Ostap::Math::Integrator::cauchy_pv_b
       */
      template <class FUNCTION1>
      inline result cauchy_pv_infinity_err 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const 
      { return cauchy_pv_infinity_ 
          ( std::cref ( f1 ) , c , 
            m_workspace , tag , rescale , 
            0 < aprecision ? aprecision : m_abs_precision_cpvi , 
            0 < rprecision ? rprecision : m_rel_precision_cpvi , width ) ; }
      // ======================================================================
    public :
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi}
       *     \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin low integration edge 
       *  @param n   number of subtractions 
       *  @see Ostap::Math::Integrtator::cauchy_pv_to_infinity 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_KK is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_KK is used)
       *  @param width the width 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double kramers_kronig 
      ( FUNCTION1            f1             , 
        const double         s              , 
        const double         xmin           , 
        const unsigned short n          = 0 , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const
      { return kramers_kronig_ 
          ( std::cref ( f1 ) , s , xmin , n , 
            m_workspace , tag , rescale ,
            0 < aprecision ? aprecision : m_abs_precision_kk , 
            0 < rprecision ? rprecision : m_rel_precision_kk , width ).first  ; }
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi}
       *     \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin low integration edge 
       *  @param n   number of subtractions 
       *  @see Ostap::Math::Integrtator::cauchy_pv_to_infinity 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_KK is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_KK is used)
       *  @param width the width 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result kramers_kronig_err  
      ( FUNCTION1            f1             , 
        const double         s              , 
        const double         xmin           , 
        const unsigned short n          = 0 , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) const
      { return kramers_kronig_ 
          ( std::cref ( f1 ) , s , xmin , n , 
            m_workspace , tag , rescale ,
            0 < aprecision ? aprecision : m_abs_precision_kk , 
            0 < rprecision ? rprecision : m_rel_precision_kk , width ) ; }
      // ======================================================================
    public :
      // ======================================================================      
      /** integration with known singular points 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
       *  @param  f1 the   function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       *  @param points known singular points 
       *  @param tag unqye tag/label 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAGP is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAGP is used)
       *  
       *  - Only singular poins between \f$ x_{min} \f$  and \f$ x_{max} \f$ are considered 
       *  - \f$ x_{min} \f$  and \f$ x_{max \f$ are also  considered as singular points 
       */
      template <class FUNCTION1>
      inline double integrate_singular
      ( FUNCTION1                  f1             , 
        const double               xmin           , 
        const double               xmax           ,
        const std::vector<double>& points         ,
        const std::size_t          tag        = 0 ,
        const double               aprecision = 0 , 
        const double               rprecision = 0 ) const
      { return integrate_singular_
          ( std::cref ( f1 ) , 
            xmin , xmax , points , 
            m_workspace , tag    , 
            0 < aprecision ? aprecision : m_abs_precision_qagp , 
            0 < rprecision ? rprecision : m_rel_precision_qagp ).first ; }
      // ======================================================================      
      /** integration with known singular points 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
       *  @param  f1 the   function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       *  @param points known singular points 
       *  @param tag unqye tag/label 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAGP is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAGP is used)
       *  
       *  - Only singular poins between \f$ x_{min} \f$  and \f$ x_{max} \f$ are considered 
       *  - \f$ x_{min} \f$  and \f$ x_{max \f$ are also  considered as singular points 
       */
      template <class FUNCTION1>
      inline result integrate_singular_err
      ( FUNCTION1                  f1             , 
        const double               xmin           , 
        const double               xmax           ,
        const std::vector<double>& points         ,
        const std::size_t          tag        = 0 ,
        const double               aprecision = 0 , 
        const double               rprecision = 0 ) const
      { return integrate_singular_
          ( std::cref ( f1 ) , 
            xmin , xmax , points , 
            m_workspace , tag    , 
            0 < aprecision ? aprecision : m_abs_precision_qagp , 
            0 < rprecision ? rprecision : m_rel_precision_qagp ) ; }      
      // ======================================================================
    public:  // specific: use double-adaptive CQUAD integrator
      // ======================================================================
      /** calculate the integral usnig double adaptive CQUAD integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CQUAD is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CQUAD is used)
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double integrate_cquad
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_cquad_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_cquad , 
            0 < rprecision ? rprecision : m_rel_precision_cquad ).first ; }

      // ======================================================================
      /** calculate the integral usnig double adaptive CQUAD integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CQUAD is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CQUAD is used)
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result integrate_cquad_err
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_cquad_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_cquad , 
            0 < rprecision ? rprecision : m_rel_precision_cquad ) ; }      
      // ======================================================================
    public:  // specific: use Romberg integrator 
      // ======================================================================
      /** calculate the integral using Romberg integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_ROMBERG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_ROMBERG is used)
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline double integrate_romberg
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_romberg_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_romberg , 
            0 < rprecision ? rprecision : m_rel_precision_romberg ).first ; }
      // ======================================================================
      /** calculate the integral using Romberg integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_ROMBERG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_ROMBERG is used)
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      inline result integrate_romberg_err
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_romberg_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_romberg , 
            0 < rprecision ? rprecision : m_rel_precision_romberg ) ; }      
      // ======================================================================
    public : // couple of explicit specializations
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower  integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used)
       *  @param rule 
       *  @return the value of the integral 
       */
      inline double integrate
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const 
      { return integrate_
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 < rule       ? rule       : m_qag_rule          ).first ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param  f1 the function 
       *  @param  xmin lower  integration edge 
       *  @param  xmax upper  integration edge
       *  @param  tag  uqniue label/tag  
       *  @param  rescale rescale function for better numerical precision  
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_QAG is used)
       *  @param  rule
       *  @return the value of the integral 
       */
      inline result integrate_err
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const 
      { return integrate_
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 < rule       ? rule       : m_qag_rule          ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param tag unique tag  (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      inline double integrate_infinity
      ( function1         f1             , 
        const std::size_t tag        = 0 , 
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_infinity_ 
          ( std::cref    ( f1 ) , 
            m_workspace , tag   , 
            0 < aprecision ? aprecision : m_abs_precision_qagi , 
            0 < rprecision ? rprecision : m_rel_precision_qagi ).first  ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param tag unique tag  (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      inline result integrate_infinity_err
      ( function1         f1             , 
        const std::size_t tag        = 0 , 
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_infinity_ 
          ( std::cref    ( f1 ) , 
            m_workspace , tag   , 
            0 < aprecision ? aprecision : m_abs_precision_qagi , 
            0 < rprecision ? rprecision : m_rel_precision_qagi ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param xmin lower integration edge
       *  @param tag unique tag (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      inline double integrate_to_infinity
      ( function1         f1             , 
        const double      xmin           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_to_infinity_
          ( std::cref ( f1 ) , xmin , 
            m_workspace      , tag  ,
            0 < aprecision ? aprecision : m_abs_precision_qagiu , 
            0 < rprecision ? rprecision : m_rel_precision_qagiu ).first ; } 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function
       *  @param xmin lower integration edge
       *  @param tag unique tag (for cache)
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      inline result integrate_to_infinity_err 
      ( function1         f1             , 
        const double      xmin           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_to_infinity_
          ( std::cref ( f1 ) , xmin , 
            m_workspace      , tag  ,
            0 < aprecision ? aprecision : m_abs_precision_qagiu , 
            0 < rprecision ? rprecision : m_rel_precision_qagiu ) ; } 
      // ======================================================================
    public:
      // ======================================================================            
      /** calculate the integral 
       *  \f[ r = \mathcal{P}\int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @param tag     unique label/tag 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      inline double integrate_from_infinity
      ( function1         f1             , 
        const double      xmax           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const
      { return integrate_from_infinity_ 
          ( std::cref ( f1 ) , xmax , 
            m_workspace , tag       , 
            0 < aprecision ? aprecision : m_abs_precision_qagil , 
            0 < rprecision ? rprecision : m_rel_precision_qagil ).first ; }
      // ========================================================================
      /** calculate the integral 
       *  \f[ r = \mathcal{P}\int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @param tag     unique label/tag 
       *  @param  aprecision absolute precision 
       *  @param  rprecision relative precision        
       *  @return the value of the integral 
       */
      inline result integrate_from_infinity_err 
      ( function1         f1             , 
        const double      xmax           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const
      { return integrate_from_infinity_ 
          ( std::cref ( f1 ) , xmax , 
            m_workspace , tag       , 
            0 < aprecision ? aprecision : m_abs_precision_qagil , 
            0 < rprecision ? rprecision : m_rel_precision_qagil ) ; }
      // ========================================================================
    public:
      // ======================================================================            
      /** calculate the integral using double adaptive CQUAD integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CQUAD is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CQUAD is used)
       *  @return the value of the integral 
       */
      inline double integrate_cquad
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_cquad_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_cquad , 
            0 < rprecision ? rprecision : m_rel_precision_cquad ).first ; }
      // ======================================================================      
      /** calculate the integral using double adaptive CQUAD integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CQUAD is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CQUAD is used)
       *  @return the value of the integral 
       */
      inline result integrate_cquad_err 
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_cquad_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_cquad , 
            0 < rprecision ? rprecision : m_rel_precision_cquad ) ; }      
      // ======================================================================
    public :
      // ======================================================================      
      /** calculate the integral using Romberg integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_ROMBERG is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_ROMBERG is used)
       *  @return the value of the integral 
       */
      inline double integrate_romberg
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_romberg_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_romberg , 
            0 < rprecision ? rprecision : m_rel_precision_romberg ).first ; }
      // =======================================================================
      /** calculate the integral using Romberg integrator  
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_ROMBERG is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_ROMBERG is used)
       *  @return the value of the integral 
       */
      inline result integrate_romberg_err 
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_romberg_ 
          ( std::cref ( f1 )    , 
            xmin , xmax ,
            m_workspace , 
            tag  , rescale      , 
            0 < aprecision ? aprecision : m_abs_precision_romberg , 
            0 < rprecision ? rprecision : m_rel_precision_romberg ) ; }
      // =======================================================================
    public: // Integrations for 2D-functions 
      // ======================================================================
      // 2D integration 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param  xmin lower integration edge in x 
       *  @param  xmax upper integration edge in x 
       *  @param  ymin lower integration edge in y 
       *  @param  ymax upper integration edge in y
       *  @param  tag unique ta g (for cache)
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CUBE2 is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CUBE2 is used)       
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      inline double integrate2
      ( FUNCTION2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 ,
        const double      rprecision = 0 ) const 
      { return integrate2_
          ( std::cref ( f2 ) , 
            xmin , xmax , 
            ymin , ymax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube2 , 
            0 < rprecision ? rprecision : m_rel_precision_cube2 ).first ; }
      // ==================================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param  xmin lower integration edge in x 
       *  @param  xmax upper integration edge in x 
       *  @param  ymin lower integration edge in y 
       *  @param  ymax upper integration edge in y
       *  @param  tag unique ta g (for cache)
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CUBE2 is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CUBE2 is used)       
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      inline result integrate2_err 
      ( FUNCTION2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 ,
        const double      rprecision = 0 ) const 
      { return integrate2_
          ( std::cref ( f2 ) , 
            xmin , xmax , 
            ymin , ymax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube2 , 
            0 < rprecision ? rprecision : m_rel_precision_cube2 ) ; }
      // ======================================================================
    public : 
      // ======================================================================
      // partial integration for 2D-functions 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
       *  @param f2 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      double integrate2X
      ( FUNCTION2            f2             , 
        const double         y              , 
        const double         xmin           ,
        const double         xmax           ,  
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,        
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const
      { return integrate2X_ 
          ( std::cref ( f2 )     , 
            y , xmin , xmax      , 
            m_workspace , tag    , rescale , 
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 != rule ? rule : m_qag_rule     ).first  ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
       *  @param f2 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      double integrate2Y
      ( FUNCTION2            f2             , 
        const double         x              , 
        const double         ymin           ,
        const double         ymax           ,
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,        
        const double         aprecision = 0 , 
        const double         rprecision = 0 ,
        const int            rule       = 0 ) const
      { return integrate2Y_ 
          ( std::cref ( f2 )     , 
            x , ymin , ymax      , 
            m_workspace , tag    , rescale , 
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 != rule ? rule : m_qag_rule     ).first  ; }
      // ======================================================================
    public: // 3D integration 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}
       *          \int_{y_{min}}^{y_{max}}
       *          \int_{z_{min}}^{z_{max}}f_3(x,y,z) dx dy dz\f]
       *  @param f3 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE3 is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE3 is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      template <class FUNCTION3>
      inline double integrate3
      ( FUNCTION3         f3             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const double      zmin           , 
        const double      zmax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate3_
          ( std::cref ( f3 ) , 
            xmin , xmax ,
            ymin , ymax , 
            zmin , zmax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube3 , 
            0 < rprecision ? rprecision : m_rel_precision_cube3 ).first ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}
       *          \int_{y_{min}}^{y_{max}}
       *          \int_{z_{min}}^{z_{max}}f_3(x,y,z) dx dy dz\f]
       *  @param f3 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE3 is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE3 is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      template <class FUNCTION3>
      inline result integrate3_err 
      ( FUNCTION3         f3             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const double      zmin           , 
        const double      zmax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate3_
          ( std::cref ( f3 ) , 
            xmin , xmax ,
            ymin , ymax , 
            zmin , zmax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube3 , 
            0 < rprecision ? rprecision : m_rel_precision_cube3 ) ; }
      // ======================================================================      
    public:
      // ======================================================================
      // partical integration for 3D-functions 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(z) = \int_{x_{min}}^{x_{max}}
       *             \int_{y_{min}}^{y_{max}} f_3(x,y,z) dxdy \f]
       *  @param f3 the function 
       *  @param z parameter z
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
       *  @return the value of the integral and the estimate of the error 
       *  @see Ostap::Math::Integrator::integrate2_ 
       */
      template <class FUNCTION3>
      double integrate3XY
      ( FUNCTION3            f3             ,
        const double         z              , 
        const double         xmin           ,
        const double         xmax           ,
        const double         ymin           ,
        const double         ymax           ,
        const std::size_t    tag        = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 )  const 
      { return integrate3XY_ 
          ( std::cref ( f3 ) , 
            z    , 
            xmin , xmax ,
            ymin , ymax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube2 , 
            0 < rprecision ? rprecision : m_rel_precision_cube2 ).first ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(y) = \int_{x_{min}}^{x_{max}}
       *             \int_{z_{min}}^{z_{max}} f_3(x,y,z) dxdz \f]
       *  @param f3 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
       *  @return the value of the integral and the estimate of the error 
       *  @see Ostap::Math::Integrator::integrate2_ 
       */
      template <class FUNCTION3>
      double integrate3XZ
      ( FUNCTION3            f3             ,
        const double         y              , 
        const double         xmin           ,
        const double         xmax           ,
        const double         zmin           ,
        const double         zmax           ,
        const std::size_t    tag        = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate3XZ_ 
          ( std::cref ( f3 ) , 
            y    , 
            xmin , xmax ,
            zmin , zmax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube2 , 
            0 < rprecision ? rprecision : m_rel_precision_cube2 ).first ; } 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(x) = \int_{y_{min}}^{y_{max}}
       *             \int_{z_{min}}^{z_{max}} f_3(x,y,z) dydz \f]
       *  @param f3 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
       *  @return the value of the integral and the estimate of the error 
       *  @see Ostap::Math::Integrator::integrate2_ 
       */
      template <class FUNCTION3>
      double integrate3YZ
      ( FUNCTION3            f3             ,
        const double         x              , 
        const double         ymin           ,
        const double         ymax           ,
        const double         zmin           ,
        const double         zmax           ,
        const std::size_t    tag        = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const
      { return integrate3YZ_ 
          ( std::cref ( f3 ) , 
            x    , 
            ymin , ymax ,
            zmin , zmax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube2 , 
            0 < rprecision ? rprecision : m_rel_precision_cube2 ).first ; } 
      // ==========================================================================
      /** calculate the integral 
       *  \f[ r(y,z) = \int_{x_{min}}^{x_{max}}f_3(x,y,z) dx  \f]
       *  @param f3 the function 
       *  @param y parameter y
       *  @param z parameter z
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      template <class FUNCTION3>
      double integrate3X
      ( FUNCTION3            f3             ,
        const double         y              , 
        const double         z              , 
        const double         xmin           ,
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const
      { return integrate3X_ 
          ( std::cref ( f3 ) , 
            y     , z    , 
            xmin  , xmax ,
            ws () , tag  , rescale ,  
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 != rule ? rule : m_qag_rule     ).first  ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(x,z) = \int_{y_{min}}^{y_{max}}f_3(x,y,z) dy  \f]
       *  @param f3 the function 
       *  @param x parameter x
       *  @param z parameter z
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      template <class FUNCTION3>
      double integrate3Y
      ( FUNCTION3            f3             ,
        const double         x              , 
        const double         z              , 
        const double         ymin           ,
        const double         ymax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const 
      { return integrate3Y_ 
          ( std::cref ( f3 ) , 
            x     , z    , 
            ymin  , ymax ,
            ws () , tag  , rescale ,  
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 != rule ? rule : m_qag_rule     ).first  ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(x,y) = \int_{z_{min}}^{z_{max}}f_3(x,y,z) dyz \f]
       *  @param f3 the function 
       *  @param x parameter x
       *  @param y parameter y
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      template <class FUNCTION3>
      double integrate3Z
      ( FUNCTION3            f3             ,
        const double         x              , 
        const double         y              , 
        const double         zmin           ,
        const double         zmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) const 
      { return integrate3Z_ 
          ( std::cref ( f3 ) , 
            x     , y    , 
            zmin  , zmax ,
            ws () , tag  , rescale ,  
            0 < aprecision ? aprecision : m_abs_precision_qag , 
            0 < rprecision ? rprecision : m_rel_precision_qag , 
            0 != rule ? rule : m_qag_rule     ).first  ; }
      // ======================================================================
    public: // set of speciailzations 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param  xmin lower integration edge in x 
       *  @param  xmax upper integration edge in x 
       *  @param  ymin lower integration edge in y 
       *  @param  ymax upper integration edge in y
       *  @param  tag unique ta g (for cache)
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CUBE2 is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CUBE2 is used)       
       *  @return the value of the integral 
       */
      inline double integrate2
      ( function2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 ,
        const double      rprecision = 0 ) const 
      { return integrate2_
          ( std::cref ( f2 ) , 
            xmin , xmax , 
            ymin , ymax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube2 , 
            0 < rprecision ? rprecision : m_rel_precision_cube2 ).first ; }
      // ======================================================================      
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param  xmin lower integration edge in x 
       *  @param  xmax upper integration edge in x 
       *  @param  ymin lower integration edge in y 
       *  @param  ymax upper integration edge in y
       *  @param  tag unique ta g (for cache)
       *  @param  aprecision absolute precision  (if non-positive s_APRECISION_CUBE2 is used) 
       *  @param  aprecision relative precision  (if non-positive s_RPRECISION_CUBE2 is used)       
       *  @return the value of the integral 
       */
      inline result integrate2_err 
      ( function2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 ,
        const double      rprecision = 0 ) const 
      { return integrate2_
          ( std::cref ( f2 ) , 
            xmin , xmax , 
            ymin , ymax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube2 , 
            0 < rprecision ? rprecision : m_rel_precision_cube2 ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}
       *          \int_{y_{min}}^{y_{max}}
       *          \int_{z_{min}}^{z_{max}}f_3(x,y,z) dx dy dz\f]
       *  @param f3 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE3 is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE3 is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      inline double integrate3
      ( function3         f3             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const double      zmin           , 
        const double      zmax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate3_
          ( std::cref ( f3 ) , 
            xmin , xmax ,
            ymin , ymax , 
            zmin , zmax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube3 , 
            0 < rprecision ? rprecision : m_rel_precision_cube3 ).first ; }      
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}
       *          \int_{y_{min}}^{y_{max}}
       *          \int_{z_{min}}^{z_{max}}f_3(x,y,z) dx dy dz\f]
       *  @param f3 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE3 is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE3 is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      inline result integrate3_err 
      ( function3         f3             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const double      zmin           , 
        const double      zmax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate3_
          ( std::cref ( f3 ) , 
            xmin , xmax ,
            ymin , ymax , 
            zmin , zmax , 
            tag  , 
            0 < aprecision ? aprecision : m_abs_precision_cube3 , 
            0 < rprecision ? rprecision : m_rel_precision_cube3 ) ; }
      // ======================================================================      
    public: // SET OF STATIC PUBLIC METHODS With explicit WorksSpace)
      // =====================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param ws   integration workspace 
       *  @param tag   unique tag/label 
       *  @param rescale rescale function for better numerical precision  
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      static inline double integrate
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 )
      { return integrate_
          ( std::cref ( f1 ) , 
            xmin , xmax ,  
            ws   , tag  , rescale , 
            aprecision  , 
            rprecision  , 
            rule        ).first ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param ws integration workspace 
       *  @param tag   unique tag/label 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      static inline double integrate_infinity
      ( FUNCTION1         f1             , 
        const WorkSpace&  ws             , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 )
      { return integrate_infinity_ 
          ( std::cref ( f1 ) , 
            ws , tag   , 
            aprecision , 
            rprecision ).first ; }        
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      static inline double integrate_to_infinity
      ( FUNCTION1         f1             , 
        const double      xmin           , 
        const WorkSpace&  ws             ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 )
      { return integrate_to_infinity_ 
          ( std::cref ( f1 ) , xmin , 
            ws , tag   , 
            aprecision , 
            rprecision ).first ; }        
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      static inline double integrate_from_infinity
      ( function1        f1              , 
        const double     xmax            , 
        const WorkSpace& ws              ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 )
      { return integrate_from_infinity_ 
          ( std::cref ( f1 ) , xmax , 
            ws , tag   , 
            aprecision , 
            rprecision ).first ; }        
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1      the function 
       *  @param c       the parameter 
       *  @param xmin    lower integration edge 
       *  @param xmax    upper  integration edge 
       *  @param ws      integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       */
      template <class FUNCTION1>
      static inline double cauchy_pv 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) 
      { return cauchy_pv_
          ( std::cref ( f1 ) , c , xmin , xmax , 
            ws , tag , rescale , 
            aprecision , 
            rprecision ).first ; }   
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param ws      integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       */
      template <class FUNCTION1>
      static inline double cauchy_pv_to_infinity 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           , 
        const WorkSpace&     ws             ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 )
      { return cauchy_pv_to_infinity_ 
          ( std::cref ( f1 ) , c , xmin , 
            ws , tag , rescale ,
            aprecision , 
            rprecision , width ).first ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       *  @param integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       */
      template <class FUNCTION1>
      static inline double cauchy_pv_from_infinity 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ,
        const double         width      = 0 )
      { return cauchy_pv_from_infinity_ 
          ( std::cref ( f1 ) , c , xmax , 
            ws , tag , rescale , 
            aprecision , 
            rprecision , width ).first ; }
      // ======================================================================
      /** get Cauchy principal value integral for the infinite range 
       *  \f[ g(c) = \mathcal{P} \int_{-\infty}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *
       *  Integral is calculated as sum of three components
       *   - \f$ \int_{-\infty}^{a}\frac{f(x)}{x-c}dx \f$
       *   - \f$ \mathcal{P} \int_{a}^{b}\frac{f(x)}{x-c}dx \f$
       *   - \f$ int_{b}^{+\infty}\frac{f(x)}{x-c}dx \f$
       *
       *  where a and b are chosen such \f$ a < c < b \f$ and 
       *  - for \f$ 0<w \f$ one has 
       *   \f[ \left( \begin{array}{c}a   \\ b  \end{array}\right) = 
       *       \left( \begin{array}{c}c-w \\ c+w\end{array}\right) \f]
       *  - for \f$ w \le 0 \f$ and \f$ \left| c \right| < 1 \f$  one has
       *   \f[ \left( \begin{array}{c}a   \\ b  \end{array}\right) = 
       *       \left( \begin{array}{c}-2  \\ 2 \end{array}\right) \f]
       *  - for \f$ w \le 0 \f$ one has
       *  - for \f$ w \le 0 \f$ and \f$ \left| c \right| < 1 \f$  one has
       *   \f[ \left( \begin{array}{c}a   \\ b  \end{array}\right) = 
       *       \left( \begin{array}{c}c-1 \\ c+1 \end{array}\right) \f]       
       *
       *  @param f1      the function 
       *  @param c       the parameter 
       *  @param ws      integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAWC is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAWC is used) 
       *  @param width   width parameter 
       *  @return value of the integral and the estimate of the uncertainty
       */
      template <class FUNCTION1>
      static inline double cauchy_pv_infinity 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ,
        const double         width      = 0 ) 
      { return cauchy_pv_infinity_ 
          ( std::cref ( f1 ) , c , 
            ws , tag , rescale , 
            aprecision , 
            rprecision , width ).first ; }
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi} 
       *   \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin  lower integration range 
       *  @param n     number of subtracion
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @return value of the dispersion integral 
       *  @see Ostap::Math::Integrator::cauchy_pv_to_infinity 
       */
      template <class FUNCTION1>
      static inline double kramers_kronig 
      ( FUNCTION1            f1             , 
        const double         s              , 
        const double         xmin           , 
        const unsigned short n              , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) 
      { return kramers_kronig_
          ( std::cref ( f1 ) , s , xmin , n , 
            ws , tag , rescale , 
            aprecision , 
            rprecision , width ).first ; }
      // ======================================================================
      /** integration with known singular points 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
       *  @param  f1 the   function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       *  @param points known singular points 
       *  @param ws integration workspace 
       *  @param tag unqye tag/label 
       *  
       *  - Only singular poins between \f$ x_{min} \f$  and \f$ x_{max} \f$ are considered 
       *  - \f$ x_{min} \f$  and \f$ x_{max \f$ are also  considered as singular points 
       */
      template <class FUNCTION1>
      static inline double integrate_singular
      ( FUNCTION1                  f1             , 
        const double               xmin           , 
        const double               xmax           ,
        const std::vector<double>& points         ,
        const WorkSpace&           ws             , 
        const std::size_t          tag        = 0 ,
        const double               aprecision = 0 , 
        const double               rprecision = 0 )
      { return integrate_singular_
          ( std::cref ( f1 ) , xmin , xmax , points , 
            ws , tag   ,  
            aprecision , 
            rprecision ).first ; }
      // ======================================================================
    public:
      // ======================================================================
      // Integration with doubly-adaptive CQUAD integrator  
      // ======================================================================
      /** calculate the integral using double adaptive  CQUAD integrator 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param ws   integration workspace 
       *  @param tag  unique tag/label 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      template <class FUNCTION1>
      static inline double integrate_cquad
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) 
      { return integrate_cquad_ 
          ( std::cref ( f1 ) , 
            xmin , xmax , 
            ws   , tag  , rescale , 
            aprecision  , 
            rprecision  ).first   ; }
      // ======================================================================
      // Integration using Romberg integrator 
      // ======================================================================
      /** calculate the integral using Romberg integratorr 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge
       *  @param ws   integration workspace 
       *  @param tag  unique tag/label 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      template <class FUNCTION1>
      static inline double integrate_romberg
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) 
      { return integrate_romberg_ 
          ( std::cref ( f1 ) , 
            xmin , xmax , 
            ws   , tag  , rescale , 
            aprecision  , 
            rprecision  ).first ; }
      // ======================================================================
    public: // the actual static methods to perform the integration 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param ws   integration workspace 
       *  @param tag   unique tag/label 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @param rule       the actual Gauss-Kronrod integration rule 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param ws integration workspace 
       *  @param tag   unique tag/label 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_infinity_
      ( function1         f1             , 
        const WorkSpace&  ws             , 
        const std::size_t tag        = 0 , 
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param integration workspace 
       *  @param tag   unique tag/label 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAGIU is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAGIU is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_to_infinity_
      ( function1         f1             , 
        const double      xmin           , 
        const WorkSpace&  ws             ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @param tag   unique tag/label 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAGIL is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAGIL is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_from_infinity_
      ( function1        f1              , 
        const double     xmax            , 
        const WorkSpace& ws              ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1      the function 
       *  @param c       the parameter 
       *  @param xmin    lower integration edge 
       *  @param xmax    upper  integration edge 
       *  @param ws      integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAWC is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAWC is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result cauchy_pv_
      ( function1            f1             , 
        const double         c              , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *
       *  Integral is calculated as sum of two components
       *   - \f$ \mathcal{P} \int_{x_{min}}^{b}\frac{f(x)}{x-c}dx \f$
       *   - \f$ int_{b}^{+\infty}\frac{f(x)}{x-c}dx \f$
       *
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param ws      integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision
       *  @param aprecision relative precision
       *  @param width width parameter 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result cauchy_pv_to_infinity_
      ( function1            f1             , 
        const double         c              , 
        const double         xmin           , 
        const WorkSpace&     ws             ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       *  @param integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision 
       *  @param aprecision relative precision 
       *  @param width width parameter 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result cauchy_pv_from_infinity_
      ( function1            f1             , 
        const double         c              , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 ,
        const double         rprecision = 0 ,
        const double         width      = 0 ) ;
      // ======================================================================
      /** get Cauchy principal value integral for infinite range 
       *  \f[ g(c) = \mathcal{P} \int_{-\infty}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *
       *  Integral is calculated as sum of three components
       *   - \f$ \int_{-\infty}^{a}\frac{f(x)}{x-c}dx \f$
       *   - \f$ \mathcal{P} \int_{a}^{b}\frac{f(x)}{x-c}dx \f$
       *   - \f$ int_{b}^{+\infty}\frac{f(x)}{x-c}dx \f$
       *
       *  where a and b are chosen such \f$ a < c < b \f$ and 
       *  - for \f$ 0<w \f$ one has 
       *   \f[ \left( \begin{array}{c}a   \\ b  \end{array}\right) = 
       *       \left( \begin{array}{c}c-w \\ c+w\end{array}\right) \f]
       *  - for \f$ w \le 0 \f$ and \f$ \left| c \right| < 1 \f$  one has
       *   \f[ \left( \begin{array}{c}a   \\ b  \end{array}\right) = 
       *       \left( \begin{array}{c}-2  \\ 2 \end{array}\right) \f]
       *  - for \f$ w \le 0 \f$ one has
       *  - for \f$ w \le 0 \f$ and \f$ \left| c \right| < 1 \f$  one has
       *   \f[ \left( \begin{array}{c}a   \\ b  \end{array}\right) = 
       *       \left( \begin{array}{c}c-1 \\ c+1 \end{array}\right) \f]
       *
       *  @param f1      the function 
       *  @param c       the parameter 
       *  @param ws      integration workspace 
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAWC is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAWC is used) 
       *  @param width width parameter \f$ w \f$
       *  @return value of the integral and the estimate of the uncertainty
       *
       */
      static result cauchy_pv_infinity_
      ( function1            f1             , 
        const double         c              , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const double         width      = 0 ) ;
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi} 
       *   \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin  lower integration range 
       *  @param n     number of subtracion
       *  @param tag     unique label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision 
       *  @param aprecision relative precision 
       *  @param width width parameter \f$ w \f$
       *  @return value of the integral and the estimate of the uncertainty
       *  @see Ostap::Math::Integrator::cauchy_pv_to_infinity 
       */
      static result kramers_kronig_
      ( function1            f1             , 
        const double         s              , 
        const double         xmin           , 
        const unsigned short n              , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 ,
        const double         rprecision = 0 ,
        const double         width      = 0 ) ;
      // ======================================================================
      /** integration with known singular points 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
       *  @param  f1 the   function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       *  @param points known singular points 
       *  @param ws integration workspace 
       *  @param tag unqye tag/label 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAGP is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAGP is used) 
       *  @return value of the integral and the estimate of the uncertainty
       *  
       *  - Only singular poins between \f$ x_{min} \f$  and \f$ x_{max} \f$ are considered 
       *  - \f$ x_{min} \f$  and \f$ x_{max \f$ are also  considered as singular points 
       */
      static result integrate_singular_
      ( function1                  f1             , 
        const double               xmin           , 
        const double               xmax           ,
        const std::vector<double>& points         ,
        const WorkSpace&           ws             , 
        const std::size_t          tag        = 0 ,
        const double               aprecision = 0 ,
        const double               rprecision = 0 ) ;
      // ======================================================================
    public :   // Integtaion using double-adaptive CQUAD integrator 
      // ======================================================================
      /** calculate the integral using double adaptive  CQUAD integrator 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param ws   integration workspace 
       *  @param tag  unique tag/label 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_cquad_
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // ======================================================================
    public :   // Integtaion using Romberg integrator 
      // ======================================================================
      /** calculate the integral using Romberg integratorr 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge
       *  @param ws   integration workspace 
       *  @param tag  unique tag/label 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_romberg_
      ( function1            f1             , 
        const double         xmin           , 
        const double         xmax           , 
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // ======================================================================
    public:  // integrals for 2D functions 
      // ======================================================================
      /** calculate the 2D integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate2_
      ( function2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) ;
      // =======================================================================
      // partial integration for 2D functions 
      // =======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
       *  @param f2 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param integration workspace 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate2X_
      ( function2            f2             ,   
        const double         y              , 
        const double         xmin           ,
        const double         xmax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
       *  @param f2 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate2Y_
      ( function2            f2             , 
        const double         x              , 
        const double         ymin           ,
        const double         ymax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ,
        const int            rule       = 0 ) ;
      // =======================================================================
    public: // integrals for 3D-function 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}
       *          \int_{y_{min}}^{y_{max}}
       *          \int_{z_{min}}^{z_{max}}f_3(x,y,z) dx dy dz \f]
       *  @param f3 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE3D is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE3D is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate3_
      ( function3         f3             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const double      zmin           , 
        const double      zmax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) ;
      // ======================================================================
      // partial integration for 3D functions 
      // =======================================================================
      /** calculate the integral 
       *  \f[ r(z) = \int_{x_{min}}^{x_{max}}
       *             \int_{y_{min}}^{y_{max}} f_3(x,y,z) dxdy \f]
       *  @param f3 the function 
       *  @param z parameter z
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
       *  @return the value of the integral and the estimate of the error 
       *  @see Ostap::Math::Integrator::integrate2_ 
       */
      static result integrate3XY_
      ( function3            f3             ,
        const double         z              , 
        const double         xmin           ,
        const double         xmax           ,
        const double         ymin           ,
        const double         ymax           ,
        const std::size_t    tag        = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(y) = \int_{x_{min}}^{x_{max}}
       *             \int_{z_{min}}^{z_{max}} f_3(x,y,z) dxdz \f]
       *  @param f3 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
       *  @return the value of the integral and the estimate of the error 
       *  @see Ostap::Math::Integrator::integrate2_ 
       */
      static result integrate3XZ_
      ( function3            f3             ,   
        const double         y              , 
        const double         xmin           ,
        const double         xmax           ,
        const double         zmin           ,
        const double         zmax           ,
        const std::size_t    tag        = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(x) = \int_{y_{min}}^{y_{max}}
       *             \int_{z_{min}}^{z_{max}} f_3(x,y,z) dydz \f]
       *  @param f3 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
       *  @return the value of the integral and the estimate of the error 
       *  @see Ostap::Math::Integrator::integrate2_ 
       */
      static result integrate3YZ_
      ( function3            f3             ,   
        const double         x              , 
        const double         ymin           ,
        const double         ymax           ,
        const double         zmin           ,
        const double         zmax           ,
        const std::size_t    tag        = 0 , 
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // =====================================================================
      /** calculate the integral 
       *  \f[ r(y,z) = \int_{x_{min}}^{x_{max}}f_3(x,y,z) dx  \f]
       *  @param f3 the function 
       *  @param y parameter y
       *  @param z parameter z
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param integration workspace 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate3X_
      ( function3            f3             ,   
        const double         y              , 
        const double         z              , 
        const double         xmin           ,
        const double         xmax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(x,z) = \int_{y_{min}}^{y_{max}}f_3(x,y,z) dy  \f]
       *  @param f3 the function 
       *  @param x parameter x
       *  @param z parameter z
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate3Y_
      ( function3            f3             ,   
        const double         x              , 
        const double         z              , 
        const double         ymin           ,
        const double         ymax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r(x,y) = \int_{z_{min}}^{z_{max}}f_3(x,y,z) dyz \f]
       *  @param f3 the function 
       *  @param x parameter x
       *  @param y parameter y
       *  @param zmin lower integration edge in z 
       *  @param zmax upper integration edge in z 
       *  @param integration workspace 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate3Z_
      ( function3            f3             ,   
        const double         x              , 
        const double         z              , 
        const double         ymin           ,
        const double         ymax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 , 
        const int            rule       = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// set the QAG integration rule 
      inline int qag_rule () const { return m_qag_rule ; }
      /// set absolute/relative precision for GAG
      inline double abs_precision_qag     () const { return m_abs_precision_qag     ; }
      /// set absolute/relative precision for GAG
      inline double rel_precision_qag     () const { return m_rel_precision_qag     ; }
      /// set absolute/relative precision for GAGI
      inline double abs_precision_qagi    () const { return m_abs_precision_qagi    ; }
      /// set absolute/relative precision for GAGI
      inline double rel_precision_qagi    () const { return m_rel_precision_qagi    ; }
      /// set absolute/relative precision for GAGIU
      inline double abs_precision_qagiu   () const { return m_abs_precision_qagiu   ; }
      /// set absolute/relative precision for GAGIU
      inline double rel_precision_qagiu   () const { return m_rel_precision_qagiu   ; }
      /// set absolute/relative precision for GAGIL
      inline double abs_precision_qagil   () const { return m_abs_precision_qagil   ; }
      /// set absolute/relative precision for GAGIL
      inline double rel_precision_qagil   () const { return m_rel_precision_qagil   ; }
      /// set absolute/relative precision for GAGP
      inline double abs_precision_qagp    () const { return m_abs_precision_qagp    ; }
      /// set absolute/relative precision for GAGP
      inline double rel_precision_qagp    () const { return m_rel_precision_qagp    ; }
      /// set absolute/relative precision for GAWC
      inline double abs_precision_qawc    () const { return m_abs_precision_qawc    ; }
      /// set absolute/relative precision for GAWC
      inline double rel_precision_qawc    () const { return m_rel_precision_qawc    ; }
      /// set absolute/relative precision for CPV
      inline double abs_precision_cpv     () const { return m_abs_precision_cpv     ; }
      /// set absolute/relative precision for CPV
      inline double rel_precision_cpv     () const { return m_rel_precision_cpv     ; }      
      /// set absolute/relative precision for KK
      inline double abs_precision_kk      () const { return m_abs_precision_kk      ; }
      /// set absolute/relative precision for KK
      inline double rel_precision_kk      () const { return m_rel_precision_kk      ; }
      /// set absolute/relative precision for CQUAD
      inline double abs_precision_cquad   () const { return m_abs_precision_cquad   ; }
      /// set absolute/relative precision for CQUAD
      inline double rel_precision_cquad   () const { return m_rel_precision_cquad   ; }      
      /// set absolute/relative precision for Romberg
      inline double abs_precision_romberg () const { return m_abs_precision_romberg ; }
      /// set absolute/relative precision for Romberg
      inline double rel_precision_romberg () const { return m_rel_precision_romberg ; }
      /// set absolute/relative precision for CUBE2
      inline double abs_precision_cube2   () const { return m_abs_precision_cube2   ; }
      /// set absolute/relative precision for CUBE2
      inline double rel_precision_cube2   () const { return m_rel_precision_cube2   ; }
      /// set absolute/relative precision for CUBE3
      inline double abs_precision_cube3   () const { return m_abs_precision_cube3   ; }
      /// set absolute/relative precision for CUBE3
      inline double rel_precision_cube3   () const { return m_rel_precision_cube3   ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// set the QAG integration rule 
      void set_qag_rule          ( const int rule ) ;
      /// set absolute/relatibe precision for GAG
      void set_precision_qag     ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for GAGI
      void set_precision_qagi    ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for GAGIL
      void set_precision_qagil   ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for GAGIU
      void set_precision_qagiu   ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for QAGP 
      void set_precision_qagp    ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for QAWC
      void set_precision_qawc    ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for Cauchy PV inetegration 
      void set_precision_cpv     ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for Cauchy PV/inf inetegration 
      void set_precision_cpvi    ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision for Kramers-Kronig integrtaion 
      void set_precision_kk      ( const double aprec , const double rprec ) ;
      /// set absolute/relative precvision for CQUAD   
      void set_precision_cquad   ( const double aprec , const double rprec ) ;
      /// set absolute/relative precvision for Romberg 
      void set_precision_romberg ( const double aprec , const double rprec ) ;
      /// set absolute/relative precvision for 2D-cubature 
      void set_precision_cube2   ( const double aprec , const double rprec ) ;
      /// set absolute/relative precvision for 3D-cubature 
      void set_precision_cube3   ( const double aprec , const double rprec ) ;
      // ======================================================================
    public: // get the workspace 
      // ======================================================================
      /// get the workspace 
      const Ostap::Math::WorkSpace& ws   () const { return m_workspace ; }
      /// get the name 
      const std::string&            name () const { return m_name      ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get low limit/cutoff for Cauchy-PV integration
      static double cauchy_pv_a  
      ( const double c     ,
        const double width ) ;
      // =======================================================================
      /// get high limit/cutoff  for Cauchy-PV integration
      static double cauchy_pv_b 
      ( const double c     ,
        const double width ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// integrator name 
      std::string m_name             ; // integrator name
      /// QAG integration rule 
      int    m_qag_rule              ;
      /// absolute precision for QAG integration 
      double m_abs_precision_qag     ; // absolute precision for QAG   integration 
      /// relative precision for QAG integration
      double m_rel_precision_qag     ; // relative precision for QAG   integration
      /// absolute precision for QAGI integration 
      double m_abs_precision_qagi    ; // absolute precision for QAGI  integration 
      /// relative precision for QAGI integration
      double m_rel_precision_qagi    ; // relative precision for QAGI  integration
      /// absolute precision for QAGIU integration 
      double m_abs_precision_qagiu   ; // absolute precision for QAGIU integration 
      /// relative precision for QAGIU integration
      double m_rel_precision_qagiu   ; // relative precision for QAGIU integration
      /// absolute precision for QAGIL integration 
      double m_abs_precision_qagil   ; // absolute precision for QAGIL integration 
      /// relative precision for QAGIL integration
      double m_rel_precision_qagil   ; // relative precision for QAGIL integration
      /// absolute precision for QAGP integration 
      double m_abs_precision_qagp    ; // absolute precision for QAGP  integration 
      /// relative precision for QAGP integration
      double m_rel_precision_qagp    ; // relative precision for QAGP  integration      
      /// absolute precision for QAWC integration 
      double m_abs_precision_qawc    ; // absolute precision for QAWC  integration 
      /// relative precision for QAWC integration
      double m_rel_precision_qawc    ; // relative precision for QAWC  integration
      /// absolute precision for Cauchy PV integration 
      double m_abs_precision_cpv     ; // absolute precion for Cauchy integration      
      /// relative precision for Cauchy PV integration 
      double m_rel_precision_cpv     ; // relative precion for Cauchy integrtaion
      /// absolute precision for Cauchy PV/inf integration 
      double m_abs_precision_cpvi    ; // absolute precion for Cauchy/inf integrtaion      
      /// relative precision for Cauchy PV/inf integration 
      double m_rel_precision_cpvi    ; // relative precion for Cauchy/inf integrtaion
      /// absolute precision for Kramers-Kronig integration 
      double m_abs_precision_kk      ; // absolute precion for Kramers-Kronig
      /// relative precision for Kramers-Kronig integration 
      double m_rel_precision_kk      ; // relative precion for Kramers-Kronig
      /// absolute precision for CQUAD integration 
      double m_abs_precision_cquad   ; // absolute precion for CQUAD 
      /// relative precision for CQUAD integration 
      double m_rel_precision_cquad   ; // relative precion for CQUAD
      /// absolute precision for Romberg integration 
      double m_abs_precision_romberg ; // absolute precion for Romberg
      /// relative precision for Romberg integration 
      double m_rel_precision_romberg ; // relative precion for Romberg
      /// absolute precision for 2d-cubature integration
      double m_abs_precision_cube2   ; // absolute precision for cubature integration 
      /// relative precision for 2d-cubature integration
      double m_rel_precision_cube2   ; // relative precision for cubature integration
      /// absolute precision for 3d-cubature integration
      double m_abs_precision_cube3   ; // absolute precision for cubature integration 
      /// relative precision for 3d-cubature integration
      double m_rel_precision_cube3   ; // relative precision for cubature integration
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace 
      Ostap::Math::WorkSpace m_workspace {}  ; // integration workspace 
      // ======================================================================
    };
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_INTEGRATOR_H
// ============================================================================
