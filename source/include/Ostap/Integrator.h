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
     *  simple numerical integrator for 1D&2D-cases 
     */
    class Integrator 
    {
    public:
      // ======================================================================
      typedef std::function<double(double)>        function1 ;
      typedef std::function<double(double,double)> function2 ;
      typedef std::pair<double,double>             result    ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with integration workspace size
      Integrator
      ( const std::size_t    size         = 0  , 
        const unsigned short size_cquad   = 0  ,
        const unsigned short size_romberg = 0  ) ;
      // ======================================================================
      /// constructor with integration workspace
      Integrator
      ( const Ostap::Math::WorkSpace& ws        ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param tag  uqniue label/tag  
       *  @param rescale rescale function for better numerical precision  
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_ ( std::cref ( f1 )    , 
                            xmin , xmax ,
                            m_workspace , 
                            tag  , rescale      , 
                            0 < aprecision ? aprecision : m_abs_precision_gaq , 
                            0 < rprecision ? rprecision : m_rel_precision_gaq , 
                            m_gaq_rule ).first ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_infinity
      ( FUNCTION1         f1             , 
        const std::size_t tag        = 0 , 
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_infinity_  ( std::cref    ( f1 ) , 
                                      m_workspace , tag   , 
                                      0 < aprecision ? aprecision : m_abs_precision_gaqi , 
                                      0 < rprecision ? rprecision : m_rel_precision_gaqi ).first  ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_to_infinity
      ( FUNCTION1         f1             , 
        const double      xmin           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const 
      { return integrate_to_infinity_ ( std::cref ( f1 ) , xmin , 
                                        m_workspace      , tag  ,
                                        0 < aprecision ? aprecision : m_abs_precision_gaqiu , 
                                        0 < rprecision ? rprecision : m_rel_precision_gaqiu ).first ; }
      // ========================================================================
      /** calculate the integral 
       *  \f[ r = \mathcal{P}\int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @param tag     unuque label/tag 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_from_infinity
      ( FUNCTION1         f1             , 
        const double      xmax           , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const
      { return integrate_from_infinity_ ( std::cref ( f1 ) , xmax , 
                                          m_workspace , tag       , 
                                          0 < aprecision ? aprecision : m_abs_precision_gaqil , 
                                          0 < rprecision ? rprecision : m_rel_precision_gaqil ).first ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1      the function 
       *  @param c       the parameter 
       *  @param xmin    lower integration edge 
       *  @param xmax    upper integration edge 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double cauchy_pv 
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           , 
        const double         xmax           , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const
      { return cauchy_pv_ ( std::cref ( f1 ) , c , xmin , xmax , 
                            m_workspace , tag , rescale  , 
                            0 < aprecision ? aprecision : m_abs_precision_cpv , 
                            0 < rprecision ? rprecision : m_rel_precision_cpv ).first  ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double cauchy_pv_to_infinity  
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmin           ,  
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const
      { return cauchy_pv_to_infinity_ ( std::cref ( f1 ) , c , xmin , 
                                        m_workspace , tag , rescale ,
                                        0 < aprecision ? aprecision : m_abs_precision_cpvi , 
                                        0 < rprecision ? rprecision : m_rel_precision_cpvi ).first ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double cauchy_pv_from_infinity  
      ( FUNCTION1            f1             , 
        const double         c              , 
        const double         xmax           , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const
      { return cauchy_pv_from_infinity_ ( std::cref ( f1 ) , c , xmax , 
                                          m_workspace , tag , rescale , 
                                          0 < aprecision ? aprecision : m_abs_precision_cpvi , 
                                          0 < rprecision ? rprecision : m_rel_precision_cpvi ).first ; }
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi}
       *     \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin low integration edge 
       *  @param n   number of subtractions 
       *  @see Ostap::Math::Integrtator::cauchy_pv_to_infinity 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical treatmenrt 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double kramers_kronig 
      ( FUNCTION1            f1             , 
        const double         s              , 
        const double         xmin           , 
        const unsigned short n          = 0 , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const
      { return kramers_kronig_ ( std::cref ( f1 ) , s , xmin , n , 
                                 m_workspace , tag , rescale ,
                                 0 < aprecision ? aprecision : m_abs_precision_kk , 
                                 0 < rprecision ? rprecision : m_rel_precision_kk ).first  ; }
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
      double integrateX
      ( FUNCTION2         f2             , 
        const double      y              , 
        const double      xmin           ,
        const double      xmax           ,  
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const
      { return integrate_ ( std::cref ( f2 )     , 
                            y , xmin , xmax      , 
                            m_workspace , tag    , 
                            0 < aprecision ? aprecision : m_abs_precision_gaq , 
                            0 < rprecision ? rprecision : m_rel_precision_gaq , 
                            m_gaq_rule     ).first  ; }
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
      double integrateY
      ( FUNCTION2         f2             , 
        const double      x              , 
        const double      ymin           ,
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) const
      { return integrate_ ( std::cref ( f2 )     , 
                            x , ymin , ymax      , 
                            m_workspace , tag    , 
                            0 < aprecision ? aprecision : m_abs_precision_gaq , 
                            0 < rprecision ? rprecision : m_rel_precision_gaq , 
                            m_gaq_rule ).first  ; }
      // ======================================================================
      /** integration with known singular points 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
       *  @param  f1 the   function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       *  @param points known singular points 
       *  @param tag unqye tag/label 
       *  
       *  - Only singular poins between \f$ x_{min} \f$  and \f$ x_{max} \f$ are considered 
       *  - \f$ x_{min} \f$  and \f$ x_{max \f$ are also  considered as singular points 
       */
      template <class FUNCTION1>
      double integrate_singular
      ( FUNCTION1                  f1             , 
        const double               xmin           , 
        const double               xmax           ,
        const std::vector<double>& points         ,
        const std::size_t          tag        = 0 ,
        const double               aprecision = 0 , 
        const double               rprecision = 0 ) const
      { return integrate_singular_ ( std::cref ( f1 ) , 
                                     xmin , xmax , points , 
                                     m_workspace , tag    , 
                                     0 < aprecision ? aprecision : m_abs_precision_gaqp , 
                                     0 < rprecision ? rprecision : m_rel_precision_gaqp ).first ; }
      
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
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_cquad
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_cquad_ ( std::cref ( f1 )    , 
                                  xmin , xmax ,
                                  m_workspace , 
                                  tag  , rescale      , 
                                  0 < aprecision ? aprecision : m_abs_precision_cquad , 
                                  0 < rprecision ? rprecision : m_rel_precision_cquad ).first ; }
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
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_romberg
      ( FUNCTION1            f1             , 
        const double         xmin           , 
        const double         xmax           ,
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) const 
      { return integrate_romberg_ ( std::cref ( f1 )    , 
                                    xmin , xmax ,
                                    m_workspace , 
                                    tag  , rescale      , 
                                    0 < aprecision ? aprecision : m_abs_precision_romberg , 
                                    0 < rprecision ? rprecision : m_rel_precision_romberg ).first ; }
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
      { return integrate_ ( std::cref ( f1 ) , 
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
      { return integrate_infinity_ ( std::cref ( f1 ) , 
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
      { return integrate_to_infinity_ ( std::cref ( f1 ) , xmin , 
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
      { return integrate_from_infinity_ ( std::cref ( f1 ) , xmax , 
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
       *  @param tag     unuque label/tag 
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
      { return cauchy_pv_ ( std::cref ( f1 ) , c , xmin , xmax , 
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
       *  @param tag     unuque label/tag 
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
        const double         rprecision = 0 )
      { return cauchy_pv_to_infinity_ ( std::cref ( f1 ) , c , xmin , 
                                        ws , tag , rescale ,
                                        aprecision , 
                                        rprecision ).first ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       *  @param integration workspace 
       *  @param tag     unuque label/tag 
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
        const double         rprecision = 0 )
      { return cauchy_pv_from_infinity_ ( std::cref ( f1 ) , c , xmax , 
                                          ws , tag , rescale , 
                                          aprecision , 
                                          rprecision ).first ; }
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi} 
       *   \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin  lower integration range 
       *  @param n     number of subtracion
       *  @param tag     unuque label/tag 
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
        const double         rprecision = 0 )        
      { return kramers_kronig_ ( std::cref ( f1 ) , s , xmin , n , 
                                 ws , tag , rescale , 
                                 aprecision , 
                                 rprecision ).first ; }
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
      { return integrate_singular_ ( std::cref ( f1 ) , xmin , xmax , points , 
                                     ws , tag   ,  
                                     aprecision , 
                                     rprecision ).first ; }
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
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQ is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQ is used) 
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
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQ is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQ is used) 
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
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQIU is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQIU is used) 
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
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQIL is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQIL is used) 
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
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAWC is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAWC is used) 
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
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param ws      integration workspace 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision
       *  @param aprecision relative precision
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
        const double         rprecision = 0 ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       *  @param integration workspace 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision 
       *  @param aprecision relative precision 
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
        const double         rprecision = 0 ) ;
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi} 
       *   \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin  lower integration range 
       *  @param n     number of subtracion
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @param aprecision absolute precision 
       *  @param aprecision relative precision 
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
        const double         rprecision = 0 ) ;
      // ======================================================================
      /** integration with known singular points 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
       *  @param  f1 the   function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       *  @param points known singular points 
       *  @param ws integration workspace 
       *  @param tag unqye tag/label 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQP is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQP is used) 
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
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQ is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQ is used) 
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
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQ is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQ is used) 
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
    public:
      // ======================================================================
      // 2D integration 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      static inline double integrate
      ( FUNCTION2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const WorkSpace&  ws             , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 ,
        const double      rprecision = 0 ) 
      { return integrate_ ( std::cref ( f2 ) , 
                            xmin , xmax ,
                            ymin , ymax , 
                            ws   , tag  , 
                            aprecision  , 
                            rprecision  ).first ; }
      // ======================================================================= 
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      static inline double integrate
      ( FUNCTION2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 ,
        const double      rprecision = 0 ) 
      { return integrate_ ( std::cref ( f2 ) , 
                            xmin , xmax , 
                            ymin , ymax , 
                            tag  , 
                            aprecision  , 
                            rprecision  ).first ; }
      // =======================================================================
      // 1D integration for 2D functions 
      // =======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
       *  @param f2 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ws integration workspace 
       *  @param tag unique label/tag 
       *  @param rescale rescale function for better numerical precision 
       *
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      static inline double integrateX
      ( FUNCTION2            f2             , 
        const double         y              , 
        const double         xmin           ,
        const double         xmax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 ,
        const double         rprecision = 0 ) 
      { return integrateX_ ( std::cref ( f2 ) , y , xmin , xmax , 
                             ws , tag , rescale , 
                             aprecision  , 
                             rprecision  ).first ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
       *  @param  f2 the function 
       *  @param  x parameter x
       *  @param  ymin lower integration edge in y 
       *  @param  ymax upper integration edge in y 
       *  @param  ws   integration workspace
       *  @param  tag unique label/tag 
       *  @param  rescale rescale function for better numerical precision 
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      static double integrateY
      ( FUNCTION2            f2             , 
        const double         x              , 
        const double         ymin           ,
        const double         ymax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 , 
        const double         aprecision = 0 ,
        const double         rprecision = 0 ) 
      { return integrateY_ ( std::cref ( f2 ) , x , ymin , ymax , 
                             ws , tag , rescale ,
                             aprecision  , 
                             rprecision  ).first ; }
      // ======================================================================
    public:
      // ======================================================================
      // 2D integration 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace (not used)
       *  @param  tag unique label/tag 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate_
      ( function2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const WorkSpace& /* ws */        , 
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) ;
      // =======================================================================
      /** calculate the integral 
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
      static result integrate_
      ( function2         f2             , 
        const double      xmin           , 
        const double      xmax           ,
        const double      ymin           , 
        const double      ymax           ,
        const std::size_t tag        = 0 ,
        const double      aprecision = 0 , 
        const double      rprecision = 0 ) ;
      // =======================================================================
      // 1D integration for 2D functions 
      // =======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
       *  @param f2 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param integration workspace 
       *  @param  tag unique label/tag 
       *  @param  rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQ is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQ is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrateX_
      ( function2            f2             ,   
        const double         y              , 
        const double         xmin           ,
        const double         xmax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 , 
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
       *  @param f2 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace 
       *  @param  tag unique label/tag 
       *  @param  rescale rescale function for better numerical precision 
       *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAQ is used) 
       *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAQ is used) 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrateY_
      ( function2            f2             , 
        const double         x              , 
        const double         ymin           ,
        const double         ymax           ,
        const WorkSpace&     ws             , 
        const std::size_t    tag        = 0 ,
        const unsigned short rescale    = 0 ,
        const double         aprecision = 0 , 
        const double         rprecision = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// set the GAQ integration rule 
      void set_gaq_rule          ( const int rule ) ;
      /// set absolute/relatibe precision fore GAG
      void set_precision_gaq     ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision fore GAGI
      void set_precision_gaqi    ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision fore GAGIL
      void set_precision_gaqil   ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision fore GAGIU
      void set_precision_gaqiu   ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision fore GAQP 
      void set_precision_gaqp    ( const double aprec , const double rprec ) ;
      /// set absolute/relative precision fore QAWC
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
      /// set absolute/relative precvision for cubature 
      void set_precision_cube    ( const double aprec , const double rprec ) ;
      // ======================================================================
    public: // get the workspace 
      // ======================================================================
      /// get the workspace 
      const Ostap::Math::WorkSpace& ws   () const { return m_workspace ; }
      /// get the name 
      const std::string&            name () const { return m_name      ; }
      // ======================================================================
    private:
      // ======================================================================
      /// integrator name 
      std::string m_name             ; // integrator name
      /// GAQ integration rule 
      int    m_gaq_rule              ;
      /// absolute precision for GAQ integration 
      double m_abs_precision_gaq     ; // absolute precision for GAQ   integration 
      /// relative precision for GAQ integration
      double m_rel_precision_gaq     ; // relative precision for GAQ   integration
      /// absolute precision for GAQI integration 
      double m_abs_precision_gaqi    ; // absolute precision for GAQI  integration 
      /// relative precision for GAQI integration
      double m_rel_precision_gaqi    ; // relative precision for GAQI  integration
      /// absolute precision for GAQIU integration 
      double m_abs_precision_gaqiu   ; // absolute precision for GAQIU integration 
      /// relative precision for GAQIU integration
      double m_rel_precision_gaqiu   ; // relative precision for GAQIU integration
      /// absolute precision for GAQIL integration 
      double m_abs_precision_gaqil   ; // absolute precision for GAQIL integration 
      /// relative precision for GAQIL integration
      double m_rel_precision_gaqil   ; // relative precision for GAQIL integration
      /// absolute precision for GAQP integration 
      double m_abs_precision_gaqp    ; // absolute precision for GAQP  integration 
      /// relative precision for GAQP integration
      double m_rel_precision_gaqp    ; // relative precision for GAQP  integration      
      /// absolute precision for QAWC integration 
      double m_abs_precision_qawc    ; // absolute precision for QAWC  integration 
      /// relative precision for QAWC integration
      double m_rel_precision_qawc    ; // relative precision for QAWC  integration
      /// absolute precision for Cauchy PV integration 
      double m_abs_precision_cpv     ; // absolute precion for Cauchy integrtaion      
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
      /// absolute precision for cubature integration
      double m_abs_precision_cube    ; // absolute precision for cubature integration 
      /// relative precision for cubature integration
      double m_rel_precision_cube    ; // relative precision for cubature integration
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
