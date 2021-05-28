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
      Integrator ( const std::size_t size = 0 ) ;
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
      ( FUNCTION1            f1          , 
        const double         xmin        , 
        const double         xmax        ,
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) const 
      { return integrate_ ( std::cref ( f1 ) , xmin , xmax , m_workspace , tag , rescale ).first  ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_infinity
      ( FUNCTION1         f1       , 
        const std::size_t tag  = 0 ) const 
      { return integrate_infinity  ( std::cref ( f1 ) , m_workspace , tag ) ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_to_infinity
      ( FUNCTION1         f1       , 
        const double      xmin     , 
        const std::size_t tag  = 0 ) const 
      { return integrate_to_infinity  ( std::cref ( f1 ) , xmin , m_workspace , tag ) ; }
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
      ( FUNCTION1         f1       , 
        const double      xmax     , 
        const std::size_t tag  = 0 ) const 
      { return integrate_from_infinity  ( std::cref ( f1 ) , xmax , m_workspace , tag ) ; }
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
      ( FUNCTION1            f1          , 
        const double         c           , 
        const double         xmin        , 
        const double         xmax        , 
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) const 
      { return cauchy_pv ( std::cref ( f1 ) , c , xmin , xmax , m_workspace , tag , rescale ) ; }
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
      ( FUNCTION1            f1          , 
        const double         c           , 
        const double         xmin        ,  
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) const 
      { return cauchy_pv_to_infinity ( std::cref ( f1 ) , c , xmin , m_workspace , tag , rescale ) ; }
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
      ( FUNCTION1            f1          , 
        const double         c           , 
        const double         xmax        , 
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) const 
      { return cauchy_pv_from_infinity ( std::cref ( f1 ) , c , xmax , m_workspace , tag , rescale ) ; }
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
      ( FUNCTION1            f1          , 
        const double         s           , 
        const double         xmin        , 
        const unsigned short n       = 0 , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) const 
      { return kramers_kronig ( std::cref ( f1 ) , s , xmin , n , m_workspace , tag , rescale ) ; }
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
      ( FUNCTION2         f2      , 
        const double      y       , 
        const double      xmin    ,
        const double      xmax    ,  
        const std::size_t tag = 0 ) const
      { return integrate ( std::cref ( f2 ) , y , xmin , xmax , m_workspace , tag ) ; }
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
      ( FUNCTION2         f2      , 
        const double      x       , 
        const double      ymin    ,
        const double      ymax    ,
        const std::size_t tag = 0 ) const 
      { return integrate ( std::cref ( f2 ) , x , ymin , ymax , m_workspace , tag ) ; }
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
      ( FUNCTION1                  f1       , 
        const double               xmin     , 
        const double               xmax     ,
        const std::vector<double>& points   ,
        const std::size_t          tag  = 0 ) const 
      { return integrate_singular ( std::cref ( f1 ) , xmin , xmax , points , m_workspace , tag ) ; }
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
      ( FUNCTION1            f1          , 
        const double         xmin        , 
        const double         xmax        , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) 
      { return integrate_ ( std::cref ( f1 ) , xmin , xmax , ws , tag , rescale ).first ; }
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
      ( FUNCTION1         f1      , 
        const WorkSpace&  ws      , 
        const std::size_t tag = 0 ) 
      { return integrate_infinity_ ( std::cref ( f1 ) , ws , tag ).first ; }        
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
      ( FUNCTION1         f1      , 
        const double      xmin    , 
        const WorkSpace&  ws      ,
        const std::size_t tag = 0 ) 
      { return integrate_to_infinity_ ( std::cref ( f1 ) , xmin , ws , tag ).first ; }        
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
      ( function1        f1       , 
        const double     xmax     , 
        const WorkSpace& ws       ,
        const std::size_t tag = 0 ) 
      { return integrate_from_infinity_ ( std::cref ( f1 ) , xmax , ws , tag ).first ; }        
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
      ( FUNCTION1            f1          , 
        const double         c           , 
        const double         xmin        , 
        const double         xmax        , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) 
      { return cauchy_pv_ ( std::cref ( f1 ) , c , xmin , xmax , ws , tag , rescale ).first ; }   
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
      ( FUNCTION1            f1          , 
        const double         c           , 
        const double         xmin        , 
        const WorkSpace&     ws          ,
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) 
      { return cauchy_pv_to_infinity_ ( std::cref ( f1 ) , c , xmin , ws , tag , rescale ).first ; }
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
      ( FUNCTION1            f1          , 
        const double         c           , 
        const double         xmax        , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) 
      { return cauchy_pv_from_infinity_ ( std::cref ( f1 ) , c , xmax , ws , tag , rescale ).first ;}
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
      ( FUNCTION1            f1          , 
        const double         s           , 
        const double         xmin        , 
        const unsigned short n           , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) 
      { return kramers_kronig_ ( std::cref ( f1 ) , s , xmin , n , ws , tag , rescale ).first ; }
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
      ( FUNCTION1                  f1       , 
        const double               xmin     , 
        const double               xmax     ,
        const std::vector<double>& points   ,
        const WorkSpace&           ws       , 
        const std::size_t          tag  = 0 ) 
      { return integrate_singular_ ( std::cref ( f1 ) , xmin , xmax , points , ws , tag ).first ; }
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
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_
      ( function1            f1          , 
        const double         xmin        , 
        const double         xmax        , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param ws integration workspace 
       *  @param tag   unique tag/label 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_infinity_
      ( function1         f1      , 
        const WorkSpace&  ws      , 
        const std::size_t tag = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param integration workspace 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_to_infinity_
      ( function1         f1      , 
        const double      xmin    , 
        const WorkSpace&  ws      ,
        const std::size_t tag = 0 ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result integrate_from_infinity_
      ( function1        f1       , 
        const double     xmax     , 
        const WorkSpace& ws       ,
        const std::size_t tag = 0 ) ;
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
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result cauchy_pv_
      ( function1            f1          , 
        const double         c           , 
        const double         xmin        , 
        const double         xmax        , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param ws      integration workspace 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result cauchy_pv_to_infinity_
      ( function1            f1          , 
        const double         c           , 
        const double         xmin        , 
        const WorkSpace&     ws          ,
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       *  @param integration workspace 
       *  @param tag     unuque label/tag 
       *  @param rescale rescale function for better numerical precision  
       *  @return value of the integral and the estimate of the uncertainty
       */
      static result cauchy_pv_from_infinity_
      ( function1            f1          , 
        const double         c           , 
        const double         xmax        , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) ;
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
       *  @return value of the integral and the estimate of the uncertainty
       *  @see Ostap::Math::Integrator::cauchy_pv_to_infinity 
       */
      static result kramers_kronig_
      ( function1            f1          , 
        const double         s           , 
        const double         xmin        , 
        const unsigned short n           , 
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) ;
      // ======================================================================
      /** integration with known singular points 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
       *  @param  f1 the   function 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       *  @param points known singular points 
       *  @param ws integration workspace 
       *  @param tag unqye tag/label 
       *  @return value of the integral and the estimate of the uncertainty
       *  
       *  - Only singular poins between \f$ x_{min} \f$  and \f$ x_{max} \f$ are considered 
       *  - \f$ x_{min} \f$  and \f$ x_{max \f$ are also  considered as singular points 
       */
      static result integrate_singular_
      ( function1                  f1       , 
        const double               xmin     , 
        const double               xmax     ,
        const std::vector<double>& points   ,
        const WorkSpace&           ws       , 
        const std::size_t          tag  = 0 ) ;
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
      ( FUNCTION2         f2      , 
        const double      xmin    , 
        const double      xmax    ,
        const double      ymin    , 
        const double      ymax    ,
        const WorkSpace&  ws      , 
        const std::size_t tag = 0 ) 
      { return integrate_ ( std::cref ( f2 ) , xmin , xmax , ymin , ymax , ws , tag ).first ; }
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
      ( FUNCTION2         f2      , 
        const double      xmin    , 
        const double      xmax    ,
        const double      ymin    , 
        const double      ymax    ,
        const std::size_t tag = 0 ) 
      { return integrate_ ( std::cref ( f2 ) , xmin , xmax , ymin , ymax , tag ).first ; }
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
      ( FUNCTION2            f2          , 
        const double         y           , 
        const double         xmin        ,
        const double         xmax        ,
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) 
      { return integrateX_ ( std::cref ( f2 ) , y , xmin , xmax , ws , tag , rescale ).first ; }
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
      ( FUNCTION2            f2          , 
        const double         x           , 
        const double         ymin        ,
        const double         ymax        ,
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) 
      { return integrateY_ ( std::cref ( f2 ) , x , ymin , ymax , ws , tag , rescale ).first ; }
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
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate_
      ( function2         f2      , 
        const double      xmin    , 
        const double      xmax    ,
        const double      ymin    , 
        const double      ymax    ,
        const WorkSpace& /* ws */ , 
        const std::size_t tag = 0 ) ;
      // =======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param f2 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrate_
      ( function2         f2      , 
        const double      xmin    , 
        const double      xmax    ,
        const double      ymin    , 
        const double      ymax    ,
        const std::size_t tag = 0 ) ;
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
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrateX_
      ( function2            f2          ,   
        const double         y           , 
        const double         xmin        ,
        const double         xmax        ,
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 , 
        const unsigned short rescale = 0 ) ;
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
       *  @return the value of the integral and the estimate of the error 
       */
      static result integrateY_
      ( function2            f2          , 
        const double         x           , 
        const double         ymin        ,
        const double         ymax        ,
        const WorkSpace&     ws          , 
        const std::size_t    tag     = 0 ,
        const unsigned short rescale = 0 ) ;
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
