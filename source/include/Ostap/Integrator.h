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
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with integrtaion workspace size 
      Integrator ( const std::size_t size = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax uppr  integration edge
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate
      ( FUNCTION1        f1   , 
        const double     xmin , 
        const double     xmax ) const
      { return integrate ( std::cref ( f1 ) , xmin , xmax , m_workspace ) ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_infinity
      ( FUNCTION1        f1   )  const 
      { return integrate_infinity  ( std::cref ( f1 ) , m_workspace ) ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_to_infinity
      ( FUNCTION1        f1   , 
        const double     xmin )  const 
      { return integrate_to_infinity  ( std::cref ( f1 ) , xmin , m_workspace ) ; }
      // ========================================================================
      /** calculate the integral 
       *  \f[ r = \mathcal{P}\int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      template <class FUNCTION1>
      double integrate_from_infinity
      ( FUNCTION1        f1   , 
        const double     xmax )  const 
      { return integrate_from_infinity  ( std::cref ( f1 ) , xmax , m_workspace ) ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param xmax upper integration edge 
       */
      template <class FUNCTION1>
      double cauchy_pv 
      ( FUNCTION1        f1   , 
        const double     c    , 
        const double     xmin , 
        const double     xmax ) const 
      { return cauchy_pv ( std::cref ( f1 ) , c , xmin , xmax , m_workspace ) ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       */
      template <class FUNCTION1>
      double cauchy_pv_to_infinity  
      ( FUNCTION1        f1   , 
        const double     c    , 
        const double     xmin ) const 
      { return cauchy_pv_to_infinity ( std::cref ( f1 ) , c , xmin , m_workspace ) ; }
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       */
      template <class FUNCTION1>
      double cauchy_pv_from_infinity  
      ( FUNCTION1        f1   , 
        const double     c    , 
        const double     xmax ) const 
      { return cauchy_pv_from_infinity ( std::cref ( f1 ) , c , xmax , m_workspace ) ; }
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi}
       *     \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin low integrtaion edge 
       *  @param n   number of subtractions 
       *  @see Ostap::Math::Integrtator::cauchy_pv_to_infinity 
       */
      template <class FUNCTION1>
      double kramers_kronig 
      ( FUNCTION1            f1       , 
        const double         s        , 
        const double         xmin     , 
        const unsigned short n    = 0 ) const 
      { return kramers_kronig ( std::cref ( f1 ) , s , xmin , n , m_workspace ) ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
       *  @param f2 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      double integrateX
      ( FUNCTION2        f2   , 
        const double     y    , 
        const double     xmin ,
        const double     xmax ) const 
      { return integrate ( std::cref ( f2 ) , y , xmin , xmax , m_workspace ) ; }
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
       *  @param f2 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      template <class FUNCTION2>
      double integrateY
      ( FUNCTION2        f2   , 
        const double     x    , 
        const double     ymin ,
        const double     ymax ) const 
      { return integrate ( std::cref ( f2 ) , x , ymin , ymax , m_workspace ) ; }
      // ======================================================================
    public: // the actual static methods to perform the integration 
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      static double integrate
      ( function1        f1   , 
        const double     xmin , 
        const double     xmax , 
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      static double integrate_infinity
      ( function1        f1   , 
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      static double integrate_to_infinity
      ( function1        f1   , 
        const double     xmin , 
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{-\infty}^{x_{max}} f_1(x) dx \f]
       *  @param f1 the function 
       *  @param xmax upper  integration edge
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      static double integrate_from_infinity
      ( function1        f1   , 
        const double     xmax , 
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param xmax upper  integration edge 
       *  @param integration workspace 
       */
      static double cauchy_pv 
      ( function1        f1   , 
        const double     c    , 
        const double     xmin , 
        const double     xmax , 
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmin lower integration edge 
       *  @param integration workspace 
       */
      static double cauchy_pv_to_infinity 
      ( function1        f1   , 
        const double     c    , 
        const double     xmin , 
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** get Cauchy principal value integral 
       *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
       *  @param f1 the function 
       *  @param c  the parameter 
       *  @param xmax upper integration edge 
       *  @param integration workspace 
       */
      static double cauchy_pv_from_infinity 
      ( function1        f1   , 
        const double     c    , 
        const double     xmax , 
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** Kramers-Kronig dispersion relation with n-subtractions 
       *  \f[ g(s) = \frac{s^n}{\pi} 
       *   \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
       *  @param f1 the function 
       *  @param s  s-parameter 
       *  @param xmin  lower integration range 
       *  @param n     number of subtracion
       *  @return value of the dispersion integral 
       *  @see Ostap::Math::Integrator::cauchy_pv_to_infinity 
       */
      static double kramers_kronig 
      ( function1            f1   , 
        const double         s    , 
        const double         xmin , 
        const unsigned short n    , 
        const WorkSpace&     ws   ) ;
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
      static double integrate
      ( function2        f2   , 
        const double     xmin , 
        const double     xmax ,
        const double     ymin , 
        const double     ymax ,
        const WorkSpace& /* ws */ ) ;
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
      static double integrate
      ( function2        f2   , 
        const double     xmin , 
        const double     xmax ,
        const double     ymin , 
        const double     ymax ) ;
      // =======================================================================
      // 1D integrtaion for 2D functions 
      // =======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
       *  @param f2 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      static double integrateX
      ( function2        f2   , 
        const double     y    , 
        const double     xmin ,
        const double     xmax ,
        const WorkSpace& ws   ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
       *  @param f2 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      static double integrateY
      ( function2        f2   , 
        const double     x    , 
        const double     ymin ,
        const double     ymax ,
        const WorkSpace& ws   ) ;
      // ======================================================================
    public: // integration with cache 
      // ======================================================================      
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
       *  @param tag unique tag 
       *  @param f1 the function 
       *  @param xmin lower integration edge 
       *  @param xmax uppr  integration edge
       *  @param integration workspace 
       *  @return the value of the integral 
       */
      static double  integrate
      ( const std::size_t tag  , 
        function1         f1   , 
        const double      xmin , 
        const double      xmax , 
        const WorkSpace&  ws   ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
       *  @param tag unique tag 
       *  @param f2 the function 
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      static double integrate
      ( const std::size_t tag  , 
        function2         f2   , 
        const double      xmin , 
        const double      xmax ,
        const double      ymin , 
        const double      ymax ,
        const WorkSpace& /* ws */ ) ;
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
      static double  integrate
      ( const std::size_t tag  , 
        function2         f2   , 
        const double      xmin , 
        const double      xmax ,
        const double      ymin , 
        const double      ymax ) ;
      // =======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
       *  @param tag unique tag 
       *  @param f2 the function 
       *  @param y parameter y
       *  @param xmin lower integration edge in x 
       *  @param xmax upper integration edge in x 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      static double integrateX
      ( const std::size_t tag  , 
        function2         f2   , 
        const double      y    , 
        const double      xmin ,
        const double      xmax ,
        const WorkSpace&  ws   ) ;
      // ======================================================================
      /** calculate the integral 
       *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
       *  @param tag unique tag 
       *  @param f2 the function 
       *  @param x parameter x
       *  @param ymin lower integration edge in y 
       *  @param ymax upper integration edge in y 
       *  @param integration workspace (not used)
       *  @return the value of the integral 
       */
      static double  integrateY
      ( const std::size_t tag  , 
        function2         f2   , 
        const double      x    , 
        const double      ymin ,
        const double      ymax ,
        const WorkSpace&  ws   ) ;
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
