// ============================================================================
#ifndef OSTAP_INTEGRATOR_H 
#define OSTAP_INTEGRATOR_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
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
       *  @param xmax uppr  integration edge
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
