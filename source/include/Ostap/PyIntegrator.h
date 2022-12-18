#ifndef OSTAP_PYINTEGRATOR_H 
#define OSTAP_PYINTEGRATOR_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
// ============================================================================
// forward declarations 
// ============================================================================
namespace Ostap { namespace Math      { class Integrator  ; } }
namespace Ostap { namespace Math      { class Workspace   ; } }
namespace Ostap { namespace Functions { class PyCallable  ; } }
namespace Ostap { namespace Functions { class PyCallable2 ; } }
namespace Ostap { namespace Functions { class PyCallable3 ; } }
// ============================================================================
/** @file Ostap/PyIntegrator.h
 *  Helper function for back-propagation of Ostap::Math::Integrator for 
 *  the  old verisons of ROOT 
 *  @see Ostap:Math::Integrator 
 *  @see Ostap:Math::PyCallable
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
#if ROOT_VERSION_CODE<ROOT_VERSION(6,18,0)
    // ========================================================================
    /** calculate the integral 
     *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
     *  @param f1 the function 
     *  @param xmin lower integration edge 
     *  @param xmax upper  integration edge
     *  @param tag  uqniue label/tag  
     *  @param rescale rescale function for better numerical precision  
     *  @return the value of the integral 
     *  @see Ostap::Math::Integrator::integrate 
     */
    double py_integrate
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,           
      const double                        xmin           , 
      const double                        xmax           ,
      const std::size_t                   tag        = 0 , 
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 , 
      const int                           rule       = 0 ) ;
    // ========================================================================
    /** calculate the integral 
     *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
     *  @param f1 the function 
     *  @return the value of the integral 
     */
    double py_integrate_infinity
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,           
      const std::size_t                   tag        = 0 , 
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 ) ;
    // ========================================================================
    /** calculate the integral 
     *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
     *  @param f1 the function 
     *  @param xmin lower integration edge 
     *  @return the value of the integral 
     */
    double py_integrate_to_infinity
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,           
      const double                        xmin           , 
      const std::size_t                   tag        = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 ) ;
    // ========================================================================
    /** calculate the integral 
     *  \f[ r = \mathcal{P}\int_{-\infty}^{x_{max}} f_1(x) dx \f]
     *  @param f1 the function 
     *  @param xmax upper  integration edge
     *  @param integration workspace 
     *  @param tag     unuque label/tag 
     *  @return the value of the integral 
     */
    double py_integrate_from_infinity
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,           
      const double                        xmax           , 
      const std::size_t                   tag        = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 ) ;
    // ========================================================================
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
    double py_cauchy_pv 
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,           
      const double                        c              , 
      const double                        xmin           , 
      const double                        xmax           , 
      const std::size_t                   tag        = 0 , 
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 ) ;
    // ========================================================================
    /** get Cauchy principal value integral 
     *  \f[ g(c) =  \mathcal{P} \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
     *  @param f1 the function 
     *  @param c  the parameter 
     *  @param xmin lower integration edge 
     *  @param tag     unuque label/tag 
     *  @param rescale rescale function for better numerical treatmenrt 
     *  @return the value of the integral 
     */
    double py_cauchy_pv_to_infinity  
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,
      const double                        c              , 
      const double                        xmin           ,  
      const std::size_t                   tag        = 0 , 
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 , 
      const double                        width      = 0 ) ;
    // ========================================================================
    /** get Cauchy principal value integral 
     *  \f[ g(c) =  \mathcal{P} \int_{-\infty}^{x_{max}}\frac{f(x)}{x-c}dx \f]
     *
     *  @param f1 the function 
     *  @param c  the parameter 
     *  @param xmax upper integration edge 
     *  @param tag     unuque label/tag 
     *  @param rescale rescale function for better numerical treatmenrt 
     *  @return the value of the integral 
     */
    double cauchy_pv_from_infinity  
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,      
      const double                        c              , 
      const double                        xmax           , 
      const std::size_t                   tag        = 0 , 
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 , 
      const double                        width      = 0 ) ;
    // ========================================================================
    /** get Cauchy principal value integral for the infinite range 
     *  \f[ g(c) = \mathcal{P} \int_{-\infty}^{+\infty}\frac{f(x)}{x-c}dx \f]
     *
     *  Integral is calculated as sum of three components
     *   - \f$ \int_{-\infty}^{a}\frac{f(x)}{x-c}dx \f$
     *   - \f$ \mathcal{P} \int_{a}^{b}\frac{f(x)}{x-c}dx \f$
     *   - \f$ int_{b}^{+\infty}\frac{f(x)}{x-c}dx \f$
     *
     *  @param f1      the function 
     *  @param c       the parameter 
     *  @param ws      integration workspace 
     *  @param tag     unuque label/tag 
     *  @param rescale rescale function for better numerical precision  
     *  @param aprecision absolute precision  (if non-positive s_APRECISION_GAWC is used) 
     *  @param aprecision relative precision  (if non-positive s_RPRECISION_GAWC is used) 
     *  @param width   width parameter
     *  @return value of the integral and the estimate of the uncertainty
     *
     *  @see Ostap::Math::Integrator::cauchy_pv
     *  @see Ostap::Math::Integrator::integrate_to_infinity 
     *  @see Ostap::Math::Integrator::integrate_from_infinity 
     *  @see Ostap::Math::Integrator::cauchy_pv_b
     *  @see Ostap::Math::Integrator::cauchy_pv_b
     */
    double py_cauchy_pv_infinity 
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,
      const double                        c              , 
      const std::size_t                   tag        = 0 ,
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 , 
      const double                        width      = 0 ) ;
    // ======================================================================== 
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
    double py_kramers_kronig 
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,
      const double                        s              , 
      const double                        xmin           , 
      const unsigned short                n          = 0 , 
      const std::size_t                   tag        = 0 ,
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 , 
      const double                        width      = 0 ) ;
    // ========================================================================
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
    double py_integrate_singular
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,
      const double                        xmin           , 
      const double                        xmax           ,
      const std::vector<double>&          points         ,
      const std::size_t                   tag        = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 ) ;
    // ========================================================================
    /** calculate the integral using double adaptive CQUAD integrator  
     *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
     *  @param f1 the function 
     *  @param xmin lower integration edge 
     *  @param xmax upper  integration edge
     *  @param tag  uqniue label/tag  
     *  @param rescale rescale function for better numerical precision  
     *  @return the value of the integral 
     */
    double py_integrate_cquad
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,    
      const double                        xmin           , 
      const double                        xmax           ,
      const std::size_t                   tag        = 0 , 
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 ) ;
    // ========================================================================      
    /** calculate the integral using Romberg integrator  
     *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
     *  @param f1 the function 
     *  @param xmin lower integration edge 
     *  @param xmax upper  integration edge
     *  @param tag  uqniue label/tag  
     *  @param rescale rescale function for better numerical precision  
     *  @return the value of the integral 
     */
    double py_integrate_romberg
    ( const Ostap::Math::Integrator&      I              , 
      const Ostap::Functions::PyCallable& f1             ,      
      const double                        xmin           , 
      const double                        xmax           ,
      const std::size_t                   tag        = 0 , 
      const unsigned short                rescale    = 0 ,
      const double                        aprecision = 0 , 
      const double                        rprecision = 0 ) ;
    // =========================================================================
    // Integration for 2D-functions 
    // =========================================================================
    /** calculate the integral 
     *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
     *  @param f2 the function 
     *  @param xmin lower integration edge in x 
     *  @param xmax upper integration edge in x 
     *  @param ymin lower integration edge in y 
     *  @param ymax upper integration edge in y 
     *  @return the value of the integral 
     */
    double py_integrate2
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable2& f2             ,      
      const double                         xmin           , 
      const double                         xmax           ,
      const double                         ymin           , 
      const double                         ymax           ,
      const std::size_t                    tag        = 0 ,
      const double                         aprecision = 0 ,
      const double                         rprecision = 0 ) ;
    // =========================================================================
    /** calculate the integral 
     *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
     *  @param f2 the function 
     *  @param y parameter y
     *  @param xmin lower integration edge in x 
     *  @param xmax upper integration edge in x 
     *  @return the value of the integral 
     */
    double py_integrate2X
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable2& f2             ,      
      const double                         y              , 
      const double                         xmin           ,
      const double                         xmax           ,  
      const std::size_t                    tag        = 0 ,
      const unsigned short                 rescale    = 0 ,        
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 , 
      const int                            rule       = 0 ) ;
    // ========================================================================
    /** calculate the integral 
     *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
     *  @param f2 the function 
     *  @param x parameter x
     *  @param ymin lower integration edge in y 
     *  @param ymax upper integration edge in y 
     *  @return the value of the integral 
     */
    double py_integrate2Y
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable2& f2             ,
      const double                         x              , 
      const double                         ymin           ,
      const double                         ymax           ,
      const std::size_t                    tag        = 0 ,
      const unsigned short                 rescale    = 0 ,        
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 ,
      const int                            rule       = 0 ) ;
    // ========================================================================
    // 3D integration 
    // ========================================================================
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
     *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE is used) 
     *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE is used) 
     *  @return the value of the integral and the estimate of the error 
     */
    double py_integrate3
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable3& f3             ,
      const double                         xmin           , 
      const double                         xmax           ,
      const double                         ymin           , 
      const double                         ymax           ,
      const double                         zmin           , 
      const double                         zmax           ,
      const std::size_t                    tag        = 0 ,
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 ) ;
    // ========================================================================
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
    double py_integrate3XY
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable3& f3             ,
      const double                         z              , 
      const double                         xmin           ,
      const double                         xmax           ,
      const double                         ymin           ,
      const double                         ymax           ,
      const std::size_t                    tag        = 0 , 
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 ) ;
    // ========================================================================    
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
    double py_integrate3XZ
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable3& f3             ,
      const double                         y              , 
      const double                         xmin           ,
      const double                         xmax           ,
      const double                         zmin           ,
      const double                         zmax           ,
      const std::size_t                    tag        = 0 , 
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 ) ;
    // ========================================================================  
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
    double py_integrate3YZ
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable3& f3             ,
      const double                         x              , 
      const double                         ymin           ,
      const double                         ymax           ,
      const double                         zmin           ,
      const double                         zmax           ,
      const std::size_t                    tag        = 0 , 
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 ) ;
    // ======================================================================== 
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
    double py_integrate3X
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable3& f3             ,
      const double                         y              , 
      const double                         z              , 
      const double                         xmin           ,
      const double                         xmax           ,
      const std::size_t                    tag        = 0 , 
      const unsigned short                 rescale    = 0 ,
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 , 
      const int                            rule       = 0 ) ;
    // ========================================================================
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
    double py_integrate3Y
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable3& f3             ,
      const double                         x              , 
      const double                         z              , 
      const double                         ymin           ,
      const double                         ymax           ,
      const std::size_t                    tag        = 0 , 
      const unsigned short                 rescale    = 0 ,
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 , 
      const int                            rule       = 0 ) ;
    // ========================================================================    
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
    double py_integrate3Z
    ( const Ostap::Math::Integrator&       I              , 
      const Ostap::Functions::PyCallable3& f3             ,
      const double                         x              , 
      const double                         y              , 
      const double                         zmin           ,
      const double                         zmax           ,
      const std::size_t                    tag        = 0 , 
      const unsigned short                 rescale    = 0 ,
      const double                         aprecision = 0 , 
      const double                         rprecision = 0 , 
      const int                            rule       = 0 ) ;
    // ========================================================================
#endif
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                The end osf namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYINTEGRATOR_H
// ============================================================================
