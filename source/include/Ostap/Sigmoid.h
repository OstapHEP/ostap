// ============================================================================
#ifndef OSTAP_SIGMOID_H 
#define OSTAP_SIGMOID_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <cmath>
// ============================================================================
// Ostap
// =============================================================================
#include "Ostap/Math.h"
// ============================================================================
/** @file Ostap/Sigmoid.h 
 *  Sigmoid functions 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    
    // ========================================================================
    // Sigmoid/kink functions & friends  
    // ========================================================================
    
    // ========================================================================
    /** smooth transition function 
     *   \f[ \phix() = left\{ 
     *   \begin{array}{ll}
     *     0 &  x\le a \			\
     *     1 &  x\ge b \\ 
     *    smooth & 
     *    \end{array} \right. \f] 
     */
    double smooth_transition 
    ( const double x     , 
      const double a = 0 ,
      const double b = 1 ) ;
    // ========================================================================
    
    // ========================================================================
    /** smooth (polynomial) step function
     *  @see https://en.wikipedia.org/wiki/Smoothstep
     *  Transition function for \f$ 0 \le x \le 1\f$ 
     *  - \f$ f(x)=0\f$ for  \f$ x \le 0 \f$
     *  - \f$ f(x)=1\f$ for  \f$ x \ge 1 \f$
     *  - \f$ f(x)\f$ if a \f$ 2n+1 \f$ polynomial fuction inbetween      
     *  @param x variable
     *  @param n index, polynomial of order \f$ 2n+1 \f$     
     *  @see Ostap::Math::smooth_transtion 
     *  @see Ostap::Math::clamp 
     */
    double smooth_step
    ( const double         x     , 
      const unsigned short n = 1 ) ;
    // ========================================================================


    /// the Sigmoid type 
    enum class SigmoidType
    {
      //
      Logistic              , // based on logistic function 
      Hyperbolic            , // based on tanh 
      Trigonometric         , // based on atan 
      Error                 , // Based on error function 
      Gudermannian          , // Based on Gudermannian function
      Algebraic             , // 0.5 * ( 1 + 2 * x / hypot ( 1 , 2*x ) )
      SmoothTransition      , // Based on "smooth transition" function
      //
      Polynomial_n0         , // Based on "smooth step" with n=0
      Polynomial_n1         , // Based on "smooth step" with n=1
      Polynomial_n2         , // Based on "smooth step" with n=2
      Polynomial_n3         , // Based on "smooth step" with n=3
      Polynomial_n4         , // Based on "smooth step" with n=4
      Polynomial_n5         , // Based on "smooth step" with n=5
      Polynomial_n6         , // Based on "smooth step" with n=6
      //
      Sine                  , // based on sine-function
      AbsAlgebraic          , // 0.5  + x / ( 1 + abs  ( 2 * x ) ) 
      //
      First = Logistic      , 
      Last  = AbsAlgebraic  ,  
    } ;      	
    // ========================================================================
    /** sigmoid type
     *  @param  name the case-insensitive name of of sigmoid function
     *  @return ID   of sigmoid function
     */
    SigmoidType sigmoid_type
    ( const std::string& name = "Hyperbolic" ) ;
    // ========================================================================
    /*  sigmoid type
     *  @param  sigmoid ID 
     *  @return the name of of sigmoid function
     */
    // ========================================================================
    std::string sigmoid_name 
    ( const SigmoidType stype ) ;
    // ======================================================================= 
    /** Sigmoid functions
     *  All sigmoid fuctions \f$ \sigma(z) \f$ are normalized & scaled such
     *  - \f$ \sigma(-\infty) =0\f$ 
     *  - \f$ \sigma(+\infty) =1\f$ 
     *  - \f$ \sigma^\prime(0)=1\f$
     * @see Ostap::Math::SigmoidType
     */
    double sigmoid
    ( const double      x                           ,
      const SigmoidType t = SigmoidType::Hyperbolic ) ;
    // ========================================================================

    // ========================================================================
    /** Gompertz' curve/sigmoid
     *  \f$ f(x;a,c,x_0) = \left|a|right|\mathrm{e}^{ - \mathrm{e}^{ - c ( x - x_0) } }\f$
     *  @see https://en.wikipedia.org/wiki/Gompertz_function
     *  - \f$ b = \mathrm{e}^{cx_0}\f$ 
     */
    double gompertz
    ( const double x      ,
      const double a  = 1 ,
      const double c  = 1 ,
      const double x0 = 0 ) ;
    // ========================================================================
    
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_SIGMOID_H   
// ============================================================================
