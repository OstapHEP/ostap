// ============================================================================
#ifndef OSTAP_ELLIPTIC_H 
#define OSTAP_ELLIPTIC_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <cstdint>
#include <complex>
#include <cmath>
// ============================================================================
// Ostap
// =============================================================================
/** @file Ostap/Elliptis.h
 *  Elliptic functions and integrals 
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
    // Elliptic functions and integrals 
    // ========================================================================

    // ========================================================================
    // Complete Elliptic integrals 
    // ========================================================================

    // ========================================================================
    /** Complete elliptic integral \f$ K(k) \f$  
     *  \[ K(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     */
    double elliptic_K 
    ( const double k   ) ;
    // ========================================================================
    /** Complete elliptic integral \f$ E(k) \f$  
     *  \[ E(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     */
    double elliptic_E 
    ( const double k   ) ;    
    // ========================================================================
    /** Complete elliptic integral \f$ K[m] \f$  as function of parameter m 
     *  \f[ K[m] = K(k) = 
     *     \frac{\pi}{2 \mathrm{agm} (1, \sqrt{1-k^2}} = 
     *     \frac{\pi}{2 \mathrm{agm} (1, \sqrt{1-m}  } = 
     *     \frac{\pi}{2 \mathrm{agm} (1, \sqrt{ m^{\prime}}} \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     */
    double elliptic_Km 
    ( const double m   ) ;        
    // ========================================================================
    /** Complete elliptic integral \f$ E[m] \f$ as function of parameter m  
     *  \[ E(,) \equiv F ( \frac{\pi}{2}, m ) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  \f[ E(m) \equiv 2 R_{G}(0, 1 - m , 1 ) \f]
     *  @see Eq. (55) in arXiv:math/9409227
     */
    double elliptic_Em 
      ( const double m   ) ;    
    // ========================================================================
    /** Complete elliptic integral \f$ K(k) \f$  
     *  \[ K(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @attention GSL is used for calculation 
     */
    double elliptic_K_gsl 
    ( const double k   ) ;
    // ========================================================================
    /** Complete elliptic integral \f$ E(k) \f$  
     *  \[ E(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @attention GSL is used for calculation 
     */
    double elliptic_E_gsl 
    ( const double k   ) ;
    // ========================================================================

    // ========================================================================
    // Incomplete Elliptic integrals 
    // ========================================================================

    // ========================================================================
    /** Trigonometric form of incomplete elliptic integral \f$ F(\phi,k) \f$
     *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \dfrac{ d \psi }{\sqrt{1-k^2 \sin^2 \psi }}\f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Eq. (59) in arXiv:math/9409227
     */
    double elliptic_F
    ( const double phi , 
      const double k   ) ;
    // ========================================================================
    /** Trigonometric form of incomplete elliptic integral \f$ E(\phi,k) \f$
     *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \sqrt{1-k^2 \sin^2 \psi } d \psi \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Eq. (60) in arXiv:math/9409227
     */
    double elliptic_E
    ( const double phi , 
      const double k   ) ;
    // ========================================================================
    /** Trigonometric form of incomplete elliptic integral \f$ F(\phi,m) \f$
     *  \f[ F(\phi,m) \equiv \int_{0}^{\phi} \dfrac{ d \psi }{\sqrt{1-m \sin^2 \psi }}\f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Eq. (59) in arXiv:math/9409227
     */
    double elliptic_Fm
    ( const double phi , 
      const double m   ) ;    
    // ========================================================================
    /** Trigonometric form of incomplete elliptic integral \f$ E(\phi,m) \f$
     *  \f[ F(\phi,m) \equiv \int_{0}^{\phi} \sqrt{1-m \sin^2 \psi } d \psi \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Eq. (60) in arXiv:math/9409227
     */
    double elliptic_Em
    ( const double phi , 
      const double m   ) ;    
    // ========================================================================
    /** Trigonometric form of incomplete elliptic integral \f$ F(\phi,k) \f$
     *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \dfrac{ d \psi }{\sqrt{1-k^2 \sin^2 \psi }}\f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     */
    double elliptic_F_gsl 
    ( const double phi , 
      const double k   ) ;
    // ========================================================================
    /** Trigonometric form of incomplete elliptic integral \f$ E(\phi,k) \f$
     *  \f[ F(\phi,k) \equiv \int_{0}^{\phi} \sqrt{1-k^2 \sin^2 \psi } d \psi \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     */
    double elliptic_E_gsl 
    ( const double phi , 
      const double k   ) ;    
    // ========================================================================
    /** difference in complete elliptic integrals  \f$ K(k) \f$ and \f$ E(k) \f$
     *  \f[ K(k) - E(k) = \frac{k^2}{3}R_D\left(0,1-k^2,1\right)\f],
     *  where \f$ R_D(x,y,z)\f$ is a symmetric Carlson form 
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double elliptic_KmE 
    ( const double k   ) ;    
    // ========================================================================
    /** Jacobi zeta function
     *  \f[ K(k) Z( \beta , k ) = K(k) E(\beta, k ) - E(k) F(\beta,k) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  http://functions.wolfram.com/EllipticIntegrals/JacobiZeta/introductions/IncompleteEllipticIntegrals/ShowAll.html
     */
    double elliptic_Z  
    ( const double beta , 
      const double k    ) ;
    // ========================================================================
    /** Product of Jacobi zeta function \f$ Z(\beta,k) \f$
     *  and complete elliptic integral \f$ K(k) \f$
     *  \f[ K(k) Z( \beta , k ) = \frac{k^2}{3} \sin \beta \cos \beta 
     *   \sqrt{ 1 - k^2\sin^2\beta } R_J\left(0,1-k^2, 1 , 1-k^2\sin^2\beta\right)\f], 
     *  where \f$ R_J(x,y,z,t)\f$ is a symmetric Carlson form  
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double elliptic_KZ 
    ( const double beta  , 
      const double k   ) ;
    // ========================================================================
    /** elliptic \f$ \Pi(\alpha^2,k)\f$ function 
     *  - \f$ alpha^2 < 1 \f$ 
     *  - \f$ k      < 1 \f$ 
     *  \f[ \Pi(\alpha^2, k) - K(k) = 
     *   \frac{1}{3}\alpha^2 R_J( 0, 1-k^2, 1 , 1 - \alpha^2) \f] 
     */ 
    double elliptic_PI
    ( const double alpha2 , 
      const double k      ) ;
    // ========================================================================
    /** elliptic \f$ \Pi(\alpha^2,k) - K(k) \f$ function 
     *  \f[ \Pi(\alpha^2, k) - K(k) \equiv  
     *   \frac{1}{3}\alpha^2 R_J( 0, 1-k^2, 1 , 1 - \alpha^2) \f] 
     *  - \f$ alpha^2 < 1 \f$ 
     *  - \f$ k      < 1 \f$ 
     */ 
    double elliptic_PImK  
    ( const double alpha2 , 
      const double k      ) ;
    // ========================================================================

    // ========================================================================
    // Symmetric Carlson forms 
    // ========================================================================
    
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_F(x,y,z) = \int_0^{+\infty} 
     *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2}dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RF 
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_F(x,y,z) = \int_0^{+\infty} 
     *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2}dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RF_gsl  
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_F(x,y,z) = \int_0^{+\infty} 
     *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2}dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RF_int 
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_J(x,y,z,p) = \int_0^{+\infty} 
     *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2} (t+p) dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RJ
    ( const double x , 
      const double y , 
      const double z , 
      const double p ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_J(x,y,z,p) = \int_0^{+\infty} 
     *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2} (t+p) dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RJ_gsl
    ( const double x , 
      const double y , 
      const double z , 
      const double p ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_J(x,y,z,p) = \int_0^{+\infty} 
     *   \left[  (t+x)(t+y)(t+z) \right]^{-1/2} (t+p) dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RJ_int
    ( const double x , 
      const double y , 
      const double z , 
      const double p ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_C(x,y) = R_F(x,y,y) = \int_0^{+\infty} 
     *   (t+x)^{-1/2}(t+y)^{-1}dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     *  For negative y, Cauchy principal value it returned 
     *  \f[ R_C(x,-y) = \left(\frac{x}{x+y}\right)R_C(x+y,y), 0 < y \f] 
     */
    double carlson_RC
    ( const double x , 
      const double y ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_C(x,y) = R_F(x,y,y) = \int_0^{+\infty} 
     *   (t+x)^{-1/2}(t+y)^{-1}dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     *  For negative y, Cauchy principal value it returned 
     *  \f[ R_C(x,-y) = \left(\frac{x}{x+y}\right)R_C(x+y,y), 0 < y \f] 
     */
    double carlson_RC_gsl 
    ( const double x , 
      const double y ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_C(x,y) = R_F(x,y,y) = \int_0^{+\infty} 
     *   (t+x)^{-1/2}(t+y)^{-1}dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     *  For negative y, Cauchy principal value it returned 
     *  \f[ R_C(x,-y) = \left(\frac{x}{x+y}\right)R_C(x+y,y), 0 < y \f] 
     */
    double carlson_RC_int
    ( const double x , 
      const double y ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_D (x,y,z) = R_J(x,y,z,z) = \int_0^{+\infty} 
     *  \left[  (t+x)(t+y)\right]^{-1/2} (t+z)^{-3/2} dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RD 
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_D (x,y,z) = R_J(x,y,z,z) = \int_0^{+\infty} 
     *  \left[  (t+x)(t+y)\right]^{-1/2} (t+z)^{-3/2} dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RD_gsl  
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_D (x,y,z) = R_J(x,y,z,z) = \int_0^{+\infty} 
     *  \left[  (t+x)(t+y)\right]^{-1/2} (t+z)^{-3/2} dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RD_int  
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_G (x,y,z) = \rac{1}{4} \int_0^{+\infty} 
     *  \left[  (t+x)(t+y)\right]^{-1/2} 
     *   \left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right)  dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RG 
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_G (x,y,z) = \rac{1}{4} \int_0^{+\infty} 
     *  \left[  (t+x)(t+y)\right]^{-1/2} 
     *   \left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right)  dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RG_gsl 
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_G (x,y,z) = \rac{1}{4} \int_0^{+\infty} 
     *  \left[  (t+x)(t+y)\right]^{-1/2} 
     *   \left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right)  dt \f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RG_int 
    ( const double x , 
      const double y , 
      const double z ) ;
    // ========================================================================
    // specific cases of symmetric forms 
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_F(x,y) = R_F(x,y,0)\f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RF
    ( const double x , 
      const double y ) ;
    // ========================================================================
    /** Symmetric Carlson form 
     *  \f[ R_G(x,y) = R_G(x,y,0)\f]
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     *  @see Carlson, B.C., "Numerical computation of real or complex elliptic integrals", 
     *                Numerical Algorithms, 10, 1995,  13
     *  @see https://doi.org/10.1007/BF02198293
     *  @see https://arxiv.org/abs/math/9409227
     */
    double carlson_RG 
    ( const double x , 
      const double y ) ;
    // ========================================================================
    
    // ========================================================================
    // Jacobi elliptic functions (real argument)
    // ========================================================================

    // ========================================================================
    /** Elliptic amplitude \f$ \mathrm{am}(u,m)=\phi\f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    double am
    ( const double u ,
      const double m ) ;
    // ========================================================================
    /** Elliptic delta amplitude \f$ \mathrm{sn} (u,m)=\frac{d}{du} \mathrm{am} ( u, m ) \f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    double dn
    ( const double u ,
      const double m ) ;
    // ========================================================================
    /** Elliptic sine amplitude \f$ \mathrm{sn} (u,m)=\sin \mathrm{am} ( u, m ) \f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    double sn
    ( const double u ,
      const double m ) ;
    // ========================================================================
    /** Elliptic sine amplitude \f$ \mathrm{sn} (u,m)=\sin \mathrm{am} ( u, m ) \f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    double sn_
    ( const double u ,
      const double m ) ;
    // ========================================================================
    /** Elliptic cosine amplitude \f$ \mathrm{sn} (u,m)=\cos \mathrm{am} ( u, m ) \f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    double cn
    ( const double u ,
      const double m ) ;
    // ========================================================================
    /** elliptic function 
     *  \f[ \matmrm{sc}\,(u,m) = \frac{ \mathrm{sn} ( u, m) } { \mathrm{cn} ( u , m ) } \f] 
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */ 
    double sc
    ( const double u ,
      const double m ) ;
    // ========================================================================

    // ========================================================================
    // Jacobi elliptic functions (complex argument)
    // ========================================================================
    
    // ========================================================================
    /** Elliptic sine amplitude \f$ \mathrm{sn} (u,m)=\sin \mathrm{am} ( u, m ) \f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    std::complex<double> sn
    ( const std::complex<double>& u ,
      const double                m ) ;
    // ========================================================================
    /** Elliptic cosine amplitude \f$ \mathrm{sn} (u,m)=\cos \mathrm{am} ( u, m ) \f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    std::complex<double> cn
    ( const std::complex<double>& u ,
      const double                m ) ;    
    // ========================================================================
    /** Elliptic delta amplitude \f$ \mathrm{sn} (u,m)=\frac{d}{du} \mathrm{am} ( u, m ) \f$, where 
     *  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */
    std::complex<double> dn
    ( const std::complex<double>& u ,
      const double                m ) ;
    // ========================================================================
    /** elliptic function 
     *  \f[ \matmrm{sc}\,(u,m) = \frac{ \mathrm{sn} ( u, m) } { \mathrm{cn} ( u , m ) } \f] 
     *  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
     */ 
    std::complex<double> sc
    ( const std::complex<double>& u ,
      const double                m ) ;
    // ========================================================================
    
    // ========================================================================
    // Dixon elliptic functions 
    // ========================================================================

    // ========================================================================
    /** Dixon (or Dixonian) elliptic function cm for real argument 
     *  @see https://en.wikipedia.org/wiki/Dixon_elliptic_functions
     *  @param x argument 
     *  @return value of Dixon elliptic function cm 
     */
    double cm ( const double x ) ;
    // ========================================================================
    /** Dixon (or Dixonian) elliptic function sm for real arguemnt 
     *  @see https://en.wikipedia.org/wiki/Dixon_elliptic_functions     
     *  @param x argument
     *  @return value of Dixon elliptic function sm 
     */
    double sm ( const double x ) ;
    // ========================================================================
    /** Dixon (or Dixonian) elliptic function cm for complex argument 
     *  @see https://en.wikipedia.org/wiki/Dixon_elliptic_functions     
     *  @param z argument
     *  @return value of Dixon elliptic function sm 
     */
    std::complex<double> cm ( const std::complex<double>& z ) ;
    // ========================================================================
    /** Dixon (or Dixonian) elliptic function sm for complex argument 
     *  @see https://en.wikipedia.org/wiki/Dixon_elliptic_functions     
     *  @param z argument
     *  @return value of Dixon elliptic function sm 
     */
    std::complex<double> sm ( const std::complex<double>& z ) ;
    // ========================================================================

    // ========================================================================
    /** Lemniscate elliptic function cl for real argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param x the argument
     *  @return the value of lemniscate elliptic functinon cl
     */
    double cl ( const double x ) ;
    // ========================================================================
    /** Lemniscate elliptic function sl for real argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param x the argument
     *  @return the value of lemniscate elliptic function sl
     */
    double sl ( const double x ) ;
    // ========================================================================
    /** Lemniscate elliptic function cl for complex argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param x the argument
     *  @return the value of lemniscate elliptic function cl
     */
    std::complex<double> cl ( const std::complex<double>& z ) ;  
    // ========================================================================
    /** Lemniscate elliptic function sl for complex argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param  z the argument
     *  @return the value of lemniscate elliptic function sl
     */
    std::complex<double> sl ( const std::complex<double>& z ) ;  
    // ========================================================================

    // ========================================================================
    /** Lemniscate elliptic function hyperbolic cosine clh for real argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param x the argument
     *  @return the value of lemniscate elliptic functinon clh
     */
    double clh ( const double x ) ;
    // ========================================================================
    /** Lemniscate elliptic function hyperbolic sine slh for real argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param x the argument
     *  @return the value of lemniscate elliptic functinon slh
     */
    double slh ( const double x ) ;
    // ========================================================================
    /** Lemniscate elliptic hyperboilic cosine function clh for complex argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param x the argument
     *  @return the value of lemniscate elliptic function clh
     */
    std::complex<double> clh ( const std::complex<double>& z ) ;  
    // ========================================================================
    /** Lemniscate elliptic hyperbolic sine function slh for complex argument 
     *  @see https://en.wikipedia.org/wiki/Lemniscate_elliptic_functions
     *  @param  z the argument
     *  @return the value of lemniscate elliptic function slh
     */
    std::complex<double> slh ( const std::complex<double>& z ) ;  
    // ========================================================================
    
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_ELLIPTIC_H
// ============================================================================
