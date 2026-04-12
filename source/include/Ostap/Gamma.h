// ============================================================================
#ifndef OSTAP_GAMMA_H 
#define OSTAP_GAMMA_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <cstdint>
#include <complex>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    // gamma & friends 
    // ========================================================================

// ========================================================================
    /** Gamma function
     *  \f$ \Gamma ( x )\f$ 
     *  @see std::tgamma 
     */
    inline double gamma   
    ( const double x ) { return std::tgamma ( x ) ; }
    // ========================================================================
    /** Gamma function
     *  \f$ \Gamma ( x )\f$ 
     *  @see Ostapo::Math::tgamma
     *  @see std::tgamma  
     */
    inline long double gamma   
    ( const long double x ) { return std::tgamma ( x ) ; }
    // ========================================================================
    /** Gamma function of complex argument 
     *  \f$ \Gamma ( x ) \f$ 
     * @see Ostap::Math::lgamma 
     */
    std::complex<double> gamma ( const std::complex<double>& z ) ;
    // ========================================================================
    /** Gamma function of complex argument 
     *  \f$ \Gamma ( x ) \f$ 
     * @see Ostap::Math::lgamma 
     */
    std::complex<long double> gamma ( const std::complex<long double>& z ) ;
    // ========================================================================
   
    // ========================================================================
    /** logarithm of gamma function
     *  \f$ \log \Gamma ( x ) \f$ 
     */
    inline double lgamma  
    ( const double x ) { return std::lgamma ( x ) ; }
    // ========================================================================
    /** logarithm of gamma function
     *  \f$ \log \Gamma ( x ) \f$ 
     */
    inline long double lgamma  
    ( const long double x ) { return std::lgamma ( x ) ; }
    // ========================================================================
    /** Logarithm of gamma function of complex argument 
     *  \f$ \log \Gamma ( x ) \f$  \f%  -\pi < arg < +\pi \f$ 
     * 
     *   Note that the imaginary part (arg) is not well-determined 
     *   when |z| is very large, due to inevitable roundoff in restricting 
     *   to (-\pi,\pi]. This will result in a GSL_ELOSS error when it occurs. 
     *   The real part (lnr), however, never suffers from loss of precision.
     */
    std::complex<double>      lgamma ( const std::complex<double>&      z ) ;
    // =========================================================================
    /** Logarithm of gamma function of complex argument 
     *  \f$ \log \Gamma ( x ) \f$  \f%  -\pi < arg < +\pi \f$ 
     * 
     *   Note that the imaginary part (arg) is not well-determined 
     *   when |z| is very large, due to inevitable roundoff in restricting 
     *   to (-\pi,\pi]. This will result in a GSL_ELOSS error when it occurs. 
     *   The real part (lnr), however, never suffers from loss of precision.
     */
    std::complex<long double> lgamma ( const std::complex<long double>& z ) ;
    // ========================================================================
  
    // ========================================================================
    /** Compute inverse Gamma function 
     *  \f[ f(x) = \frac{1}{\Gamma(x)} \f]
     *  @return the value of inverse Gamma functions 
     */
    double               igamma ( const double x ) ;
    // =======================================================================
    /** Compute inverse Gamma function 
     *  \f[ f(x) = \frac{1}{\Gamma(x)} \f]
     *  @return the value of inverse Gamma functions 
     */
    long double          igamma ( const long double x ) ;
    // =======================================================================
    /** Compute inverse Gamma function 
     *  \f[ f(n) = \frac{n}{\Gamma(n)} \f]
     *  @return the value of inverse Gamma functions 
     */
    double               igamma ( const int n ) ;
    // ========================================================================
    /** Compute inverse Gamma function 
     *  \f[ f(n) = \frac{n}{\Gamma(n)} \f]
     *  @return the value of inverse Gamma functions 
     */
    double               igamma ( const unsigned int n ) ;
    // ========================================================================
    /** Compute inverse Gamma function 
     *  \f[ f(x) = \frac{1}{\Gamma(x)} \f]
     *  @return the value of inverse Gamma functions 
     */
    std::complex<double> igamma ( const std::complex<double>& z ) ;
    // ========================================================================
    /** Compute inverse Gamma function 
     *  \f[ f(x) = \frac{1}{\Gamma(x)} \f]
     *  @return the value of inverse Gamma functions 
     */
    std::complex<long double> igamma ( const std::complex<long double>& z ) ;
    // ========================================================================

    // ========================================================================
    /** factorial function
     *   \f[ n! = \Gamma ( n + 1 ) \f]
     *  - Factorials are precomputed up to the max
     *  @attention \f$ +\infty\f$  is returned for n > 170 \f$ 
     */
    double factorial ( const std::uint8_t n ) ;
    // ========================================================================

    // ========================================================================
    /** Get the n-th derivative of gamma function at x=1 
     *   \f[  f(n) = \left. \frac{d^{n} \Gamma(x) }{dx^n} \right|_{x=1} \f]  
     *  @see https://www.researchgate.net/publication/228557691_Evaluation_of_higher-order_derivatives_of_the_Gamma_function
     *  @see Choi, Junesang & Srivastava, Hari. (2000). 
     *       Evaluation of higher-order derivatives of the Gamma function. Publ. Elektrotehn. Fak. Ser. Mat. 11. 9-18. 
     */
    double dgamma_at_1 ( const unsigned short n ) ; 
    // ========================================================================


    // ========================================================================
    /** Compute psi function 
     *  \f[ f(x) = \frac{d}{dx}\ln \Gamma(x)\f]
     * @return the value of psi function 
     * @see Ostap::Math::digamma
     * @see Ostap::Math::trigamma
     */
    double psi ( const double x ) ;
    // =======================================================================
    /** Compute psi function 
     *  \f[ f(x) = \frac{d}{dx}\ln \Gamma(x)\f]
     * @return the value of psi function 
     * @see Ostap::Math::digamma
     * @see Ostap::Math::trigamma
     */
    long double psi ( const long double x ) ;
    // =======================================================================
  
    // =======================================================================    
    /** compute psi/polygamma function 
     * \f$ \psi^{(n)}(x) = \left(\frac{d}{dx}\right)^{(n)}\psi(x) = 
     * = \left(\frac{d}{dx}\right)^{(n+1)} \log \Gamma (x) \f$ 
     * @return value of polygamma function
     * @see Ostap::Math::polygamma
     * @see Ostap::Math::digamma
     * @see Ostap::Math::trigamma
     */
    double psi  
    ( const double         x , 
      const unsigned short n ) ;
    // =======================================================================    
    /** compute psi/polygamma function 
     * \f$ \psi^{(n)}(x) = \left(\frac{d}{dx}\right)^{(n)}\psi(x) = 
     * = \left(\frac{d}{dx}\right)^{(n+1)} \log \Gamma (x) \f$ 
     * @return value of polygamma function
     * @see Ostap::Math::polygamma
     * @see Ostap::Math::digamma
     * @see Ostap::Math::trigamma
     */
    long double psi  
    ( const long double    x , 
      const unsigned short n ) ;
    // =======================================================================

    // =======================================================================
    /** Compute digamma function 
     *  \f[ f(x) = \psi^{(0)}(x) = \frac{d}{dx}\ln \Gamma(x)\f]
     * @return the value of psi function 
     * @see Ostap::Math::psi
     * @see Ostap::Math::digamma
     * @see Ostap::Math::polygamma
     */
    inline double digamma 
    ( const double x ) 
    { return psi ( x ) ; }
    // ========================================================================
    /** Compute trigamma function 
     * @return the value of psi function 
     * @see Ostap::Math::psi
     * @see Ostap::Math::digamma
     * @see Ostap::Math::polygamma
     */
    inline double trigamma 
    ( const double x ) 
    { return psi ( x , 1 ) ; }
    // ========================================================================
    /** compute polygamma function 
     * \f$ \psi^{(n)}(x) = \left(\frac{d}{dx}\right)^{(n)}\psi(x) = 
     * = \left(\frac{d}{dx}\right)^{(n+1)} \log \Gamma (x) \f$ 
     * @return value of polygamma function
     * @see Ostap::Math::psi
     */
    inline double polygamma 
    ( const double         x , 
      const unsigned short n ) 
    { return psi ( x , n ) ; }
    // ========================================================================
    
    // ========================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    double beta
    ( const double x ,
      const double y ) ;
    // =======================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    long double beta
    ( const long double x ,
      const long double y ) ;
    // ========================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    double beta
    ( const unsigned short x ,
      const unsigned short y ) ;
   // =======================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    double beta
    ( const double         x ,
      const unsigned short y ) ;
  // =======================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    double beta
    ( const unsigned short y , 
      const double         x ) ; 
    // ========================================================================


    // ========================================================================
    /** natural logarithm of beta function 
     *  \f$ \log \Beta(x,y) = \log \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of logarithm of beta function 
     */
    double lnbeta
    ( const double x ,
      const double y ) ;
    // ========================================================================
    /** natural logarithm of beta function 
     *  \f$ \log \Beta(x,y) = \log \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of logarithm of beta function 
     */
    double lnbeta
    ( const unsigned short x ,
      const unsigned short y ) ;
    // ========================================================================
    /** natural logarithm of beta function 
     *  \f$ \log \Beta(x,y) = \log \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of logarithm of beta function 
     */
    double lnbeta
    ( const unsigned short x ,
      const double         y ) ;
    // ========================================================================
    /** natural logarithm of beta function 
     *  \f$ \log \Beta(x,y) = \log \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of logarithm of beta function 
     */
    double lnbeta
    ( const double         x , 
      const unsigned short y ) ;
    // ========================================================================

    // ========================================================================
    /** reciprocal beta function for 
     *  \f$ f(x,y) = \frac{1}{B(x,y)} = \frac{\Gamma(x+y)}{\Gamma(x)\Gamma(y)}\f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of reciprocal beta function 
     */
    // ========================================================================
    double ibeta
    ( const unsigned short x ,
      const unsigned short y ) ;
    // ========================================================================
    /** reciprocal beta function for 
     *  \f$ f(x,y) = \frac{1}{B(x,y)} = \frac{\Gamma(x+y)}{\Gamma(x)\Gamma(y)}\f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of reciprocal beta function 
     */
    // ========================================================================
    double ibeta
    ( const unsigned short x ,
      const double         y ) ;
    // ========================================================================
    /** reciprocal beta function for 
     *  \f$ f(x,y) = \frac{1}{B(x,y)} = \frac{\Gamma(x+y)}{\Gamma(x)\Gamma(y)}\f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of reciprocal beta function 
     */
    // ========================================================================
    double ibeta
    ( const double         x , 
      const unsigned short y ) ;
    // ========================================================================
    /** reciprocal beta function for 
     *  \f$ f(x,y) = \frac{1}{B(x,y)} = \frac{\Gamma(x+y)}{\Gamma(x)\Gamma(y)}\f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of reciprocal beta function 
     */
    // ========================================================================
    double ibeta
    ( const double         x , 
      const double         y ) ;
    // ========================================================================

    // ========================================================================    
    /** Normalized incomplete Beta function
     *  - CDF for beta-distribution
     *
     *  \f$ f ( \alpha_1,\alpha_2, z ) = 
     *      I_z( \alpha_1, \alpha_2 ) = 
     *      \frac{\Beta_z(\alpha_1,\alpha_2}}
     *           {\Beta  (\alpha_1,\alpha_2}
     *
     *  - \f$ 0 \le z \le 1\f$
     *  - \f$ 0<\alpha_1\f$
     *  - \f$ 0<\alpha_2\f$
     *  @see Ostap::Math::I
     */
    double beta_inc 
    ( const double alpha1 , 
      const double alpha2 , 
      const double z      ) ;
    // =======================================================================
    /** Normalized incomplete Beta function
     *  - CDF for beta-distribution
     *
     *  \f$ f ( \alpha_1,\alpha_2, z ) = 
     *      I_z( \alpha_1, \alpha_2 ) = 
     *      \frac{\Beta_z(\alpha_1,\alpha_2}}
     *           {\Beta  (\alpha_1,\alpha_2}
     *
     *  - \f$ 0 \le z \le 1\f$
     *  - \f$ 0<\alpha_1\f$
     *  - \f$ 0<\alpha_2\f$
     *  @see Ostap::Math::beta_inc
     */
    inline double I
    ( const double alpha1 , 
      const double alpha2 , 
      const double z      ) 
    { return beta_inc ( alpha1 , alpha2 , z ) ; }
    // ========================================================================


    // ========================================================================
    /** Pochhammer symbol, aka "rising factorial"
     *  \f[ P(a,n) = a ( a + 1) ( a + 1 ) ... ( a + n - 1 ) = \Pi^{k-1}_{k=0} (a + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::rising_factorial
     */
    double pochhammer
    ( const double         a , 
      const unsigned short n ) ;    
    // ========================================================================
    /** Pochhammer symbol, aka "rising factorial"
     *  \f[ P(a,x) = \frac{ \Gamma ( a + x ) } { \Gamma ( a ) } \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     */
    double pochhammer
    ( const double a , 
      const double x ) ;   
    // ========================================================================== 
    /*  Pochhammer symbol, aka rising factorial 
     *  \f[ P(x,n) = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::rising_factorial
     */
    std::complex<double> pochhammer
    ( const std::complex<double>&  a ,
      const unsigned short         n ) ;
    // =============================================================================
    /*  Pochhammer symbol, aka rising factorial 
     *  \f[ P(x,n) = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::rising_factorial
     */
    std::complex<double> pochhammer
    ( const std::complex<double>&  a ,
      const double                 x ) ;  
    // =============================================================================
    /*  Pochhammer symbol, aka rising factorial 
     *  \f[ P(x,n) = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::rising_factorial
     */
    std::complex<double> pochhammer
    ( const std::complex<double>&  a ,
      const std::complex<double>&  x ) ; 
    // =========================================================================

    // ========================================================================
    /** Rising  factorial, aka "Pochhammer's symbol"
     *  \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::pochhammer 
     *  @see Ostap::Math::falling_factorial
     */
    double rising_factorial
    ( const double         x ,
      const unsigned short n ) ;
    // =========================================================================
    /** Falling factorial
     *  \f[ (x)_n = \Pi^{k-1}_{k=0}  (x - k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::rising_factorial
     */
    double falling_factorial
    ( const double         x ,
      const unsigned short n ) ;
    // =========================================================================
    /** Pochhammer symbol, aka "rising factorial" and its derivative 
     *  \f[ P(x,n) = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::rising_factorial
     *  @see Ostap::Math::pochhammer 
     */
    std::pair<double,double>
    pochhammer_with_derivative 
    ( const double         x , 
      const unsigned short n ) ;    
    // ========================================================================
   
    



    // ========================================================================
   } //                                        The end of namespace Ostap::Math
  // ==========================================================================
} // The end of namesapce Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_GAMMA_H 
// ============================================================================