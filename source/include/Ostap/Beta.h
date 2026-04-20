// ============================================================================
#ifndef OSTAP_BETA_H 
#define OSTAP_BETA_H 1
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
    // Beta & friends 
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
    std::complex<double> beta
    ( const std::complex<double>& x ,
      const std::complex<double>& y ) ;
    // ========================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    std::complex<long double> beta
    ( const std::complex<long double>& x ,
      const std::complex<long double>& y ) ;
    // ========================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    std::complex<double> beta
    ( const std::complex<double>& x ,
      const double                y ) ;
    // =========================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    inline std::complex<double> beta
    ( const double                x , 
      const std::complex<double>& y ) { return beta ( y , x ) ; }
    // =========================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    std::complex<long double> beta
    ( const std::complex<long double>& x ,
      const long double                y ) ;
    // =========================================================================
    /** beta function for 
     *  \f$ \Beta(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)} \f$ 
     *  - \f$ 0<x\f$
     *  - \f$ 0<y\f$ 
     *  @return value of beta function 
     */
    inline std::complex<long double> beta
    ( const long double                x , 
      const std::complex<long double>& y ) { return beta ( y , x ) ; }
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
     */
    double beta_inc 
    ( const unsigned short alpha1 , 
      const unsigned short alpha2 , 
      const double         z      ) ;
    // ========================================================================
    /** Derivative of the normalized incomplete Beta function
     *  - PDF for beta-distribition
     *
     *  \f$ f ( \alpha_1,\alpha_2, z ) = 
     *      I_z( \alpha_1, \alpha_2 ) = 
     *      \frac{\Beta_z(\alpha_1,\alpha_2}}
     *           {\Beta  (\alpha_1,\alpha_2}
     *
     *  - \f$ 0<z<1\f$
     *  - \f$ 0<\alpha_1\f$
     *  - \f$ 0<\alpha_2\f$ 
     */
    double dbeta_inc 
    ( const double alpha1 , 
      const double alpha2 , 
      const double z      ) ;    
    // ========================================================================
    /** Derivative of the normalized incomplete Beta function
     *  - PDF for beta-distribition
     *
     *  \f$ f ( \alpha_1,\alpha_2, z ) = 
     *      I_z( \alpha_1, \alpha_2 ) = 
     *      \frac{\Beta_z(\alpha_1,\alpha_2}}
     *           {\Beta  (\alpha_1,\alpha_2}
     *
     *  - \f$ 0<z<1\f$
     *  - \f$ 0<\alpha_1\f$
     *  - \f$ 0<\alpha_2\f$ 
     */
    double dbeta_inc 
    ( const unsigned short alpha1 , 
      const unsigned short alpha2 , 
      const double         z      ) ;    
    // ========================================================================

    // ========================================================================
    /** Get PDF for beta distribution 
     *  \f[ f(x,\alpha, \beta ) =  
     *   \frac{ x^(\alpha-1) (1-x)^{\beta-1}} { B(\alpha,\beta} \f] 
     *  - \f$ 0 < x < 1 \f$ 
     *  - \f$ 0 < alpha \f$ 
     *  - \f$ 0 < beta   \f$
     *  @see Ostap::Math::dbeta_inc 
     *  @attention Note the order of arguments!
     */
    double beta_pdf
    ( const double x     ,
      const double alpha ,
      const double beta  ) ;
    // ========================================================================
    /** Get CDF for beta distribution 
     *  \f[ F(x,\alpha, \beta ) = I_x(\alpha,\beta)\f] 
     *  - \f$ 0 < x < 1 \f$ 
     *  - \f$ 0 < alpha \f$ 
     *  - \f$ 0 < beta   \f$
     *  @see Ostap::Math::beta_inc 
     *  @attention Note the order of arguments!     
     */
    double beta_cdf
    ( const double x     ,
      const double alpha ,
      const double beta  ) ;
    // ========================================================================
    /**  Quantile function CDF for beta distribution 
     *  - \f$ 0 \le  p \le 1 \f$ 
     *  - \f$ 0 < alpha \f$ 
     *  - \f$ 0 < beta   \f$ 
     */
    double beta_quantile 
    ( const double p     ,
      const double alpha ,
      const double beta  ) ;
    // ========================================================================

    
    // ========================================================================       
   } //                                        The end of namespace Ostap::Math
  // ==========================================================================
} // The end of namesapce Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_BETA_H 
// ============================================================================
