// ============================================================================
#ifndef OSTAP_ZETA_H 
#define OSTAP_ZETA_H 1
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
    // zeta, eta &friends 
    // ========================================================================
    /** Riemann's Zeta function \f$ n \ne 1\f$:
     *  \f$ \zeta ( n ) = \sum_k k^{-n}\f$ 
     */
    double zeta ( const int    n ) ;
    // ========================================================================
    /** Riemann's Zeta function \f$ s \ne 1\f$:
     *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
     */
    double zeta ( const double s ) ;
    // ========================================================================
    /** Riemann's Zeta function \f$ s \ne 1\f$:
     *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
     */
    long double zeta ( const long double s ) ;
    // =======================================================================
    /** Riemann's Zeta function \f$ s \ne 1\f$:
     *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
     */
    std::complex<double> zeta ( const std::complex<double>& z  ) ;
    // =======================================================================
    /** Riemann's Zeta function \f$ s \ne 1\f$:
     *  \f$ \zeta ( s ) = \sum_k k^{-s}\f$ 
     */
    std::complex<long double> zeta ( const std::complex<long double>& z  ) ;

    // ========================================================================
    /** helper chi-function:
     * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{s} \Gamma (1-s), \f]
     * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s ) \f$  
     */
    double chi ( const double s ) ;
    // ========================================================================
    /** helper chi-function:
     * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{s} \Gamma (1-s), \f]
     * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s ) \f$  
     */
    long double chi ( const long double s ) ;
    // ========================================================================
   /** helper chi-function:
     * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{s} \Gamma (1-s), \f]
     * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s ) \f$  
     */
    std::complex<double> chi ( const std::complex<double>&  s ) ;
    // ========================================================================
    /** helper chi-function:
     * \f[ \chi(s) = \frac{(2\pi)^s}{\pi} \sin \frac{\pi s}{s} \Gamma (1-s), \f]
     * such as \f$  \zeta(s) = \chi(s) \zeta ( 1 - s ) \f$  
     */
    std::complex<long double> chi ( const std::complex<long double>&  s ) ;
    // ========================================================================


    // ========================================================================
    /** Riemann's Zeta function minus 1 \f$ n\ne 1\f$:
     *  \f$ f ( n ) = \zeta ( n ) - 1 = -1 + \sum_k k^{-n}\f$ 
     */
    double zetam1 ( const int    n ) ;
    // ========================================================================
    /** Riemann's Zeta function minus 1 \f$ s \ne 1\f$:
     *  \f$ f ( s ) = \zeta ( s ) - 1 = -1 + \sum_k k^{-s}\f$ 
     */
    double zetam1 ( const double s ) ;
    // ========================================================================
    /** Riemann's Zeta function minus 1 \f$ s \ne 1\f$:
     *  \f$ f ( s ) = \zeta ( s ) - 1 = -1 + \sum_k k^{-s}\f$ 
     */
    long double zetam1 ( const long double s ) ;
    // ========================================================================

    // ========================================================================
    /** Dirichlet's Eta function 
     *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
     */
    double eta ( const int    n ) ;
    // ========================================================================
    /** Dirichlet's Eta function 
     *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
     */
    double eta ( const double s ) ;
    // ========================================================================
    /** Dirichlet's Eta function 
     *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
     */
    long double eta ( const long double s ) ;
    // ========================================================================

    // ========================================================================
    /** Dirichlet's Eta function 
     *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
     */
    inline 
    double dirichlet_eta ( const int    n ) { return eta ( n ) ; }
    // ========================================================================
    /** Dirichlet's Eta function 
     *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
     */
    inline 
    double dirichlet_eta ( const double s ) { return eta ( s ) ; }
    // ========================================================================
    /** Dirichlet's Eta function 
     *  \f$ \eta ( z ) = ( 1 - 2 ^{1-s} ) \zeta ( s ) 
     */
    inline 
    long double dirichlet_eta ( const long double s ) { return eta ( s ) ; }
    // ========================================================================

    // ========================================================================
    /** Dirichlet's beta function
     *  \f[ \beta ( s ) 
     *   \equiv \sum_{n=0}^{+\infty}\frac{(-1)^n}{(2n+1)^s} \]
     *    \f[ \beta ( s ) = \frac{1}{\Gamma(s)}
     * \int_0^{+\infty}\frac{x^{s-1}\mathrm{e}^{-x}}{1+\mathrm{e}^{-2x}}dx \f]
     *   \f[ \beta ( s ) = 4^{-1}\left( \zeta ( z , \frac{1}{4} - \zeta ( s , frac{3}{4} ) \right) \f] 
     *  where \f$\zeta(a,b)\f$ is Hurwitz's zeta function
     */
    double dirichlet_beta
    ( const int x ) ;
    // ========================================================================
    /** Dirichlet's beta function
     *  \f[ \beta ( s ) 
     *   \equiv \sum_{n=0}^{+\innfty}\frac{(-1)^n}{(2n+1)^s} \]
     *    \f[ \beta ( s ) = \frac{1}{\Gamma(s)}
     * \int_0^{+\infty}\frac{x^{s-1}\mathrm{e}^{-x}}{1+\mathrm{e}^{-2x}}dx \f]
     *   \f[ \beta ( s ) = 4^{-1}\left( \zeta ( z , \frac{1}{4} - \zeta ( s , frac{3}{4} ) \right) \f] 
     *  where \f$\zeta(a,b)\f$ is Hurwitz's zeta function
     */
    double dirichlet_beta
    ( const double x ) ;
    // ========================================================================
    /** Dirichlet's beta function
     *  \f[ \beta ( s ) 
     *   \equiv \sum_{n=0}^{+\innfty}\frac{(-1)^n}{(2n+1)^s} \]
     *    \f[ \beta ( s ) = \frac{1}{\Gamma(s)}
     * \int_0^{+\infty}\frac{x^{s-1}\mathrm{e}^{-x}}{1+\mathrm{e}^{-2x}}dx \f]
     *   \f[ \beta ( s ) = 4^{-1}\left( \zeta ( z , \frac{1}{4} - \zeta ( s , frac{3}{4} ) \right) \f] 
     *  where \f$\zeta(a,b)\f$ is Hurwitz's zeta function
     */
    long  double dirichlet_beta
    ( const long double x ) ;
    // ========================================================================

    // ========================================================================
    /** Hurwitz Zeta function 
     *  \f$ zeta ( s , q ) = \sum_k  ( k + q )^{-s}\f$
     *  - \f$ 1 < s \f$
     *  - \f$ 0 < q \f$
     */
    double hurwitz
    ( const double s     ,
      const double q = 1 ) ;
    // ========================================================================
    /** Hurwitz Zeta function 
     *  \f$ zeta ( s , q ) = \sum_k  ( k + q )^{-s}\f$
     *  - \f$ 1 < s \f$
     *  - \f$ 0 < q \f$
     */
    long double hurwitz
    ( const long double s     ,
      const long double q = 1 ) ;
    // ========================================================================
    /** Hurwitz Zeta function 
     *  \f$ zeta ( s , q ) = \sum_k  ( k + q )^{-s}\f$
     *  - \f$ 1 < s \f$
     *  - \f$ 0 < q \f$
     */
    inline  
    double hurwitz_zeta
    ( const double s     ,
      const double q = 1 ) { return hurwitz ( s , q ) ; }
    // ========================================================================

    // ========================================================================
   } //                                        The end of namespace Ostap::Math
  // ==========================================================================
} // The end of namesapce Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ZETA_H 
// ============================================================================