// ============================================================================
#ifndef OSTAP_POLYLOG_H 
#define OSTAP_POLYLOG_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <complex>
// ============================================================================
/** @file Ostap/PolyLog.h
 *  Polylogarithm & friends
 */
 // ===========================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {    
    // ========================================================================
    // Polylogaritm & friends    
    // ========================================================================
    /** Imaginary part of polylogarithm function \f$ Im Li_n(x) f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @see Eq (3.1) in Wood
     *  @see Ostap::Math::Li 
     *  @param n parameter
     *  @param x argument
     *  @return Imaginary part of polylogarithm function \f$ Li_n(x) \f$
     */
    double ImLi
    ( const short  n , 
      const double x ) ;
    // ========================================================================
    /** Imaginary part of polylogarithm function \f$ Im Li_n(x) f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @see Eq (3.1) in Wood
     *  @see Ostap::Math::Li 
     *  @param s parameter
     *  @param x argument
     *  @return Imaginary part of polylogarithm function \f$ Li_s(x) \f$
     */
    double ImLi
    ( const double s , 
      const double x ) ;    
    // ========================================================================
    /** Polylogarithm function  \f$ Li_n(x)  = \sum \frac{x^k}{k^n} f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @param n parameter
     *  @param x argument
     *  @return value of polylogarithm function \f$ Li_n(x) \f$
     *  @attention For \$ 1 < x \f$ the real part is returned, for imaginary part see      
     *  @see  Ostap::Math::ImLi
     # 
     *  The algorithm:
     *
     *  Treat the explicit cases:
     *  - \f$ Li_n(0)  = 0 \f$
     *  - \f$ Li_n(1)  = \zeta(n) \f$
     *  - \f$ Li_n(-1) = \eta(n) \f$
     *  - Explicit case   \f$ Li_1(x)\f$
     *  - Explicit case   \f$ Li_0(x)\f$
     *  - Explicit case   \f$ Li_{-1}(x)\f$
     *  - Explicit case   \f$ Li_{-2}(x)\f$
     *  - Rational case   \f$ Li_{n}(x)\f$ for negative \f$ n \f$  larger then N_EULERIAN_MAX 
     *  - For dilogarithm \f$ Li_2(x) \f$ first reduce the arguumetn to \f$ \left| x \right| \le \frac{1}{2} \f$
     *
     * Generic actions:
     *  - Use the power series around \f$ x=0\f$ for small \f$ x \f$
     *  - Use the power series around \f$ x=1\f$ if \f$ \left| \log x \right| \f$ is not too large  
     *  - Transform negative arguments to positive using square formula
     *  - Reduce large argumens using reciprocal formula 
     */
    double Li
    ( const short  n ,
      const double x ) ;
    // ========================================================================
    /** Polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @param s parameter
     *  @param x argument
     *  @return value of polylogarithm function \f$ Li_s(x) \f$
     *  @attention For \$ 1 < x \f$ the real part is returned, for imaginary part see      
     *  @see  Ostap::Math::ImLi
     */
    double Li
    ( const double s ,
      const double x ) ;
    // =======================================================================
    /** Polylogarithm function  \f$ Li_n(x)  = \sum \frac{x^k}{k^n} f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @param s parameter
     *  @param x argument
     */
    std::complex<double> Li
    ( const short                 n ,
      const std::complex<double>& z ) ;
    // =======================================================================
    /** Polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @param s parameter
     *  @param x argument
     *  @return value of polylogarithm function \f$ Li_s(x) \f$
     */
    std::complex<double> Li
    ( const double                s ,
      const std::complex<double>& z ) ;
    // =======================================================================

    // =======================================================================
    // Related functions 
    // =======================================================================

    // =======================================================================
    // Inverse tangent integral
    // =======================================================================
    
    // =======================================================================
    /** (Real part of) inverse tangent integral
     *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
     *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
     */ 
    double Ti
    ( const short  n ,
      const double x ) ;
    // =======================================================================
    /** (Real part of) inverse tangent integral
     *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
     *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
     */ 
    double Ti
    ( const double s ,
      const double x ) ;
    // ======================================================================
    /** Inverse tangent integral
     *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
     *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
     */ 
    std::complex<double> Ti
    ( const short                 n ,
      const std::complex<double>& z ) ;
    // =======================================================================
    /** Inverse tangent integral
     *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
     *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
     */ 
    std::complex<double> Ti
    ( const double                s ,
      const std::complex<double>& z ) ;
    // =======================================================================

    // =======================================================================
    // Legendre chi-function
    // =======================================================================
    
    // =======================================================================
    /** (Real part of) Legendre chi-function
     *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
     *  https://en.wikipedia.org/wiki/Legendre_chi_function
     */
    double legendre_chi
    ( const short  n ,
      const double x ) ;
    // =======================================================================
    /** (Real part of) Legendre chi-function
     *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
     *  https://en.wikipedia.org/wiki/Legendre_chi_function
     */
    double legendre_chi
    ( const double s ,
      const double x ) ;
    // =======================================================================
    /** Legendre chi-function
     *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
     *  https://en.wikipedia.org/wiki/Legendre_chi_function
     */
    std::complex<double> legendre_chi
    ( const short                 n ,
      const std::complex<double>& z ) ;
    // =======================================================================
    /** Legendre chi-function
     *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
     *  https://en.wikipedia.org/wiki/Legendre_chi_function
     */
    std::complex<double> legendre_chi
    ( const double                s ,
      const std::complex<double>& z ) ;
    // =======================================================================
    
    // =======================================================================    
  } //                                         The end of namesace Ostap::Math
  // =========================================================================
} //                                                The end of namespace Ostap
// ===========================================================================
#endif // OSTAP_POLYLOG_H 
// ===========================================================================
//                                                                     The END
// =========================================================================== 
    
