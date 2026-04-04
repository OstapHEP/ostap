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
    /** polylogarithm function  \f$ Li_n(x)  = \sum \frac{x^k}{k^n} f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @param n parameter
     *  @param x argument
     *  @return value of polylogarithm function \f$ Li_n(x) \f$
     *  @attention for \$ 1 < x \f$ the real part is returned, for imaginary part see      
     *  @see  Ostap::Math::ImLi
     */
    double Li
    ( const short  n ,
      const double x ) ;
    // ========================================================================
    /** polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
     *  @see https://en.wikipedia.org/wiki/Polylogarithm
     *  @see David Wood, "The computation of polylogarithms",
     *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
     *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
     *  @param s parameter
     *  @param x argument
     *  @return value of polylogarithm function \f$ Li_n(x) \f$
     *  @attention for \$ 1 < x \f$ the real part is returned, for imaginary part see      
     *  @see  Ostap::Math::ImLi
     */
    double Li
    ( const double s ,
      const double x ) ;
    // =======================================================================
  } //                                         The end of namesace Ostap::Math
  // =========================================================================
} //                                                The end of namespace Ostap
// ===========================================================================
#endif // OSTAP_POLYLOG_H 
// ===========================================================================
//                                                                     The END
// =========================================================================== 
    