// ============================================================================
#ifndef OSTAP_CESARO_H 
#define OSTAP_CESARO_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <numeric>
#include <algorithm>
// ============================================================================
/** @file Ostap/Cesaro.h
 *  Utilities to calualte Cesarp's sums 
 *  @see https://en.wikipedia.org/wiki/Ces%C3%A0ro_summation
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // =========================================================================
    /** Calculate "corrected" coefficeincy for Cesaro's sum
     *  @see https://encyclopediaofmath.org/index.php?title=Ces%C3%A0ro_summation_methods
     *  @param k      (INPUT) summationorder, k=1 correspond sot regular sum 
     *  @param first  (INPUT) start of the sequnce of coefficients 
     *  @param last   (INPUT) end   of the sequnce of coefficients 
     *  @param output (OUTPUT) iterator for the sequnce of updated coefficients 
     *  @author Vanya BELYAEV Ivan.Belyave@cern.ch
     *  @date   2025-04-29 
     */
    template <class INPUT,
             class OUTPUT>
    inline OUTPUT cesaro
    ( const unsigned short K      ,
      INPUT                first  ,
      INPUT                last   ,
      OUTPUT               output )
      {
        if ( first == last || !K ) { return std::copy ( first , last , output ) ; } 
        // 
        const std::size_t N = std::distance ( first , last ) ;
        long double       a = ( K + N + 1 ) * 1.0L / ( N + 1 ) ; 
        for ( std::size_t k = 0 ; first != last ; ++k , ++first , ++output )
        {
          a *= ( N - k + 1 ) * 1.0L / ( K + N - k + 1 ) ;
          *output = a * (*first) ;
        }
        return output ;
      }   
    // =========================================================================
    /** Calculate the Cesaro sum 
     *  @see https://encyclopediaofmath.org/index.php?title=Ces%C3%A0ro_summation_methods
     *  @param k      (INPUT) summationorder, k=1 correspond sot regular sum 
     *  @param first  (INPUT) start of the sequnce of coefficients 
     *  @param last   (INPUT) end   of the sequnce of coefficients 
     *  @param return Cesaro sum of order k 
     *  @author Vanya BELYAEV Ivan.Belyave@cern.ch
     *  @date   2025-04-29 
     */
    template <class ITERATOR>
    inline double cesaro_sum
    ( const unsigned short K     ,
      ITERATOR             first ,
      ITERATOR             last  )
      {
        if ( first == last ) { return 0 ; ; }
        if ( !K            ) { return std::accumulate ( first , last , 0.0L ) ; }
        //
        const std::size_t N = std::distance ( first , last ) ;
        //
        long double a      = ( K + N + 1 ) * 1.0L / ( N + 1 ) ;
        long double result = 0 ;
        for ( std::size_t k = 0 ;  first != last ; ++k , ++first )
        {
          a      *= ( N - k + 1 ) * 1.0L / ( K + N - k + 1 ) ;
          result += a * ( *first ) ;
        }
        return result ;
      } 
    // =========================================================================
  }
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CESARO_H
// ============================================================================
