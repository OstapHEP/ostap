// ============================================================================
#ifndef CHOOSE_UTILS_H 
#define CHOOSE_UTILS_H 1
// ============================================================================
// Include files
// ============================================================================
//  STD&STL
// ============================================================================
#include <climits>
#include <cmath>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    namespace Utils 
    {
      // ======================================================================
      /** calculate the binomial coefficient C(k,n) = n!/((n-k)!*k!)
       *  In case of overflow std::numeric_limits<unsigned long long>::max is returned 
       */
      inline unsigned long long choose 
      ( unsigned short n , 
        unsigned short k ) 
      {
        //
        // ====================================================================
        typedef  std::numeric_limits<unsigned long long> ULLTYPE ;
        static_assert ( ULLTYPE::is_specialized, 
                        "numeric_limits<unsigned long long> is not specialzaed!" ) ;
        static const unsigned long long s_ullmax = ULLTYPE::max () ;
        // ===================================================================
        //
        if      ( k > n                ) { return 0 ; }
        else if ( 0 == k || n == k     ) { return 1 ; }
        else if ( 1 == k || n == k + 1 ) { return n ; }
        else if ( 2 == k || n == k + 2 ) { return 1ull * n * ( n - 1 ) / 2 ; }
        //
        k = std::min ( k , (unsigned short) ( n - k ) ) ;
        unsigned long long r = 1  ;
        for ( unsigned short d = 1 ; d <= k  ; ++d , --n ) 
        {
          if ( r > s_ullmax / n * d ) { return s_ullmax ; }  //  RETURN
          r = ( r / d ) * n + ( r % d ) * n / d;
        }
        return r ;
      }
      // ======================================================================
      /// evaluate the binomila  coefficient as long double C(k,n) = n!/((n-k)!*k!)
      inline long double choose_long_double 
      ( const unsigned short n ,
        const unsigned short k ) 
      {
        //
        // ====================================================================
        typedef  std::numeric_limits<unsigned long long> ULLTYPE ;
        static_assert ( ULLTYPE::is_specialized, 
                        "numeric_limits<unsigned long long> is not specialzaed!" ) ;
        static const unsigned long long s_ullmax = ULLTYPE::max ()              ;
        static const long double        s_emax   = std::log ( 0.2L * s_ullmax ) ;
        static const unsigned short     s_digits = ULLTYPE::digits - 2          ;
        // ===================================================================
        //
        if      ( k > n            ) { return 0 ; }
        else if ( 0 == k || n == k ) { return 1 ; }
        else if ( n < s_digits     ) { return choose ( n , k ) ; }
        else if ( n <= 67          ) { return choose ( n , k ) ; }
        //
        const unsigned  short k1 =  2 * k < n ? k : n - k ;
        if ( k1 * std::log2 ( M_E * n / k1 ) < 63 ) { return choose ( n , k ) ; }
        //
        //
        long double a = std::lgamma ( ( long double )     n + 1 )  ;
        if ( a < s_emax ) { return choose ( n , k ) ; }
        a            -= std::lgamma ( ( long double ) n - k + 1 ) ;
        if ( a < s_emax ) { return choose ( n , k ) ; }
        a            -= std::lgamma ( ( long double )     k + 1 ) ;
        if ( a < s_emax ) { return choose ( n , k ) ; }
        //
        return std::exp ( a ) ;
      }
      // =====================================================================
      /** calculate the inverse binomial coefficient 
       */
      inline long double ichoose 
      ( const unsigned short n , 
        const unsigned short k ) 
      {
        //
        if      ( k > n                ) { return 0        ; }
        else if ( 0 == k || n == k     ) { return 1        ; }
        else if ( 1 == k || n == k + 1 ) { return 1.0L / n ; }
        else if ( 2 == k || n == k + 2 ) { return 2.0L / ( 1ull * n * ( n - 1 ) ) ; }
        //
        const unsigned short kk = std::min ( k , (unsigned short) ( n - k ) ) ;
        //
        long double r = 1  ;
        for ( unsigned short i = 1 ; i <= kk  ; ++i ) 
        {
          r *= i ;
          r /= ( n + 1 - i ) ;  
        }
        return r ;
      }
      // ======================================================================
      /*  calculate the logarithm of binomial coefficient
       *  \f$ \log C^n_k \f$
       */
      // ======================================================================
      inline long double log_choose
      ( const unsigned short n ,
        const unsigned short k ) 
      {
        if      ( k <= 1 || k >= n ) { return 0 ; } //
        else if ( n <= 67 )          { return std::log( (long double) choose ( n , k ) ) ; }
        //
        const unsigned  short k1 = 2*k < n ? k : n - k ;
        if ( k1 * std::log2 ( M_E * n / k1 ) < 63 ) 
        { return std::log( (long double) choose ( n , k ) ) ; }
        // 
        return 
          std::lgamma ( (long double) ( n     + 1 ) ) - 
          std::lgamma ( (long double) ( k     + 1 ) ) - 
          std::lgamma ( (long double) ( n - k + 1 ) ) ;
      }
      // ======================================================================
    } //                                The end of namespace Ostap::Math::Utils 
    // ========================================================================
  } //                                         The end of namepsace Ostap::Math 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // CHOOSE_UTILS_H
// ============================================================================
