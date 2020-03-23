#ifndef BERNSTEIN_UTILS_H 
#define BERNSTEIN_UTILS_H 1
// ============================================================================
// Inclde files 
// ============================================================================
#include <iterator>
// ============================================================================
// local
// ============================================================================
#include "choose_utils.h"
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
      /// multiply two Bernstein polynomials 
      template <class OUTPUT , 
                class INPUT1 , 
                class INPUT2 >
      OUTPUT b_multiply ( INPUT1 a_begin , 
                          INPUT1 a_end   , 
                          INPUT2 b_begin , 
                          INPUT2 b_end   ,
                          OUTPUT output  )
      {
        const unsigned int M = 
          std::max ( typename std::iterator_traits<INPUT1>::difference_type ( 0 )  ,
                     std::distance ( a_begin , a_end ) ) ;
        const unsigned int N = 
          std::max ( typename std::iterator_traits<INPUT2>::difference_type ( 0 )  ,
                     std::distance ( b_begin , b_end ) ) ;
        //
        if      ( 1 == M ) 
        {
          const long double a = *a_begin ;
          for ( unsigned short k = 0 ; k < N ; ++k , ++output ) 
          {
            const long double bk =  *( b_begin + k ) ;
            *output = a * bk ;
          }
          return output ;
        }
        else if ( 1 == N ) 
        {
          const long double b = *b_begin ;
          for ( unsigned short k = 0 ; k < M ; ++k , ++output ) 
          {
            const long double ak =  *( a_begin + k ) ;
            *output = ak * b ;
          }
          return output ;      
        }
        //
        const unsigned int  m = M - 1 ;
        const unsigned int  n = N - 1 ;
        //
        const unsigned long K = m + n ;
        for ( unsigned long k = 0 ; k <= K ; ++k , ++output ) 
        {
          const unsigned int jmin = ( n <= k ) ? k - n : 0 ;
          const unsigned int jmax = ( k <= m ) ? k     : m ;
          //
          long double ck = 0 ;
          for ( unsigned int j = jmin ; j <= jmax ; ++j ) 
          {
            const long double a  = * ( a_begin +       j   ) ;
            const long double b  = * ( b_begin + ( k - j ) ) ;
            const long double ab = a * b ;
            if ( 0 != a  && 0 != b && 0 != ab && !s_zero ( ab ) ) 
            {
              ck += ab * 
                Ostap::Math::Utils::choose_long_double ( m     ,     j ) * 
                Ostap::Math::Utils::choose_long_double ( n     , k - j ) *
                Ostap::Math::Utils::ichoose            ( m + n , k     ) ;
            }
          }
          //
          *output = ck ;
        }
        return output ;
      } // 
      // ======================================================================
      /// (recursive) De Casteljau's algorithm
      template <class ITERATOR>
      long double casteljau
      ( ITERATOR          first ,
        ITERATOR          last  ,
        const long double t0    ,
        const long double t1    )
      {
        // the trivial cases
        if      ( first == last    ) { return 0       ; }
        //
        const std::size_t len  = std::distance ( first , last  ) ;
        //
        if      ( 1 == len ) { return       *first                        ; }
        else if ( 2 == len ) { return t1 * (*first) + t0 * ( *(first+1) ) ; }
        //
        ITERATOR second = --last ;
        //
        // prepare recursion
        for ( ITERATOR it = first ; it != second ; ++it )
        { *it = t1 * ( *it )  + t0 * ( *( it+1 ) ) ; }
        //
        // recursion
        return casteljau ( first , second , t0 , t1 ) ;
      }
      // ======================================================================
    } //                                The end of namespace Ostap::Math::Utils 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // BERNSTEIN_UTILS_H 1
// ============================================================================




