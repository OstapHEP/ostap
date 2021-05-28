#ifndef BERNSTEIN_UTILS_H 
#define BERNSTEIN_UTILS_H 1
// ============================================================================
// Inclde files 
// ============================================================================
#include <iterator>
#include <algorithm>
#include <array>
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
        if      ( 0 == M ) { return std::copy ( b_begin , b_end , output ) ; }
        else if ( 0 == N ) { return std::copy ( a_begin , a_end , output ) ; }
        else if ( 1 == M ) 
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
      /// multiply two Bernstein polynomials 
      template <class OUTPUT   , 
                class ITERATOR , 
                class TYPE     , 
                std::size_t   N>
      OUTPUT b_multiply 
      ( ITERATOR                  a_begin , 
        ITERATOR                  a_end   , 
        const std::array<TYPE,N>& b       , 
        OUTPUT                    output  )
      { return b_multiply ( a_begin , a_end , b.begin() , b.end () , output ) ; }
      // =========================================================================
      /// multiply two Bernstein polynomials 
      template <class OUTPUT   , 
                class TYPE     , 
                std::size_t   N,
                std::size_t   K>
      OUTPUT b_multiply 
      ( const std::array<TYPE,N>& a       , 
        const std::array<TYPE,K>& b       , 
        OUTPUT                    output  )
      { return b_multiply ( a.begin () , a.end () , 
                            b.begin () , b.end () , output ) ; }
      // ====================================================================== 
      /** Create Bernstein coefficients for the linear polynomial 
       *  \f$ p(x) = x-x_0 = \alpha_0 (1-x) + \alpha_1 x \f$  
       *  @param x0 root of linear polynomial \f$ p(x) = x - x0 \f$
       *  @return array of coefficients \f$ \alpha_0,  \alpha_1 \f$
       */
      inline std::array<double,2> 
      bernstein1_from_roots 
      ( const long double x0 ) 
      { return {{ double ( -x0 ) , double ( 1 - x0 ) }} ; }
      // ========================================================================
      /** Create Bernstein coefficients for the linear polynomial 
       *  \f$ p(x) = x-x_0 = \alpha_0 (1-x) + \alpha_1 x \f$  
       *  @param x0 root of linear polynomial \f$ p(x) = x - x0 \f$
       */
      template <class TYPE>
      void bernstein1_from_roots 
      ( const long double   x0 , 
        std::array<TYPE,2>& b  )
      {
        b [ 0 ] = x0 < 0.5 ?   - x0 : x0     ; 
        b [ 1 ] = x0 < 0.5 ? 1 - x0 : x0 - 1 ;
      }
      // ======================================================================
      /** Create Bernstein coefficients for the quadratic polynomial with two real roots 
       *  \f$ p(x) = (x-x_0)(x-x_1) = \alpha_0 (1-x)^2 + \alpha_1 2x(1-x) + \alpha_2 x^2\f$
       *  @param x0 root of quadratic polynomial \f$ p(x) = (x - x_0)(x-x_1) \f$
       *  @param x1 root of quadratic polynomial \f$ p(x) = (x - x_0)(x-x_1) \f$
       *  @return array of coefficients \f$ \alpha_0,  \alpha_1 , \alpha_2 \f$
       */ 
      inline std::array<double,3> 
      bernstein2_from_roots 
      ( const long double x0 ,
        const long double x1 ) 
      { 
        const long double s = x0 + x1 ;
        const long double p = x0 * x1 ;
        return {{ double ( p ) , double ( p - 0.5 * s ) , double ( 1 + p - s ) }}; 
      }
      // ======================================================================
      /** Create Bernstein coefficients for the quadratic polynomial with two real roots 
       *  \f$ p(x) = (x-x_0)(x-x_1) = \alpha_0 (1-x)^2 + \alpha_1 2x(1-x) + \alpha_2 x^2\f$
       *  @param x0 root of quadratic polynomial \f$ p(x) = (x - x_0)(x-x_1) \f$
       *  @param x1 root of quadratic polynomial \f$ p(x) = (x - x_0)(x-x_1) \f$
       *  @return array of coefficients \f$ \alpha_0,  \alpha_1 , \alpha_2 \f$
       */ 
      template <class TYPE>
      void bernstein2_from_roots 
      ( const long double     x0 ,
        const long double     x1 , 
        std::array<TYPE,3>&   b  )
      { 
        const long double s = x0 + x1 ;
        const long double p = x0 * x1 ;
        //
        b [ 0 ] = p           ;
        b [ 1 ] = p - 0.5 * s ;
        b [ 2 ] = 1 + p - s   ;
      }
      // ======================================================================
      /** Create Bernstein coefficients for the quadratic polynomial with complex root
       *  \f$ p(x) = (x-x_0)(x-x_0^*) = \alpha_0 (1-x)^2 + \alpha_1 2x(1-x) + \alpha_2 x^2\f$
       *  @param x0 complex root of quadratic polynomial \f$ p(x) = (x - x_0)(x-x_0^*) \f$
       *  @return array of coefficients \f$ \alpha_0,  \alpha_1 , \alpha_2 \f$
       */ 
      template <class TYPE>
      void bernstein2_from_roots 
      ( const std::complex<double> x0 ,
        std::array<TYPE,3>&        b  )
      { 
        const long double s = 2 * x0.real ()   ;
        const long double p = std::norm ( x0 ) ;
        //
        b [ 0 ] = p           ;
        b [ 1 ] = p - 0.5 * s ;
        b [ 2 ] = 1 + p - s   ;
      }
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




