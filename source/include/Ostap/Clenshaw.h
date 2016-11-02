// $Id$ 
// ============================================================================
#ifndef OSTAP_CLENSHAW_H 
#define OSTAP_CLENSHAW_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <utility>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** Namespace with collection of Clenshaw summation algorithms 
     *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2015-08-06
     */
    namespace Clenshaw
    {
      // ======================================================================
      /** Clenshaw algorithm for summation of monomial series (aka "Horner's rule")
       *  \f$  f_1(x) = \sum_{i=0}^{n} a_i x^i     \f$ and 
       *  \f$  f_2(x) = \sum_{i=0}^{n} a_i x^{n-1} \f$
       *
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @see https://en.wikipedia.org/wiki/Horner%27s_method
       *
       *  @code
       *  VECTOR_LIKE_TYPE a { { 1 , 2 , 3 } } ;
       *  double x     = 0 ;
       *  
       *  // f1 = a[0]*x*x +a[1]*x + a[2] 
       *  double f1 = monomial_sum ( a.begin()  , a.end()  , x ).first ;
       *  
       *  // f2 = a[0]+a[1]*x+a[2]*x*x 
       *  double f2 = monomial_sum ( a.rbegin() , a.rend() , x ).first ;
       *
       *  @endcode 
       *
       *  With this algorithm one also gets a first derivative:
       *
       *  @code
       *  VECTOR_LIKE_TYPE a { { 1 , 2 , 3 } } ;
       *  double x     = 0 ;
       *  
       *  // f1 = a[0]*x*x +a[1]*x + a[2] 
       *  std::pair<double,double> r  = monomial_sum ( a.begin()  , a.end()  , x ) ;
       *  double f1 = r.first  ; // polynomial
       *  double d1 = r.second ; // derivative 
       *
       *  @endcode        
       *  
       *  @param first iterator for start of the sequence 
       *  @param last  iterator for end   of the sequence 
       *  @param x    x-value 
       *  @return value of polynom and the first derivative at point x 
       *
       *  @see Ostap::Math::clenshaw_polynom
       *  @see Ostap::Math::horner_a0
       *  @see Ostap::Math::horner_aN
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      std::pair<long double,long double>
      monomial_sum 
      ( ITERATOR          first  , 
        ITERATOR          last   , 
        const long double x      ) 
      {
        if ( first == last ) { return std::make_pair(0,0) ; }
        //
        long double p = *first ;
        long double q = 0      ;
        while ( ++first != last ) 
        {
          q = std::fma ( x , q ,  p     ) ; // x * q + p       ;
          p = std::fma ( x , p , *first ) ; // x * p + *first  ;
        }
        //
        return std::make_pair ( p , q ) ;
      }
      // ======================================================================
      /** Clenshaw algorithm for summation of Legendre series 
       *  \f$ f(x) = \sum_{i=0}^{n} a_i L_i(x) \f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @see Ostap::Math::LegendreSum 
       *  @see Ostap::Math::clenshaw_legendre
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double 
      legendre_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        // 
        if ( first == last ) { return 0 ; }
        //
        long double b2 = 0 ;
        long double b1 = 0 ;
        long double b0 = 0 ;
        while ( first != last ) 
        {
          --last ;
          const int j          =  last - first ;
          const long double pj = *last ;
          //
          b2 = b1 ;
          b1 = b0 ;
          b0 = pj + ( 2 * j + 1 ) * x * b1 / ( j + 1 ) - ( j + 1 ) *  b2 / ( j + 2 ) ;
        }
        //
        return b0 ;
      }
      // ======================================================================
      /** Clenshaw algorithm for summation of Chebyshev series 
       *  \f$ f(x) = \sum_{i=0}^{n} a_i C_i(x) \f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @see Ostap::Math::ChebyshevSum 
       *  @see Ostap::Math::clenshaw_chebyshev
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      chebyshev_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2 = 0 ;
        long double b1 = 0 ;
        long double b0 = 0 ;
        while ( first != last ) 
        {
          --last  ;
          b2 = b1 ;
          b1 = b0 ;
          // b0 = (*ia) + 2 * x * b1 - b2 ;
          b0 = std::fma ( 2 * x , b1 , (*last) - b2 ) ;
        }
        //
        b0 += *first ;
        //
        return 0.5 * ( b0 - b2) ;
      }
      // ======================================================================
      /** Clenshaw algorithm for summation of Hermite's series 
       *  \f$ f(x) = \sum_{i=0}^{n} a_i He_i(x) \f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @see Ostap::Math::clenshaw_hermite
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      hermite_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2 = 0 ;
        long double b1 = 0 ;
        long double b0 = 0 ;
        //
        unsigned long k = std::distance ( first , last ) ;
        while ( first != last ) 
        {
          --last  ;
          b2 = b1 ;
          b1 = b0 ;
          b0 = std::fma ( x , b1 , (*last) - k * b2 ) ;
          --k     ;
        }
        //
        return std::fma ( x , b1 , (*first) - b2  ) ;
      }
      // ======================================================================      
      // Trigonometric sums 
      // ======================================================================      
      /** Clenshaw algorithm for summation of cosine-series 
       *  \f$ f(x) = \frac{a_0}{2} + \sum_{i=k}^{n} a_k \cos( k x) \f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      cosine_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2 = 0 ;
        long double b1 = 0 ;
        long double b0 = 0 ;
        //
        const long double cosx = std::cos ( x ) ;
        while ( first != last ) 
        {
          --last  ;
          b2 = b1 ;
          b1 = b0 ;
          b0 = std::fma ( 2 * cosx  , b1 , *last - b2 ) ;
        }
        //
        return std::fma ( cosx ,  b1 , 0.5L * (*first) - b2 ) ;
      }
      // ======================================================================      
      /** Clenshaw algorithm for summation of sine-series 
       *  \f$ f(x) = \sum_{i=k}^{n} a_k \sin( k x) \f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      sine_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2   = 0 ;
        long double b1   = 0 ;
        long double b0   = 0 ;
        //
        long double sinx = std::sin ( x ) ;
        long double cosx = std::cos ( x ) ;
        //
        while ( 1 < 2 )  
        {
          b2 = b1 ;
          b1 = b0 ;
          if ( first == last ) { break ; }  // BREAK 
          --last  ;   // advace 
          b0 = std::fma ( 2 * cosx  , b1 , *last - b2 ) ;
        }
        //
        return b1 * sinx ;
      }
      // ======================================================================      
      /** Clenshaw algorithm for summation of Fourier-series 
       *  \f$ f(x) = \frac{a_0}{2} + \sum_{i=k}^{n} a_{2k-1}\sin(kx)+a_{2k}\cos(kx) \f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      fourier_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2s  = 0 ;
        long double b1s  = 0 ;
        long double b0s  = 0 ;
        long double b2c  = 0 ;
        long double b1c  = 0 ;
        long double b0c  = 0 ;
        //
        long double sinx = std::sin ( x ) ;
        long double cosx = std::cos ( x ) ;
        //
        while ( first != last ) 
        {
          //
          // cosine 
          //
          --last    ;   // advance 
          b2c = b1c ;
          b1c = b0c ;
          b0c = std::fma ( 2 * cosx   , b1c , *last - b2c ) ;
          //
          // sine 
          //
          b2s = b1s ;
          b1s = b0s ;
          // 
          if ( last == first )  { break ; }
          //
          --last    ;   // advance 
          b0s = std::fma ( 2 * cosx   , b1s , *last - b2s ) ;
          //
        }
        //
        return std::fma ( cosx ,  b1c , 0.5 * (*first) - b2c + b1s * sinx ) ;
      }
      // ======================================================================      
      /** Clenshaw algorithm for Fejer sums for cosine-series 
       *  For the series of partial sums 
       *  \f$ f_n(x) = \frac{a_0}{2} + \sum_{i=k}^{n} a_k \cos( k x) \f$
       *  Fejer sums are defiend as 
       *  \f$ F_n(s) \equiv \frac{1}{N+1}\sum_{i=0}^{N} f_i(x)\f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      fejer_cosine_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2 = 0 ;
        long double b1 = 0 ;
        long double b0 = 0 ;
        //
        const long double cosx = std::cos ( x ) ;
        //
        const unsigned long N = std::distance ( first , last ) ;
        const long double   d = 1.0L/N ;
        unsigned long k = 0 ;
        while ( first  != last ) 
        {
          ++k     ;
          --last  ;
          b2 = b1 ;
          b1 = b0 ;
          b0 = std::fma ( 2 * cosx  , b1 , ( *last ) * k * d  - b2 ) ;
        }
        //
        return std::fma ( cosx ,  b1 , 0.5L * (*first) - b2 ) ;
      }
      // ======================================================================      
      /** Clenshaw algorithm for Fejer sums for sine-series 
       *  For the series of partial sums 
       *  \f$ f_n(x) = \frac{a_0}{2} + \sum_{i=k}^{n} a_k \sin( k x) \f$
       *  Fejer sums are defiend as 
       *  \f$ F_n(s) \equiv \frac{1}{N+1}\sum_{i=0}^{N} f_i(x)\f$
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      fejer_sine_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2   = 0 ;
        long double b1   = 0 ;
        long double b0   = 0 ;
        //
        long double sinx = std::sin ( x ) ;
        long double cosx = std::cos ( x ) ;
        //
        const unsigned long N = std::distance ( first , last ) ;
        const long double   d = 1.0L/N ;
        unsigned long       k = 0 ;        
        //
        while ( first != last ) 
        {
          ++k     ;
          --last  ;
          b2 = b1 ;
          b1 = b0 ;
          b0 = std::fma ( 2 * cosx  , b1 , ( *last ) * k * d - b2 ) ;
        }
        //
        return b1 * sinx ;
      }
      // ======================================================================      
      /** Clenshaw algorithm for Fejer sums for Fourier-series 
       *  \f$ f_n(x) = \frac{a_0}{2} + \sum_{i=k}^{n} a_{2k-1}\sin(kx)+a_{2k}\cos(kx) \f$
       *  \f$ F_n(x) = \frac{1}{n}\sum_{k=0}{n} f_n(x)\f$ 
       *  @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2015-02-10
       */
      template <class ITERATOR>
      inline 
      long double
      fejer_sum 
      ( ITERATOR           first , 
        ITERATOR           last  ,
        const long double  x     ) 
      {
        if ( first == last ) { return 0 ; }
        //
        long double b2s  = 0 ;
        long double b1s  = 0 ;
        long double b0s  = 0 ;
        long double b2c  = 0 ;
        long double b1c  = 0 ;
        long double b0c  = 0 ;
        //
        long double sinx = std::sin ( x ) ;
        long double cosx = std::cos ( x ) ;
        //
        const unsigned long N = std::distance ( first , last ) ;
        const long double   d = 1.0L/( N + 1 ) ;
        unsigned long   k     = 0 ;
        while ( first != last ) 
        {
          ++k  ;
          //
          // cosine 
          //
          --last    ;   // advance 
          b2c = b1c ;
          b1c = b0c ;
          b0c = std::fma ( 2 * cosx , b1c , ( *last ) * 2 * k * d  - b2c ) ;
          //
          // sine 
          //
          b2s = b1s ;
          b1s = b0s ;
          // 
          if ( last == first )  { break ; }
          //
          --last    ;   // advance 
          b0s = std::fma ( 2 * cosx , b1s ,  ( *last ) * 2 * k * d  - b2s ) ;
          //
        }
        //
        return std::fma ( cosx ,  b1c , 0.5 * (*first) - b2c + b1s * sinx ) ;
      }
      // ======================================================================
    } //                             The end of namespace Ostap::Math::Clenshaw     
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Gaudi
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CLENSHAW_H
// ============================================================================
