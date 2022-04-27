// ============================================================================
#ifndef OSTAP_DIFFERENCES_H 
#define OSTAP_DIFFERENCES_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
#include <iterator>
// ============================================================================
/** @file Ostap/Differences.h
 *  Collection of classes and functions to deal with the finite differences 
 */
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @namespace Ostap::Math::Differences
     *  Collection of classes and functions to deal with the finite differences 
     *  @see https://en.wikipedia.org/wiki/Divided_differences 
     *  @author Vanya Belyaev
     *  @date   2018-07-30
     */
    namespace Differences 
    {
      // ======================================================================
      // Functions for divided differences :
      // ======================================================================
      /** Divided Forward differences of order-0 from the function 
       *  @code
       *  auto  fun = []  ( double x ) { return std::sin(x) ; }
       *  double d0 = divided ( fun , 0.5 ) ;
       *  @endcode
       *  @see https://en.wikipedia.org/wiki/Divided_differences 
       */ 
      template <class FUNCTION>
      inline double divided 
      ( FUNCTION     fun , 
        const double x ) { return fun  ( x ) ; }
      // ======================================================================
      /** Divided Forward differences of high order from the function 
       *  @code
       *  auto  fun = []  ( double x ) { return std::sin(x) ; }
       *  double d0 = divided ( fun , 0.5 ) ;
       *  double d1 = divided ( fun , 0.5 , 0.6 ) 
       *  double d2 = divided ( fun , 0.5 , 0.6 , 0.7 )
       *  double d3 = divided ( fun , 0.5 , 0.6 , 0.7 , 0.8 )       
       *  double d4 = divided ( fun , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 )       
       *  @endcode
       *  @see https://en.wikipedia.org/wiki/Divided_differences 
       */ 
      template <class FUNCTION, typename... Args>
      inline double divided
      ( FUNCTION       fun  , 
        const double   x0   , 
        const Args...  args , 
        const double   xn   ) 
      {
        const auto cfunr = std::cref ( fun ) ;
        return ( divided ( cfunr , args... , xn ) - 
                 divided ( cfunr , x0 , args... ) ) /  ( xn -  x0 ) ;
      }
      // ======================================================================
      /** Divided forward differences from two sequnces
       *  @param xbegin the start of sequence of abscissas 
       *  @param xend   the end   of sequence of abscissas 
       *  @param ybegin the start of sequence of function values 
       *  @return divided  difference clacualetd from these sequences
       */
      template <class XITERATOR, 
                class YITERATOR> 
      inline double divided
      ( XITERATOR xbegin , 
        XITERATOR xend   , 
        XITERATOR ybegin )
      {
        const auto N = std::distance ( xbegin , xend ) ;
        //
        if       ( 0 >= N ) { return 0.0     ; } // RETURN
        else  if ( 1 == N ) { return *ybegin ; } // RETURN
        //
        const long double x0 = * ( xbegin          ) ;
        const long double xn = * ( xbegin +  N - 1 ) ;
        //
        return 
          ( divided ( xbegin + 1 , xend           , ybegin + 1 ) * 1.0L - 
            divided ( xbegin     , xbegin + N - 1 , ybegin     ) ) / ( xn - x0 ) ;
      }
      // ======================================================================
      /** Divided forward differences from two sequnces
       *  @param xbegin the start of sequence of abscissas 
       *  @param xend   the end   of sequence of abscissas 
       *  @param ybegin the start of sequence of function values 
       *  @param xvalue adapter to get x-value from dereferenced x-iterator 
       *  @param yvalue adapter to get y-value from dereferenced y-iterator 
       *  @return divided  difference clacualetd from these sequences
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTER ,
                class YADAPTER >
      inline double divided 
      ( XITERATOR xbegin , 
        XITERATOR xend   , 
        YITERATOR ybegin , 
        XADAPTER  xvalue , 
        YADAPTER  yvalue ) 
      {
        const auto N = std::distance ( xbegin , xend ) ;
        //
        if       ( 0 >= N ) { return 0.0                ; } // RETURN
        else  if ( 1 == N ) { return yvalue ( *ybegin ) ; } // RETURN
        //
        const long double x0 = xvalue ( * ( xbegin          ) ) ;
        const long double xn = xvalue ( * ( xbegin +  N - 1 ) ) ;
        //
        const auto cxv = std::cref ( xvalue ) ;
        const auto cyv = std::cref ( yvalue ) ;
        //
        return 
          ( divided ( xbegin + 1 , xend            , ybegin + 1 , cxv , cyv ) * 1.0L - 
            divided ( xbegin     , xbegin + N - 1  , ybegin     , cxv , cyv ) ) 
          / ( xn - x0 ) ;
      }
      // ======================================================================
      // Finite differences
      // ======================================================================
      /** @class Forward
       *  simple evaluator of the Nth forward difference
       *  @see https://en.wikipedia.org/wiki/Finite_difference
       */
      template <unsigned int N> 
      class Forward
      {
      public:
        // ====================================================================
        /** get Nth difference
         *  @param fun the fnuction
         *  @param x   x-value 
         *  @param h   step 
         *  @return Nth difference 
         */
        template <class FUNCTION>
        inline 
        static double evaluate 
        ( FUNCTION     fun , 
          const double x   , 
          const double h   ) ;
        // ====================================================================
      } ;
      // ======================================================================
      /// specialization for the 0th difference (stop recursion)
      template <> 
      class Forward<0>
      {
      public:
        // ====================================================================
        template <class FUNCTION>
        inline 
        static double evaluate 
        ( FUNCTION        fun  , 
          const double    x    , 
          const double /* h */ ) { return fun ( x ) ; }
        // ====================================================================
      };
      // ====================================================================== 
      /// recursive method to get the Nth forward difference 
      template <unsigned int N> 
      template <class FUNCTION> 
      inline 
      double Forward<N>::evaluate 
      ( FUNCTION fun      , 
        const    double x , 
        const    double h ) 
      {
        const auto cfunr = std::cref ( fun ) ;
        return 
          Forward<N-1>::evaluate ( cfunr , x + h , h ) - 
          Forward<N-1>::evaluate ( cfunr , x     , h ) ;
      }
      // ======================================================================
      /** @class Backward
       *  simple evaluator of the Nth backward difference
       *  @see https://en.wikipedia.org/wiki/Finite_difference
       */
      template <unsigned int N> 
      class Backward
      {
      public:
        // ====================================================================
        /** get Nth difference
         *  @param fun the fnuction
         *  @param x   x-value 
         *  @param h   step 
         *  @return Nth difference 
         */
        template <class FUNCTION>
        inline 
        static double evaluate
        ( FUNCTION     fun , 
          const double x   , 
          const double h   ) ;
        // ====================================================================
      } ;
      // ======================================================================
      /// specialization for the 0th difference (stop recursion)
      template <> 
      class Backward<0>
      {
      public:
        // ====================================================================
        template <class FUNCTION>
        inline 
        static double evaluate 
        ( FUNCTION        fun  , 
          const double    x    , 
          const double /* h */ ) { return fun ( x ) ; }
        // ====================================================================
      };
      // ====================================================================== 
      /// recursive method to get the Nth forward difference 
      template <unsigned int N> 
      template <class FUNCTION> 
      inline 
      double Backward<N>::evaluate
      ( FUNCTION fun      , 
        const    double x , 
        const    double h ) 
      {
        const auto cfunr = std::cref ( fun ) ;
        return 
          Forward<N-1>::evaluate ( cfunr , x     , h ) - 
          Forward<N-1>::evaluate ( cfunr , x - h , h ) ;
      }
      // ======================================================================
      /** @class Central
       *  simple evaluator of the Nth central difference
       *  @see https://en.wikipedia.org/wiki/Finite_difference
       */
      template <unsigned int N> 
      class Central 
      {
      public:
        // ====================================================================
        /** get Nth difference
         *  @param fun the fnuction
         *  @param x   x-value 
         *  @param h   step 
         *  @return Nth difference 
         */
        template <class FUNCTION>
        inline 
        static double evaluate
        ( FUNCTION     fun , 
          const double x   , 
          const double h   ) ;
        // ====================================================================
      } ;
      // ======================================================================
      /// specialization for the 0th difference (stop recursion)
      template <> 
      class Central <0>
      {
      public:
        // ====================================================================
        template <class FUNCTION>
        inline 
        static double evaluate 
        ( FUNCTION        fun  , 
          const double    x    , 
          const double /* h */ ) { return fun ( x ) ; }
        // ====================================================================
      };
      // ====================================================================== 
      /// recursive method to get the Nth forward difference 
      template <unsigned int N> 
      template <class FUNCTION> 
      inline 
      double Central<N>::evaluate 
      ( FUNCTION fun      , 
        const    double x , 
        const    double h ) 
      {
        const auto cfunr = std::cref ( fun ) ;
        return 
          Central<N-1>::evaluate ( cfunr , x     , h ) - 
          Central<N-1>::evaluate ( cfunr , x - h , h ) ;
      }
      // =======================================================================
      /** recursively evaluate N-th forward difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th forward dirrerence 
       */
      template <class FUNCTION>
      inline double forward 
      ( FUNCTION             fun , 
        const unsigned short N   , 
        const double         x   , 
        const double         h   ) 
      {
        const auto cfunr = std::cref ( fun );
        return 
          0 == N ? fun ( x ) :
          forward ( cfunr , N - 1 , x + h , h ) - 
          forward ( cfunr , N - 1 , x     , h ) ;  
      }
      // ======================================================================
      /** recursively evaluate N-th backward difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th backward dirrerence 
       */
      template <class FUNCTION>
      inline double backward 
      ( FUNCTION             fun , 
        const unsigned short N   , 
        const double         x   , 
        const double         h   ) 
      {
        const auto cfunr = std::cref ( fun );
        return 
          0 == N ? fun ( x ) :
          backward ( cfunr , N - 1 , x     , h ) - 
          backward ( cfunr , N - 1 , x - h , h ) ;  
      }
      // ======================================================================
      /** recursivelyu evaluate N-th central difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th central dirrerence 
       */
      template <class FUNCTION>
      inline double central
      ( FUNCTION             fun , 
        const unsigned short N   , 
        const double         x   , 
        const double         h   ) 
      {
        const auto cfunr = std::cref ( fun );
        return 
          0 == N ? fun ( x ) :
          central ( cfunr , N - 1 , x + 0.5 * h , h ) - 
          central ( cfunr , N - 1 , x - 0.5 * h , h ) ;  
      }
      // ======================================================================
    } //                          The end of namespace Ostap::Math::Differences  
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DIFFERENCES_H
// ============================================================================
