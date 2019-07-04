// ============================================================================
#ifndef OSTAP_POLYNOMIALS_H 
#define OSTAP_POLYNOMIALS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <functional>
#include <vector>
#include <cmath>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Clenshaw.h"
// ============================================================================
/** @file Ostap/Polynomials.h
 *  various polinomials 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    namespace detail 
    {
      // ======================================================================
      template <typename F, std::size_t ... Is>
      auto make_array(F f, std::index_sequence<Is...>)
        -> std::array<std::decay_t<decltype(f(0u))>, sizeof...(Is)>
      { return {{f(Is)...}}; }
      // ======================================================================
    }
    // ========================================================================
    // Chebyshev 
    // ========================================================================
    namespace detail
    {
      // ======================================================================
      /** Evaluate Chebyshev polynom using the recurrence relation, 
       *  based on the fictive summation of the Chebyshev series 
       *  (0,0,...,0,1) using Clenshaw algorithm 
       *  @param N (input) the order of Chebyshev polynomial
       *  @param x (input) the point
       *  @return the value of Chebyshev polynomial of order <code>N</code> at point <code>x</code>
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2019-06-24
       *  @see Ostap::Math::detail::chebyshev_value
       */
      template <unsigned int N, typename TYPE>
      struct   Chebyshev_eval_ ;
      template <unsigned int N, typename TYPE>
      struct   Chebyshev_eval_
      {
        /// start recursion 
        static TYPE evaluate ( const TYPE x  , const TYPE b1 , const TYPE b2 ) 
        { return Chebyshev_eval_<N-1,TYPE>::evaluate 
            ( x , std::fma ( 2 * x ,  b1  , - b2 ) , b1 ) ; }
      } ;
      /// stop  recursion 
      template <typename TYPE>
      struct   Chebyshev_eval_<0,TYPE>
      {
        static TYPE evaluate ( const TYPE x  , const TYPE b1 , const TYPE b2 ) 
        { return x * b1 - b2 ; }
      } ;
      // ======================================================================
      constexpr inline long double chebyshev_eval_
      ( const unsigned int N  , 
        const long double  x  , 
        const long double  b1 , 
        const long double  b2 )
      {
        return 
          0 == N ? x * b1 - b2 : chebyshev_eval_ 
          ( N - 1 , x , std::fma ( 2 * x , b1 , -b2 ) , b1 ) ; 
      }
      // =============================================================================
    }
    // ===============================================================================
    /** Evaluate Chebyshev polynom using the recurrence relation, 
     *  based on the fictive summation of the Chebyshev series 
     *  (0,0,...,0,1) using Clenshaw algorithm 
     *  @param N (input) the order of Chebyshev polynomial
     *  @param x (input) the point
     *  @return the value of Chebyshev polynomial of order <code>N</code> at point <code>x</code>
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-06-24
     *  @see Ostap::Math::detail::chebyshev_eval_ 
     */
    inline double chebyshev_value ( const unsigned int N , const double x ) 
    {
      return 
        0 == N ? 1 : 
        1 == N ? x :
        detail::chebyshev_eval_ ( N - 1 , x , 1 , 0 ) ; 
    }
    // ========================================================================
    //  Chebyshev 1st kind 
    // ========================================================================
    template <unsigned int N> class  Chebyshev_  ;
    // ========================================================================
    /** @class Chebychev_
     *  Efficient evaluator of Chebyshev polynomial
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    template <unsigned int N>
    class  Chebyshev_
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double    x    ) const 
      { return evaluate ( x ) ; }
      /// evaluate the polynomial 
      static inline double evaluate ( const double x ) ;
      // ======================================================================      
    public:
      // ======================================================================      
      /// get the array of roots 
      static inline const std::array<double,N>&   roots   () ;
      /// get the array of extrema (the endpoints are not included)
      static inline const std::array<double,N-1>& extrema () ;
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for N=0
    template <>
    class  Chebyshev_<0>
    {
    public:
      // ======================================================================
      /// evaluate it!
      inline double operator() ( const double /* x */ ) const { return    1 ; }
      /// evaluate it!
      static inline double evaluate ( const double /* x */ ) { return 1 ; }
      // ======================================================================      
    public:
      // ======================================================================      
      /// get roots 
      static inline std::array<double,0> roots   () { return {{}} ; }
      /// get extrema
      static inline std::array<double,0> extrema () { return {{}} ; }
      // ======================================================================      
    } ;
    // ========================================================================
    /// specialization for N=1
    template <>
    class  Chebyshev_<1>
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double    x    ) const { return   x ; }
      /// the only one important method
      static inline double evaluate ( const double x  ) { return x ; }
      // ======================================================================      
    public: 
      // ======================================================================      
      /// get roots 
      static inline std::array<double,1> roots   () { return {{ 0.0 }} ; }
      ///     extrema 
      static inline std::array<double,0> extrema () { return {{}}      ; }
      // ======================================================================      
    } ;
    // ========================================================================
    /// the basic recursive method 
    template <unsigned int N>
    inline double Chebyshev_<N>::evaluate ( const double x ) 
    {
      // partly optmized recursion :
      // return 
      //   0 == N % 2 ?  
      //   2 * std::pow ( Chebyshev_<N/2>::evaluate ( x ) , 2 )                    - 1 :
      //   2 * Chebyshev_<N/2>::evaluate ( x ) * Chebyshev_<N/2+1>::evaluate ( x ) - x ;
      // optimized recursion 
      return detail::Chebyshev_eval_<N-1,long double>::evaluate ( x , 1 , 0 ) ;
    }
    // ========================================================================
    /// get the array of roots 
    template <unsigned int N>
    inline const std::array<double,N>& Chebyshev_<N>::roots ()
    {
      auto root = []( unsigned int k ) -> double  
        { return -std::cos ( ( 2 * k + 1 ) * M_PIl / ( 2 * N ) ) ; } ;
      static const std::array<double,N> s_roots = 
        detail::make_array ( root , std::make_index_sequence<N>() ) ;
      return s_roots ;
    }
    // ========================================================================
    /// get the array of extrema (the endpoints are not included)
    template <unsigned int N>
    inline const std::array<double,N-1>& Chebyshev_<N>::extrema ()
    {
      auto extremum = []( unsigned int k ) -> double 
        { return -std::cos ( ( k + 1 ) * M_PIl / N ) ; } ;
      static const std::array<double,N-1> s_extrema = 
        detail::make_array ( extremum , std::make_index_sequence<N-1>() ) ;
      return s_extrema ;
    }
    // ========================================================================
    //  Chebyshev 2nd kind 
    // ========================================================================
    template <unsigned int N> class  ChebyshevU_ ;
    // ========================================================================
    /** @class ChebychevU_
     *  Efficient evaluator of Chebyshev polynomial of the secon kind 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    template <unsigned int N>
    class  ChebyshevU_
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double    x    ) const { return evaluate ( x ) ; }
      /// the only one important method
      static inline double evaluate ( const double x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the array of roots 
      static inline const std::array<double,N>&   roots   () ;
      // ======================================================================
    } ;
    // ========================================================================    
    /// specialization for N=0
    template <>
    class  ChebyshevU_<0>
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double /* x */ ) const { return 1 ; }
      /// the only one important method
      static inline double evaluate ( const double /* x */ )  { return 1 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// roots 
      static inline std::array<double,0> roots   () { return {{}} ; }
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for N=1
    template <>
    class  ChebyshevU_<1> 
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double x ) const { return 2 * x ; }
      // ======================================================================
      static inline double evaluate ( const double x )  { return 2 * x ; }
      // ======================================================================
    public:
      // ======================================================================
      /// roots 
      static inline std::array<double,1> roots   () { return {{ 0.0 }} ; }
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for N=2
    template <>
    class  ChebyshevU_<2> 
    {
      // ======================================================================
    public: // evaluate 
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double x ) const { return 4 * x * x - 1 ; }
      /// the only one important method
      static inline double evaluate ( const double x )  { return 4 * x * x - 1 ; }
      // ======================================================================
    public: // roots & extrema 
      // ======================================================================
      /// roots 
      // ======================================================================
      static inline std::array<double,2> roots   () { return {{ -0.5 , 0.5 }} ; }
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for N=3
    template <>
    class  ChebyshevU_<3> 
    {
      // ======================================================================
    public: // evaluate 
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double x ) const 
      { return 4 * x * ( 2 * x *  x - 1 ) ; }
      // ======================================================================
      static inline double evaluate ( const double x ) 
      { return 4 * x * ( 2 * x * x - 1 ) ; }
      // ======================================================================
    public: // roots & extrema 
      // ======================================================================
      /// roots 
      static inline const std::array<double,3>& roots   () ;
      // ======================================================================
    } ;
    // ========================================================================
    /// the basic recurrence 
    template <unsigned int N>
    inline double ChebyshevU_<N>::evaluate ( const double x ) 
    {
      return 
        ChebyshevU_<N-2>::evaluate ( x ) * ( ChebyshevU_<2>::evaluate ( x ) - 1 )  
        - ChebyshevU_<N-4>::evaluate ( x ) ;
    }
    // ========================================================================
    /// get the array of roots 
    template <unsigned int N>
    inline const std::array<double,N>& ChebyshevU_<N>::roots ()
    {
      auto root = []( unsigned int k ) -> double 
        { return - std::cos ( ( k + 1 ) * M_PIl / ( N + 1 ) ) ; } ;
      static const std::array<double,N> s_roots = 
        detail::make_array ( root , std::make_index_sequence<N>() ) ;
      return s_roots ;
    }
    // ========================================================================

    // ========================================================================
    //  Chebyshev 3rd kind 
    // ========================================================================
    template <unsigned int N> class  Chebyshev3_ ;
    // ========================================================================
    /** @class Chebychev3_
     *  Efficient evaluator of Chebyshev polynomial of the third kind: 
     *  \f$ V_n^{(3)} = \frac{ \cos \left( n+\frac{1}{2}\right) \theta}
     *                       { \cos \frac{1}{2} \theta} \f$, where 
     *  \f$ x = \cos \theta\f$
     *  Also known as "Air-flow or airfoil polynomials"
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    template <unsigned int N>
    class  Chebyshev3_
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double    x    ) const { return evaluate ( x ) ; }
      /// the only one important method
      static inline double evaluate ( const double x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// roots 
      static inline const std::array<double,N>& roots () ; // roots 
      // ======================================================================
    } ;
    // ========================================================================    
    /// specialization for N=0
    template <>
    class  Chebyshev3_<0>
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double /* x */ ) const { return 1 ; }
      /// the only one important method
      static inline double evaluate ( const double /* x */ )  { return 1 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// roots 
      static inline std::array<double,0> roots   () { return {{}} ; }
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for N=1
    template <>
    class  Chebyshev3_<1> 
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double x ) const { return 2 * x - 1 ; }
      // ======================================================================
      static inline double evaluate ( const double x )  { return 2 * x - 1 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// roots 
      static inline std::array<double,1> roots   () { return {{ 0.5 }} ; }
      // ======================================================================
    } ;
    // ========================================================================
    /// the basic recursive method 
    template <unsigned int N>
    inline double Chebyshev3_<N>::evaluate ( const double x ) 
    {
      return 
        2 * x * Chebyshev3_<N-1>::evaluate ( x ) 
        -       Chebyshev3_<N-2>::evaluate ( x ) ;  
    } ;
    // ========================================================================
    /// get the array of roots 
    template <unsigned int N>
    inline const std::array<double,N>& Chebyshev3_<N>::roots ()
    {
      auto root = []( unsigned int k ) -> double 
        { return std::cos ( ( 2 * N - 2 * k - 1 ) * M_PIl / ( 2 * N + 1 ) ) ; } ;
      static const std::array<double,N> s_roots = 
        detail::make_array ( root , std::make_index_sequence<N>() ) ;
      return s_roots ;
    }
    // ========================================================================

    // ========================================================================
    //  Chebyshev 4th kind 
    // ========================================================================
    template <unsigned int N> class  Chebyshev4_ ;
    // ========================================================================
    /** @class Chebychev4_
     *  Efficient evaluator of Chebyshev polynomial of the third kind: 
     *  \f$ W_n^{(4)} = \frac{ \sin \left( n+\frac{1}{2}\right) \theta}
     *                       { \sin \frac{1}{2} \theta} \f$, where 
     *  \f$ x = \cos \theta\f$
     *  Also known as "Air-flow or airfoil polynomials"
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    template <unsigned int N>
    class  Chebyshev4_
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// the only one important method
      static inline double evaluate ( const double x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// roots 
      static inline const std::array<double,N>& roots () ; // roots 
      // ======================================================================
    } ;
    // ========================================================================
    /// the basic evaluaton method 
    template <unsigned int N>
    inline double Chebyshev4_<N>::evaluate ( const double x ) 
    { return ( N % 2 ? -1 : 1 ) * Chebyshev3_<N>::evaluate ( -x ) ; } ;
    // ========================================================================
    /// get the array of roots 
    template <unsigned int N>
    inline const std::array<double,N>& Chebyshev4_<N>::roots ()
    {
      auto root = []( unsigned int k ) -> double 
        { return std::cos ( ( 2 * N - 2 * k ) * M_PIl / ( 2 * N + 1 ) ) ; } ;
      static const std::array<double,N> s_roots = 
        detail::make_array ( root , std::make_index_sequence<N>() ) ;
      return s_roots ;
    }
    // ========================================================================
    /** Calculate the k-th root of Legendre polynomial of order n
     *  @param k root number
     *  @param n legendre polynomial order 
     *  @return k-th root of Legendre polynomial of order n
     */
    double legendre_root ( const unsigned short k , 
                           const unsigned short n ) ;    
    // ========================================================================
    namespace detail
    {
      // ======================================================================
      /** Evaluate Legendre polynom using the recurrence relation, 
       *  based on the fictive summation of the Legendre series 
       *  (0,0,...,0,1) using Clenshaw algorithm 
       *  @param N (input) the order of Legendere polynomial
       *  @param x (input) the point
       *  @return the value of Legendre polynomial of order <code>N</code> at point <code>x</code>
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2019-06-24
       *  @see Ostap::Math::detail::legendre_value
       */
      template <unsigned int N, typename TYPE>
      struct   Legendre_eval_ ;
      template <unsigned int N, typename TYPE>
      struct   Legendre_eval_
      {
        /// start recursion 
        static TYPE evaluate ( const TYPE x  , const TYPE b1 , const TYPE b2 ) 
        { return Legendre_eval_<N-1,TYPE>::evaluate 
            ( x , ( 2 * N + 1 ) * x * b1 / ( N + 1 ) - ( N + 1 ) * b2 / ( N + 2 ) , b1 ) ; }
      };
      /// stop  recursion 
      template <typename TYPE>
      struct   Legendre_eval_<0,TYPE>
      {
        static TYPE evaluate ( const TYPE x  , const TYPE b1 , const TYPE b2 ) 
        { return x * b1 - b2 / 2 ; }
      };
      // ============================================================================
      constexpr inline long double legendre_eval_
      ( const unsigned int N  , 
        const long double  x  , 
        const long double  b1 , 
        const long double  b2 )
      {
        return 
          0 == N ? x * b1 - b2 / 2 : legendre_eval_ 
          ( N - 1 , x , ( 2 * N + 1 ) * x * b1 / ( N + 1 ) - ( N + 1 ) * b2 / ( N + 2 ) , b1 ) ; 
      }
      // ======================================================================
    }
    // ========================================================================
    /** Evaluate Legendre polynom using the recurrence relation, 
     *  based on the fictive summation of the Legendre series 
     *  (0,0,...,0,1) using Clenshaw algorithm 
     *  @param N (input) the order of Legendere polynomial
     *  @param x (input) the point
     *  @return the value of Legendre polynomial of order <code>N</code> at point <code>x</code>
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-06-24
     *  @see Ostap::Math::detail::Legendre_eval_
     */
    inline double legendre_value ( const unsigned int N , const double x ) 
    {
      return 
        0 == N ? 1 : 
        1 == N ? x :
        detail::legendre_eval_ ( N - 1 , x , 1 , 0 ) ; 
    }
    // ========================================================================
    /** calculate sequence of Legendre polynomials \f$ a_i = P_i(x) \f$  
     *  @param begin start  of the sequece
     *  @param end   end  of the sequence
     *  @param x     x-value 
     */
    template <class ITERATOR>
    inline void legendre_values
    ( ITERATOR          begin, 
      ITERATOR          end  ,
      const long double x    )
    {
      if ( begin == end ) { return ; }
      long double p_0 = 1 ;
      *begin = p_0 ; ++begin;
      if ( begin == end ) { return ; }
      long double p_1 = x ;
      *begin = p_1 ; ++begin;
      if ( begin == end ) { return ; }
      long double p_i = 0 ;
      unsigned int  i = 2 ;
      while ( begin != end ) 
      {
        //
        p_i    =  ( ( 2 * i - 1 ) * x * p_1  - ( i - 1 ) * p_0 ) / i ;
        p_0    = p_1 ;
        p_1    = p_i ;
        //
        *begin = p_i ;
        ++begin ;
        ++i     ;
      }
    } 
    // ========================================================================
    /** calculate the integral for Legendre polynomial
     *  \f$ \int_{x_{low}}^{x_{high}}P_N(x)\deriv x \f$ 
     *  @param N the order/degree of Legendre polynomial
     *  @param xlow  the low  edge 
     *  @param xhigh the high edge 
     *  @return the integral
     */
    inline long double legendre_integral 
    ( const unsigned int N     , 
      const long double  xlow  , 
      const long double  xhigh ) 
    {
      return 
        0 == N ?       xhigh - xlow                      :
        0 == 1 ? 0.5*( xhigh - xlow ) * ( xhigh + xlow ) :
        ( detail::legendre_eval_( N - 2 , xhigh , -1 + ( 2 * N - 1 ) * xhigh * xhigh / N , xhigh ) -
          detail::legendre_eval_( N - 2 , xlow  , -1 + ( 2 * N - 1 ) * xlow  * xlow  / N , xlow  ) ) / ( N + 1 ) ;
    }
    // ========================================================================
    //  Legendre 
    // ========================================================================
    template <unsigned int N> class  Legendre_   ;
    // ========================================================================
    /** @class Legendre_
     *  Efficient evaluator of Legendre polynomial
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    template <unsigned int N>
    class  Legendre_ 
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator()        ( const double x ) const { return evaluate ( x ) ; }
      /// calculate the polynomial
      static inline double evaluate   ( const double x ) ;
      /// calculate the derivative 
      static inline double derivative ( const double x ) ;
      // ======================================================================
      /// get the roots of Legendre polynomial 
      static inline const std::array<double,N>& roots() ;
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for 0
    template <>
    class  Legendre_<0> 
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator()        ( const double /* x */ ) const { return 1 ; }
      /// calculate the polynomial
      static inline double evaluate   ( const double /* x */ )       { return 1 ; }
      /// calculate the derivative 
      static inline double derivative ( const double /* x */ )       { return 0 ; }
      /// get the roots of Legendre polymonial 
      static inline std::array<double,0> roots() { return {{}} ; }
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for N=1
    template <>
    class  Legendre_<1> 
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator()        ( const double x ) const { return x ; }
      /// calculaet the polynomial
      static inline double evaluate   ( const double x )       { return x ; }
      /// calculate the derivative 
      static inline double derivative ( const double /* x */ ) { return 1 ; }
      /// get the roots of Legendre polymonial 
      static inline std::array<double,1> roots() { return {{ 0.0 }} ; }
      // ======================================================================
    } ;
    // ========================================================================
    /// should one use here legendre_eval ? 
    template <unsigned int N>
    inline double Legendre_<N>::evaluate ( const double x ) 
    {
      // naive recursive :
      // return ( ( 2 * N - 1 ) * x * Legendre_<N-1>::evaluate ( x )  - 
      // (     N - 1 )     * Legendre_<N-2>::evaluate ( x ) ) / N ;
      // optimized recursive 
      return detail::Legendre_eval_<N-1,long double>::evaluate ( x , 1 , 0 ) ;
    }
    // ========================================================================
    /// get the array of roots 
    template <unsigned int N>
    inline const std::array<double,N>& Legendre_<N>::roots ()
    {
      auto root = []( unsigned int k ) -> double { return legendre_root ( k , N ) ; } ;
      static const std::array<double,N> s_roots = 
        detail::make_array ( root , std::make_index_sequence<N>() ) ;
      return s_roots ;
    }
    // =======================================================================
    /// calculate the derivative 
    template <unsigned int N>
    inline double Legendre_<N>::derivative ( const double x ) 
    {
      // /// 1) naive recursive algorithm, not safe when |x| is  close to 1 
      // static const Ostap::Math::Equal_To<double> s_equal {} ;
      // return 
      //   x >  0.999 && s_equal ( x ,  1 ) ? 0.5 * N * ( N + 1 )                      :
      //   x < -0.999 && s_equal ( x , -1 ) ? 0.5 * N * ( N + 1 ) * ( N % 2 ? 1 : -1 ) :
      //   N * ( x * Legendre_<N>::evaluate ( x ) - Legendre_<N-1>::evaluate ) / ( x * x - 1 ) ;
      //
      // /// 2) the recursive algorithm that is safe for |x| close to 1 
      // return 
      //   N * Legendre_<N-1>::evaluate   ( x ) + 
      //   x * Legendre_<N-1>::derivative ( x ) ;
      //
      /// 3) and this algorithm is much faster (linear):
      auto ak = []  ( const unsigned int k )  -> unsigned int 
        { return (k+N)%2 ? 2*k+1 : 0 ; };
      auto alpha = [] ( const unsigned int k , const long double    y    ) -> long double 
        { return (2*k+1)*y/(k+1) ; } ;
      auto beta  = [] ( const unsigned int k , const long double /* y */ ) -> long double 
        { return -1.0L*k*1.0L/(k+1) ; } ;
      auto phi0  = [] ( const long double /* y */ ) -> long double { return 1 ; } ;
      auto phi1  = [] ( const long double    y    ) -> long double { return y ; } ;
      //
      return Ostap::Math::Clenshaw::sum 
        ( x  , N-1  , ak , alpha , beta , phi0 , phi1 );
    }
    // ========================================================================
    // Hermite
    // ========================================================================
    namespace detail
    {
      // ======================================================================
      /** Evaluate Hermite  polynom using the recurrence relation, 
       *  based on the fictive summation of the Hermite series 
       *  (0,0,...,0,1) using Clenshaw algorithm 
       *  @param N (input) the order of Hermite polynomial
       *  @param x (input) the point
       *  @return the value of Chebyshev polynomial of order <code>N</code> at point <code>x</code>
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2019-06-24
       *  @see Ostap::Math::detail::hermite_value
       */
      template <unsigned int N, typename TYPE>
      struct   Hermite_eval_ ;
      template <unsigned int N, typename TYPE>
      struct   Hermite_eval_
      {
        /// start recursion 
        static TYPE evaluate ( const TYPE x  , const TYPE b1 , const TYPE b2 ) 
        { return Hermite_eval_<N-1,TYPE>::evaluate 
            ( x , x * b1 - ( N + 1 ) * b2 ,  b1 ) ; }
      } ;
      /// stop  recursion 
      template <typename TYPE>
      struct   Hermite_eval_<0,TYPE>
      {
        static TYPE evaluate ( const TYPE x  , const TYPE b1 , const TYPE b2 ) 
        { return x * b1 - b2 ; }
      } ;
      // ======================================================================
      constexpr inline long double hermite_eval_
      ( const unsigned int N  , 
        const long double  x  , 
        const long double  b1 , 
        const long double  b2 )
      {
        return 
          0 == N ? x * b1 - b2 : hermite_eval_ 
          ( N - 1 , x , x * b1 - ( N + 1 ) * b2 , b1 ) ; 
      }
      // =============================================================================
    }
    // ===============================================================================
    /** Evaluate Hermite polynom using the recurrence relation, 
     *  based on the fictive summation of the Hermite series 
     *  (0,0,...,0,1) using Clenshaw algorithm 
     *  @param N (input) the order of Hermite polynomial
     *  @param x (input) the point
     *  @return the value of Hermite polynomial of order <code>N</code> at point <code>x</code>
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-06-24
     *  @see Ostap::Math::detail::hermite_eval_ 
     */
    inline double hermite_value ( const unsigned int N , const double x ) 
    {
      return 
        0 == N ? 1 : 
        1 == N ? x :
        detail::hermite_eval_ ( N - 1 , x , 1 , 0 ) ; 
    }
    // ========================================================================
    // Hermite
    // ========================================================================
    template <unsigned int N> class  Hermite_ ;
    // ========================================================================
    /** @class Hermite_
     *  Efficienct evaluator of Hermite polynomial
     *  These are "probabilistic" polinomials,
     *  \f$He(x)\f$ 
     *  such as coefficienst at maximar degree is always equal to 1 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    template <unsigned int N>
    class  Hermite_ 
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double    x    ) const
      { return evaluate ( x ) ; }
      // ======================================================================
      static inline double evaluate ( const double x ) ;
      // ======================================================================
    } ;
    // ========================================================================
    /// specialization for 0
    template <>
    class  Hermite_<0> 
    {
    public:
      // ======================================================================
      inline double operator() ( const double /* x */ ) const { return 1 ; }
      // ======================================================================
      static inline double evaluate ( const double /* x */ ) { return  1 ; }
      // ======================================================================     
    } ;
    // ========================================================================
    /// specialization for N=1
    template <>
    class  Hermite_<1> 
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double    x    ) const { return   x ; }
      // ======================================================================
      static inline double evaluate ( const double x ) { return  x ; }
      // ======================================================================     
    } ;
    // ========================================================================
    template <unsigned int N>
    inline double Hermite_<N>::evaluate ( const double x ) 
    { 
      // naive recursion 
      // return x * Hermite_<N-1>::evaluate ( x ) - ( N - 1 ) * Hermite_<N-2>::evaluate ( x ) ; 
      // optimized recursion 
      return detail::Hermite_eval_<N-1,long double>::evaluate ( x , 1 , 0 );
    }
    // ========================================================================
    // Associated Legendre functions 
    // ========================================================================
    namespace detail
    {
      // ======================================================================
      /// evaluate normalized assiciated Legendre polynomials/functions
      template <unsigned int L, unsigned int M , typename = void>
      class PLegendre_helper_ ;
      // 
      template <unsigned int L, unsigned int M>
      class PLegendre_helper_<L,M,typename std::enable_if<L==0&&M==0>::type>
      {
      public : 
        inline double operator()      ( const double    x ) const { return evaluate ( x ) ; }
        static inline double evaluate ( const double /* x */ )    
        {
          static const long double s_P00 = std::sqrt ( 1.0L / ( 4 * M_PIl ) ) ;
          return s_P00 ; 
        }
      };
      //
      template <unsigned int L, unsigned int M>
      class PLegendre_helper_<L,M,typename std::enable_if<L==1&&M==1>::type>
      {
      public : 
        inline double operator()      ( const double x ) const { return evaluate ( x ) ; }
        static inline double evaluate ( const double x )    
        { 
          static const long double s_n = - std::sqrt( 3.0L / 2.0L ) * 
            PLegendre_helper_<0,0>::evaluate ( 0 ) ;
          return s_n * std::sqrt ( 1.0L - x * x ) ;
        }
      };
      //
      template <unsigned int L, unsigned int M>
      class PLegendre_helper_<L,M,typename std::enable_if<L==M&&L!=0&&L!=1>::type>
      {
      public : 
        inline double operator()      ( const double x ) const { return evaluate ( x ) ; }
        static inline double evaluate ( const double x )    
        { 
          static const long double s_n = 
            std::sqrt (  ( 2 * M + 1 ) * ( 2 * M - 1 ) * 0.25L / ( M * ( M - 1 ) ) ) ;
          return s_n * ( 1.0L -  x * x ) * PLegendre_helper_<L-2,L-2>::evaluate ( x ) ;
        }
      };
      //
      template <unsigned int L, unsigned int M>
      class PLegendre_helper_<L,M,typename std::enable_if<L==M+1>::type>
      {
      public : 
        inline double operator()      ( const double x ) const { return evaluate ( x ) ; }
        static inline double evaluate ( const double x )    
        { 
          static const long double s_n = std::sqrt ( 2 * M + 3.0L ) ;
          return s_n * x * PLegendre_helper_<M,M>::evaluate ( x ) ;
        }
      };
      //
      template <unsigned int L, unsigned int M>
      class PLegendre_helper_<L,M,typename std::enable_if<M+2<=L>::type>
      {
      public : 
        inline double operator()      ( const double x ) const { return evaluate ( x ) ; }
        static inline double evaluate ( const double x )    
        {
          //
          long double p0 = PLegendre_helper_<M  ,M>::evaluate ( x ) ;
          long double p1 = PLegendre_helper_<M+1,M>::evaluate ( x ) ;
          long double pN = 0 ;
          //
          unsigned int N = M + 2 ;
          while ( L >= N ) 
          {
            pN = a ( N ) * x * p1 - b ( N ) * p0 ;
            p0 = p1 ;
            p1 = pN ;
            ++N ;
          }
          //
          return pN ;
        }
      private :
        // ====================================================================
        static inline long double a ( const unsigned int J ) 
        {
          auto   afun = [] ( const unsigned int K ) -> long  double 
            {
              const unsigned int I = K + M + 1 ;
              return std::sqrt ( ( 2 * I - 1 ) * 1.0L * ( 2 * I + 1 ) / ( I * I - M * M ) ) ;
            } ;
          static const std::array<long double, L-M> s_a = 
            detail::make_array ( afun , std::make_index_sequence<L-M>() ) ;
          return s_a [ J - ( M + 1 ) ] ;
        }
        static inline long double b ( const unsigned int J ) 
        { 
          auto   bfun = [] ( const unsigned int K ) -> long  double 
            { 
              const unsigned int I = K + M + 2 ;
              return  PLegendre_helper_<L,M>::a(I)/PLegendre_helper_<L,M>::a(I-1) ;
            } ;
          static const std::array<long double, L-M-1> s_b = 
            detail::make_array ( bfun , std::make_index_sequence<L-M-1>() ) ;
          return s_b [ J - ( M + 2 ) ] ;
        }
        // ====================================================================
      };
      // =====================================================================
      /** evaluate the normalized associated legendre polynomials/function 
       *  \$ P^{m}_{l}(z)\$ with normalization suitbale for spherical harmonics :
       *  \f$ \int_{-1}{+1}P^m_l(x)P^m_{l}(x)\deriv x = \frac{1}{2\pi}\f$
       */
      inline long double plegendre_eval_ 
      ( const unsigned int L , 
        const unsigned int M , 
        const long double  x ) 
      {
        //
        if ( M  >  L ){ return 0 ; }
        else if ( L == 0 && M == 0   ) { return PLegendre_helper_<0,0>::evaluate ( x ) ; }
        else if ( L == 1 && M == 1   ) { return PLegendre_helper_<1,1>::evaluate ( x ) ; }
        else if ( L == M ) 
        {
          long double result = ( 0 == L%2 ) ? 
            PLegendre_helper_<0,0>::evaluate ( x ) :
            PLegendre_helper_<1,1>::evaluate ( x ) ;
          const unsigned int L0 = ( 0 == L%2 ) ? 2 : 3 ;
          for ( unsigned int l  = L0 ; l <= M ; l += 2 ) 
          {
            result *= ( 1.0L - x * x ) *   
              std::sqrt ( ( l + 0.5L ) * ( l - 0.5L ) / ( ( l - 1.0L ) * l ) ) ;
          } 
          return result ;
        }
        else if ( L == M && 0 == L%2 ) 
        {
          long double result = PLegendre_helper_<0,0>::evaluate ( x ) ;
          for ( unsigned int l = 2;  l <= M ; l += 2 ) 
          {
            result *= ( 1.0L - x * x ) *   
              std::sqrt ( ( l + 0.5L ) * ( l - 0.5L ) / ( ( l - 1.0L ) * l ) ) ;
          } 
          return result ;
        }
        else if ( L == M + 1 ) 
        { return std::sqrt ( 2 * M + 3.0L ) * x * plegendre_eval_ ( M , M , x ) ; }
        //
        /// regular case 
        //
        long double p0 = plegendre_eval_ ( M     , M , x ) ;
        long double p1 = plegendre_eval_ ( M + 1 , M , x ) ;
        // long double p1 = std::sqrt ( 2 * M + 3.0L ) * x * p0 ; 
        long double pN = 0 ;
        //
        auto afun = []      ( const unsigned int J , const unsigned int M ) -> long double 
          { return std::sqrt ( ( 2 * J - 1 ) * 1.0L * ( 2 * J + 1 ) / ( J * J - M * M ) ) ; } ;
        auto bfun = [&afun] ( const unsigned int J , const unsigned int M ) -> long double 
          { return afun ( J , M ) / afun ( J - 1 , M ) ; } ;
        //
        unsigned int N = M + 2 ;
        while ( L >= N ) 
        {
          pN = afun ( N , M ) * x * p1 - bfun ( N , M ) * p0 ;
          p0 = p1 ;
          p1 = pN ;
          ++N ;
        } 
        //
        return pN ;
      } ;
      // ======================================================================
    }
    // ========================================================================
    /** @class PLegendre_
     *  The normalized associated Legendre polynomials/functions 
     *  Normalization is suitable for usage of them for the spherical harmonics.
     *  @see https://arxiv.org/abs/1410.1748
     *  \f$ \int _{-1}^{1} P_l^m(x)P_l^{m}(x) \deriv x = \frac{1}{2\pi}\f$  
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-06-27
     */
    template  <unsigned int L , unsigned int M>
    class PLegendre_ : public detail::PLegendre_helper_<L,M>
    { static_assert ( M <= L , "PLegendre_ : M must me  M <= L" ); }; 
    // =====================================================================
    /** evaluate the normalized associated legendre polynomials/function 
     *  \$ P^{m}_{l}(z)\$ with normalization suitbale for spherical harmonics :
     *  \f$ \int_{-1}{+1}P^m_l(x)P^m_{l}(x)\deriv x = \frac{1}{2\pi}\f$
     */
    inline double plegendre_value
    ( const unsigned int L , 
      const unsigned int M , 
      const double       x ) 
    { return  detail::plegendre_eval_ ( L  , M , x ) ; }
    // ========================================================================
    // Non-templated 
    // ========================================================================
    /** @class Chebyshev
     *  evaluate the chebyshev polynomials
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class Chebyshev 
    {
    public :
      // ======================================================================
      /// constructor
      Chebyshev ( const unsigned int N = 0 ) : m_N ( N ) {}
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the polynomial
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate the polynomial
      inline double evaluate   ( const double x ) const 
      { return chebyshev_value ( m_N , x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      unsigned int degree () const { return m_N ; }
      // ======================================================================
    public:
      // ======================================================================
      /// derivative 
      double derivative ( const double x    ) const ;
      // ======================================================================
      /// get integral between low and high 
      double integral   ( const double low  , 
                          const double high ) const ;
      // ======================================================================
    public: // roots & extrema 
      // ======================================================================
      /// get all roots   of the polynomial 
      std::vector<double> roots   () const ;
      /// get all extrema of the polynomial 
      std::vector<double> extrema () const ;
      // ======================================================================
    private:
      // ======================================================================
      unsigned int m_N ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ChebyshevU
     *  evaluate the chebyshev polynomials of the second kind 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class ChebyshevU 
    {
    public :
      // ======================================================================
      /// constructor
      ChebyshevU ( const unsigned int N = 0 ) : m_N ( N ) {}
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the polynomial
      double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      unsigned int degree () const { return m_N ; }
      // ======================================================================
    public:
      // ======================================================================
      /// derivative 
      double derivative ( const double x ) const ;
      // ======================================================================
      /// get integral between low and high 
      double integral   ( const double low  , 
                          const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      unsigned int m_N ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Hermite
     *  evaluate the Hermite polynomials
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class Hermite
    {
    public :
      // ======================================================================
      /// constructor
      Hermite ( const unsigned int N = 0  ) ; 
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the polynomial
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      inline double evaluate   ( const double x ) const 
      { return hermite_value ( m_N , x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      unsigned int degree () const { return m_N ; }
      // ======================================================================
    public:
      // ======================================================================
      /// derivative 
      double derivative  ( const double x ) const 
      { return 0 == m_N ? 0 : m_N * hermite_value ( m_N - 1 , x ) ; }
      // ======================================================================
      /// get integral between low and high 
      double integral    ( const double low  , const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      unsigned int m_N ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Legendre
     *  evaluate the Legendre polynomials
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class Legendre 
    {
    public :
      // ======================================================================
      /// constructor
      Legendre ( const unsigned int N = 0  ) : m_N ( N ) {}
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the polynomial
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      inline double evaluate   ( const double x ) const 
      { return legendre_value ( m_N , x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      unsigned int degree () const { return m_N ; }
      // ======================================================================
    public:
      // ======================================================================
      /// derivative 
      double derivative  ( const double x ) const ;
      // ======================================================================
      /// get integral between low and high 
      double integral    ( const double low  , 
                           const double high ) const ;
      // ======================================================================
    public: // roots 
      // ======================================================================
      /// get the root of the Legendre polynomial
      double                     root  ( const unsigned short i ) const ;
      /// get all roots of the Legendre polynomial
      const std::vector<double>& roots () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the root of the Legendre polynomial
      double calculate_root ( const unsigned short i ) const ;
      // ======================================================================
    private:
      // ======================================================================
      unsigned int m_N ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PLegendre
     *  evaluate the associative Legendre polynomials/functions 
     *  \f$P^{m}_{l}(x)\f$
     *  Normalization is sutable for spherical harmonics 
     *  \f$ \int_{-1}^{+1} P^{m}_{l}(x)P^{m}_{l}(x) \deriv x = \frac{1}{2\pi} \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class PLegendre
    {
    public :
      // ======================================================================
      /// constructor
      PLegendre ( const unsigned int L = 0 , 
                  const unsigned int M = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the polynomial
      inline double operator() ( const double x ) const { return   evaluate ( x ) ; }
      inline double evaluate   ( const double x ) const 
      { return plegendre_value ( m_L , m_M , x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      unsigned int L () const { return m_L ; }
      unsigned int M () const { return m_M ; }
      unsigned int l () const { return m_L ; }
      unsigned int m () const { return m_M ; }
      // ======================================================================
    public:
      // ======================================================================
      unsigned int m_L ;
      unsigned int m_M ;
      // ======================================================================
    } ;
    // ========================================================================
    /** affine transformation of polynomial
     *  \f$ x ^{\prime} = \alpha x + \beta \f$
     *  @param input  (INPUT)  input polynomial coefficients 
     *  @param result (UPDATE) coefficients of transformed polynomial 
     *  @param alpha  (INPUT)  alpha
     *  @param beta   (INPUT)  beta
     *  @return true for valid transformations
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-09
     */
    bool affine_transform
    ( const std::vector<double>& input      , 
      std::vector<double>&       result     ,  
      const double               alpha  = 1 , 
      const double               beta   = 0 ) ;              
    // ========================================================================
    /// forward declarations 
    class Bernstein ; // forward declarations 
    // ========================================================================
    /** @class Parameters 
     *  Holder for parameters 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-10
     */
    class Parameters 
    {
    public:
      // =======================================================================
      /// constructor from number of parameters 
      Parameters ( const unsigned int          np = 1 ) ; 
      /// constructor from  the list of parameters 
      Parameters ( const std::vector<double>&  pars   ) ;
      /// constructor from  the list of parameters 
      Parameters (       std::vector<double>&& pars   ) ;
      /// templated constructor from the sequnce of parameters 
      template <class ITERATOR>
      Parameters ( ITERATOR begin , 
                   ITERATOR end   )
        : m_pars ( begin , end )
      {}
      /// copy constructor  
      /// Parameters ( const Parameters&  ) = default ;
      /// move constructor  
      /// Parameters (       Parameters&& ) = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of parameters 
      unsigned short npars  () const { return m_pars.size()     ; }
      /// all parameters are zero ?
      bool           zero   () const ;
      /** set k-parameter
       *  @param k index
       *  @param value new value 
       *  @return true if parameter is actually changed 
       */
      bool setPar          ( const unsigned short k , const double value ) ;
      /** set k-parameter
       *  @param k index
       *  @param value new value 
       *  @return true iof parameter is actually changed 
       */
      bool setParameter    ( const unsigned short k , const double value )
      { return setPar      ( k , value ) ; }
      /// get the parameter value
      double  par          ( const unsigned short k ) const
      { return ( k < m_pars.size() ) ? m_pars[k] : 0.0 ; }
      /// get the parameter value
      double  parameter    ( const unsigned short k ) const { return par ( k ) ; }
      /// get all parameters:
      const std::vector<double>& pars () const { return m_pars ; }
      // ======================================================================
    public: // simple  manipulations with parameters 
      // ======================================================================
      /// simple  manipulations with parameters: scale it! 
      /// Parameters& operator *= ( const double a ) ;     // scale it! 
      /// simple  manipulations with parameters: scale it! 
      /// Parameters& operator /= ( const double a ) ;     // scale it! 
      // ======================================================================
    protected:
      // ======================================================================
      /// swap two parameter sets 
      void swap ( Parameters& right ) ;
      // ======================================================================
    protected :
      // ======================================================================
      /// parameters 
      std::vector<double> m_pars ; //  vector of parameters 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PolySum
     *  Base class for polynomial sums 
     *  \f$ f(x) = \sum_i \alpha_i P_i(x) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-10
     */
    class PolySum : public Parameters 
    {
    public:
      // ======================================================================
      /// constructor from polynomial degree 
      PolySum ( const unsigned short degree = 0  ) ;
      /// constructor from vector of parameters 
      PolySum ( const std::vector<double>&  pars ) ;
      /// constructor from vector of parameters 
      PolySum (       std::vector<double>&& pars ) ;
      /// constructor from sequence of parameters 
      template <class ITERATOR>
        PolySum ( ITERATOR begin , 
                  ITERATOR end   )
        : Parameters ( begin , end )
      { if ( m_pars.empty() ) { m_pars.push_back ( 0 ) ; } }
      // ======================================================================
    public:
      // ======================================================================
      /// degree  of polynomial 
      unsigned short degree () const { return m_pars.size() - 1 ; }
      // ======================================================================
    } ;  
    // ========================================================================    
    /// forward declarations 
    class Bernstein    ; // forward declarations 
    class Polynomial   ; // forward declarations 
    class LegendreSum  ; // forward declarations 
    class ChebyshevSum ; // forward declarations 
    class HermiteSum   ; // forward declarations 
    // ========================================================================
    // Polynomial sums
    // ========================================================================
    /** @class Polynomial
     *  Trivial polynomial
     *  \f$ f(x) = \sum_i p_i x^i\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-22
     */
    class Polynomial : public PolySum 
    {
    public:
      // =====================================================================
      /// constructor from the degree  
      Polynomial ( const unsigned short       degree =   0  , 
                   const double               xmin   =  -1  , 
                   const double               xmax   =   1  ) ;
      // ======================================================================
      /// constructor from the parameter list 
      Polynomial ( const std::vector<double>& pars          , 
                   const double               low    =  -1  , 
                   const double               high   =   1  ) ;
      /// template constructor from sequence of parameters 
      template <class ITERATOR>
        Polynomial ( ITERATOR                 first , 
                     ITERATOR                 last  , 
                     const double             xmin  , 
                     const double             xmax  ) 
        : Ostap::Math::PolySum ( first , last ) 
        , m_xmin ( std::min ( xmin, xmax ) )
        , m_xmax ( std::max ( xmin, xmax ) )
      {}
      // ======================================================================
      /// copy 
      Polynomial ( const Polynomial&  ) = default ;
      /// move 
      Polynomial (       Polynomial&& ) = default ;
      // ======================================================================
      ///  constructor from Bernstein polinomial (efficient) 
      explicit Polynomial ( const Bernstein&     poly ) ;
      ///  constructor from Legendre polinomial  (efficient)
      explicit Polynomial ( const LegendreSum&   poly ) ;
      ///  constructor from Chebyshev polinomial (delegation) 
      explicit Polynomial ( const ChebyshevSum&  poly ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double evaluate    ( const double x ) const ;
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const 
      { return x < m_xmin ? 0 : x > m_xmax ? 0 : evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin ; }
      /// get upper edge
      double xmax  () const { return m_xmax ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const
      { return  0.5 * ( t * ( m_xmax - m_xmin ) +   m_xmax + m_xmin ) ; }
      double t ( const double x ) const
      { return (  2 *   x   - m_xmax - m_xmin ) / ( m_xmax - m_xmin ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const ;
      /// get the integral between low and high 
      double integral   ( const double low , const double high ) const ;
      /// get the derivative at point "x" 
      double derivative ( const double x     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get indefinte integral 
      Polynomial indefinite_integral ( const double C = 0 ) const ;
      /// get the derivative
      Polynomial derivative          () const ;
      // ======================================================================
     public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it! 
      Polynomial& operator += ( const double a ) ; 
      /// simple  manipulations with polynoms: shift it! 
      Polynomial& operator -= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it 
      Polynomial& operator *= ( const double a ) ; 
      /// simple  manipulations with polynoms: scale it  
      Polynomial& operator /= ( const double a ) ;
      /// negate it! 
      Polynomial  operator-() const ; // negate it! 
      // ======================================================================
    public:
      // ======================================================================
      /// Add       polynomials (with the same domain!)
      Polynomial sum      ( const Polynomial& other ) const ;
      /// Subtract  polynomials (with the same domain!)
      Polynomial subtract ( const Polynomial& other ) const ;
      // ======================================================================      
    public:
      // ======================================================================      
      /// Add       polynomials (with the same domain!)
      Polynomial __add__   ( const Polynomial& other ) const ;
      /// Subtract  polynomials (with the same domain!)
      Polynomial __sub__   ( const Polynomial& other ) const ;
      // ======================================================================      
    public:
      // ======================================================================      
      Polynomial& __iadd__      ( const double a )       ;
      Polynomial& __isub__      ( const double a )       ;
      Polynomial& __imul__      ( const double a )       ;
      Polynomial& __itruediv__  ( const double a )       ;
      Polynomial& __idiv__      ( const double a )       { return __itruediv__ ( a ) ; }
      // ======================================================================
    public:
      // ======================================================================
      Polynomial  __add__       ( const double a ) const ;
      Polynomial  __sub__       ( const double a ) const ;
      Polynomial  __mul__       ( const double a ) const ;
      Polynomial  __truediv__   ( const double a ) const ;
      Polynomial  __div__       ( const double a ) const { return __truediv__  ( a ) ; }      
      // ======================================================================
    public:
      // ======================================================================
      Polynomial  __radd__  ( const double a ) const ;
      Polynomial  __rsub__  ( const double a ) const ;
      Polynomial  __rmul__  ( const double a ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// Negate it! 
      Polynomial  __neg__   () const ; // Negate it! 
      // ======================================================================      
    private:
      // ======================================================================
      /// x-min 
      double              m_xmin ; // x-min
      /// x-max 
      double              m_xmax ; // x-max
      // ======================================================================      
    } ;    
    // ========================================================================
    inline Polynomial operator+( const Polynomial& a , const Polynomial& b ) 
    { return a.sum       ( b ) ; }    
    inline Polynomial operator-( const Polynomial& a , const Polynomial& b ) 
    { return a.subtract  ( b ) ; }
    inline Polynomial operator+( const Polynomial& a , const double      b ) 
    { return a.__add__   ( b ) ; }    
    inline Polynomial operator+( const double      b , const Polynomial& a )
    { return a.__add__   ( b ) ; }
    inline Polynomial operator-( const Polynomial& a , const double      b ) 
    { return a.__sub__   ( b ) ; }    
    inline Polynomial operator-( const double      b , const Polynomial& a )
    { return a.__rsub__  ( b ) ; }
    inline Polynomial operator*( const Polynomial& a , const double      b ) 
    { return a.__mul__   ( b ) ; }    
    inline Polynomial operator*( const double      b , const Polynomial& a )
    { return a.__mul__   ( b ) ; }
    inline Polynomial operator/( const Polynomial& a , const double      b ) 
    { return a.__truediv__   ( b ) ; }    
    // ========================================================================
    /** @class ChebyshevSum 
     *  Sum of chebychev polinomials 
     *  \f$ f(x) = \sum_i p_i T_i(x)\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-22
     */
    class ChebyshevSum : public PolySum 
    {
    public:
      // =====================================================================
      /// constructor from the degree  
      ChebyshevSum ( const unsigned short       degree =  0 , 
                     const double               xmin   = -1 , 
                     const double               xmax   =  1 )  ;
      // ======================================================================
      /// constructor from the parameter list 
      ChebyshevSum ( const std::vector<double>& pars       , 
                     const double               xmin  = -1 , 
                     const double               xmax  =  1 )  ;
      /// template constructor from sequence of parameters 
      template <class ITERATOR>
      ChebyshevSum ( ITERATOR                 first , 
                     ITERATOR                 last  , 
                     const double             xmin  , 
                     const double             xmax  ) 
        : Ostap::Math::PolySum ( first , last ) 
        , m_xmin ( std::min ( xmin, xmax ) )
        , m_xmax ( std::max ( xmin, xmax ) )
      {}
      // ======================================================================
      /// copy 
      ChebyshevSum ( const ChebyshevSum&  ) = default ;
      /// copy 
      ChebyshevSum (       ChebyshevSum&& ) = default ;
      // ======================================================================
      ///  constructor from Polinomial           (efficient)
      explicit ChebyshevSum ( const Polynomial&  poly ) ;
      ///  constructor from Bernstein            (delegation) 
      explicit ChebyshevSum ( const Bernstein&   poly ) ;
      ///  constructor from Legendre             (delegation) 
      explicit ChebyshevSum ( const LegendreSum& poly ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double evaluate    ( const double x ) const ;
      /// get the value
      double operator () ( const double x ) const 
      { return x < m_xmin ? 0 : x > m_xmax ? 0 : evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin () const { return m_xmin ; }
      /// get upper edge
      double xmax () const { return m_xmax ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const
      { return  0.5 * ( t * ( m_xmax - m_xmin ) +   m_xmax + m_xmin ) ; }
      double t ( const double x ) const
      { return (  2 *   x   - m_xmax - m_xmin ) / ( m_xmax - m_xmin ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const ;
      /// get the integral between low and high 
      double integral   ( const double low , const double high ) const ;
      /// get the derivative at point "x" 
      double derivative ( const double x     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get indefinte integral 
      ChebyshevSum indefinite_integral ( const double C = 0 ) const ;
      /// get the derivative 
      ChebyshevSum derivative          () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it! 
      ChebyshevSum& operator += ( const double a ) ; 
      /// simple  manipulations with polynoms: shift it! 
      ChebyshevSum& operator -= ( const double a ) ; 
      /// simple  manipulations with polynoms: scale it!
      ChebyshevSum& operator *= ( const double a ) ; 
      /// simple  manipulations with polynoms: scale it! 
      ChebyshevSum& operator /= ( const double a ) ; 
      // ======================================================================
    public:
      // ======================================================================
      /// negate it! 
      ChebyshevSum  operator-() const ; // negate it! 
      // ======================================================================
    public:
      // ======================================================================
      /// add      chebyshev sum (with the same domain)
      ChebyshevSum sum      ( const ChebyshevSum& other ) const ;
      /// subtract chebyshev sum (with the same domain)
      ChebyshevSum subtract ( const ChebyshevSum& other ) const ;
      // ======================================================================
    public:
      // ======================================================================
      ChebyshevSum& __iadd__      ( const double a ) ;
      ChebyshevSum& __isub__      ( const double a ) ;
      ChebyshevSum& __imul__      ( const double a ) ;
      ChebyshevSum& __itruediv__  ( const double a ) ;
      ChebyshevSum& __idiv__      ( const double a ) { return __itruediv__ ( a ) ; }
      // ======================================================================
    public:
      // ======================================================================
      ChebyshevSum  __add__       ( const double a ) const ;
      ChebyshevSum  __sub__       ( const double a ) const ;
      ChebyshevSum  __mul__       ( const double a ) const ;
      ChebyshevSum  __truediv__   ( const double a ) const ;
      ChebyshevSum  __div__       ( const double a ) { return __truediv__  ( a ) ; }
      // ======================================================================
    public:
      // ======================================================================
      ChebyshevSum  __radd__  ( const double a ) const ;
      ChebyshevSum  __rsub__  ( const double a ) const ;
      ChebyshevSum  __rmul__  ( const double a ) const ;
      // ======================================================================
    public:
      // ======================================================================
      ChebyshevSum  __add__   ( const ChebyshevSum& a ) const ;
      ChebyshevSum  __sub__   ( const ChebyshevSum& a ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // negate it! 
      ChebyshevSum __neg__ () const ; // negate it!
      // ======================================================================
    private:
      // ======================================================================
      /// x-min 
      double              m_xmin ; // x-min
      /// x-max 
      double              m_xmax ; // x-max
      // ======================================================================      
    } ;
    // ========================================================================
    inline ChebyshevSum operator+( const ChebyshevSum& a , const ChebyshevSum& b ) 
    { return a.sum       ( b ) ; }    
    inline ChebyshevSum operator-( const ChebyshevSum& a , const ChebyshevSum& b ) 
    { return a.subtract  ( b ) ; }
    inline ChebyshevSum operator+( const ChebyshevSum& a , const double        b ) 
    { return a.__add__   ( b ) ; }    
    inline ChebyshevSum operator+( const double        b , const ChebyshevSum& a )
    { return a.__add__   ( b ) ; }
    inline ChebyshevSum operator-( const ChebyshevSum& a , const double        b ) 
    { return a.__sub__   ( b ) ; }    
    inline ChebyshevSum operator-( const double        b , const ChebyshevSum& a )
    { return a.__rsub__  ( b ) ; }
    inline ChebyshevSum operator*( const ChebyshevSum& a , const double        b ) 
    { return a.__mul__   ( b ) ; }    
    inline ChebyshevSum operator*( const double        b , const ChebyshevSum& a )
    { return a.__mul__   ( b ) ; }
    inline ChebyshevSum operator/( const ChebyshevSum& a , const double        b ) 
    { return a.__truediv__   ( b ) ; }    
    // ========================================================================
    /** @class LegendreSum 
     *  Sum of Legendre polinomials 
     *  \f$ f(x) = \sum_i p_i P_i(x)\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-22
     */
    class LegendreSum : public PolySum
    {
    public:
      // =====================================================================
      /// constructor from the degree 
      LegendreSum ( const unsigned short        degree =  0 , 
                    const double                xmin   = -1 , 
                    const double                xmax   =  1 )  ;
      // ======================================================================
      /// constructor from the parameter list 
      LegendreSum ( const std::vector<double>&  pars       , 
                    const double                xmin   = -1 , 
                    const double                xmax   =  1 )  ;
      /// template constructor from sequence of parameters 
      template <class ITERATOR>
        LegendreSum ( ITERATOR                  first , 
                      ITERATOR                  last  , 
                      const double              xmin  , 
                      const double              xmax  ) 
        : Ostap::Math::PolySum ( first , last ) 
        , m_xmin ( std::min ( xmin, xmax ) )
        , m_xmax ( std::max ( xmin, xmax ) )
      {}
      // ======================================================================
      /// copy 
      LegendreSum  ( const LegendreSum&  ) = default ;
      /// move
      LegendreSum  (       LegendreSum&& ) = default ;
      // ======================================================================
      /**  constructor from Bernstein polinomial (efficient)
       *  @see http://www.sciencedirect.com/science/article/pii/S0377042700003769 eq.21
       */
      explicit LegendreSum ( const Bernstein&     poly ) ;
      /// constructor from polynoimial            (delegation)
      explicit LegendreSum ( const Polynomial&    poly ) ;
      /// constructor from Chebyshev              (delegation)
      explicit LegendreSum ( const ChebyshevSum&  poly ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double evaluate    ( const double x ) const ;
      /// get the value
      double operator () ( const double x ) const 
      { return x < m_xmin ? 0 : x > m_xmax ? 0 : evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin () const { return m_xmin ; }
      /// get upper edge
      double xmax () const { return m_xmax ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const
      { return  0.5 * ( t * ( m_xmax - m_xmin ) +   m_xmax + m_xmin ) ; }
      double t ( const double x ) const
      { return (  2 *   x   - m_xmax - m_xmin ) / ( m_xmax - m_xmin ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const ;
      /// get the integral between low and high 
      double integral   ( const double low , const double high ) const ;
      /// get the derivative at point "x" 
      double derivative ( const double x     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get indefinte integral 
      LegendreSum indefinite_integral ( const double C = 0 ) const ;
      /// get the derivative 
      LegendreSum derivative          () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it! 
      LegendreSum& operator += ( const double a ) ;
      /// simple  manipulations with polynoms: shift it! 
      LegendreSum& operator -= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it  
      LegendreSum& operator *= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it 
      LegendreSum& operator /= ( const double a ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// negate it! 
      LegendreSum  operator-() const ; // negate it! 
      // ======================================================================
    public:
      // ======================================================================
      /// add      legendre sum (with the same domain)
      LegendreSum sum      ( const LegendreSum& other ) const ;
      /// subtract legendre sum (with the same domain)
      LegendreSum subtract ( const LegendreSum& other ) const ;
      // ======================================================================
    public:
      // ======================================================================
      LegendreSum& __iadd__      ( const double a ) ;
      LegendreSum& __isub__      ( const double a ) ;
      LegendreSum& __imul__      ( const double a ) ;
      LegendreSum& __itruediv__  ( const double a ) ;
      LegendreSum& __idiv__      ( const double a ) { return __itruediv__ ( a ) ; }
      // ======================================================================
    public:
      // ======================================================================
      LegendreSum  __add__       ( const double a ) const ;
      LegendreSum  __sub__       ( const double a ) const ;
      LegendreSum  __mul__       ( const double a ) const ;
      LegendreSum  __truediv__   ( const double a ) const ;
      LegendreSum  __div__       ( const double a ) { return __truediv__  ( a ) ; }
      // ======================================================================
    public:
      // ======================================================================
      LegendreSum  __radd__      ( const double a ) const ;
      LegendreSum  __rsub__      ( const double a ) const ;
      LegendreSum  __rmul__      ( const double a ) const ;
      // ======================================================================
    public:
      // ======================================================================
      LegendreSum  __add__       ( const LegendreSum& a ) const ;
      LegendreSum  __sub__       ( const LegendreSum& a ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // negate it! 
      LegendreSum __neg__ () const ; // negate it!
      // ======================================================================
    public:
      // ======================================================================
      /** update  the Legendre expansion by addition of one "event" with 
       *  the given weight
       *  @code
       *  LegendreSum sum = ... ;
       *  for ( auto x : .... ) { sum.fill ( x ) ; }
       *  @endcode
       *  This is a useful function to make an unbinned parameterization 
       *  of certain distribution and/or efficiency 
       *  @parameter x      the event content 
       *  @parameter weight the weight 
       */
      bool fill ( const double x , const double weight = 1 ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// x-min 
      double              m_xmin ; // x-min
      /// x-max 
      double              m_xmax ; // x-max
      // ======================================================================      
    } ;
    // ========================================================================
    inline LegendreSum operator+( const LegendreSum& a , const LegendreSum& b ) 
    { return a.sum       ( b ) ; }    
    inline LegendreSum operator-( const LegendreSum& a , const LegendreSum& b ) 
    { return a.subtract  ( b ) ; }
    inline LegendreSum operator+( const LegendreSum& a , const double       b ) 
    { return a.__add__   ( b ) ; }    
    inline LegendreSum operator+( const double       b , const LegendreSum& a )
    { return a.__add__   ( b ) ; }
    inline LegendreSum operator-( const LegendreSum& a , const double       b ) 
    { return a.__sub__   ( b ) ; }    
    inline LegendreSum operator-( const double       b , const LegendreSum& a )
    { return a.__rsub__  ( b ) ; }
    inline LegendreSum operator*( const LegendreSum& a , const double       b ) 
    { return a.__mul__   ( b ) ; }    
    inline LegendreSum operator*( const double       b , const LegendreSum& a )
    { return a.__mul__   ( b ) ; }
    inline LegendreSum operator/( const LegendreSum& a , const double       b ) 
    { return a.__truediv__   ( b ) ; }    
    // ========================================================================
    /** @class HermiteSum 
     *  Sum of Hermite polinomials 
     *  \f$ f(x) = \sum_i p_i He_i(x)\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-08-08
     */
    class HermiteSum : public PolySum
    {
    public:
      // =====================================================================
      /// constructor from the degree 
      HermiteSum ( const unsigned short       degree =   0  ,
                   const double               xmin   =  -1  , 
                   const double               xmax   =   1  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double evaluate    ( const double x ) const ;
      /// get the value
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin () const { return m_xmin ; }
      /// get upper edge
      double xmax () const { return m_xmax ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const 
      { return 0.5 * ( t / m_scale + m_xmin + m_xmax ) ; }
      double t ( const double x ) const 
      { return m_scale * ( 2 * x   - m_xmin - m_xmax ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high 
      double integral   ( const double low , const double high ) const ;
      /// get the derivative at point "x" 
      double derivative ( const double x     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get indefinte integral 
      HermiteSum indefinite_integral ( const double C = 0 ) const ;
      /// get the derivative 
      HermiteSum derivative          () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it! 
      HermiteSum& operator += ( const double a ) ;
      /// simple  manipulations with polynoms: shift it! 
      HermiteSum& operator -= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it!
      HermiteSum& operator *= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it! 
      HermiteSum& operator /= ( const double a ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// negate it! 
      HermiteSum  operator-() const ; // negate it! 
      // ======================================================================
    public:
      // ======================================================================
      /// add      legendre sum (with the same domain)
      HermiteSum sum      ( const HermiteSum& other ) const ;
      /// subtract legendre sum (with the same domain)
      HermiteSum subtract ( const HermiteSum& other ) const ;
      // ======================================================================
    public:
      // ======================================================================
      HermiteSum& __iadd__      ( const double a ) ;
      HermiteSum& __isub__      ( const double a ) ;
      HermiteSum& __imul__      ( const double a ) ;
      HermiteSum& __itruediv__  ( const double a ) ;
      HermiteSum& __idiv__      ( const double a ) { return __itruediv__ ( a ) ; }
      // ======================================================================
    public:
      // ======================================================================
      HermiteSum  __add__       ( const double a ) const ; 
      HermiteSum  __sub__       ( const double a ) const ;
      HermiteSum  __mul__       ( const double a ) const ;
      HermiteSum  __truediv__   ( const double a ) const ;
      HermiteSum  __div__       ( const double a ) { return __truediv__  ( a ) ; }
      // ======================================================================
    public:
      // ======================================================================
      HermiteSum  __radd__  ( const double a ) const ;
      HermiteSum  __rsub__  ( const double a ) const ;
      HermiteSum  __rmul__  ( const double a ) const ;
      // ======================================================================
    public:
      // ======================================================================
      HermiteSum  __add__   ( const HermiteSum& a ) const ;
      HermiteSum  __sub__   ( const HermiteSum& a ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // negate it! 
      HermiteSum __neg__    () const ; // negate it!
      // ======================================================================
    private:
      // ======================================================================
      /// low  edge 
      double m_xmin  ; // low  edge 
      /// high edge 
      double m_xmax  ; // high edge 
      /// scale 
      double m_scale ; // scale 
      // ======================================================================      
    } ;
    // ========================================================================
    inline HermiteSum operator+( const HermiteSum& a , const HermiteSum& b ) 
    { return a.sum       ( b ) ; }    
    inline HermiteSum operator-( const HermiteSum& a , const HermiteSum& b ) 
    { return a.subtract  ( b ) ; }
    inline HermiteSum operator+( const HermiteSum& a , const double      b ) 
    { return a.__add__   ( b ) ; }    
    inline HermiteSum operator+( const double      b , const HermiteSum& a )
    { return a.__add__   ( b ) ; }
    inline HermiteSum operator-( const HermiteSum& a , const double      b ) 
    { return a.__sub__   ( b ) ; }    
    inline HermiteSum operator-( const double      b , const HermiteSum& a )
    { return a.__rsub__  ( b ) ; }
    inline HermiteSum operator*( const HermiteSum& a , const double      b ) 
    { return a.__mul__   ( b ) ; }    
    inline HermiteSum operator*( const double      b , const HermiteSum& a )
    { return a.__mul__   ( b ) ; }
    inline HermiteSum operator/( const HermiteSum& a , const double      b ) 
    { return a.__truediv__   ( b ) ; }    
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    // helper utilities for integration of product of polynomial and an exponent
    // ========================================================================
    /** get the integral between low and high for a product of Bernstein
     *  polynom and the exponential function with the exponent tau
     *  \f[  \int_{a}^{b} \mathcal{B} e^{\tau x } \mathrm{d}x \f] 
     *  @param poly  bernstein polynomial
     *  @param tau   slope parameter for exponential 
     *  @param a     low  integration range 
     *  @param b     high integration range 
     */
    double integrate 
    ( const Ostap::Math::Bernstein& poly ,
      const double                  tau  ,
      const double                  a    , 
      const double                  b    ) ;
    // ========================================================================    
    /** get the integral between low and high for a product of
     *  polynom and the exponential function with the exponent tau
     *  \f[ r = \int_{a}^{b} \mathcal{P} e^{\tau x } \mathrm{d}x \f] 
     *  @param poly  polynomial
     *  @param tau   slope parameter for exponential 
     *  @param a     low  integration range 
     *  @param b     high integration range 
     */
    double integrate 
    ( const Ostap::Math::Polynomial& poly ,
      const double                   tau  ,
      const double                   a    , 
      const double                   b    ) ;
    // ========================================================================    
    /** get the integral between low and high for a product of
     *  Chebyshev polynom and the exponential function with the exponent tau
     *  \f[ r = \int_{a}^{b} \mathcal{T} e^{\tau x } \mathrm{d}x \f] 
     *  @param poly  chebyshev polynomial
     *  @param tau   slope parameter for exponential 
     *  @param a     low  integration range 
     *  @param b     high integration range 
     */
    double integrate 
    ( const Ostap::Math::ChebyshevSum& poly ,
      const double                     tau  ,
      const double                     a    , 
      const double                     b    ) ;
    // ========================================================================    
    /** get the integral between low and high for a product of
     *  Legendre polynom and the exponential function with the exponent tau
     *  \f[ r = \int_{a}^{b} \mathcal{L} e^{\tau x } \mathrm{d}x \f] 
     *  @param poly  Legendre polynomial
     *  @param tau   slope parameter for exponential 
     *  @param a     low  integration range 
     *  @param b     high integration range 
     */
    double integrate 
    ( const Ostap::Math::LegendreSum& poly ,
      const double                    tau  ,
      const double                    a    , 
      const double                    b    ) ;
    // ========================================================================    

    // ========================================================================    
    // special cases:
    // ========================================================================    
    
    // ========================================================================    
    /** get the integral between low and high for a product of
     *  polynom and the exponential function with the exponent tau
     *  \f[  r = \int_{a}^{b} \mathcal{P} e^{\tau x } \mathrm{d}x \f] 
     *  @param poly  polynomial
     *  @param tau   slope parameter for exponential 
     */
    double integrate 
    ( const Ostap::Math::Polynomial& poly ,
      const double                   tau  ) ;
    // ========================================================================    
    
    // ========================================================================
    /** construct chebyshev approximation for arbitrary function 
     *  @param func the function
     *  @param x_min low edge
     *  @param x_max high edge 
     *  @return Chebyshev approximation 
     *  @see ChebyshevSum 
     *  @code 
     *  FUNC func = ...
     *  ChebyshevSum a = chebyshev_sum<6> ( func , x_min , x_max ) ;
     *  @endcode 
     */
    template <unsigned short N, class FUNCTION>
    inline ChebyshevSum 
    chebyshev_sum ( FUNCTION     func  , 
                    const double x_min , 
                    const double x_max ) 
    { 
      // array of precomputed function values 
      std::array<double,N>   fv ;
      // 
      const double xmin = std::min ( x_min , x_max ) ;
      const double xmax = std::max ( x_min , x_max ) ;
      //
      const double xhs = 0.5 * ( xmin + xmax ) ;
      const double xhd = 0.5 * ( xmax - xmin ) ;
      const long double pi_N = M_PIl / N ;
      auto _xi_ = [xhs,xhd,pi_N] ( const unsigned short k ) 
        { return std::cos ( pi_N * ( k + 0.5 ) ) * xhd + xhs ; } ;
      for ( unsigned short i = 0 ; i < N ; ++i ) { fv[i] = func ( _xi_ ( i ) ) ; }
      //
      ChebyshevSum cs ( N , xmin , xmax ) ;
      for ( unsigned short i = 0 ; i < N + 1 ; ++i ) 
      {
        double c_i = 0 ;
        if ( 0 == i ) 
        { for ( unsigned short k = 0 ; k < N ; ++k ) { c_i += fv[k] ; } }
        else 
        {
          for ( unsigned short k = 0 ; k < N ; ++k ) 
          { c_i += fv[k] * std::cos ( pi_N * i * ( k + 0.5 ) ) ; }
        }
        c_i *= 2.0 / N ;
        if ( 0 == i ) { c_i *= 0.5 ;}
        cs.setPar ( i, c_i ) ;
      }
      return cs ;
    }              
    // ========================================================================
    /*  construct chebyshev approximation for arbitrary function 
     *  @param func the function
     *  @param N    degree of polynomial 
     *  @param x_min low edge
     *  @param x_max high edge 
     *  @return Chebyshev approximation 
     *  @see ChebyshevSum 
     *  @code 
     *  FUNC func = ...
     *  ChebyshevSum a = chebyshev_sum ( func , 10 ,  xmin , xmax ) ;
     *  @endcode 
     */
    ChebyshevSum 
    chebyshev_sum ( std::function<double(double)> func  , 
                    const unsigned short          N     , 
                    const double                  x_min , 
                    const double                  x_max ) ;    
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_POLYNOMIALS_H
// ============================================================================
