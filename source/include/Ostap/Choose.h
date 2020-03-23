// ============================================================================ 
#ifndef OSTAP_CHOOSE_H 
#define OSTAP_CHOOSE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <array>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Clenshaw.h"
// ============================================================================
/** @file Ostap/Choose.h
 *  Binomial coefficients, Pochhammer's symbols and Stirling numbers 
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {    
    // ========================================================================
    /** calculate the binomial coefficient C(n,k) = n!/((n-k)!*k!)
     *  the result is exact for all n,k<=67
     *  @warning In case of overflow std::numeric_limits<unsigned long long>::max is returned 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    unsigned long long  choose 
    ( const unsigned short n ,
      const unsigned short k ) ;
    // ========================================================================
    /** calculate the inverse binomial conefficient 
     *  \f$ a = C(n,k)^{-1} = \frac{ (n-k)!k!}{n!}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2020-01-31
     */
    double ichoose 
    ( const unsigned short n , 
      const unsigned short k ) ;
    // ========================================================================
    /** calculate the logarithm of binomial coefficient
     *  \f$ \log C^n_k \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double log_choose
    ( const unsigned short n ,
      const unsigned short k ) ;
    // ========================================================================
    /** calculate the binomial coefficient C(k,n) = n!/((n-k)!*k!)
     *  @author Vanya BELYAEV Ivan.Belyaev@irep.ru
     *  @date 2015-03-08
     */
    double choose_double 
    ( const unsigned short n , 
      const unsigned short k ) ;
    // ========================================================================
    /** calculate the generalized binomial coefficient C(a,k) 
     *  \f$C(\alpha,k) = \frac{\alpha}{k}\frac{\alpha-1}{k-1}...\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double gen_choose
    ( const double         a ,
      const unsigned short k ) ;
    // ========================================================================
    /** calculate the generalized binomial coefficient C(n/2,k) 
     *  \f$C(n,k) = \frac{n/2}{k}\frac{n/2-1}{k-1}...\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double choose_half
    ( const int            n ,
      const unsigned short k ) ;
    // ========================================================================
    
    // ========================================================================
    /** @class Choose 
     *  Binomial coefficents    \f$ C^n_k = \frac{n!}{k!(n-k)!}\f$
     *  computed recursively as \f$ C^{n-1}_{k-1} + C^{n-1}_{k}\f$ with
     *  with initial/boundary conditions 
     *  \f$ C^n_0 = 1 \f$  and \f$C^n_n=1\f$
     */
    template <unsigned short N, unsigned short K>
    struct Choose
    {
      enum _ : unsigned long long { value =  
          N <  K           ? 0 :
          0 == K || K == N ? 1 :
          Choose<N-1,K-1>::value + 
          Choose<N-1,K  >::value } ;
    };
    // ========================================================================
    /** @class Stirling1
     *  unsigned Stirling numbers of 1st kind 
     *  @see https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
     */
    template <unsigned short N, unsigned short K>
    struct Stirling1 ;
    // ========================================================================
    /// stop recursion: \f$ S^0_0 = 1 \f$
    template <>
    struct Stirling1<0,0> { enum _ :  unsigned long long { value = 1 } ; } ;    
    // ========================================================================
    /// stop recursion \f$ S^n_0 = 0 \f$
    template <unsigned short N>
    struct Stirling1<N,0> { enum _ : unsigned long long { value = 0 } ; } ;    
    // ========================================================================
    /// stop recursion \f$ S^0_s = 0 \f$
    template <unsigned short K>
    struct Stirling1<0,K> { enum _ : unsigned long long { value = 0 } ; } ;
    /// start recursion  \f$ S^{n+1}_k = n S^n_k + S^n_{k-1}\f$ 
    template <unsigned short N, unsigned short K>
    struct Stirling1 
    {
      enum _ : unsigned long long { value =  
          Stirling1<N-1,K  >::value * ( N - 1 ) +
          Stirling1<N-1,K-1>::value } ;
    } ;
    // ========================================================================
    /** Compile-time generation of the sequence of Stirling numbers of 1st kind 
     *  useful for implementatinoof Pochhammer symbols 
     *  @code
     *  std::array<double,7> p7 = stirling_array<double,6>() ;
     *  @endcode 
     */
    template<class TYPE,unsigned short N,size_t... i>
    constexpr auto stirling1_array ( std::index_sequence<i...> )
    { return std::array<TYPE,N+1>{{Stirling1<N,N-i>::value...}}; }
    template<class TYPE,unsigned short N>
    constexpr auto stirling1_array() 
    { return stirling1_array<TYPE,N>(std::make_index_sequence<N+1>{}); }
    // ========================================================================
    /** @class Pochhammer_ 
     *  Pochhammer symbols as polynomials 
     *  \f[ P(x,n) = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @see Ostap::Math::rising_factorial
     *  @see Ostap::Math::pochhammer
     *  @see Ostap::Math::stirling1_array
     */
    template <unsigned short N>
    class Pochhammer_
    {
    public:
      // ======================================================================
      /// the only one important method
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================      
    public:
      // ======================================================================      
      /// evaluate the polynomial 
      static inline double evaluate   ( const double x ) 
      { return Ostap::Math::Clenshaw::monomial_sum ( s_coeff , x ).first  ; }
      /// get the derivative 
      static inline double derivative ( const double x ) 
      { return Ostap::Math::Clenshaw::monomial_sum ( s_coeff , x ).second ; }      
      /// get the value and the derivative 
      static inline std::pair<double,double> value_with_derivative  ( const double x ) 
      { return Ostap::Math::Clenshaw::monomial_sum ( s_coeff , x )        ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// the polynomial coefficients 
      static const std::array<unsigned long long,N+1> s_coeff ;
      // ======================================================================
    } ;
    // ========================================================================
    template<unsigned short N>  
    constexpr std::array<unsigned long long,N+1>
    Pochhammer_<N>::s_coeff { Ostap::Math::stirling1_array<unsigned long long,N>() } ;
    // ========================================================================
    /** calculate unsigned Stirling number of 1st kind 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    unsigned long long stirling1 ( const unsigned short n ,
                                   const unsigned short k ) ;
    // =========================================================================
    /** calculate unsigned Stirling number of 1st kind 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double stirling1_double ( const unsigned short n ,
                              const unsigned short k ) ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CHOOSE_H
// ============================================================================
