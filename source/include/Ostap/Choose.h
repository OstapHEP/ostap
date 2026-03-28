// ============================================================================ 
#ifndef OSTAP_CHOOSE_H 
#define OSTAP_CHOOSE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <cstdint>
#include <array>
#include <vector>
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
    /** @var N_CHOOSE_MAX
     *  maximal valeu of n, fo which all  \f$ C^n_k \f$ fit into 64 bit integer
     */	
    inline constexpr unsigned short N_CHOOSE_MAX = 68 ;
    // ========================================================================
    /** @struct Choose_
     *  Binomial coefficient \f$ C^n_k = \left(\begin{array}{c}n \\ k \end{array}\rght) = \frac{n!}{k!(n-k)!}\f$
     */
    template <unsigned short N,
	      unsigned short K>
    struct Choose_ ;
    // ========================================================================
    // various template specializations:  \f$ C^0_0 = 1\f$ 
    template <>
    struct Choose_<0,0>   { enum _ { value = 1 } ; } ;
    //
    /// \f$ C^n_0 = 1\f$ 
    template <unsigned short N>
    struct Choose_<N,0>   { enum _ { value = 1 } ; } ;
    template <unsigned short N>
    struct Choose_<N,N>   { enum _ { value = 1 } ; } ; 
    //
    template <unsigned short N,
	      unsigned short K>
    struct Choose_ : public std::enable_if<(0<K)&&(K<N)&&(N<=N_CHOOSE_MAX),void>
    {
      // ======================================================================
      static_assert( N<=N_CHOOSE_MAX , "Choose: argument is too large" ) ;
      // ======================================================================
      enum _ : std::uintmax_t { value = Choose_<N,K-1>::value * ( N + 1 - K ) / K } ;
    } ;
    // ========================================================================
    /** @class Pascal
     *  The row in Pascal triangle
     */  
    template <unsigned short N>
    struct Pascal_ : public std::enable_if<(N<=N_CHOOSE_MAX)>
    {
      // ======================================================================
      static_assert( N<=N_CHOOSE_MAX , "Pascal's triange: argument is too large" ) ;
      // ======================================================================
      /// the actual type of the row in trianlge
      typedef std::array<std::uintmax_t,N+1>  Row ;
      /// get the row of coefficients 
      template <std::size_t... i>
      inline static constexpr Row choose_array ( std::index_sequence<i...> )
      { return Row { { Choose_<N,i>::value... } } ; }
      // ======================================================================
      /// get the row of coefficients 
      inline static constexpr Row choose_array ()
      { return Pascal_<N>::choose_array ( std::make_index_sequence<N+1> {} ) ; }
      // ======================================================================
      /// get the row of coefficients 
      static constexpr const Row row { choose_array () } ;
      // ======================================================================
    } ;
    // ========================================================================
    /** get the row in Pascal's triangle - the row of binomial coefficients 
     *  @see Pascal_
     */
    template<unsigned short N>
    inline constexpr const auto& choose_array() { return Ostap::Math::Pascal_<N>::row ; }
    // ========================================================================
    /** calculate the binomial coefficient C(n,k) = n!/((n-k)!*k!)
     *  the result is exact for all n,k<=67
     *  @warning In case of overflow std::numeric_limits<std::uintmax_t>::max is returned 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    std::uintmax_t choose 
    ( const unsigned short n ,
      const unsigned short k ) ;
    // ========================================================================
    /** calculate the inverse binomial coefficient 
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
    long double
    choose_double 
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
    /** @class Stirling1
     *  unsigned Stirling numbers of 1st kind 
     *  @see https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
     */
    template <unsigned short N, unsigned short K>
    struct Stirling1 ;
    // ========================================================================
    /// stop recursion: \f$ S^0_0 = 1 \f$
    template <>
    struct Stirling1<0,0> { enum _ : std::uintmax_t { value = 1 } ; } ;    
    // ========================================================================
    /// stop recursion \f$ S^n_0 = 0 \f$
    template <unsigned short N>
    struct Stirling1<N,0> { enum _ : std::uintmax_t { value = 0 } ; } ;    
    // ========================================================================
    /// stop recursion \f$ S^0_s = 0 \f$
    template <unsigned short K>
    struct Stirling1<0,K> { enum _ : std::uintmax_t { value = 0 } ; } ;
    /// start recursion  \f$ S^{n+1}_k = n S^n_k + S^n_{k-1}\f$ 
    template <unsigned short N, unsigned short K>
    struct Stirling1 
    {
      enum _ : std::uintmax_t { value =  
          Stirling1<N-1,K  >::value * ( N - 1 ) +
          Stirling1<N-1,K-1>::value } ;
    } ;
    // ========================================================================
    /** Compile-time generation of the sequence of Stirling numbers of 1st kind 
     *  useful for implementation of Pochhammer symbols 
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
      static const std::array<std::uintmax_t,N+1> s_coeff ;
      // ======================================================================
    } ;
    // ========================================================================
    template<unsigned short N>  
    constexpr std::array<std::uintmax_t,N+1>
    Pochhammer_<N>::s_coeff { Ostap::Math::stirling1_array<std::uintmax_t,N>() } ;
    // ========================================================================
    /** calculate unsigned Stirling number of 1st kind 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    std::uintmax_t 
    stirling1
    ( const unsigned short n ,
      const unsigned short k ) ;
    // =========================================================================
    /** calculate unsigned Stirling number of 1st kind 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double stirling1_double
    ( const unsigned short n ,
      const unsigned short k ) ;
    // ========================================================================
    /** @struct Eulerian_
     *  Eulerian numbers A(n,k)
     *  @see https://en.wikipedia.org/wiki/Eulerian_number     
     */
    template <unsigned short N,
              unsigned short K>
    struct Eulerian_
    {
      // enum : unsigned long long 
    };
    // ========================================================================
    /** Eulerian number A(n,k)
     *  @see https://en.wikipedia.org/wiki/Eulerian_number
     *  @param n   \f$ 0 \le n \f$
     *  @param k   \f$ 0 \le k \le n \f$ 
     *  @return euleria number A(n,k)  
     */
    std::uintmax_t 
    eulerian 
    ( const unsigned short n , 
      const unsigned short k ) ;
    // ========================================================================
    /** Eulerian number A(n,k)
     *  @see https://en.wikipedia.org/wiki/Eulerian_number
     *  @param n   \f$ 0 \le n \f$
     *  @param k   \f$ 0 \le k \le n \f$ 
     *  @return euleria number A(n,k)  
     */
     double 
     eulerian_double 
     ( const unsigned short n , 
       const unsigned short k ) ; 
    // ========================================================================
    /** Get a row of Eulerian numbers for given N
     *  @see https://en.wikipedia.org/wiki/Eulerian_number
     *  @param n   \f$ 0 \le n \f$
     *  @param k   \f$ 0 \le k \le n \f$ 
     *  @return euleria number A(n,k)  
     */

    
    // const std::vector<double>&
    // eulerian
    // ( const unsigned short n ) ;

    
    // ========================================================================    
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CHOOSE_H
// ============================================================================
