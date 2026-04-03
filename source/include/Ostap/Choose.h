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
     *  maximal value of n scuh that for all k the \f$ C^n_k \f$ fits into 64 bit integer
     */	
    inline constexpr unsigned short N_CHOOSE_MAX = 67 ;
    /// Use this type for Binomial coeffients 
    template <unsigned int N>
    using Choose_t = typename std::conditional<(N<=N_CHOOSE_MAX),std::uintmax_t,long double>::type ;    
    // ========================================================================
    template <unsigned short N,
	      unsigned short K>
    struct Choose_ ;
    // ========================================================================
    namespace details
    {      
      // ======================================================================
      template <unsigned short N ,
		unsigned short K ,
		bool           A ,  // N is not *too* large 
		bool           B >  //  K <= N 
      struct _Choose ;
      // ======================================================================
      template <unsigned short N,
		unsigned short K>
      struct _Choose<N,K,true,true>
      {
	typedef std::uintmax_t type ;
	enum _ : std::uintmax_t { value =
				  std::uintmax_t ( Ostap::Math::Choose_<N-1,K-1>::value ) +
				  std::uintmax_t ( Ostap::Math::Choose_<N-1,K>::value   ) } ; 
      } ;
      // ======================================================================
      template <unsigned short N,
		unsigned short K>
      struct _Choose<N,K,false,true>
      {
	typedef Ostap::Math::Choose_t<N> type ;
	inline static constexpr type value { 1.0L * type ( Ostap::Math::Choose_<N,K-1>::value ) * ( N + 1 - K ) / K } ; 	  
      } ;
    // ========================================================================
    }
    // ========================================================================
    /** @struct Choose_
     *  Binomial coefficient \f$ C^n_k = \left(\begin{array}{c}n \\ k \end{array}\rght) = \frac{n!}{k!(n-k)!}\f$
     */
    template <unsigned short N,
	      unsigned short K>
    struct Choose_ : public details::_Choose<N,K,N<=N_CHOOSE_MAX,(K<=N)>
    {} ; 
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
    // ========================================================================
    /** @class Pascal
     *  The row in Pascal triangle
     */  
    template <unsigned short N>
    struct Pascal_
    {
      /// the actual type fof elements in the row 
      // typedef typename std::conditional<(N<=N_CHOOSE_MAX),std::uintmax_t,long double>::type the_type ;
      typedef Choose_t<N> the_type ;
      // ======================================================================
      /// the actual type of the row in trianlge
      typedef std::array<the_type,N+1>  Row ;
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
    long double
    log_choose
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

    // ========================================================================
    // Stirling numbers 
    // ========================================================================    
    /** @class Stirling1
     *  Unsigned Stirling numbers of 1st kind 
     *  @see https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
     */
    template <unsigned short N, unsigned short K>
    struct Stirling1_ ;
    // ========================================================================
    /// stop recursion: \f$ S^0_0 = 1 \f$
    template <>
    struct Stirling1_<0,0> { enum _ : std::uintmax_t { value = 1 } ; } ;    
    // ========================================================================
    /// stop recursion \f$ S^n_0 = 0 \f$
    template <unsigned short N>
    struct Stirling1_<N,0> { enum _ : std::uintmax_t { value = 0 } ; } ;    
    // ========================================================================
    /// stop recursion \f$ S^0_s = 0 \f$
    template <unsigned short K>
    struct Stirling1_<0,K> { enum _ : std::uintmax_t { value = 0 } ; } ;
    /// start recursion  \f$ S^{n+1}_k = n S^n_k + S^n_{k-1}\f$ 
    template <unsigned short N,
	      unsigned short K>
    struct Stirling1_ 
    {
      enum _ : std::uintmax_t { value =  
				Stirling1_<N-1,K  >::value * ( N - 1 ) +
				Stirling1_<N-1,K-1>::value } ;
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
    { return std::array<TYPE,N+1> { { Stirling1_<N,N-i>::value... } } ; }
    template<class TYPE,unsigned short N>
    constexpr auto stirling1_array() 
    { return stirling1_array<TYPE,N> ( std::make_index_sequence<N+1> { } ) ; }
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
    /** @class Pochhammer
     *  Pochhammer symbol as polynomial
     *  @see Ostap::Math::rising_factorial
     *  @see Ostap::Math::pochhammer
     *  @see Ostap::Math::Pochhammer_
     */
    class Pochhammer
    {
    public :
      // =====================================================================
      /// constructor 
      Pochhammer ( const unsigned short N = 0 ) ;
      // =====================================================================
    public :
      // =====================================================================
      /// evaluate Pochhammer symbol/polynomial 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Pochhammer symbol/polynomial 
      double        evaluate   ( const double x ) const ;
      /// get the derivative 
      double        derivative ( const double x ) const ;
      // =====================================================================
    public :
      // =====================================================================
      /// degree of polynomial 
      inline unsigned short degree () const { return m_N ; }
      /// degree of polynomial 
      inline unsigned short N      () const { return m_N ; }
      // =====================================================================
    public :
      // =====================================================================
      /// vector of roots 
      inline std::vector<double> roots () const { return roots ( m_N ) ; }
      // =====================================================================
      /// vector of roots 
      static std::vector<double> roots ( const unsigned short N ) ;
      // =====================================================================
    private: 
      // =====================================================================
      /// degree order 
      unsigned short m_N { 0 } ; // degree order
      // =====================================================================
    } ;
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
    /** @class Stirling2
     *  Unsigned Stirling numbers of 2nd kind
     *  @see https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind
     */
    template <unsigned short N,
	      unsigned short K>
    struct Stirling2_ ;
    // ========================================================================
    /// stop recursion: \f$ S^0_0 = 1 \f$
    template <>
    struct Stirling2_<0,0> { enum _ { value = 1 } ; } ;    
    // ========================================================================
    /// stop recursion: \f$ S^n_n = 1 \f$
    template <unsigned short N>
    struct Stirling2_<N,N> { enum _ { value = 1 } ; } ;    
    // ========================================================================
    /// stop recursion \f$ S^n_0 = 0 \f$
    template <unsigned short N>
    struct Stirling2_<N,0> { enum _ { value = 0 } ; } ;    
    // ========================================================================
    /// start recursion  \f$ S^{n+1}_k = n S^n_k + S^n_{k-1}\f$
    template <unsigned short N,
	      unsigned short K>
    struct Stirling2_
    {
      enum _ : std::uintmax_t { value =  
				Stirling2_<N-1,K  >::value * K +
				Stirling2_<N-1,K-1>::value     } ;
    } ;
    // ========================================================================
    
    // ========================================================================
    /** @var N_EULERIAN_MAX
     *  maximal value of n such that for all k the \f$ A(n,k) \f$ fits into 64 bit integer
     */	
    inline constexpr unsigned short N_EULERIAN_MAX = 20 ;    
    // ========================================================================
    /// Use this type for Eulerian numbers 
    template <unsigned int N>
    using Eulerian_t = typename std::conditional<(N<=N_EULERIAN_MAX),std::uintmax_t,long double>::type ;    
    // ========================================================================
    /** @struct Eulerian_
     *  Eulerian numbers A(n,k)
     *  @see https://en.wikipedia.org/wiki/Eulerian_number
     *  - we need it for polylogarithms 
     */
    template <unsigned short N,
              unsigned short K>
    struct Eulerian_ ;
    // =========================================================================
    namespace details
    {      
      // ======================================================================
      template <unsigned short N ,
		unsigned short K ,
		bool           A ,  // N is not *too* large 
		bool           B >  //  K <= N 
      struct _Eulerian ;
      // ======================================================================
      template <unsigned short N,
		unsigned short K>
      struct _Eulerian<N,K,true,true>
      {
	typedef std::uintmax_t type ;
	enum _ : std::uintmax_t { value =
				  std::uintmax_t ( Ostap::Math::Eulerian_<N-1,K-1>::value ) * ( N - K ) +
				  std::uintmax_t ( Ostap::Math::Eulerian_<N-1,K>::value   ) * ( K + 1 ) } ; 
      } ;
      // ======================================================================
      template <unsigned short N,
		unsigned short K>
      struct _Eulerian<N,K,false,true>
      {
	typedef Ostap::Math::Eulerian_t<N> type ;
	inline static constexpr type value { type ( Ostap::Math::Eulerian_<N-1,K-1>::value ) * ( N - K ) +
					     type ( Ostap::Math::Eulerian_<N-1,K>::value   ) * ( K + 1 ) } ;
      } ;
      // ========================================================================
    }
    // =========================================================================    
    template <unsigned short N,
	      unsigned short K>
    struct Eulerian_ : public details::_Eulerian<N,K,N<=N_EULERIAN_MAX,(K<=N)>
    {} ; 
    template <>
    struct Eulerian_<0,0>   { enum _ { value = 1 } ; } ;    
    template <unsigned short N>
    struct Eulerian_<N,0>   { enum _ { value = 1 } ; } ;
    // 
    // ========================================================================
    template <unsigned short N>
    struct Eulerian_<N,N>   { enum _ { value = 0 } ; } ;
    // ========================================================================
    /** @class EulerianRaw
     *  The row in Euler Triangle triangle
     */  
    template <unsigned short N>
    struct EulerianRow_
    {
      // ======================================================================
      /// the actual type fof elements in the row 
      typedef Eulerian_t<N> the_type ;
      // ======================================================================
      /// the actual type of the row in trianlge
      typedef std::array<the_type,N>  Row ;
      // ======================================================================
      /// get the row of coefficients 
      template <std::size_t... i>
      inline static constexpr Row eulerian_array ( std::index_sequence<i...> )
      { return Row { { Eulerian_<N,i>::value... } } ; }
      // ======================================================================
      /// get the row of coefficients 
      inline static constexpr Row eulerian_array ()
      { return EulerianRow_<N>::eulerian_array ( std::make_index_sequence<N> {} ) ; }
      // ======================================================================
      /// get the row of coefficients 
      static constexpr const Row row { eulerian_array () } ;
      // ======================================================================
    } ;
    // ========================================================================
    /** get the row in Euler's triangle - the row of Eulerian number coefficients 
     *  @see EulerianRow_
     */
    template<unsigned short N>
    inline constexpr const auto& eulerian_array() { return Ostap::Math::EulerianRow_<N>::row ; }
    // ========================================================================    
    /** Eulerian number A(n,k)
     *  @see https://en.wikipedia.org/wiki/Eulerian_number
     *  @param n   \f$ 0 \le n \f$
     *  @param k   \f$ 0 \le k \le n \f$ 
     *  @return Eulerian number A(n,k)  
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
     *  @return Eulerian number A(n,k)  
     */
     double 
     eulerian_double 
     ( const unsigned short n , 
       const unsigned short k ) ; 
    // ========================================================================
    /** Eulerian polynomials
     *  \f$  A_n(t) = \sum_k  A(k,k) t^k \f$ 
     *  @see https://en.wikipedia.org/wiki/Eulerian_number
     */
    class Eulerian
    {
      // ======================================================================
    public :
      // ======================================================================
      /// constructor
      Eulerian ( const unsigned short N ) ;      
      // ======================================================================
    public: 
      // ======================================================================
      /// evaluate Eulerian polynomial 
      inline double operator () ( const double x ) const {  return evaluate ( x ) ; } 
      /// evaluate Eulerian polynomial 
      double evaluate   ( const double x ) const ;
      /// get the derivative 
      double derivative ( const double x ) const ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get N : index 
      inline unsigned short N      () const { return m_N ; }
      inline unsigned short index  () const { return m_N ; }
      /// ATTENTION the actual degree is N-1 
      inline unsigned short degree () const     // ATTENTION! 
      { return !m_N ? 0 : ( m_N - 1 ) ; }       // ATTENTION!
      // ======================================================================
    private: 
      // ======================================================================
      /// polynomial index, the actual degree is oneunit less 
      unsigned short m_N { 0 } ;              // polynomial index 
      // ======================================================================      
    } ;
    // ========================================================================            
    /** Get a row of Eulerian numbers for given N
     *  @see https://en.wikipedia.org/wiki/Eulerian_number
     *  @param n   \f$ 0 \le n \f$
     *  @param k   \f$ 0 \le k \le n \f$ 
     *  @return eulerian number A(n,k)  
     */
    const std::vector<double>&
    eulerian
    ( const unsigned short n ) ;    
    // ========================================================================   
    /** @struct Stieltjes_
     * Stiltjes's constants
    *  @see https://en.wikipedia.org/wiki/Stieltjes_constants
     */
    template  <unsigned short N>
    struct Stieltjes_ ;
    // ========================================================================
    template  <>
    struct Stieltjes_<0>
    {
	     inline static constexpr long double value {
        +0.5772156649015328606065120900824024310421593359L } ; 
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<1>
    {
	     inline static constexpr long double value {
       -0.0728158454836767248605863758749013191377363383L } ;
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<2>
    {
	     inline static constexpr long double value {
        -0.0096903631928723184845303860352125293590658061L } ;
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<3>
    {
	     inline static constexpr long double value {
        +0.0020538344203033458661600465427533842857158044L } ; 
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<4>
    {
	     inline static constexpr long double value {
        +0.0023253700654673000574681701775260680009044694L } ; 
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<5>
    {
	     inline static constexpr long double value {
        +0.0007933238173010627017533348774444448307315394L } ;
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<6>
    {
	     inline static constexpr long double value {
        -0.0002387693454301996098724218419080042777837151L } ; 
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<7>
    {
	     inline static constexpr long double value {
        -0.0005272895670577510460740975054788582819962534L } ; 
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<8>
    {
	     inline static constexpr long double value {
        -0.0003521233538030395096020521650012087417291805L} ; 
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<9>
    {
	     inline static constexpr long double value {
        -0.0000343947744180880481779146237982273906207895L } ;
    }; 
    // ========================================================================
    template  <>
    struct Stieltjes_<10>
    {
	     inline static constexpr long double value {
        +0.0002053328149090647946837222892370653029598537L } ;
    }; 
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CHOOSE_H
// ============================================================================
