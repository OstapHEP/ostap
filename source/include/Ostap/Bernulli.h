// ============================================================================
#ifndef OSTAP_BERNULLI_H 
#define OSTAP_BERNULLI_H 1
// ============================================================================
// STD & STL
// ============================================================================
#include <ratio>
#include <complex>
// ============================================================================
/** @file Ostap/Bernulli.h
 *  Bernulli numbers and polynomials 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @var N_BERNULLI_MAX
     * number  of exact Bernully numbers
     */
    constexpr unsigned short N_BERNULLI_MAX = 40 ;
    // ========================================================================
    /// Bernulli's numbers 
    template <unsigned short N>
    struct Bernulli_ ;    
    // ========================================================================
    template <>
    struct Bernulli_<0>
    { static constexpr std::ratio<+1,1>        ratio {}  ;  } ;       
    // ========================================================================
    template <>
    struct Bernulli_<1>
    { static constexpr std::ratio<-1,2>        ratio {}  ; } ;
    //  
    template <>
    struct Bernulli_<2>
    { static constexpr std::ratio<+1,6>        ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<4>
    { static constexpr std::ratio<-1,30>       ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<6>
    { static constexpr std::ratio<+1,42>       ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<8>
    { static constexpr std::ratio<-1,30>       ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<10>
    { static constexpr std::ratio<+5,66>       ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<12>
    { static constexpr std::ratio<-691,2730>   ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<14>
    { static constexpr std::ratio<+7,6>        ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<16>
    { static constexpr std::ratio<-3617,510>   ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<18>
    { static constexpr std::ratio<+43867,798>  ratio {}  ;  } ;
    //
    template <>
    struct Bernulli_<20>
    { static constexpr std::ratio<-174611,330> ratio {}  ;  } ;
    //
    template <>
    struct Bernulli_<22>
    { static constexpr std::ratio<+854513,138> ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<24>
    { static constexpr std::ratio<-236364091,2730> ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<26>
    { static constexpr std::ratio<+8553103,6> ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<28>
    { static constexpr std::ratio<-23749461029,870> ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<30>
    { static constexpr std::ratio<+8615841276005,14322> ratio {} ; } ;
    //
    template <>
    struct Bernulli_<32>
    { static constexpr std::ratio<-7709321041217,510> ratio {}  ; } ;
    //
    template <>
    struct Bernulli_<34>
    { static constexpr std::ratio<+2577687858367,6> ratio {} ; } ;
    //
    template <unsigned short N>
    struct Bernulli_ : public std::enable_if<(3<=N)&&(1== N%2)>
    {
      static constexpr std::ratio<+0,1>  ratio {}  ; 
      static constexpr const long double value { 0.0 } ;
    } ;
    // ========================================================================
    /// get sequence of Bernulli' numbers from 0 to N (inclusive)
    template <unsigned short N, std::size_t...i>
    constexpr std::array<long double,N+1> bernulli ( std::index_sequence<i...> )
    { return std::array<long double,N+1> {{ Bernulli_<i>::value ...}} ; } ;  
    /// get sequence of Bernulli' numbers fmor 0 to N (inclusive)    
    template <unsigned short N>
    constexpr std::array<long double,N+1> bernulli ()
    { return bernulli<N> ( std::make_index_sequence<N+1> {} ) ; }
    // ========================================================================
    /** Get Bernulli number 
     *  - \f$ B_0 = 1 \f$ 
     *  - \f$ B_1 = -\frac{1}{2} \f$ 
     *  - \f$ B_{2k+1) = 0 \f$ 
     */
    double bernulli ( const unsigned short k ) ;
    // =========================================================================
    /** @class Bernulli
     *  evaluate the Bernulli polynomials
     *  @see https://en.wikipedia.org/wiki/Bernoulli_polynomials
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2025-04-05
     */
    class Bernulli
    {
    public:
      // ======================================================================
      /// constructor from the order 
      Bernulli ( const unsigned short N = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the polynomial for real values 
      inline               double  operator()
      ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate the polynomial for complex values 
      inline  std::complex<double> operator()
      ( const std::complex<double>& z ) const { return evaluate ( z ) ; }
      /// evaluate Bernulli polynomial for real values 
      double evaluate
      ( const double x ) const ;
      /// evaluate Bernulli  polynomial for complex values 
      std::complex<double>
      evaluate 
      ( const std::complex<double>& z ) const ;
      // ======================================================================
    public:
      // ======================================================================
      inline unsigned short N      () const { return m_N ; }
      inline unsigned short degree () const { return m_N ; }
      inline double         xmin   () const { return 0   ; }
      inline double         xmax   () const { return 1   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// derivative at real point x 
      double derivative
      ( const double x ) const ;
      /// derivative at complex point z 
      std::complex<double> derivative 
      ( const std::complex<double>& z ) const ;
      /// integral
      double integral
      ( const double xmin ,
        const double xmax ) const ;
      // ======================================================================
    private:
      // ======================================================================
      unsigned short       m_N { 0 } ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_BERNULLI_H
// ============================================================================
