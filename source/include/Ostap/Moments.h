// ============================================================================
#ifndef OSTAP_MOMENTS_H 
#define OSTAP_MOMENTS_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <array>
#include <type_traits>
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Choose.h"
#include "Ostap/ValueWithError.h"
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Statistic
     *  Helper abstract base class for statistic counters 
     */
    class Statistic
    {
    public :
      // ======================================================================
      virtual ~Statistic() ;
      /// add new value to the counter 
      virtual void update ( const double x ) = 0 ;
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Moment
     *  Helper (empty) base class for weighted statistic counters 
     *  - it is not really neeeded for C++, but it simplifies python decorations 
     */
    class Moment : public Statistic
    {
    public :
      // ======================================================================
      virtual ~Moment() ;
      // ======================================================================
    } ;
    // ========================================================================
    /// forward declaration 
    template <unsigned short N> class Moment_   ;
    /// template spccialization for N=0
    template <>                 class Moment_<0>;
    /// template spccialization for N=1
    template <>                 class Moment_<1>;
    // ========================================================================
    /** @class Moment_  Ostap/Moments.h 
     *  Simple class to keep/calculate 
     *  the high-order central momentts
     *  \f[  \mu_n \equiv \frac{1}{N} \sum_{i}  \left( x_i - \bar{x} \right)^n \f] 
     *  It implements the (single-pass) algorithm 
     *  @see  Pebay, P., Terriberry, T.B., Kolla, H. et al. 
     *        "Numerically stable, scalable formulas for parallel and online 
     *        computation of higher-order multivariate central moments with 
     *        arbitrary weights". Comput Stat 31, 1305–1325 (2016). 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    template <unsigned short N>
    class Moment_ : public Moment 
    {
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = N } ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor 
      Moment_ () = default ;
      /// constructor from the value and previous moment 
      Moment_
      ( const double        mom  ,
	const Moment_<N-1>& prev )
	: m_prev ( prev              )
	, m_M    ( mom * prev.size() )
      {} 
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of the Nth moment 
       *  \f[ \mu_N \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      inline double moment () const { return this->ok () ? this->M ( N ) / this->size () : 0 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const
      { return 
          N <  k       ? 0 :
          0 == k       ? 1 :
          1 == k       ? 0 :
          !this->ok () ? 0 : ( this->M ( k ) / this->size() ) ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { 
        return 
          N <  k    ? 0 : 
          0 == k    ? 1 : 
          1 == k    ? 0 :
          2 == k    ? 1 :
          !this->ok()  ? 0 : this->moment ( k ) / std::pow ( this->moment ( 2 ) , 0.5 * k ) ;
      }
      // ======================================================================
      /// get number of entries
      inline unsigned long long size  () const { return m_prev.size  () ; }
      /// get effective number of entries 
      inline unsigned long long nEff  () const { return m_prev.nEff  () ; }
      /// get the mean value (if \f$ 1 \le N \$4)
      inline long double        mu    () const { return m_prev.mu    () ; }
      /// empty ?
      inline bool               empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool               ok    () const { return m_prev.ok    () ; }
      // ======================================================================
    public: // basic operations for the counter 
      // ======================================================================
      /// increment with some value 
      inline Moment_& operator+= ( const double   x ) { return add ( x ) ; } 
      /// increment with other moment 
      inline Moment_& operator+= ( const Moment_& x ) { return add ( x ) ; }
      // ======================================================================
    public: // add more values to the counter 
      // ======================================================================
      /// add single value 
      inline Moment_& add ( const double   x ) ;
      /// add another moment 
      inline Moment_& add ( const Moment_& x ) ;
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public: // python operations 
      // ======================================================================
      Moment_   __add__ ( const double   x ) const { Moment_ r (*this) ; r += x ; return r ; }
      Moment_  __radd__ ( const double   x ) const { return __add__ ( x ) ; }
      Moment_& __iadd__ ( const double   x )       { return   add   ( x ) ; }
      // ======================================================================
      Moment_   __add__ ( const Moment_& x ) const { Moment_ r (*this) ; r += x ; return r ; }
      Moment_  __radd__ ( const Moment_& x ) const { return __add__ ( x ) ; }
      Moment_& __iadd__ ( const Moment_& x )       { return   add   ( x ) ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// update the counter 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// get "previos" moment 
      inline const Moment_<N-1>& previous () const { return this->m_prev ; }
      // ======================================================================
    public:      
      // ======================================================================
      /** get the value of \f$ M_k = \sum_i \left( x_i - \bar{x} \right)^k \f$ 
       *  for \f$ k \le N \f$
       *  @param k the order   \f$  0 \le k \le N \f$
       *  @return the value of \$ M_k \f$  if \f$  0 \le k \le N \f$, 0 otherwise 
       */   
      inline long double M ( const unsigned short k = N ) const
      { return k >  N ? 0.0L : k == N ? this->m_M : this->m_prev.M( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( Moment_& right )
      {
        m_prev.swap ( right.m_prev ) ;
        std::swap ( m_M , right.m_M ) ;
      }
      // ======================================================================
    public:  // templated prev
      // ======================================================================
      template <unsigned  int K, typename std::enable_if<(K+1==N),int>::type = 0 >
      inline const Moment_<K>& prev_ () const { return this->m_prev ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(N>K+1),int>::type = 0 >
      inline const Moment_<K>& prev_ () const { return this->m_prev.template prev_<K> () ; }
      // =====================================================================
    public: //  templated M_ 
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==N),int>::type = 0 >
      inline long double M_ () const { return this->m_M ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(N>K),int>::type = 0 >
      inline long double M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moment
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0)&&(K<=N),int>::type = 0 >
      inline long double moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(N>=K),int>::type = 0 >
      inline long double moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(K==N),int>::type = 0 >
      inline long double moment_ () const
      { return !this->ok () ? 0 : this->m_M / this->size  () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(N<2*K)&&(N>K),int>::type = 0 >
      inline long double moment_ () const
      { return this->m_prev.template moment_<K> () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError
      moment_ () const
      {
        if  ( !this -> ok () ) { return 0 ; }
        const unsigned long long n = this->size() ;
        //
        const long double muo  = this->template M_<K>   () / n ;
        //
        if ( this -> size() < 2 * K ) { return muo ; } // error estimate is impossible 
        //
        const long double mu2o = this->template M_<2*K> () / n ;
        const long double muop = this->template M_<K+1> () / n ;
        const long double muom = this->template M_<K-1> () / n ;  
        const long double mu2  = this->template M_<2>   () / n ;  
        //
        long double cov2 = mu2o     ;
        cov2 -= 2 * K * muop * muom ;
        cov2 -=         muo  * muo  ;
        cov2 += K * K * mu2  * muom * muom ;
        cov2 /= n ;
        //
        return Ostap::Math::ValueWithError ( muo , cov2 ) ;
      }
      // =======================================================================
    public: // templated std_moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==2)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<K)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const
      {
        return
	  !this->ok () ? 0 :
          this->template moment_<K> () / std::pow ( this->moment_<2>() , 0.5 * K ) ;
      }      
      // ======================================================================
    public:
      // ======================================================================      
      /// 1st cumulant 
      template <unsigned int K, typename std::enable_if<(1==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const { return this->mu() ; }
      /// 2nd cumulant 
      template <unsigned int K, typename std::enable_if<(2==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { return !this->ok() ? 0.0 : this->template M_<K>() / this->size() ; }
      /// 3rd cumulant 
      template <unsigned int K, typename std::enable_if<(3==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { return !this->ok() ? 0.0 : this->template M_<K>() / this->size() ; }
      /// 4th cumulant 
      template <unsigned int K, typename std::enable_if<(4==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { 
        if ( !this->ok() ) { return 0.0 ; }
        //
        const long double m4 = this->template M_<4> () / this->size() ;
        const long double m2 = this->template M_<2> () / this->size() ;
        //
        return m4 - 3 * m2 * m2 ;
      }
      /// 5th cumulant 
      template <unsigned int K, typename std::enable_if<(5==K)&&(K<=N),int>::type = 0 >
        inline long double cumulant_ () const 
      { 
        if ( !this->ok() ) { return 0.0 ; }
        //
        const long double m5 = this->template M_<5> () / this->size() ;
        const long double m3 = this->template M_<3> () / this->size() ;
        const long double m2 = this->template M_<2> () / this->size() ;
        //
        return m5 - 10  *m3 * m2 ;
      }
      /// 6th cumulant 
      template <unsigned int K, typename std::enable_if<(6==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { 
        if ( !this->ok() ) { return 0.0 ; }
        //
        const long double m6 = this->template M_<6> () / this->size() ;       
        const long double m4 = this->template M_<4> () / this->size() ;
        const long double m3 = this->template M_<3> () / this->size() ;
        const long double m2 = this->template M_<2> () / this->size() ;
        //
        return m6 - 15 * m4 * m2 - 10 * m3 * m3 + 30 * m2 * m2 * m2 ;
      }
      /// 7th cumulant 
      template <unsigned int K, typename std::enable_if<(7==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { 
        if ( !this->ok() ) { return 0.0 ; }
        //
        const long double m7 = this->template M_<7> () / this->size() ;       
        const long double m5 = this->template M_<5> () / this->size() ;       
        const long double m4 = this->template M_<4> () / this->size() ;
        const long double m3 = this->template M_<3> () / this->size() ;
        const long double m2 = this->template M_<2> () / this->size() ;
        //
        return m7 - 21 * m5 * m2 - 35 * m4 * m3 + 210 *m3 * m2 * m2 ;
      }
      /// 8th cumulant 
      template <unsigned int K, typename std::enable_if<(8==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { 
        if ( !this->ok() ) { return 0.0 ; }
        //
        const long double m8 = this->template M_<8> () / this->size() ;       
        const long double m6 = this->template M_<6> () / this->size() ;       
        const long double m5 = this->template M_<5> () / this->size() ;       
        const long double m4 = this->template M_<4> () / this->size() ;
        const long double m3 = this->template M_<3> () / this->size() ;
        const long double m2 = this->template M_<2> () / this->size() ;
        //
        return m8 - 28 * m6 * m2 - 56 * m5 * m3 - 35 * m4 * m4 
        + 420 * m4 * m2 * m2 + 560 * m3 * m3 * m2 - 630 * m2 * m2 * m2 * m2 ;
        
      }
      /// 9th cumulant 
      template <unsigned int K, typename std::enable_if<(9==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { 
        if ( !this->ok() ) { return 0.0 ; }
        //
        const long double m9 = this->template M_<9> () / this->size() ;       
        const long double m7 = this->template M_<7> () / this->size() ;       
        const long double m6 = this->template M_<6> () / this->size() ;       
        const long double m5 = this->template M_<5> () / this->size() ;       
        const long double m4 = this->template M_<4> () / this->size() ;
        const long double m3 = this->template M_<3> () / this->size() ;
        const long double m2 = this->template M_<2> () / this->size() ;
        //
        return m9 - 36 * m7 * m2 - 84 * m6 * m3 - 126 * m5 * m4 
        + 756 * m5 * m2 * m2 + 2520 * m4 * m3 * m2 
        + 560 * m3 * m3 * m3 - 7560 * m3 * m2 * m2 * m2 ; // NB? typo here? 
        
      }
      /// 10th cumulant 
      template <unsigned int K, typename std::enable_if<(10==K)&&(K<=N),int>::type = 0 >
      inline long double cumulant_ () const 
      { 
        if ( !this->ok() ) { return 0.0 ; }
        //
        const long double m10 = this->template M_<10> () / this->size() ;       
        const long double m8  = this->template M_<8>  () / this->size() ;       
        const long double m7  = this->template M_<7>  () / this->size() ;       
        const long double m6  = this->template M_<6>  () / this->size() ;       
        const long double m5  = this->template M_<5>  () / this->size() ;       
        const long double m4  = this->template M_<4>  () / this->size() ;
        const long double m3  = this->template M_<3>  () / this->size() ;
        const long double m2  = this->template M_<2>  () / this->size() ;
        //
        return m10 - 45 * m8 * m2 - 120 * m7 * m3 - 210 * m6 * m4 
        + 1260  * m6 * m2 * m2       - 126   * m5 * m5
        + 5040  * m5 * m3 * m2       + 3150  * m4 * m4 * m2 
        + 4200  * m4 * m3 * m3       - 18900 * m4 * m2 * m2 * m2  
        - 37800 * m3 * m3 * m2 * m2  + 22680 * m2 * m2 * m2 * m2 * m2 ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// counter of (N-1)th order
      Moment_<N-1> m_prev {}   ;  // counter of (N-1)th order
      /// the current value of \f$ N_N \f$  
      long double   m_M   {0}  ; // the current value
      // ======================================================================
    private:
      // ======================================================================
      /// helper array of binomial coefficients 
      static const std::array<unsigned long long,N+1> s_Ck ;
      // ======================================================================
    } ;
    // ========================================================================
    /// initialize the helper array of binomial coefficients
    template <unsigned short N>
    constexpr std::array<unsigned long long,N+1> 
    Moment_<N>::s_Ck { Ostap::Math::choose_array<N>() } ;
    // ========================================================================
    /// increment with some value
    template <unsigned short N>
    inline Moment_<N>& Moment_<N>::add ( const double  x )
    {
      //
      const auto nA = this->size() ;
      const unsigned long long nB =      1 ;
      const unsigned long long nN = nA + 1 ;
      //
      const long double delta = x - this -> mu () ;
      const long double b_n   =      -1.0L / nN ;
      const long double a_n   =  nA * 1.0L / nN ;
      const long double d_n   = - delta    / nN ;
      //
      m_M += ( nA * std::pow ( b_n , N ) + std::pow ( a_n , N ) ) * std::pow ( delta , N ) ;
      long double d = 1 ;
      for ( unsigned int k = 1 ; k + 2 <= N ; ++k )
      {
        d   *= d_n ;
        m_M += s_Ck [ k ] * this-> M ( N - k ) * d   ;
      }
      /// update previous 
      this->m_prev += x ; // update previous
      //
      return *this ;
    }
    // ========================================================================
    /// increment with some other counter 
    template <unsigned short N>
    inline Moment_<N>& Moment_<N>::add ( const Moment_<N>&  x )
    {
      //
      if      ( x     . empty () ) {               return *this ; }
      else if ( this -> empty () ) { (*this) = x ; return *this ; }
      //
      const unsigned long long nA    = this->size () ;
      const unsigned long long nB    =    x. size () ;
      const unsigned long long nN    = nA + nB       ;
      const long double        delta = x.mu() - this->mu() ;
      const long double        b_n   =  ( -1.0L * nB ) / nN ;
      const long double        a_n   =  (  1.0L * nA ) / nN ;
      //
      m_M += x.m_M ;
      m_M += nA * std::pow ( b_n * delta , N ) + nB * std::pow ( a_n * delta , N ) ;
      //
      long double a = 1 ;
      long double b = 1 ;
      long double d = 1 ;
      //
      for ( unsigned short k = 1 ; k + 2 <= N ; ++k )
      {
        a   *= a_n   ;
        b   *= b_n   ;
        d   *= delta ;        
        m_M += s_Ck [ k ] * d * ( this-> M( N -k ) * b + x. M ( N - k ) * a ) ;
      }
      /// update previous 
      this->m_prev += x.m_prev ; // update previous
      //
      return *this ;
    }
    // =======================================================================
    /// specialization for \f$N=0f\$
    template <>
    class Moment_<0> : public Moment 
    {
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 0 } ;
      // ======================================================================
    public:
      // ======================================================================
      /// (default) constructor 
      Moment_
      ( const unsigned long long size = 0 )
	: m_size ( size )
      {}
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of the Nth moment 
       *  \f[ \mu_N \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      inline double moment () const { return 1 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const { return 0 < k ? 0 : 1 ; }
      // ========================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 0 == k || 2 == k ? 1 : 0 ; }
      // ======================================================================
      /// get number of entries
      inline unsigned long long size  () const { return m_size  ; }
      /// get effective number of entries 
      inline unsigned long long nEff  () const { return size () ; }
      /// empty ?
      inline bool               empty () const { return 0 == m_size ; }
      /// ok ?
      inline bool               ok    () const { return      m_size ; }
      // ======================================================================
    public: // basic operations 
      // ======================================================================
      /// increment with some value 
      inline Moment_& operator+= ( const double    x ) { add ( x ) ; return *this ; }
      /// increment with other moment 
      inline Moment_& operator+= ( const Moment_&  x ) { add ( x ) ; return *this ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add single value
      inline void add ( const double   /* x */ ) { ++m_size           ; }
      /// add single value
      inline void add ( const Moment_&    x    ) { m_size += x.m_size ; }
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { m_size += std::distance  ( begin , end ) ; }
      // ======================================================================
    public: // python operations 
      // ======================================================================
      Moment_   __add__ ( const double   x ) const { Moment_ r (*this) ; r += x ; return r ; }
      Moment_  __radd__ ( const double   x ) const { return __add__ ( x ) ; }
      Moment_& __iadd__ ( const double   x )       { (*this) += x ; return *this ; }
      // ======================================================================
      Moment_   __add__ ( const Moment_& x ) const { Moment_ r (*this) ; r += x ; return r ; }
      Moment_  __radd__ ( const Moment_& x ) const { return __add__ ( x ) ; }
      Moment_& __iadd__ ( const Moment_& x )       { (*this) += x ; return *this ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// update the counter 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of \f$ M_k = \sum_i \left( x_i - \bar{x} \right)^k \f$ 
       *  for \f$ k \le N \f$
       *  @param k the order   \f$  0 \le k \le N \f$
       *  @return the value of \$ M_k \f$  if \f$  0 \le k \le  \f$, 0 otherwise 
       */   
      inline long double M ( const unsigned short k = 0 ) const
      { return  0 < k ? 0 : m_size  ; }
      // ======================================================================
    public:
      // ======================================================================      
      void swap ( Moment_& right )
      { std::swap ( m_size , right.m_size ) ; }
      // ======================================================================
    public: //  templated M_ 
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double M_ () const { return this->m_size ; }
      // ======================================================================
    public:  // templated moment
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double moment_ () const { return 1 ; }
      // ======================================================================
    public: // templated std_moment 
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================      
    private:
      // ======================================================================
      unsigned long long m_size { 0 } ;
      // ======================================================================
    private:
      // ======================================================================
      /// helper array of binomial coefficients 
      static const std::array<unsigned long long,1> s_Ck ;
      // ======================================================================
    } ;
    // ========================================================================
    /// specializaton for \f$ N=1 \f$ 
    template <>
    class Moment_<1> : public Moment 
    {
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 1 } ;
      // ======================================================================
    public: 
      // ======================================================================
      /// default constructor from mu and size 
      Moment_
      ( const double             mu   = 0 ,
	const unsigned long long size = 0 ) 
	: m_prev ( size ) 
	, m_mu   ( mu   )
      {}
      /// constructor from mu and previous 
      Moment_
      ( const double      mu   ,
	const Moment_<0>& prev ) 
	: m_prev ( prev ) 
	, m_mu   ( mu   )
      {}
      // ======================================================================
    public: 
      // ======================================================================
      /** get the value of the Nth moment 
       *  \f[ \mu_N \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      inline double moment () const { return 0 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const { return 0 ; }
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 0 == k || 2 == k ? 1 : 0 ; }
      // ======================================================================
      /// get number of entries
      inline unsigned long long size  () const { return m_prev.size ()  ; }
      /// get effective number of entries 
      inline unsigned long long nEff  () const { return m_prev.nEff ()  ; }
      // get the mean value
      inline long double        mu    () const { return m_mu            ; } 
      /// empty ?
      inline bool               empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool               ok    () const { return m_prev.ok    () ; }
      // ======================================================================
    public: // basic operations 
      // ======================================================================
      /// increment with some value 
      inline Moment_& operator+= ( const double  x ) { return add ( x ) ; }
      /// increment with other moment 
      inline Moment_& operator+= ( const Moment_ x ) { return add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add single value 
      inline Moment_& add ( const double x )
      {
        const auto n = m_prev.size() ;
        m_mu = ( n * m_mu + x ) /  ( n + 1 );  // calculate new mean value 
        m_prev += x ;                          // updated previous  
        return *this ;
      }
      /// add the moment 
      inline Moment_& add ( const Moment_& x )
      {
        if      ( x     . empty () ) {               return *this ; }
        else if ( this -> empty () ) { (*this) = x ; return *this ; }
        //
        const unsigned long long n1 =   m_prev.size() ;
        const unsigned long long n2 = x.m_prev.size() ;
        //
        m_mu = ( n1 * m_mu + n2 * x.m_mu ) / ( n1 + n2 ) ; // update mean 
        m_prev += x.m_prev ;                               // update previous
        //
        return *this ;
      }
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public: // python operations 
      // ======================================================================
      Moment_   __add__ ( const double   x ) const { Moment_ r (*this) ; r += x ; return r ; }
      Moment_  __radd__ ( const double   x ) const { return __add__ ( x ) ; }
      Moment_& __iadd__ ( const double   x )       { return   add   ( x ) ; }
      // ======================================================================
      Moment_   __add__ ( const Moment_& x ) const { Moment_ r (*this) ; r += x ; return r ; }
      Moment_  __radd__ ( const Moment_& x ) const { return __add__ ( x ) ; }
      Moment_& __iadd__ ( const Moment_& x )       { return   add   ( x ) ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// update the counter 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// get "previos" moment 
      inline const Moment_<0>& previous () const { return this->m_prev ; }
      // ======================================================================
    public: 
      // ======================================================================
      inline long double M ( const unsigned short k = 1 ) const
      { return 1 < k  ? 0 : 1 == k ? 0 : this->m_prev. M ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( Moment_& right )
      {
        m_prev.swap ( right.m_prev ) ;
        std::swap ( m_mu , right.m_mu ) ;
      }
      // ======================================================================
    public:  // templated prev
      // ======================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      const Moment_<K>& prev_ () const { return this->m_prev ; }
      // =====================================================================
    public: //  templated M_ 
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline long double M_ () const { return 0 ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moments
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline long double moment_ () const { return 0 ; }
      // =======================================================================
    public: // templated std_moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline long double std_moment_ () const { return 0 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// 1st cumulant 
      template <unsigned int K, typename std::enable_if<(1==K),int>::type = 0 >
      inline long double cumulant_ () const { return this->mu() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// mean value 
      inline double mean () const { return mu() ; }      
      // =======================================================================
    private :
      // ======================================================================
      Moment_<0>  m_prev {   } ;
      /// mean value
      long double m_mu   { 0 } ; // mean value
      // ======================================================================
    private:
      // ======================================================================
      /// array of binomial coefficients 
      static const std::array<unsigned long long,2> s_Ck ;
      // ======================================================================
    } ;
    // ========================================================================
    // Operations with counters 
    // =========================================================================
    /// add two counters        
    template <unsigned short N>
    inline Moment_<N> operator+ ( const Moment_<N>&  a , const Moment_<N>&  b  )
    { Moment_<N> r  ( a ) ; r += b ; return r ; }
    /// add counter and the value
    template <unsigned short N>
    inline Moment_<N> operator+ ( const Moment_<N>&  a , const double       b  )
    { Moment_<N> r  ( a ) ; r += b ; return r ; }
    /// add counter and the value
    template <unsigned short N>
    inline Moment_<N> operator+ ( const double       a , const Moment_<N>&  b  )
    { return b + a ; }
    // ========================================================================
    /// swap two counters 
    template <unsigned int N>
    inline void swap ( Moment_<N>&  a , Moment_<N>&  b  )
    { a.swap( b ) ; }
    // ========================================================================
    
    // ========================================================================
    // Weighted moments 
    // ========================================================================
    
	
    // ========================================================================
    /** @class WStatistic
     *  Helper (empty) base class for weighted statistics 
     */
    class WStatistic
    {
    public :
      // ======================================================================
      virtual ~WStatistic () ;
      /// add new value to the counter 
      virtual void update ( const double x , const double w = 1 ) = 0 ;
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class WMoment
     *  Helper (empty) base class for weighted moment-counters 
     *  - it is not really neeeded for C++, but it simplifies python decorations 
     */
    class WMoment : public WStatistic
    {
    public :
      // ======================================================================
      virtual ~WMoment() ;
      // ======================================================================      
    } ;
    // ========================================================================
    /// forward declaration 
    template <unsigned short N> class WMoment_    ;
    /// template spccialization for N=0
    template <>                 class WMoment_<0> ;
    /// template spccialization for N=1
    template <>                 class WMoment_<1> ;
    // ========================================================================
    /** @class WMoment_  Ostap/Moments.h 
     *  Simple class to keep/calculate 
     *  the summable  high-order weighted central momentts
     *  \f[  \mu_n \equiv \frac{1}{\sum w_i} \sum_{i}  w_i \left( x_i - \bar{x} \right)^n \f] 
     *  It implements the (single-pass) algorithm 
     *  @see  Pebay, P., Terriberry, T.B., Kolla, H. et al. 
     *        "Numerically stable, scalable formulas for parallel and online 
     *        computation of higher-order multivariate central moments with 
     *        arbitrary weights". Comput Stat 31, 1305–1325 (2016). 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    template <unsigned short N>
    class WMoment_ : public WMoment 
    {
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = N } ;
      // ======================================================================
    public: 
      // ======================================================================
      /// default constructor 
      WMoment_ () = default ;
      /// constructor from the value and previous moment 
      WMoment_
      ( const double         mom  ,
	const WMoment_<N-1>& prev )
	: m_prev ( prev            )
	, m_M    ( mom * prev.w () )
      {} 
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of the Nth weighted  moment 
       *  \f[ \mu_N \equiv  \frac{1}{\sum w_i} \sum w_i \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      inline double moment () const { return this-> ok () ? this->M ( N ) / this-> w () : 0.0 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{\sum w_i} \sum w_i \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const
      { return
	  N <  k       ? 0  :
	  0 == k       ? 1  :
	  1 == k       ? 0  :
	  !this->ok () ? 0  : ( this->M ( k ) / this-> w () ) ; }
      // ======================================================================
      /** get value of the kth standartized moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double std_moment ( const unsigned short k ) const
      { return 
          N <  k       ? 0 : 
          0 == k       ? 1 : 
          1 == k       ? 0 :
          2 == k       ? 1 :
	  !this->ok () ? 0 : this->moment ( k ) / std::pow ( this->moment ( 2 ) , 0.5 * k ) ;
      }
      // ======================================================================
      /// get number of entries
      inline unsigned long long size  () const { return m_prev.size () ; }
      /// get effective number of entries \f$  \frac{(\sum w_i)^2}{\sum w_i^2} \f$
      inline long double        nEff  () const { return m_prev.nEff () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline long double        w     () const { return m_prev.w    () ; }
      /// get sum of weights squared 
      inline long double        w2    () const { return m_prev.w2   () ; }
      /// get the weighted mean value (if \f$ 1 \le N \$4)
      inline long double        mu    () const { return m_prev.mu   () ; } 
      /// empty ?
      inline bool               empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool               ok    () const { return m_prev.ok    () ; }
      // ======================================================================
    public: // basic operations for the counter 
      // ======================================================================
      /// increment with other moment 
      inline WMoment_& operator+= ( const WMoment_& x ) { return add  ( x ) ; }
      // ======================================================================
    public: // add more values to the counter 
      // ======================================================================
      /// add single value 
      inline WMoment_& add ( const double x , const double w = 1 ) ;
      /// add aother moment   
      inline WMoment_& add ( const WMoment_& x ) ;
      // ======================================================================
    public: // python operations 
      // ======================================================================
      WMoment_   __add__ ( const WMoment_& x ) const { WMoment_ r (*this) ; r += x ; return r ; }
      WMoment_  __radd__ ( const WMoment_& x ) const { return __add__ ( x ) ; }
      WMoment_& __iadd__ ( const WMoment_& x )       { return   add   ( x ) ; }
      // ======================================================================
   public:      
      // ======================================================================
      /// update the counter 
      void update ( const double x , const double w = 1 ) override { add ( x , w ) ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// get "previos" moment 
      const WMoment_<N-1>& previous () const { return this->m_prev ; }
      // ======================================================================
    public:      
      // ======================================================================
      /** get the value of \f$ M_k = \sum_i w_i \left( x_i - \bar{x} \right)^k \f$ 
       *  for \f$ k \le N \f$
       *  @param k the order   \f$  0 \le k \le N \f$
       *  @return the value of \$ M_k \f$  if \f$  0 \le k \le N \f$, 0 otherwise 
       */   
      inline long double M ( const unsigned short k = N ) const
      { return k >  N ? 0.0L : k == N ? this->m_M : this->m_prev.M( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( WMoment_& right )
      {
        m_prev.swap ( right.m_prev ) ;
        std::swap ( m_M , right.m_M ) ;
      }
      // ======================================================================
    public:  // templated prev
      // ======================================================================
      template <unsigned  int K, typename std::enable_if<(K+1==N),int>::type = 0 >
      inline const WMoment_<K>& prev_ () const { return this->m_prev ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(N>K+1),int>::type = 0 >
      inline const WMoment_<K>& prev_ () const { return this->m_prev.template prev_<K> () ; }
      // =====================================================================
    public: //  templated M_ 
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==N),int>::type = 0 >
      inline long double M_ () const { return this->m_M ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(N>K) ,int>::type = 0 >
      inline long double M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moment
      // ======================================================================      
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0)&&(N>=K),int>::type = 0 >
      inline long double moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(N>=K),int>::type = 0 >
      inline long double moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(K==N),int>::type = 0 >
      inline long double moment_ () const
      { return !this->ok() ?  0 : this->m_M / this->w () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(N<2*K)&&(N>K),int>::type = 0 >
      inline long double moment_ () const
      { return this->m_prev.template moment_<K> () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(2*K<=N),int>::type = 0 >
      inline
      Ostap::Math::ValueWithError
      moment_ () const
      {
        //
        if ( !this->ok() ) { return 0 ; }
        //
		const long double n = this->w() ; // ATENTION!
        //
        const long double muo  = this->template M_ <K>  () / n ;
        const long double mu2o = this->template M_<2*K> () / n ;
        const long double muop = this->template M_<K+1> () / n ;
        const long double muom = this->template M_<K-1> () / n ;  
        const long double mu2  = this->         M_<2>   () / n ;  
        //
        long double cov2 = mu2o     ;
        cov2 -= 2 * N * muop * muom ;
        cov2 -=         muo  * muo  ;
        cov2 += N * N * mu2  * muom * muom ;
        cov2 /= n ;
        //
        return Ostap::Math::ValueWithError ( muo , cov2 ) ;
      }
      // =======================================================================
    public: // templated std_moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==2)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<K)&&(N>=K),int>::type = 0 >
      inline long double std_moment_ () const
      {
        return !this->ok() ? 0 :
          this->template moment_<K> () / std::pow ( this->moment_<2>() , 0.5 * K ) ;
      }      
      // ======================================================================
    private:
      // ======================================================================
      /// counter of (N-1)th order
      WMoment_<N-1> m_prev {}   ;  // counter of (N-1)th order
      /// the current value of \f$ M_N \f$  
      long double   m_M   {0}  ; // the current value
      // ======================================================================
    private:
      // ======================================================================
      /// helper array of binomial coefficients 
      static const std::array<unsigned long long,N+1> s_Ck ;
      // ======================================================================
    } ;
    // ========================================================================
    /// initialize the helper array of binomial coefficients
    template <unsigned short N>
    constexpr std::array<unsigned long long,N+1> WMoment_<N>::s_Ck { Ostap::Math::choose_array<N>() } ;
    // ========================================================================
    /// increment with some value
    template <unsigned short N>
    inline WMoment_<N>& WMoment_<N>::add ( const double  x , const double w )
    {
      //
      if ( !w ) { return *this ; }  // ZERO weights are ignored 
      //
      const long double wA    = this->w () ;
      const long double wB    =       w    ;
      const long double wW    = wA + wB    ;
      const long double delta = x - this->mu() ;
      //
      const long double b_n =  -1.0L * wB / wW ;
      const long double a_n =   1.0L * wA / wW ;
      const long double d_n = - delta     / wW ;
      //
      m_M += ( wA * std::pow ( b_n , N ) + wB * std::pow ( a_n , N ) ) * std::pow ( delta , N ) ;
      long double d = 1 ;
      for ( unsigned int k = 1 ; k + 2 <= N ; ++k )
      {
        d   *= d_n ;
        m_M += s_Ck [ k ] * this-> M ( N - k ) * d   ;
      }
      /// update previous 
      this->m_prev.add ( x , w ) ; // update previous
      //
      return *this ;
    }
    // ========================================================================
    /// increment with some other counter 
    template <unsigned short N>
    inline WMoment_<N>& WMoment_<N>::add ( const WMoment_<N>&  x )
    {
      //
      if      ( x     . empty () ) {               return *this ; }
      else if ( this -> empty () ) { (*this) = x ; return *this ; }
      //
      const long double wB    =    x. w () ;
      if ( !wB ) { return *this ; }                    // RETURN! 
      //
      const long double wA    = this->w () ;
      const long double wW    = wA + wB    ;
      const long double delta = x.mu() - this->mu()  ;
      const long double b_n   =  ( -1.0L * wB ) / wW ;
      const long double a_n   =  (  1.0L * wA ) / wW ;
      //
      m_M += x.m_M ;
      m_M += wA * std::pow ( b_n * delta , N ) + wB * std::pow ( a_n * delta , N ) ;
      //
      long double a = 1 ;
      long double b = 1 ;
      long double d = 1 ;
      //
      for ( unsigned short k = 1 ; k + 2 <= N ; ++k )
      {
        a   *= a_n   ;
        b   *= b_n   ;
        d   *= delta ;        
        m_M += s_Ck [ k ] * d * ( this-> M( N - k ) * b + x. M ( N - k ) * a ) ;
      }
      /// update previous 
      this->m_prev += x.m_prev ; // update previous
      //
      return *this ;
    }
    // =======================================================================
    /// specialization for \f$ N=0 f\$
    template <>
    class WMoment_<0> : public WMoment 
    {
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 0 } ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor
      WMoment_ () = default ;
      /** full constructor
       *  @param size number of entries 
       *  @param sumw sum of weights 
       *  @param sumw sum of squared weights 
       */
      WMoment_
      ( const unsigned long long size ,
	const double             sumw  ,
	const double             sumw2 ) ;
      // ======================================================================
    public : 
      // ======================================================================
      /** get the value of the 0th moment 
       *  \f[ \mu_N \equiv  \frac{1}{\sum w_i } \sum w_i \left( x_i - \bar{x} \right)^0 \f]
       *  @return the value of the 0th central moment
       */
      inline double moment () const { return 1 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{\sum w_i} \sum w_i \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > 0 \f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const { return 0 < k ? 0 : 1 ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 0 == k || 2 == k ? 1 : 0 ; }
      // ======================================================================
      /// get number of entries
      inline unsigned long long size  () const { return m_size ; }
      /// get effective number of entries \f$  \frac{(\sum w_i)^2}{\sum w_i^2} \f$
      inline long double        nEff  () const { return m_w2 ? m_w * m_w / m_w2 : -1.0 ; }
      /// get sum of weights  \f$  \sum w_i \f$  
      inline long double        w     () const { return m_w    ; }
      /// get sum of weights squared  \f$  \sum w_i^2 \f$  
      inline long double        w2    () const { return m_w2   ; }
      /// empty ?
      inline bool               empty () const { return 0 == m_size ; }
      /// ok ?
      inline bool               ok    () const { return m_size && m_w && ( 0 < m_w2 ) ; }
      // ======================================================================
    public: // basic operations 
      // ======================================================================
      /// increment with other moment 
      inline WMoment_& operator+= ( const WMoment_&  x ) { return add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add a single value
      inline WMoment_& add
      ( const double /* x     */ ,
	const double    w = 1    )
      { 
        if ( !w ) { return *this ; }  // ZERO weights are ignored 
        //
        ++m_size ; 
        m_w  += w     ; 
        m_w2 += w * w ;
        //
        return *this ; 
      }
      /// add a single value
      inline WMoment_& add ( const WMoment_& x )
      {
        if ( !x.m_w ) { return *this ; }
        //
        m_size += x.m_size ;
        m_w    += x.m_w    ;
        m_w2   += x.m_w2   ;
        //
        return *this ;
      }
      // ======================================================================
    public: // python operations 
      // ======================================================================
      WMoment_   __add__ ( const WMoment_& x ) const { WMoment_ r (*this) ; r += x ; return r ; }
      WMoment_  __radd__ ( const WMoment_& x ) const { return __add__ ( x ) ; }
      WMoment_& __iadd__ ( const WMoment_& x )       { return   add   ( x ) ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// update the counter 
      void update ( const double x , const double w = 1 ) override { add ( x , w ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of \f$ M_k = \sum_i w_i \left( x_i - \bar{x} \right)^k \f$ 
       *  for \f$ k \le N \f$
       *  @param k the order   \f$  0 \le k \le N \f$
       *  @return the value of \$ M_k \f$  if \f$  0 \le k \le  \f$, 0 otherwise 
       */   
      inline long double M ( const unsigned short k = 0 ) const
      { return  0 < k ? 0 : m_w ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( WMoment_& right )
      {
        std::swap ( m_size , right.m_size ) ;
        std::swap ( m_w    , right.m_w    ) ;
        std::swap ( m_w2   , right.m_w2   ) ;
      }
      // ======================================================================
    public: //  templated M_ 
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double M_ () const { return this->m_w ; }
      // ======================================================================
    public:  // templated moment
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double moment_ () const { return 1 ; }
      // ======================================================================
    public: // templated std_moment 
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================
    private:
      // ======================================================================
      /// number of entries 
      unsigned long long m_size { 0 } ; // number of entries
      /// sum of weights \f$  \sum w_i \f$ 
      long double        m_w    { 0 } ; // sum of weights
      /// sum of weights squared \f$  \sum w_i^2 \f$ 
      long double        m_w2   { 0 } ; // sum of weights squared 
      // ======================================================================
    private:
      // ======================================================================
      /// helper array of binomial coefficients 
      static const std::array<unsigned long long,1> s_Ck ;
      // ======================================================================
    } ;
    // ========================================================================
    /// specializaton for \f$ N=1 \f$ 
    template <>
    class WMoment_<1> : public WMoment 
    {
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 1 } ;
      // ======================================================================
    public: 
      // ======================================================================
      /// default constructor 
      WMoment_ () = default ;
      /// constructor from mu and previous moment 
      WMoment_
      ( const double       mu   ,
	const WMoment_<0>& prev )
	: m_prev ( prev )
	, m_mu   ( mu   )
      {} 
      // ======================================================================
    public :
      // ======================================================================
      /** get the value of the 1st moment 
       *  \f[ \mu_N \equiv  \frac{1}{N} \sum w_i \left( x_i - \bar{x} \right)^1 \f]
       *  @return the value of the 1st central moment
       */
      inline double moment () const { return 0 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{\sum w_i } \sum w_i \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const { return 1 <= k ? 0 : 1 ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 0 == k || 2 == k ? 1 : 0 ; }
      // ======================================================================
      /// get number of entries
      inline unsigned long long size  () const { return m_prev.size ()  ; }
      /// get effective number of entries \f$  \frac{(\sum w_i)^2}{\sum w_i^2} \f$
      inline long double        nEff  () const { return m_prev.nEff () ; }
      /// get sum of weights \f$ \sum w_i \f$
      inline  long double       w     () const { return m_prev.w    () ; }
      /// get sum of weights squared 
      inline long double        w2    () const { return m_prev.w2   () ; }
      // get the mean value
      inline long double        mu    () const { return m_mu           ; } 
      /// empty ?
      inline bool               empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool               ok    () const { return m_prev.ok    () ; }
      // ======================================================================
    public: // basic operations 
      // ======================================================================
      /// increment with other moment 
      inline WMoment_& operator+= ( const WMoment_ x ) { return add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add single value 
      inline WMoment_& add ( const double x , const double w = 1 )
      {
        if ( !w ) { return *this ; }  // ZERO weights are ignored 
        //
        const long double wA = this -> w() ;
        const long double wB = w ;
        //
        m_mu = ( wA * m_mu + wB * x ) / ( wA + wB ) ; // update mean
        //
        this->m_prev.add ( x , w ) ;
        //
        return *this ;
      }
      /// add another moment 
      inline WMoment_& add ( const WMoment_& x )  
      {
        if      ( x     . empty () ) {               return *this ; }
        else if ( this -> empty () ) { (*this) = x ; return *this ; }
        //
        const long double wB = x.m_prev.w () ;
        if ( !wB ) { return *this ; }
        //
        const long double wA = m_prev.w   () ;
        m_mu = ( wA * m_mu + wB * x.m_mu ) / ( wA + wB ) ; // update mean
        //
        m_prev += x.m_prev ;                                 // update previous
        //
        return *this ;
      }
      // ====================================================================== 
   public: // python operations 
      // ======================================================================
      WMoment_   __add__ ( const WMoment_& x ) const { WMoment_ r (*this) ; r += x ; return r ; }
      WMoment_  __radd__ ( const WMoment_& x ) const { return __add__ ( x ) ; }
      WMoment_& __iadd__ ( const WMoment_& x )       { return   add   ( x ) ; }
      // ======================================================================
   public:      
      // ======================================================================
      /// update the counter 
      void update ( const double x , const double w = 1 ) override { add ( x , w ) ; }
      // ======================================================================
    public:      
      // ======================================================================
      /// get "previos" moment 
      const WMoment_<0>& previous () const { return this->m_prev ; }
      // ======================================================================
    public: 
      // ======================================================================
      inline long double M ( const unsigned short k = 1 ) const
      { return 1 < k  ? 0 : 1 == k ? 0 : this->m_prev. M ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( WMoment_& right )
      {
        m_prev.swap ( right.m_prev ) ;
        std::swap ( m_mu , right.m_mu ) ;
      }
      // ======================================================================
    public:  // templated prev
      // ======================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      const WMoment_<K>& prev_ () const { return this->m_prev ; }
      // =====================================================================
    public: //  templated M_ 
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline long double M_ () const { return 0 ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moments
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline long double moment_ () const { return 0 ; }
      // =======================================================================
    public: // templated std_moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline long double std_moment_ () const { return 0 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// mean value 
      inline double mean () const { return mu() ; }      
      // =======================================================================
    private:
      // ======================================================================
      WMoment_<0>  m_prev {   } ;
      /// mean value
      long double m_mu    { 0 } ; // mean value
      // ======================================================================
    private:
      // ======================================================================
      /// array of binomial coefficients 
      static const std::array<unsigned long long,2> s_Ck ;
      // ======================================================================
    } ;
    // ========================================================================
    /// add two counters        
    template <unsigned short N>
    inline WMoment_<N> operator+ ( const WMoment_<N>&  a , const WMoment_<N>&  b  )
    { WMoment_<N> r  ( a ) ; r += b ; return r ; }
    
    // ========================================================================
    /// swap two counters 
    template <unsigned int N>
    inline void swap ( WMoment_<N>&  a , WMoment_<N>&  b  )
    { a.swap( b ) ; }
    // ========================================================================
    
    // ========================================================================
    // Other type of counters 
    // ========================================================================

    // ========================================================================
    /** @class GeometricMean 
     *  Calculate the geometric mean 
     *  \f$ \left(x_1x_2...x_n\right)^{\frac{1}{n}} \f$
     *  @see https://en.wikipedia.org/wiki/Geometric_mean
     */
    class GeometricMean : public Statistic 
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef Moment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor 
      GeometricMean () = default ;
      /// constructor for the counters 
      GeometricMean ( const Counter& cnt ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the geometric mean 
      inline double mean  () const { return value () ; }
      inline double value () const { return std::pow ( 2 , m_log.mean() ) ; }   
	  // ======================================================================
    public:
      // ======================================================================
      inline GeometricMean& operator+=( const double         x ) { return add ( x ) ; }
      inline GeometricMean& operator+=( const GeometricMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of log2 values 
      const Counter& counter () const { return m_log ; }
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only positive entries 
      GeometricMean&        add ( const double         x ) ;
      /// sum of two counters 
      inline GeometricMean& add ( const GeometricMean& x )
      {
	m_log.add ( x.m_log ) ;
	return *this ;
      }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_log.size  () ; }
      /// empty ?
      inline bool               empty () const { return m_log.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_log.ok    () ; } 
      // ======================================================================      
    private:
      // ======================================================================
      /// get the counter of log2(x) 
      Counter m_log {} ;
      // ======================================================================
    };
    // ========================================================================
    /** @class HarmonicMean 
     *  Calcualet  the harmonic mean 
     *  \f$ \frac{n}{ \frac{1}{x_1} + ... + \frac{1}{x_n}} \f$
     *  @see https://en.wikipedia.org/wiki/Harmonic_mean
     */
    class HarmonicMean : public Statistic
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef Moment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// default construictor
      HarmonicMean () = default ;
      /// constructor for the counters 
      HarmonicMean ( const Counter& cnt ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the harmonic mean 
      inline double value () const { return 1. / m_inv.mean() ; }
      inline double mean  () const { return value () ; }
	  // ======================================================================
    public:
      // ======================================================================
      inline HarmonicMean& operator+=( const double        x ) { return add ( x ) ; }
      inline HarmonicMean& operator+=( const HarmonicMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only non-zero entries 
      HarmonicMean&        add ( const double        x ) ;
      /// add two counters togather 
      inline HarmonicMean& add ( const HarmonicMean& x )
      {
	m_inv.add ( x.m_inv ) ;
	return *this ;
      }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of 1/x values 
      const Counter counter () const { return m_inv ; }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_inv.size  () ; }
      /// empty ?
      inline bool               empty () const { return m_inv.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_inv.ok    () ; } 
      // ======================================================================      
    private:
      // ======================================================================
      /// get the counter of 1/x 
      Counter m_inv {} ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PowerMean 
     *  Calculate  the power mean 
     *  \f$ \left(\frac{1}{n}\sum x_i^p \right)^{\frac{1}{p}}\f$
     *  @see https://en.wikipedia.org/wiki/Power_mean
     */
    class PowerMean : public Statistic 
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef Moment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// defautl construictor
      PowerMean ( const double = 1 ) ;
      /// constructor for the counters 
      PowerMean ( const double p , const Counter& cnt ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the power mean 
      inline double value () const { return std::pow ( m_pow.mean() , 1 / m_p ) ; }
      inline double mean  () const { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline PowerMean& operator+=( const double     x ) { return add ( x ) ; }
      inline PowerMean& operator+=( const PowerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      const Counter counter () const { return m_pow ; }
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only positive entries 
      PowerMean& add ( const double     x ) ;
      /// add two counters togather if p is common 
      PowerMean& add ( const PowerMean& x ) ;
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_pow.size  () ; }
      /// empty ?
      inline bool               empty () const { return m_pow.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_pow.ok    () ; }
      /// power 
      inline double             p     () const { return m_p            ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the power
      double  m_p   {1} ;
      /// get the counter of x^p
      Counter m_pow {}  ;
      // ======================================================================
    };
    // ========================================================================
    /** @class LehmerMean 
     *  Calculate the Lehmer mean 
     *  \f$ \frac{ \sum_i x_i^p}{\sum_i x_i^{p-1}} \f$
     *  - \f$ p \rigtharrow - \infty\f$ : minimal value
     *  - \f$ p = 0 \f$   : harmonic mean 
     *  - \f$ p = 1/2 \f$ : geometric mean 
     *  - \f$ p = 1 \f$   : arithmetic mean 
     *  - \f$ p = 2 \f$   : contraharmonic mean 
     *  - \f$ p \rigtharrow + \infty\f$ : maximal value
     *  @see https://en.wikipedia.org/wiki/Lehmer_mean
     */
    class LehmerMean : public Statistic 
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef Moment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// defautl construictor
      LehmerMean
      ( const double = 1 ) ;
      /// constructor for the counters 
      LehmerMean
      ( const double   p    ,
	const Counter& cnt1 , 
	const Counter& cnt2 ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the power mean 
      inline double value () const { return m_lp.mean() / m_lpm1.mean () ; }
      inline double mean  () const { return value () ; }
	  // ======================================================================
    public:
      // ======================================================================
      inline LehmerMean& operator+=( const double      x ) { return add ( x ) ; }
      inline LehmerMean& operator+=( const LehmerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      const Counter counter1 () const { return m_lp   ; }
      /// the counter of x**(p-1)  values 
      const Counter counter2 () const { return m_lpm1 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only positive entries 
      LehmerMean& add ( const double      x ) ;
      /// add two counters togather if p is common 
      LehmerMean& add ( const LehmerMean& x ) ;
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_lp.size   () ; }
      /// empty ?
      inline bool               empty () const { return m_lp.empty  () ; } 
      /// ok ?
      inline bool               ok    () const { return m_lp.ok     () ; }
      /// power 
      inline double             p     () const { return m_p            ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the power
      double m_p          {1} ;
      /// get the counter of x^p
      Moment_<1> m_lp     {}  ;
      /// get the counter of x^(p-1)
      Moment_<1> m_lpm1   {}  ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class WGeometricMean 
     *  Calculate the weighted geometric mean 
     *  @see https://en.wikipedia.org/wiki/Geometric_mean
     */
    class WGeometricMean : public WStatistic
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef WMoment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor 
      WGeometricMean () = default ;
      /// constructor for the counters 
      WGeometricMean ( const Counter& cnt ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the geometric mean 
      inline double value () const { return std::pow ( 2 , m_log.mean() ) ; }
      inline double mean  () const { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline WGeometricMean& operator+=( const WGeometricMean& x ) { return add ( x ) ; }   
      inline WGeometricMean& operator*=( const WGeometricMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of log2 values 
      const Counter& counter () const { return m_log ; }
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only positive entries 
      WGeometricMean& add
      ( const double x      ,
	const double w  = 1 ) ;
      // sum of to counters 
      inline WGeometricMean& add ( const WGeometricMean& x )
      {
	m_log.add ( x.m_log ) ;
	return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a WMoment 
      void update
      ( const double x ,
	const double w ) override { add ( x , w ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_log.size  () ; }
      /// number of effective entries
      inline long double        nEff  () const { return m_log.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline long double        w     () const { return m_log.w     () ; }
      /// get sum of weights squared 
      inline long double        w2    () const { return m_log.w2    () ; }
      /// empty ?
      inline bool               empty () const { return m_log.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_log.ok    () ; } 
      // ======================================================================      
    private:
      // ======================================================================
      /// get the counter of log2(x) 
      Counter m_log {} ;
      // ======================================================================
    };
    // ========================================================================
    /** @class WHarmonicMean 
     *  Calcualet  the weighted harmonic mean 
     *  @see https://en.wikipedia.org/wiki/Harmonic_mean
     */
    class WHarmonicMean : public WStatistic 
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef WMoment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// default construictor
      WHarmonicMean () = default ;
      /// constructor for the counters 
      WHarmonicMean ( const Counter& cnt ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the harmonic mean 
      inline double value () const { return 1. / m_inv.mean() ; }
	  inline double mean  () const { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline WHarmonicMean& operator+=( const WHarmonicMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of 1/x values 
      const Counter counter () const { return m_inv ; }
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only non-zero entries 
      WHarmonicMean&
      add ( const double         x     ,
	    const double         w = 1 ) ;
      /// add two counters togather 
      inline WHarmonicMean& add ( const WHarmonicMean& x )
      {
	m_inv.add ( x.m_inv ) ;
	return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a WMoment 
      void update
      ( const double x ,
	const double w ) override { add ( x , w ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_inv.size  () ; }
      /// number of effective entries
      inline long double        nEff  () const { return m_inv.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline long double        w     () const { return m_inv.w     () ; }
      /// get sum of weights squared 
      inline long double        w2    () const { return m_inv.w2    () ; }
      /// empty ?
      inline bool               empty () const { return m_inv.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_inv.ok    () ; } 
      // ======================================================================      
    private:
      // ======================================================================
      /// get the counter of 1/x 
      Counter m_inv {} ;
      // ======================================================================
    };
    // ========================================================================
    /** @class WPowerMean 
     *  Calculate  the weighted power mean 
     *  @see https://en.wikipedia.org/wiki/Power_mean
     */
    class WPowerMean : public WStatistic 
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef WMoment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// defautl construictor
      WPowerMean ( const double = 1 ) ;
      /// constructor for the counters 
      WPowerMean ( const double p , const Counter& cnt ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the weighter power mean 
      inline double value () const { return std::pow ( m_pow.mean() , 1 / m_p ) ; }
      inline double mean  () const { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline WPowerMean& operator+=( const WPowerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      const Counter counter () const { return m_pow ; }
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only positive entries 
      WPowerMean& add
      ( const double  x     ,
	const double  w = 1 ) ;
      /// add two counters togather if p is common 
      WPowerMean& add ( const WPowerMean& x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update
      ( const double x ,
	const double w ) override { add ( x , w ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_pow.size  () ; }
      /// number of effective entries
      inline long double        nEff  () const { return m_pow.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline long double        w     () const { return m_pow.w     () ; }
      /// get sum of weights squared 
      inline long double        w2    () const { return m_pow.w2    () ; }
      /// empty ?
      inline bool               empty () const { return m_pow.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_pow.ok    () ; }
      /// power 
      inline double             p     () const { return m_p            ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the power
      double  m_p   {1} ;
      /// get the counter of x^p
      Counter m_pow {}  ;
      // ======================================================================
    };
    // ========================================================================
    /** @class WLehmerMean 
     *  Calculate the weighted Lehmer mean 
     *  @see https://en.wikipedia.org/wiki/Lehmer_mean
     */
    class WLehmerMean : public WStatistic 
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef WMoment_<1> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// defautl construictor
      WLehmerMean
      ( const double = 1 ) ;
      /// constructor for the counters 
      WLehmerMean
      ( const double   p    ,
	const Counter& cnt1 , 
	const Counter& cnt2 ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the Lehmer  mean 
      inline double value () const { return m_lp.mean () / m_lpm1.mean () ; }
      inline double mean  () const { return value () ; } 
      // ======================================================================
    public:
      // ======================================================================
      inline WLehmerMean& operator+=( const WLehmerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      const Counter counter1 () const { return m_lp   ; }
      /// the counter of x**(p-1)  values 
      const Counter counter2 () const { return m_lpm1 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only positive entries with non-zero weight       
      WLehmerMean& add
      ( const double x     ,
	const double w = 1 ) ;
      /// add two counters togather if p is common 
      WLehmerMean& add ( const WLehmerMean& x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a WMoment 
      void update
      ( const double x ,
	const double w ) override { add ( x , w ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_lp.size   () ; }
      /// number of effective entries
      inline long double        nEff  () const { return m_lp.nEff   () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline long double        w     () const { return m_lp.w      () ; }
      /// get sum of weights squared 
      inline long double        w2    () const { return m_lp.w2     () ; }
      /// empty ?
      inline bool               empty () const { return m_lp.empty  () ; } 
      /// ok ?
      inline bool               ok    () const { return m_lp.ok     () ; }
      /// power 
      inline double             p     () const { return m_p            ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the power
      double  m_p     {1} ;
      /// get the counter of x^p
      Counter m_lp    {}  ;
      /// get the counter of x^(p-1)
      Counter m_lpm1  {}  ;
      // ======================================================================
    };
    // ========================================================================
    /** @class ArithmeticMean 
     *  Calculate the arithmetic mean 
     */
    class ArithmeticMean : public  Moment_<1>
    {
    public:
      // ======================================================================
      typedef Moment_<1>  Counter ;
      // ======================================================================
    public:
      // ======================================================================
      ArithmeticMean() = default ;
      ArithmeticMean ( const Counter& cnt ) ; 
      // ======================================================================
    public:
      // ======================================================================
      inline double value() const { return this->mean () ; }
      // ======================================================================
    public:
      // ======================================================================
      Counter counter() const { return *this ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class WArithmeticMean 
     *  Calculate the weighted arithmetic mean 
     */
    class WArithmeticMean : public WMoment_<1>
    {
    public:
      // ======================================================================
      typedef WMoment_<1>  Counter ;
      // ======================================================================
    public:
      // ======================================================================
      WArithmeticMean() = default ;
      WArithmeticMean ( const Counter& cnt ) ; 
      // ======================================================================
    public:
      // ======================================================================
      inline double value() const { return this->mean () ; }
      // ======================================================================
    public:
      // ======================================================================
      Counter counter() const { return *this ; }
      // ======================================================================
    } ;       
    // ========================================================================
    /** @class MinMaxValue 
     */
    class MinMaxValue : public Statistic 
    {
      // ======================================================================
    public:
      // ======================================================================
      typedef Moment_<0> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor
      MinMaxValue () ;
      /// constructor from min/max & counter 
      MinMaxValue
      ( const double   min ,
	const double   max ,
	const Counter& cnt ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the min/max value  
      inline std::pair<double,double> value () const
      { return std::make_pair ( m_min , m_max ) ; }
      /// get the minvalue  
      inline double min  () const { return m_min ; }
      /// get the max value  
      inline double max  () const { return m_max ; }
      // ======================================================================
    public:
      // ======================================================================
      inline MinMaxValue& operator+=( const double       x ) { return add ( x ) ; }
      inline MinMaxValue& operator+=( const MinMaxValue& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only positive entries 
      inline MinMaxValue& add ( const double x )
      {
	m_min = std::min ( m_min , x ) ;
	m_max = std::max ( m_max , x ) ;
	m_cnt.add ( x ) ;
	return *this ;
      }
      inline MinMaxValue& add ( const MinMaxValue& x )
      {
	m_min = std::min ( m_min , x.m_min ) ;
	m_max = std::max ( m_max , x.m_max ) ;
	m_cnt.add ( x.m_cnt  ) ;
	return *this ;
      }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public:
      // ======================================================================
      const Counter& counter() const { return m_cnt ; }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_cnt.size  () ; }
      /// empty ?
      inline bool               empty () const { return m_cnt.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_cnt.ok    () ; } 
      // ======================================================================      
    private:
      // ======================================================================
      /// minimal value   ;
      double     m_min    ; 
      /// maximal value   ;
      double     m_max    ; 
      /// get the counter
      Counter    m_cnt {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class WMinMaxValue 
     */
    class WMinMaxValue : public WStatistic 
    {
      // ====================================================================== 
    public:
      // ======================================================================
      typedef WMoment_<0> Counter ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor
      WMinMaxValue () ;
      /// constructor from min/max & counter 
      WMinMaxValue
      ( const double   min ,
	const double   max ,
	const Counter& cnt ) ;
      // ======================================================================
   public:
      // ======================================================================
      /// get the min/max value  
      inline std::pair<double,double> value () const
      { return std::make_pair ( m_min , m_max ) ; }
      /// get the minvalue  
      inline double min  () const { return m_min ; }
      /// get the max value  
      inline double max  () const { return m_max ; }
      // ======================================================================
    public:
      // ======================================================================
      inline WMinMaxValue& operator+=( const WMinMaxValue& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// accumulate only no-zero weights  
      inline WMinMaxValue& add
      ( const double x     ,
	const double w = 1 )
      {
	const unsigned long long sbefore = m_cnt.size() ;
	m_cnt.add ( x , w ) ;
	if ( m_cnt.size() != sbefore  )
	  {
	    m_min = std::min ( m_min , x ) ;
	    m_max = std::max ( m_max , x ) ;
	  }
	return *this ;
      }
      inline WMinMaxValue& add ( const WMinMaxValue& x )
      {
	m_min = std::min ( m_min , x.m_min ) ;
	m_max = std::max ( m_max , x.m_max ) ;
	m_cnt.add ( x.m_cnt  ) ;
	return *this ;
      }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a WMoment 
      void update
      ( const double x      ,
	const double w  = 1 ) override { add ( x , w ) ; }
      // ======================================================================
    public:
      // ======================================================================
      const Counter& counter() const { return m_cnt ; }
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline unsigned long long size  () const { return m_cnt.size  () ; }
      /// empty ?
      inline bool               empty () const { return m_cnt.empty () ; } 
      /// ok ?
      inline bool               ok    () const { return m_cnt.ok    () ; }
      /// get effective number of entries \f$  \frac{(\sum w_i)^2}{\sum w_i^2} \f$
      inline long double        nEff  () const { return m_cnt.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline long double        w     () const { return m_cnt.w     () ; }
      /// get sum of weights squared 
      inline long double        w2    () const { return m_cnt.w2    () ; }
      // ======================================================================      
    private:
      // ======================================================================
      /// minimal value    ;
      double      m_min    ; 
      /// maximal value    ;
      double      m_max    ; 
      /// get the counter
      Counter     m_cnt {} ;
      // ======================================================================
    } ;
    // ========================================================================
    // More decorations 
    // ========================================================================
    /**  @class Moments 
     *   Collection of static functions dealing with moments 
     */
    class Moments
    {
      // ======================================================================
      typedef Ostap::Math::ValueWithError VE    ;
      // ======================================================================
    private: 
      // ======================================================================
      /** @var s_INVALID_MOMENT ; 
       *  the invalid value of the central moment 
       */ 
      static const double s_INVALID_MOMENT ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of invalid momment 
      static double invalid_moment () { return s_INVALID_MOMENT ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the unbiased estimator for the 2nd order moment:
       *  \f[ \hat{\mu}_2 \equiv \frac{n}{n-1} \mu_2 \f] 
       *  @param  m input counter
       *  @return the unbiased estimate (if number of entries exceeds 2) and 
       *   <code>s_INVALID_MOMENT</code> otherwise 
       */
      template <unsigned short N, typename std::enable_if<(N>=2),int>::type = 0 >
      static inline double unbiased_2nd ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size() ;
        return !m.ok() || n < 2 ? s_INVALID_MOMENT  : m.template M_<2> () / ( n - 1 ) ;  
      }
      // ======================================================================
      /** get the unbiased estimator for the 3rd order moment:
       *  \f[ \hat{\mu}_3 \equiv \frac{n^2}{(n-1)(n-2)} \mu_3 \f] 
       *  @param  m input counter
       *  @return the unbiased estimate (if number of entries exceeds 3) and 
       *   <code>s_INVALID_MOMENT</code> otherwise 
       */
      template <unsigned short N, typename std::enable_if<(N>=3),int>::type = 0 >
      static inline double unbiased_3rd ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size() ;
        return !m.ok() || n < 3 ? s_INVALID_MOMENT : m.template M_<3> * n / ( ( n - 1.0L ) * (  n - 2.0L  ) ) ;  
      }
      // ======================================================================
      /** get the unbiased estimator for the 4th  order moment:
       *  \f[ \hat{\mu}_4 \equiv  \frac{(n-1)(n^2-3n+3)}{n^3}\mu_4 + 
       *      + \frac{3(2n-3)(n-1)}{n^3}\mu_2^2 \f] 
       *  @see Ya. Dodge and V. Rousson J, "The Complications of the Fourth Central Moment",
       *            The American Statistician, 53 (1999), 276, (doi=10.1080/00031305.1999.10474471)
       *  @see https://amstat.tandfonline.com/doi/abs/10.1080/00031305.1999.10474471
       *  @see https://amstat.tandfonline.com/doi/pdf/10.1080/00031305.1999.10474471
       *  @param  m input counter
       *  @return the unbiased estimate (if number of entries exceeds 4) and 
       *   <code>s_INVALID_MOMENT</code> otherwise
       */
      template <unsigned short N, typename std::enable_if<(N>=4),int>::type = 0 > 
      static inline double unbiased_4th ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( !m.ok() ||  (n < 4 ) ) { return s_INVALID_MOMENT  ; }
        //
        const long double m4 = m.template M_<4> / n ;
        const long double m2 = m.template M_<2> / n ;
        //
        return ( n * m4  * ( 1.0L * n * n - 2.0L * n + 3 ) - 3.0L * n * ( 2.0L * n - 3 ) * m2 * m2 )
          / ( ( n - 1.0L ) * ( n - 2.0L  ) * ( n - 3.0L  ) ) ;
      }
      // ======================================================================
      /** get the unbiased estimator for the 5th  order moment:
       *  \f[ \hat{\mu}_5 \equiv \frac{(n-1)(n-2)}{n^4}\left[
       *      10(n-2)\mu_2\mu_3 + ( n^2-2n+2)\mu_5 \right] \f] 
       *  @param  m input counter
       *  @return the unbiased estimate and 
       *   <code>s_INVALID_MOMENT</code> otherwise
       */
      template <unsigned short N, typename std::enable_if<(N>=5),int>::type = 0 >
      static inline double unbiased_5th ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( !m.ok() ||  ( n < 5 ) ) { return s_INVALID_MOMENT  ; }
        //
        const long double m5 = m.template M_<5> () / n ;
        const long double m2 = m.template M_<2> () / n ;
        const long double m3 = m.template M_<3> () / n ;
        const auto n4 = std::pow ( n * 1.0L , 4 ) ;
        //
        return ( n - 1.0L ) * ( n - 2.0L  ) / n4 *
          ( 10 * ( n - 2.0L ) * m2 * m3 + ( 1.0L * n * n - 2.0L * n +2 ) * m5 ) ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// get the mean
      static inline double mean ( const Moment_<1>& m ) { return m.mu () ; }
      /// get the mean      
      template <unsigned short N, typename std::enable_if<(N>1),int>::type = 0 >
      static inline VE     mean ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size () ;
        if ( !m.ok () || n  < 2  ) { return m.mu (); }  // ATTENTION! 
        const  double mu  = m.mu ()            ;
        const  double m2  = unbiased_2nd ( m ) ;        // ATTENTION! 
        return VE ( mu , m2 / n ) ;
      }
      // ======================================================================
      /// get the variance  
      template <unsigned short N, typename std::enable_if<(1<=N)&&(4>N),int>::type = 0 >
      static inline double variance ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( !m.ok() || 2 > n ) { return s_INVALID_MOMENT ; }        
        return unbiased_2nd ( m ) ;
      }
      // ======================================================================
      /// get the unbiased sample variance with uncertainty
      template <unsigned short N, typename std::enable_if<(N>3),int>::type = 0 >
      static inline VE variance ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size () ;
        if ( !m.ok() || n  < 2  ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        //
        const double m2 = unbiased_2nd ( m ) ;
        //
        if ( n < 4              ) { return m2 ; }                  // ATTENTION! 
        //
        const double m4   = m.template M_<4> () / n ;
        //
        const double cov2 = ( m4 - m2 * m2 * ( n - 3.0L ) / ( n - 1.0L ) ) / n ;
        //
        return VE ( m2 , cov2 ) ;
      }
      // ======================================================================
      /// get the estimate for the sample skewness  \f$ \frac{m_3}{\sigma^{3/2}}\f$
      template <unsigned short N, typename std::enable_if<(N>=3),int>::type = 0 >
      static inline VE skewness ( const Moment_<N>& m ) 
      {
        const unsigned long long n = m.size() ;
        if ( !m.ok() || n < 3 ) { return VE  ( s_INVALID_MOMENT , -1 )  ; } // RETURN
        //
        const double m3   = unbiased_3rd ( m ) ;
        const double m2   = m.template M_<2>  () / n ;
        const double skew = m3 / std::pow ( m2 , 3.0/2 ) ;
        //
        const double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        return VE ( skew , cov2 ) ;
      }
      // ======================================================================
      /// get the estimate for the sample (excessive) kurtosis \f$ \frac{m_4}{\sigma^{4}}-3\f$
      template <unsigned short N, typename std::enable_if<(N>=4),int>::type = 0 >
      static inline VE kurtosis ( const Moment_<N>& m ) 
      {
        const unsigned long long n = m.size() ;
        if ( !m.ok() || n < 4 ) { return VE (  s_INVALID_MOMENT , -1 )  ;}
        const double m4   = unbiased_4th ( m ) ;
        const double m2   = m.template M_<2> () / n  ;
        //
        const double k    =  m4 / ( m2  * m2  ) - 3  ;
        double       cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        cov2 *= 4.0L * ( n * 1.0L * n -1 ) / ( ( n - 3.0L ) * ( n + 5.0L ) ) ;
        //
        return VE  ( k , cov2 ) ;
      }      
      // ======================================================================

      // ======================================================================
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==0),int>::type = 1 >
      static inline double central_moment ( const Moment_<N>& /* m */ ) { return 1 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      static inline double central_moment ( const Moment_<N>& /* m */ ) { return 0 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K>=2)&&(2*K>N)&&(K<=N),int>::type = 1 >
      static inline double central_moment ( const Moment_<N>& m )
      { return ( !m.ok() || m.size() < K ) ? s_INVALID_MOMENT : m.template moment_<K> () ; }
      /// get the central moment of order \f$ N \f$      
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(2<=K) && (N>=2*K),int>::type = 0 >
      static inline VE  central_moment ( const Moment_<N>& m )
      { return ( !m.ok() || m.size() < K ) ? VE(s_INVALID_MOMENT,-1) : m.template moment_<K> () ; }
      // ======================================================================
      
      // ======================================================================
      /// get the standartized central moment of order 1
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==0),int>::type = 1 >
      static inline double std_moment ( const Moment_<N>& /* m */ ) { return 1 ; }
      /// get the standartised central moment of order \f$ K \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      static inline double std_moment ( const Moment_<N>& /* m */ ) { return 0 ; }
      /// get the standartised central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K>=2)&&(K<=N),int>::type = 1 >
      static inline double std_moment ( const Moment_<N>& m )
      { return ( !m.ok() || m.size() < K ) ? s_INVALID_MOMENT : m.template moment_<K> () ; }
      // ======================================================================


      // ======================================================================
      /// get the first cumulant, well, it is actually the first central moment 
      template <unsigned short N,
                typename std::enable_if<(N==1),int>::type = 1 >
      static inline double cumulant_1st ( const Moment_<N>& m ) 
      { return ( !m.ok() || ( m.size() < 1 ) ) ? s_INVALID_MOMENT : m.mu () ; }
      /// get the first cumulant, well, it is actually the first central moment 
      template <unsigned short N,
                typename std::enable_if<(N>1),int>::type = 1 >
      static inline VE     cumulant_1st ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || m.size() < 1 ) { return s_INVALID_MOMENT ; }
        //
        const double m1 = m.mu()  ;
        if ( m.size () < 2 ) { return m1 ; }
        //
        const double m2 = m.template moment_<2>() ;
        const double c2 = m2 / m.size() ;
        //
        return VE ( m1 , c2 ) ;
      }
      /// get the second unbiased cumulant, well it is actually the second unbiased central moment 
      template <unsigned short N,
                typename std::enable_if<(N>=2)&&(N<=3),int>::type = 1 >
        static inline double cumulant_2nd ( const Moment_<N>& m ) 
      { return ( !m.ok() || ( m.size() < 2 ) ) ? s_INVALID_MOMENT : m.size() * m.template M_<2> () / ( m.size() - 1 ) ; }
      /// get the second unbiased cumulant, well it is actually the second unbiased central moment 
      template <unsigned short N,
                typename std::enable_if<(N>=4),int>::type = 1 >
      static inline VE     cumulant_2nd ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 2 ) ) { return s_INVALID_MOMENT ;} 
        //
        const auto n = m.size() ;
        //
        const double k2 = m.template M_<2> () / ( n - 1 ) ; // unbinased estiimate 
        if ( m.size() < 4 ) { return k2 ; }
        //
        const double m2 = m.template M_<2> () / n ; 
        const double m4 = m.template M_<4> () / n ;
        //
        const double k4 = 
          ( ( n + 1 ) * m4  - 3 * m2 * m2 * ( n - 1 ) )  * n * n 
          / ( ( n -1 ) * ( n - 2 ) * ( n - 3 ) ) ;
        // unbinased variance 
        const double c2 = ( 2 * k2 * k2 * n + ( n - 1 ) * k4 ) / ( n * ( n + 1 ) ) ;
        //
        return VE ( k2 , c2 ) ;
      }
      /// get the third unbiased cumulant, well it is actually the third unbiased central moment 
      template <unsigned short N,
                typename std::enable_if<(N>=3)&&(N<6),int>::type = 1 >
      static inline double cumulant_3rd  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 3 ) ) { return s_INVALID_MOMENT ; }
        //
        const auto n = m.size() ;
        const double m3 = m.template M_<3> () / n ; 
        //
        return m3 * n * n  / ( ( n - 1 ) * ( n - 2 ) ) ; 
      }
      /// get the third unbiased cumulant, well it is actually the third unbiased central moment 
      template <unsigned short N,
                typename std::enable_if<(N>=6),int>::type = 1 >
      static inline VE cumulant_3rd  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 3 ) ) { return s_INVALID_MOMENT ; }
        //
        const auto n = m.size() ;
        const double m3 = m.template M_<3> () / n ; 
        //
        const double k3u = m3 * n * n  / ( ( n - 1 ) * ( n - 2 ) ) ; 
        //
        const double k6 = m.template cumulant_<6> () ;
        const double k4 = m.template cumulant_<4> () ;
        const double k3 = m.template cumulant_<3> () ;
        const double k2 = m.template cumulant_<2> () ;
        //
        const double c2 = k6 / n 
          + 9 *     k4 * k2      /   ( n - 1 ) 
          + 9 *     k3 * k3      /   ( n - 1 ) 
          + 6 * n * k2 * k2 * k2 / ( ( n - 1 ) * ( n - 2 ) ) ;
        //
        return VE ( k3u , c2 ) ;
      }
      /// get the 4th unbiased cumulant 
      template <unsigned short N,
                typename std::enable_if<(N>=4)&&(N<8),int>::type = 1 >
      static inline double cumulant_4th  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 4 ) ) { return s_INVALID_MOMENT ; }
        //
        const auto n = m.size() ;
        // 
        const double m2 = m.template M_<2> () / n ; 
        const double m4 = m.template M_<4> () / n ;
        //
        const double k4 = ( ( n + 1 ) * m4 - 3 * m2  *m2 * ( n - 1 ) ) 
          / ( ( n - 1 ) * ( n - 2 ) * ( n -3 ) ) ;
        //
        return k4 ;
      }
      //
      template <unsigned short N,
                typename std::enable_if<(N>=8),int>::type = 1 >
      static inline VE cumulant_4th  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 4 ) ) { return s_INVALID_MOMENT ; }
        //
        const auto n = m.size() ;
        // 
        const double m2  = m.template M_<2> () / n ; 
        const double m4  = m.template M_<4> () / n ;
        //
        const double k4u = ( ( n + 1 ) * m4 - 3 * m2  *m2 * ( n - 1 ) ) 
          / ( ( n - 1 ) * ( n - 2 ) * ( n -3 ) ) ;
        //
        const double k8 = m.template cumulant_<8> () ;
        const double k6 = m.template cumulant_<6> () ;
        const double k5 = m.template cumulant_<5> () ;
        const double k4 = m.template cumulant_<4> () ;
        const double k3 = m.template cumulant_<3> () ;
        const double k2 = m.template cumulant_<2> () ;
        //
        const double c2 =  k8 / n 
          +  16     * k6 * k2      /   ( n - 1 ) 
          +  48     * k5 * k3      /   ( n - 1 ) 
          +  34     * k4 * k4      /   ( n - 1 ) 
          +  72 * n * k4 * k2 * k2 / ( ( n - 1 ) * ( n - 2 ) ) 
          + 144 * n * k3 * k3 * k2 / ( ( n - 1 ) * ( n - 2 ) ) 
          +  24 * n * ( n + 1 ) * std::pow ( m2 , 4 ) / 
          ( ( n -1 ) * ( n -2 ) * ( n -3 ) )  ;
        //
        return VE ( k4u , c2 ) ;
      }
      // ======================================================================

      // ======================================================================
      // Weighted
      // ======================================================================

      // ======================================================================
      /// get the mean 
      static inline double mean ( const WMoment_<1>& m ) { return m.mu () ; }
      /// get the mean      
      template <unsigned short N, typename std::enable_if<(N>1),int>::type = 0 >
      static inline VE     mean ( const WMoment_<N>& m )
      {
        if ( !m.ok() || m.size() < 2 ) { return m.mu() ; }
        //
        const auto n = m.nEff () ;
        const  double mu  = m.mu     (   ) ;
        const  double m2  = m.template moment_<2> () ;
        return VE ( mu , m2 / n ) ;
      }
      
      // ======================================================================
      /// get the variance  
      template <unsigned short N, typename std::enable_if<(1<=N)&&(4>N),int>::type = 0 >
      static inline double variance ( const WMoment_<N>& m )
      { return !m.ok() || ( m.size() < 2 ) ? s_INVALID_MOMENT : m.template moment_<2>() ; }      
      // ======================================================================
      /// get the  sample variance with uncertainty
      template <unsigned short N, typename std::enable_if<(N>3),int>::type = 0 >
      static inline VE variance ( const WMoment_<N>& m )
      {
        if ( !m.ok() || m.size() < 2  ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        //
        const double m2 = m.template M_<2> () / m.w()  ;
	    // 
        if ( 0 >  m2 ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        //
        if ( m.size() < 4  ) { return m2 ; }  // RETURN 
        //
        const double m4 = m.template M_<4> () / m.w () ;
        const auto n = m.nEff () ;
        if ( !n || !m4 || m4 <= 0 ) { return m2 ; }  // RETURN 
        //
        const double cov2 = ( m4 - m2 * m2 * ( n - 3 ) / ( n - 1 ) ) / n ;
        //
        return 0 <= cov2 ? VE ( m2 , cov2 ) : VE ( m2 , 0.0 ) ;
      }
      
      // ======================================================================
      /// get the estimate for the sample skewness  \f$ \frac{m_3}{\sigma^{3/2}}\f$
      template <unsigned short N, typename std::enable_if<(N>=3),int>::type = 0 >
      static inline VE skewness ( const WMoment_<N>& m ) 
      {
        if ( !m.ok() || ( m.size() < 3 ) ) { return VE  ( s_INVALID_MOMENT , -1 )  ; }
        //
        const auto   n    = m.nEff () ;
        const double m3   = m.template M_<3> () / m.w () ;
        const double m2   = m.template M_<2> () / m.w () ;
        const double skew = m3 / std::pow ( m2 , 3.0/2 ) ;
        const double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        //
        return 0 <= cov2 ? VE ( skew , cov2 ) : VE ( skew , 0.0 ) ;  
      }
      // ======================================================================
      /// get the estimate for the sample (excessive) kurtosis \f$ \frac{m_4}{\sigma^{4}}-3\f$
      template <unsigned short N, typename std::enable_if<(N>3),int>::type = 0 >
      static inline VE kurtosis ( const WMoment_<N>& m ) 
      {
        if ( !m.ok() || ( m.size() < 4 ) ) { return VE  ( s_INVALID_MOMENT , -1 )  ; }
        //
        const auto   n  = m.nEff () ;
        const double m4 = m.template M_<4> () / m.w() ;
        const double m2 = m.template M_<2> () / m.w() ;
        const double k  = m4 / ( m2  * m2  ) - 3  ;
        double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        cov2 *= 4.0L * ( n * 1.0L * n -1 ) / ( ( n - 3.0L ) * ( n + 5.0L ) ) ;
        //
        return 0 <= cov2 ? VE ( k , cov2 ) : VE  ( k , 0.0 ) ;
      }
      // ======================================================================

      // ======================================================================
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==0),int>::type = 1 >
      static inline double central_moment ( const WMoment_<N>& /* m */ ) { return 1 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      static inline double central_moment ( const WMoment_<N>& /* m */ ) { return 0 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K>=2)&&(2*K>N)&&(K<=N),int>::type = 1 >
      static inline double central_moment ( const WMoment_<N>& m )
      { return ( !m.ok() || m.size() < K ) ? s_INVALID_MOMENT : m.template moment_<K> () ; }
      /// get the central moment of order \f$ N \f$      
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(2<=K) && (N>=2*K),int>::type = 0 >
      static inline VE     central_moment ( const WMoment_<N>& m )
      { return ( !m.ok() || m.size() < K ) ? VE( s_INVALID_MOMENT , -1 ) : m.template moment_<K> () ; }
      
      // ======================================================================
      /// get the standartized moment of order 1
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==0),int>::type = 1 >
      static inline double std_moment ( const WMoment_<N>& /* m */ ) { return 1 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      static inline double std_moment ( const WMoment_<N>& /* m */ ) { return 0 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K>=2)&&(K<=N),int>::type = 1 >
      static inline double std_moment ( const WMoment_<N>& m )
      { return ( !m.ok() || m.size() < K ) ? s_INVALID_MOMENT : m.template moment_<K> () ; }
      
      // ======================================================================
    } ; //                                The end of class Ostap::Math::Moments 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END  
// ============================================================================
#endif // OSTAP_MOMENTS_H
// ============================================================================

