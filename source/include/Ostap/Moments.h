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
    /** @struct Moment
     *  Helper (empty) base class for moment-counters 
     *  - it is not really neeeded for C++, but it simplifies python decorations 
     */
    class Moment
    {
    public :
      // ======================================================================
      virtual ~Moment() ;
      /// add new value to the counter 
      virtual void update ( const double x ) = 0 ;
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
      /** get the value of the Nth moment 
       *  \f[ \mu_N \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      inline double moment () const { return this->M ( N ) / this->size () ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const
      { return N <  k ? 0 : 0 == k ? 1 : 1 == k ? 0 : ( this->M ( k ) / this->size() ) ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 
          N <  k ? 0 : 
          0 == k ? 1 : 
          1 == k ? 0 :
          2 == k ? 1 : this->moment ( k ) / std::pow ( this->moment ( 2 ) , 0.5 * k ) ;
      }
      // ======================================================================
      /// get number of entries
      inline unsigned long long size  () const { return m_prev.size  () ; }
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
      inline long double M ( const unsigned short k ) const
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
      { return this->empty () ? 0 : this->m_M / this->size  () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(N<2*K)&&(N>K),int>::type = 0 >
      inline long double moment_ () const
      { return this->m_prev.template moment_<K> () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError
      moment_ () const
      {
        if  ( this -> empty() ) { return 0 ; }
        const unsigned long long n = this->size() ;
        //
        const long double muo  = this->template M_<K>   () / n ;
        const long double mu2o = this->template M_<2*K> () / n ;
        const long double muop = this->template M_<K+1> () / n ;
        const long double muom = this->template M_<K-1> () / n ;  
        const long double mu2  = this->         M_<2>   () / n ;  
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
        return this->empty () ? 0 :
          this->template moment_<K> () / std::pow ( this->moment_<2>() , 0.5 * K ) ;
      }      
      // =======================================================================
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
      inline unsigned long long size  () const { return m_size ; }
      /// empty ?
      inline bool               empty () const { return 0 == m_size ; }
      /// ok ?
      inline bool               ok    () const { return 0 != m_size ; }
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
      inline long double M ( const unsigned short k ) const
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
      unsigned long long m_size ;
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
      inline unsigned long long size  () const { return m_prev.size()  ; }
      // get the mean value
      inline long double        mu    () const { return m_mu ; } 
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
      inline long double M ( const unsigned short k ) const
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
    /** @struct WMoment
     *  Helper (empty) base class for weighted moment-counters 
     *  - it is not really neeeded for C++, but it simplifies python decorations 
     */
    class WMoment
    {
    public :
      // ======================================================================
      virtual ~WMoment() ;
      /// add new value to the counter 
      virtual void update ( const double x , const double w = 1 ) = 0 ;
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
      /** get the value of the Nth weighted  moment 
       *  \f[ \mu_N \equiv  \frac{1}{\sum w_i} \sum w_i \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      inline double moment () const { return this->M ( N ) / this-> w () ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{\sum w_i} \sum w_i \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const
      { return N <  k ? 0 : 0 == k ? 1 : 1 == k ? 0 : ( this->M ( k ) / this-> w () ) ; }
      // ======================================================================
      /** get value of the kth standartized moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double std_moment ( const unsigned short k ) const
      { return 
          N <  k ? 0 : 
          0 == k ? 1 : 
          1 == k ? 0 :
          2 == k ? 1 : this->moment ( k ) / std::pow ( this->moment ( 2 ) , 0.5 * k ) ;
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
      inline long double M ( const unsigned short k ) const
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
        const auto n = this->nEff () ;
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
      if ( !w ) { return *this ; }
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
      inline bool               ok    () const { return 0 != m_size && 0 != m_w ; }
      // ======================================================================
    public: // basic operations 
      // ======================================================================
      /// increment with other moment 
      inline WMoment_& operator+= ( const WMoment_&  x ) { return add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add a single value
      inline WMoment_& add ( const double /* x */ , const double w = 1 )
      { 
        if ( !w ) { return *this ; }
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
      inline long double M ( const unsigned short k ) const
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
      unsigned long long m_size ; // number of entries
      /// sum of weights \f$  \sum w_i \f$ 
      long double        m_w    ; // sum of weights
      /// sum of weights squared \f$  \sum w_i^2 \f$ 
      long double        m_w2   ; // sum of weights squared 
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
        if ( !w ) { return *this ; }
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
      inline long double M ( const unsigned short k ) const
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
    private:
      // ======================================================================
      WMoment_<0>  m_prev { } ;
      /// mean value
      long double m_mu   {0} ; // mean value
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
      template <unsigned short N, typename std::enable_if<(2<=N),int>::type = 0 >
      static inline double unbiased_2nd ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size() ;
        return  n < 2 ? s_INVALID_MOMENT  : m.M ( 2 ) / ( n - 1 ) ;  
      }
      // ======================================================================
      /** get the unbiased estimator for the 3rd order moment:
       *  \f[ \hat{\mu}_3 \equiv \frac{n^2}{(n-1)(n-2)} \mu_3 \f] 
       *  @param  m input counter
       *  @return the unbiased estimate (if number of entries exceeds 3) and 
       *   <code>s_INVALID_MOMENT</code> otherwise 
       */
        template <unsigned short N, typename std::enable_if<(3<=N),int>::type = 0 >
        static inline double unbiased_3rd ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size() ;
        return  n < 3 ? s_INVALID_MOMENT : m.M ( 3 ) * n / ( ( n - 1.0L ) * (  n - 2.0L  ) ) ;  
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
      template <unsigned short N, typename std::enable_if<(4<=N),int>::type = 0 > 
      static inline double unbiased_4th ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( 4 > n ) { return s_INVALID_MOMENT  ; }
        const long double m4 = m.M ( 4 ) / n ;
        const long double m2 = m.M ( 2 ) / n ;
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
        template <unsigned short N, typename std::enable_if<(4<N),int>::type = 0 >
        static inline double unbiased_5th ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( 5 > n ) { return s_INVALID_MOMENT  ; }
        const long double m5 = m.M ( 5 ) / n ;
        const long double m2 = m.M ( 2 ) / n ;
        const long double m3 = m.M ( 3 ) / n ;
        const auto n4 = std::pow ( n * 1.0L , 4 ) ;
        //
        return ( n - 1.0L ) * ( n - 2.0L  ) / n4 *
                  ( 10 * ( n - 2.0L ) * m2 * m3 + ( 1.0L * n * n - 2.0L * n +2 ) * m5 ) ;
      }
      // ======================================================================
      /// get the mean
      static inline double mean ( const Moment_<1>& m ) { return m.mu () ; }
      /// get the mean      
      template <unsigned short N, typename std::enable_if<(1<N),int>::type = 0 >
      static inline VE     mean ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size () ;
        if ( n  < 2  ) { return s_INVALID_MOMENT ; }
        const  double mu  = m.mu(   )          ;
        const  double m2  = unbiased_2nd ( m ) ;
        return VE ( mu , m2 / n ) ;
      }        
      // ======================================================================
      /// get the variance  
      static inline double variance ( const Moment_<2>& m ) { return unbiased_2nd ( m ) ; }
      /// get the variance  
      static inline double variance ( const Moment_<3>& m ) { return unbiased_2nd ( m ) ; }
      // ======================================================================
      /// get the unbiased sample variance with uncertainty
      template <unsigned short N, typename std::enable_if<(3<N),int>::type = 0 >
      static inline VE variance ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size () ;
        if ( n  < 2  ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        //
        const double m2 = unbiased_2nd ( m ) ;
        if ( 0 > m2  ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        if ( 0 == m2 ) { return VE ( 0 , 0 )     ; } // ?
        //      
        const double m4 = m.moment ( 4 ) ;
        //
        const double cov2 = ( m4 - m2 * m2 * ( n - 3.0L ) / ( n - 1.0L ) ) / n ;
        //
        return VE ( m2 , cov2 ) ;
      }
      // ======================================================================
      /// get the estimate for the sample skewness  \f$ \frac{m_3}{\sigma^{3/2}}\f$
      template <unsigned short N, typename std::enable_if<(2<N),int>::type = 0 >
      static inline VE skewness ( const Moment_<N>& m ) 
      {
        const unsigned long long n = m.size() ;
        if ( n < 3 ) { return VE  ( s_INVALID_MOMENT , -1 )  ; }
        const double m3   = unbiased_3rd ( m ) ;
        const double m2   = m.moment     ( 2 ) ;
        const double skew =  m3 / std::pow ( m2 , 3.0/2 ) ;
        const double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        return VE ( skew , cov2 ) ;
      }
      // ======================================================================
      /// get the estimate for the sample (excessive) kurtosis \f$ \frac{m_4}{\sigma^{4}}-3\f$
      template <unsigned short N, typename std::enable_if<(3<N),int>::type = 0 >
      static inline VE kurtosis ( const Moment_<N>& m ) 
      {
        const unsigned long long n = m.size() ;
        if ( n < 4 ) { return VE (  s_INVALID_MOMENT , -1 )  ;}
        const double m4 = unbiased_4th ( m ) ;
        const double m2 = m.moment     ( 2 ) ;
        const double k  =  m4 / ( m2  * m2  ) - 3  ;
        double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        cov2 *= 4.0L * ( n * 1.0L * n -1 ) / ( ( n - 3.0L ) * ( n + 5.0L ) ) ;
        return VE  ( k , cov2 ) ;
      }
      // ======================================================================
      /** get the central moment of order \f$ N \f$  
       *  @aparam m counter 
       *  @return moment for non-empty counter 
       *          <code>s_INVALID_MOMENT</code> for empty counters 
       */
      template <unsigned short N, unsigned short K,
                typename std::enable_if< (1<N) && (N<=K) && (K<2*N),int>::type = 1 >
      static inline double central_moment ( const Moment_<K>& m )
      { return m.moment ( N ) ; }
      // ======================================================================
      /** get the central moment of order \f$ N \f$  with 
       *  the estimate of the uncertainty (with \f$O(n^{-2})\f$~precision
       *  - the error estimate is possible only when \f$ 2N \le K \f$!
       *  @aparam m counter 
       *  @return moment with uncertainty for non-empty counter 
       *          <code>s_INVALID_MOMENT</code> for empty counters 
       */
      template <unsigned short N, unsigned short K,
                typename std::enable_if< (1<N) && (2*N<=K),int>::type = 0 >
      static inline VE  central_moment ( const Moment_<K>& m )
      {
        //
        const unsigned long long n = m.size() ;
        if ( 0 == n    ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        //
        const long double muo  = m.M ( N     ) / n ;
        const long double mu2o = m.M ( 2 * N ) / n ;
        const long double muop = m.M ( N + 1 ) / n ;
        const long double muom = m.M ( N - 1 ) / n ;  
        const long double mu2  = m.M ( 2     ) / n ;  
        //
        long double cov2 = mu2o     ;
        cov2 -= 2 * N * muop * muom ;
        cov2 -=         muo  * muo  ;
        cov2 += N * N * mu2  * muom * muom ;
        cov2 /= n ;
        //
        return VE ( muo , cov2 ) ;
      }
      // ======================================================================
      /** get the central moment of order \f$ N \f$  with 
       *  the estimate of the uncertainty (with \f$O(n^{-2})\f$~precision
       *  - the error estimate is possible only when \f$ 2N \le K \f$!
       *  @aparam m counter 
       *  @return moment with uncertainty for non-empty counter 
       *          <code>s_INVALID_MOMENT</code> for empty counters 
       */
      template <unsigned short N, unsigned short K,
                typename std::enable_if<(1<N)&&(2*N<=K),int>::type = 0 >
      static inline VE _central_moment_2 ( const Moment_<K>& m )
      { return central_moment<N> ( m ) ; }
      // ======================================================================
      /** get the standartized moment of order 1 
       */
      template <unsigned short N, unsigned short K,
                 typename std::enable_if< (1==N) && (N<=K) ,int>::type = 1 > 
      static inline double std_moment ( const Moment_<K>& /* m */ ) { return 0 ; }
      // =======================================================
      /** get the standartized moment of order 2 
       */
      template <unsigned short N, unsigned short K,
                 typename std::enable_if< (2==N) && (N<=K) ,int>::type = 1 > 
      static inline double std_moment ( const Moment_<K>& /* m */ ) { return 1 ; }   
      /** get the standartized moment of order N
       */
      template <unsigned short N, unsigned short K,
                 typename std::enable_if< (2<N) && (N<=K) ,int>::type = 1 > 
      static inline double std_moment ( const Moment_<K>& m ) 
      { 
        const double m2 = m.moment ( 2 ) ;
        return m2 ? m.moment ( N ) / std::pow ( m2 , 0.5 * N ) : s_INVALID_MOMENT ;
      }

      // ======================================================================
      // Weighted
      // ======================================================================
        
      // ======================================================================
      /// get the mean 
      static inline double mean ( const WMoment_<1>& m ) { return m.mu () ; }
      /// get the mean      
      template <unsigned short N, typename std::enable_if<(1<N),int>::type = 0 >
      static inline VE     mean ( const WMoment_<N>& m )
      {
        const auto n = m.nEff () ;
        if ( m.size()  < 2  ) { return VE ( s_INVALID_MOMENT , -1  ) ; }
        const  double mu  = m.mu     (   ) ;
        const  double m2  = m.moment ( 2 ) ;
        return VE ( mu , m2 / n ) ;
      }        
      // ======================================================================
      /// get the variance  
      static inline double variance ( const WMoment_<2>& m )
      { return m.size() < 2  ? s_INVALID_MOMENT : m.moment ( 2 ) ; }
      /// get the variance  
      static inline double variance ( const WMoment_<3>& m )
      { return m.size() < 2  ? s_INVALID_MOMENT : m.moment ( 2 ) ; }
      // ======================================================================
      /// get the  sample variance with uncertainty
      template <unsigned short N, typename std::enable_if<(3<N),int>::type = 0 >
      static inline VE variance ( const WMoment_<N>& m )
      {
        if ( m.size()  < 2  ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        const auto n = m.nEff () ;
        //
        const double m2 = m.moment ( 2 )  ; 
        if ( 0 > m2  ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        if ( 0 == m2 ) { return VE ( 0 , 0 )     ; } // ?
        //      
        const double m4 = m.moment ( 4 ) ;
        //
        const double cov2 = ( m4 - m2 * m2 * ( n - 3 ) / ( n - 1 ) ) / n ;
        //
        return VE ( m2 , cov2 ) ;
      }
      // ======================================================================
      /// get the estimate for the sample skewness  \f$ \frac{m_3}{\sigma^{3/2}}\f$
      template <unsigned short N, typename std::enable_if<(2<N),int>::type = 0 >
      static inline VE skewness ( const WMoment_<N>& m ) 
      {
        if ( m.size() < 3 ) { return VE  ( s_INVALID_MOMENT , -1 )  ; }
        const auto   n    = m.nEff () ;
        const double m3   = m.moment ( 3 ) ;
        const double m2   = m.moment ( 2 ) ;
        const double skew =  m3 / std::pow ( m2 , 3.0/2 ) ;
        const double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        return VE ( skew , cov2 ) ;
      }
      // ======================================================================
      /// get the estimate for the sample (excessive) kurtosis \f$ \frac{m_4}{\sigma^{4}}-3\f$
      template <unsigned short N, typename std::enable_if<(3<N),int>::type = 0 >
      static inline VE kurtosis ( const WMoment_<N>& m ) 
      {
        if ( m.size() < 4 ) { return VE  ( s_INVALID_MOMENT , -1 )  ; }
        const auto   n  = m.nEff () ;
        const double m4 = m.moment ( 4 ) ;
        const double m2 = m.moment ( 2 ) ;
        const double k  = m4 / ( m2  * m2  ) - 3  ;
        double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        cov2 *= 4.0L * ( n * 1.0L * n -1 ) / ( ( n - 3.0L ) * ( n + 5.0L ) ) ;
        return VE  ( k , cov2 ) ;
      }        
      // ======================================================================
      /** get the central moment of order \f$ N \f$  
       *  @aparam m counter 
       *  @return moment for non-empty counter 
       *          <code>s_INVALID_MOMENT</code> for empty counters 
       */
      template <unsigned short N, unsigned short K,
                typename std::enable_if< (1<N) && (N<=K) && (K<2*N),int>::type = 1 >
      static inline double central_moment ( const WMoment_<K>& m )
      { return m.moment ( N ) ; }
      // ======================================================================
      /** get the central moment of order \f$ N \f$  with 
       *  the estimate of the uncertainty (with \f$O(n^{-2})\f$~precision
       *  - the error estimaet is possible only when \f$ 2N \le K \f$!
       *  @aparam m counter 
       *  @return moment with uncertainty for non-empty counter 
       *          <code>s_INVALID_MOMENT</code> for empty counters 
       */
      template <unsigned short N, unsigned short K,
                typename std::enable_if< (1<N) && (2*N<=K),int>::type = 0 >
      static inline VE  central_moment ( const WMoment_<K>& m )
      {
        //
        if ( 0 == m.size() ) { return VE ( s_INVALID_MOMENT , -1 ) ; } // RETURN
        const auto n = m.nEff () ;
        //
        const long double muo  = m.M ( N     ) / n ;
        const long double mu2o = m.M ( 2 * N ) / n ;
        const long double muop = m.M ( N + 1 ) / n ;
        const long double muom = m.M ( N - 1 ) / n ;  
        const long double mu2  = m.M ( 2     ) / n ;  
        //
        long double cov2 = mu2o     ;
        cov2 -= 2 * N * muop * muom ;
        cov2 -=         muo  * muo  ;
        cov2 += N * N * mu2  * muom * muom ;
        cov2 /= n ;
        //
        return VE ( muo , cov2 ) ;
      }
      // ======================================================================
      /** get the central moment of order \f$ N \f$  with 
       *  the estimate of the uncertainty (with \f$O(n^{-2})\f$~precision
       *  - the error estimaet is possible only when \f$ 2N \le K \f$!
       *  @aparam m counter 
       *  @return moment with uncertainty for non-empty counter 
       *          <code>s_INVALID_MOMENT</code> for empty counters 
       */
      template <unsigned short N, unsigned short K,
                typename std::enable_if<(1<N)&&(2*N<=K),int>::type = 0 >
      static inline VE _central_moment_3 ( const WMoment_<K>& m )
      { return central_moment<N> ( m ) ; }
      // ======================================================================
      /** get the standartized moment of order 1 
       */
      template <unsigned short N, unsigned short K,
                 typename std::enable_if< (1==N) && (N<=K) ,int>::type = 1 > 
      static inline double std_moment ( const WMoment_<K>& /* m */ ) { return 0 ; }
      // =======================================================
      /** get the standartized moment of order 2 
       */
      template <unsigned short N, unsigned short K,
                 typename std::enable_if< (2==N) && (N<=K) ,int>::type = 1 > 
      static inline double std_moment ( const WMoment_<K>& /* m */ ) { return 1 ; }   
      /** get the standartized moment of order N
       */
      template <unsigned short N, unsigned short K,
                 typename std::enable_if< (2<N) && (N<=K) ,int>::type = 1 > 
      static inline double std_moment ( const WMoment_<K>& m ) 
      { 
        const double m2 = m.moment ( 2 ) ;
        return m2 ? m.moment ( N ) / std::pow ( m2 , 0.5 * N ) : s_INVALID_MOMENT ;
      }
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

