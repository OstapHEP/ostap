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
     *  Helper empty base class
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
     *  the  high-order central momentts
     *  \f[  \mu_n \equiv \sum_{i}  \left( x_i - \bar{x} \right)^n \f] 
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
      double moment () const { return this->M ( N ) / this->size () ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      double moment ( const unsigned short k ) const
      { return N <  k ? 0 : 0 == k ? 1 : 1 == k ? 0 : ( this->M ( k ) / this->size() ) ; }
      // ======================================================================
      /// get number of entries
      unsigned long long size  () const { return m_prev.size () ; }
      /// get the mean value (if \f$ 1 \le N \$4)
      long double        mu    () const { return m_prev.mu   () ; } 
      // ======================================================================
    public: // basic operations for the counter 
      // ======================================================================
      /// increment with some value 
      Moment_& operator+= ( const double   x ) ;      
      /// increment with other moment 
      Moment_& operator+= ( const Moment_& x ) ;
      // ======================================================================
    public: // add more values to the counter 
      // ======================================================================
      /// add single value 
      void  add ( const double x ){ *this += x ; }
      /// add sequence of values  
      template <class ITERATOR>
      void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
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
      /// get "previos" moment 
      const Moment_<N-1>& previous () const { return this->m_prev ; }
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
    constexpr std::array<unsigned long long,N+1> Moment_<N>::s_Ck { Ostap::Math::choose_array<N>() } ;
    // ========================================================================
    /// increment with some value
    template <unsigned short N>
    Moment_<N>& Moment_<N>::operator+= ( const double  x )
    {
      const auto nA = this->size() ;
      const long double delta = x - this->mu() ;
      const long double b_n =      -1.0L / ( nA + 1 ) ;
      const long double a_n =  nA * 1.0L / ( nA + 1 ) ;
      const long double d_n = - delta    / ( nA + 1 ) ;
      //
      m_M += ( nA * std::pow ( b_n , N ) + std::pow ( a_n , N ) ) * std::pow ( delta , N ) ;
      long double d = 1 ;
      for ( unsigned int k = 1 ; k + 2 <= N ; ++k )
      {
        d   *= d_n ;
        m_M +=  s_Ck [ k ] * this-> M ( N - k ) * d   ;
      }
      /// update previous 
      this->m_prev += x ; // update previous
      //
      return *this ;
    }
    // ========================================================================
    /// increment with some other counter 
    template <unsigned short N>
    Moment_<N>& Moment_<N>::operator+= ( const Moment_<N>&  x )
    {
      const unsigned long long nA    = this->size () ;
      const unsigned long long nB    =    x. size () ;
      const unsigned long long n     = nA + nB       ;
      const long double        delta = x.mu() - this->mu() ;
      const long double        b_n   =  ( -1.0L * nB ) / ( nA + nB ) ;
      const long double        a_n   =  (  1.0L * nA ) / ( nA + nB ) ;
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
      double moment () const { return 1 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      double moment ( const unsigned short k ) const { return 0 < k ? 0 : 1 ; }
      // ======================================================================
      /// get number of entries
      unsigned long long size  () const { return m_size ; }
      // ======================================================================
    public: // basic operations 
      // ======================================================================
      /// increment with some value 
      Moment_& operator+= ( const double /* x */ ) { ++m_size           ; return *this ; }
      /// increment with other moment 
      Moment_& operator+= ( const Moment_&  x    ) { m_size += x.m_size ; return *this ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add single value
      void add ( const double /* x */ ){ ++m_size ; }
      /// add sequence of values  
      template <class ITERATOR>
      void add ( ITERATOR begin , ITERATOR end   )
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
    public :
      // ======================================================================
      /** get the value of the Nth moment 
       *  \f[ \mu_N \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      double moment () const { return 0 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      double moment ( const unsigned short k ) const { return 1 <= k ? 0 : 1 ; }
      /// get number of entries
      unsigned long long size  () const { return m_prev.size()  ; }
      // get the mean value
      long double        mu    () const { return m_mu ; } 
      // ======================================================================
    public: // basic operations 
      // ======================================================================
      /// increment with some value 
      Moment_& operator+= ( const double  x )
      {
        const auto n = m_prev.size() ;
        m_mu = ( n * m_mu + x ) /  ( n + 1 );  // calculate new mean value 
        m_prev += x ;                          // updated previous  
        return *this ;
      }
      /// increment with other moment 
      Moment_& operator+= ( const Moment_ x )
      {
        const auto n1 =   m_prev.size() ;
        const auto n2 = x.m_prev.size() ;
        m_mu = ( n1 * m_mu + n2 * x.m_mu ) / ( n1 + n2 ) ; // update mean 
        m_prev += x.m_prev ;                               // update previous   
        return *this ;  
      }
      // ======================================================================
    public:
      // ======================================================================
      /// add single value 
      void add ( const double x ){ (*this) += x ; }
      /// add sequence of values  
      template <class ITERATOR>
      void add ( ITERATOR begin , ITERATOR end   )
      { for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; } }
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
      /// get "previos" moment 
      const Moment_<0>& previous () const { return this->m_prev ; }
      // ======================================================================
    public: 
      // ======================================================================
      inline long double M ( const unsigned short k ) const
      { return 1 < k  ? 0 : 1 == k ? 0 : this->m_prev. M ( k ) ; }
      // ======================================================================
    private:    private:
      // ======================================================================
      Moment_<0>  m_prev { } ;
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
    /**  @class Moments 
     *   collectionn  of static functions dealing with moments 
     */
    class Moments
    {
      // ======================================================================
      typedef Ostap::Math::ValueWithError VE    ;
      // ======================================================================
    private: 
      // ======================================================================
      /** @var s_INVALID_CENTRAL_MOMENT ; 
       *  the invalid value of the central moment 
       */ 
      static const double s_INVALID_MOMENT ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the unbiased estimator for the 2nd order moment:
       *  \f[ \hat{\mu}_2 \equiv \frac{n}{n-1} \mu_2 \f] 
       *  @param  m input counter
       *  @return the unbiased estimate (if number of entries exceeds 2) and 
       *   <code>s_INVALID_MOMENT</code> otherwise 
       */
      template <unsigned short N, typename std::enable_if<(1<N),int>::type = 0 >
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
      template <unsigned short N, typename std::enable_if<(2<N),int>::type = 0 >
      static inline double unbiased_3rd ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size() ;
        return  n < 3 ? s_INVALID_MOMENT : m.M ( 3 ) * n / ( ( n - 1.0 ) * (  n - 2.0  ) ) ;  
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
      template <unsigned short N, typename std::enable_if<(3<N),int>::type = 0 >
      static inline double unbiased_4th ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( 0 == n ) { return s_INVALID_MOMENT  ; }
        const long double m4 = m.M ( 4 ) / n ;
        const long double m2 = m.M ( 2 ) / n ;
        //
        return n < 4 ? s_INVALID_MOMENT  :
          ( n * m4  * ( 1.0L * n * n - 2.0L * n + 3 ) - 3.0L * n * ( 2.0L * n - 3 ) * m2 * m2 )
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
        if ( 0 == n ) { return s_INVALID_MOMENT  ; }
        const long double m5 = m.M ( 5 ) / n ;
        const long double m2 = m.M ( 2 ) / n ;
        const long double m3 = m.M ( 3 ) / n ;
        const auto n4 = std::pow ( n * 1.0L , 4 ) ;
        //
        return n < 5 ? s_INVALID_MOMENT  :
          ( n - 1.0L  ) * ( n - 2.0L  ) / n4 *
                   ( 10 * ( n - 2.0L ) * m2 * m3 + ( 1.0L * n * n - 2.0L * n +2 ) * m5 ) ;
      }
      // ======================================================================
      static inline double variance ( const Moment_<2>& m ) { return unbiased_2nd ( m ) ; }
      static inline double variance ( const Moment_<3>& m ) { return unbiased_2nd ( m ) ; }
      // ======================================================================
      /// get the unbiased sample variance with uncertainty
      template <unsigned short N, typename std::enable_if<(3<N),int>::type = 0 >
      static inline VE variance ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size () ;
        if ( n  < 2  ) { return s_INVALID_MOMENT ; }
        //
        typedef Ostap::Math::ValueWithError VE    ;
        //
        const double m2 = unbiased_2nd ( m ) ;
        if ( 0 > m2  ) { return s_INVALID_MOMENT ; } // RETURN
        if ( 0 == m2 ) { return VE ( 0 , 0 )     ; } // ?
        //      
        const double m4 = m.moment ( 4 ) ;
        //
        const double cov2 = ( m4 - m2 * m2 * ( n - 3 ) / ( n -1 ) ) / n ;
        //
        return VE ( m2 , m2 * cov2 ) ;
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
      /// get the estimaet for the sample (excessive) kurtosis \f$ \frac{m_4}{\sigma^{4}}-3\f$
      template <unsigned short N, typename std::enable_if<(3<N),int>::type = 0 >
      static inline VE kurtosis ( const Moment_<N>& m ) 
      {
        const unsigned long long n = m.size() ;
        if ( n < 4 ) { return VE(  s_INVALID_MOMENT , -1 )  ;}
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
       *  - the error estimaet is possible only when \f$ 2N \le K \f$!
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
       *  - the error estimaet is possible only when \f$ 2N \le K \f$!
       *  @aparam m counter 
       *  @return moment with uncertainty for non-empty counter 
       *          <code>s_INVALID_MOMENT</code> for empty counters 
       */
        template <unsigned short N, unsigned short K,
                typename std::enable_if<(1<N)&&(2*N<=K),int>::type = 0 >
      static inline VE _central_moment_2 ( const Moment_<K>& m )
      { return central_moment<N> ( m ) ; }
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
