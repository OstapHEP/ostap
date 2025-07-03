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
#include <limits>
#include <type_traits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Choose.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/Statistic.h"
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  namespace Math
  {
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
    protected:
      // ======================================================================
      /// return the value for invalid moment 
      double invalid_moment () const ;
      // ======================================================================
    } ;
    // ========================================================================
    /// forward declaration 
    template <unsigned short N> class Moment_   ;
    /// template specialization for N=0
    template <>                 class Moment_<0>;
    /// template specialization for N=1
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
    template <unsigned short N> class Moment_ ;
    // =======================================================================
    /// specialization for \f$N=0f\$
    template <>
    class Moment_<0> : public Moment 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// # entries
      typedef unsigned long long size_type ;
      /// data type 
      typedef long double        data_type ;
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 0 } ;
      // ======================================================================
    public:
      // ======================================================================
      /// (default) constructor 
      Moment_ ( const size_type  size = 0 ) ; 
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
      inline double moment     ( const unsigned short k ) const
      { return 0 == k ? 1.0 : this->invalid_moment() ; }
      // ========================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 0 == k ? 1.0 : this->invalid_moment() ; }
      // =======================================================================
      /// get value of the kth centralzed moment for \f$  k \le N \f$
      inline double centralized_moment 
      ( const unsigned short k       , 
        const double    /* center */ ) const 
      { return 0 == k ? 1.0 : this->invalid_moment() ; }
      // =======================================================================
    public:
      // ======================================================================
      /// get number of entries
      inline size_type size  () const { return m_size  ; }
      /// get effective number of entries 
      inline size_type nEff  () const { return size () ; }
      /// empty ?
      inline bool      empty () const { return 0 == m_size ; }
      /// ok ?
      inline bool      ok    () const { return      m_size ; }
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
      inline Moment_& add ( const double   x )
      { if ( std::isfinite ( x ) ) { ++m_size ; } ; return *this ; }
      /// add single value
      inline Moment_& add ( const Moment_& x    ) { m_size += x.m_size ; return *this ; }
      /// add sequence of values  
      template <class ITERATOR>
      inline
      Moment_&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
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
      /// Reset the content
      void reset  () override { m_size = 0 ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      inline bool isfinite () const  { return true ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of \f$ M_k = \sum_i \left( x_i - \bar{x} \right)^k \f$ 
       *  for \f$ k \le N \f$
       *  @param k the order   \f$  0 \le k \le N \f$
       *  @return the value of \$ M_k \f$  if \f$  0 \le k \le  \f$, 0 otherwise 
       */   
      inline long double M ( const unsigned short k = 0 ) const
      { return 0  == k ? 1.0 * m_size : this->invalid_moment() ; }
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
    public: // templated centralized moment 
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline long double centralized_moment_ 
      ( const double /* center */ ) const { return 1 ; }
      // ====================================================================== 
    private:
      // ======================================================================
      size_type  m_size { 0 } ;
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
      // # of entries 
      typedef typename Moment_<0>::size_type size_type ; 
      /// data type 
      typedef typename Moment_<0>::data_type data_type ; 
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 1 } ;
      // ======================================================================
    public: 
      // ======================================================================
      /// default constructor 
      Moment_ () = default ;
      /// constructor from mu and previous 
      Moment_
      ( const Moment_<0>& prev ,	
	      const double      mu   ,
	      const double      xmin ,
	      const double      xmax ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** get the value of the Nth moment 
       *  \f[ \mu_N \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^N \f]
       *  @return the value of the Nth central moment
       */
      inline double moment      () const { return 0 ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment     ( const unsigned short k ) const
      { return ( 0 == k ) ? 1.0 : ( 1 == k ) ? 0.0 : this->invalid_moment () ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return ( 0 == k ) ? 1.0 : ( 1 == k ) ? 0.0 : this->invalid_moment () ; }
      // ======================================================================
      /// get value of the kth centralized moment for \f$  k \le N \f$
      inline double centralized_moment 
      ( const unsigned short k      , 
        const double         center ) const
      { return ( 0 == k ) ? 1.0 : ( 1 == k ) ? ( m_mu - center ) : this->invalid_moment () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of entries
      inline size_type  size  () const { return m_prev.size ()  ; }
      /// get effective number of entries 
      inline size_type  nEff  () const { return m_prev.nEff ()  ; }
      // get the mean value
      inline data_type  mu    () const { return m_mu            ; } 
      /// empty ?
      inline bool       empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool       ok    () const { return m_prev.ok    () ; }
      /// minimal value
      inline double     min   () const { return m_min ; } 
      /// maximal value
      inline double     max   () const { return m_max ; } 
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
        if ( !std::isfinite ( x ) ) { return *this ; }
        //
        const size_type n = m_prev.size() ;
        m_mu = ( n * m_mu + x ) /  ( n + 1 );  // calculate new mean value 
        m_prev += x ;                          // updated previous
	      m_min   = std::min ( m_min , x ) ;
	      m_max   = std::max ( m_max , x ) ;
        return *this ;
      }
      /// add the moment 
      inline Moment_& add ( const Moment_& x )
      {
        if      ( x     . empty () ) {               return *this ; }
        else if ( this -> empty () ) { (*this) = x ; return *this ; }
        //
        const size_type n1 =   m_prev.size() ;
        const size_type n2 = x.m_prev.size() ;
        //
        m_mu = ( n1 * m_mu + n2 * x.m_mu ) / ( n1 + n2 ) ; // update mean 
        m_prev += x.m_prev ;                               // update previous
	      //
	      m_min = std::min ( m_min , x.m_min ) ;
	      m_max = std::max ( m_max , x.m_max ) ;	
        //
        return *this ;
      }
      /// add sequence of values  
      template <class ITERATOR>
      inline
      Moment_&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
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
      /// Reset the content
      void reset  () override
      {
      	m_mu  = 0 ;
	      m_min =   std::numeric_limits<double>::max () ;
	      m_max = - std::numeric_limits<double>::max () ;
	      m_prev.reset() ;
      }
    public :
      // ======================================================================      
      /// is finite ?
      inline bool isfinite () const
      { return
	      std::isfinite ( m_mu  ) &&
	      std::isfinite ( m_min ) && 
	      std::isfinite ( m_max ) && m_prev.isfinite () ;
      }
      // ======================================================================
    public:      
      // ======================================================================
      /// get "previos" moment 
      inline const Moment_<0>& previous () const { return this->m_prev ; }
      // ======================================================================
    public: 
      // ======================================================================
      inline long double M ( const unsigned short k = 1 ) const
      { return 0 == k ? this->m_prev. M ( k ) : 1 == k ? 0.0 : this->invalid_moment () ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( Moment_& right )
      {
        m_prev.swap ( right.m_prev ) ;
        std::swap ( m_mu  , right.m_mu  ) ;
        std::swap ( m_min , right.m_min ) ;
        std::swap ( m_max , right.m_max ) ;
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
      inline data_type M_ () const { return 0 ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moments
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline data_type moment_ () const { return 0 ; }
      // =======================================================================
    public: // templated std_moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline data_type std_moment_ () const { return 0 ; }
      // ======================================================================
    public: // templated centralized moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double /* center */ ) const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center ) const { return mu() - center ; }
      // ======================================================================
    public:
      // ======================================================================
      /// 1st cumulant 
      template <unsigned int K, typename std::enable_if<(1==K),int>::type = 0 >
      inline data_type cumulant_ () const { return this->mu() ; }
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
      data_type   m_mu   { 0 } ; // mean value
      /// minimal value
      double      m_min  {   std::numeric_limits<double>::max () } ;
      /// maximal value
      double      m_max  { - std::numeric_limits<double>::max () } ;
      // ======================================================================
    private:
      // ======================================================================
      /// array of binomial coefficients 
      static const std::array<unsigned long long,2> s_Ck ;
      // ======================================================================
    } ;
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
      /// # of entries 
      typedef typename Moment_<0>::size_type size_type ;
      /// data type 
      typedef typename Moment_<0>::data_type data_type ; 
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
      ( const Moment_<N-1>& prev ,	
	      const double        mom  ) 
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
       *  @param k  the central moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const
      { return
          0 == k        ? 1.0 :
          1 == k        ? 0.0 :
          N >  k        ? m_prev.moment     ( k ) :
          N <  k        ? this->invalid_moment () :
          !this->ok ()  ? this->invalid_moment () :
          this->M ( k ) / this->size() ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { 
        return
          0 == k       ? 1.0 : 
          1 == k       ? 0.0 :
          2 == k       ? 1.0 :
          N >  k       ? m_prev.std_moment ( k ) : 
          N <  k       ? this->invalid_moment () :  
          !this->ok () ? this->invalid_moment () :
          this->moment ( k ) / std::pow ( this->moment ( 2 ) , 0.5 * k ) ;
      } 
      /// get value of the kth centralized moment for \f$  k \le N \f$
      inline double centralized_moment 
      ( const unsigned short k      , 
        const double         center ) const
      { 
        if      ( 0 == k       ) { return 1              ; } 
        else if ( 1 == k       ) { return mu () - center ; }
        else if ( N >  k       ) { return m_prev.centralized_moment ( k , center ) ; }
        else if ( N <  k       ) { return this->invalid_moment () ; }
        else if ( !this->ok () ) { return this->invalid_moment () ; }
        // 
        // here we have k == N
        //
        const data_type delta  = mu () - center ; 
        data_type result = 0 ;
        data_type deltai = 1 ;  // == delta**i 
        for ( unsigned short i = 0 ; i <= N ; ++i )
        {
          result += s_Ck [ i ] * deltai * this -> moment ( N - i ) ;
          deltai *= delta ; 
        }
        return result ; 
      }
      // ======================================================================
    public: 
      // ======================================================================
      /// get number of entries
      inline size_type size  () const { return m_prev.size  () ; }
      /// get effective number of entries 
      inline size_type nEff  () const { return m_prev.nEff  () ; }
      /// get the mean value (if \f$ 1 \le N \f$)
      inline data_type mu    () const { return m_prev.mu    () ; }
      /// empty ?
      inline bool      empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool      ok    () const { return m_prev.ok    () ; }
      // ======================================================================
      /// minimal value
      inline double    min   () const { return m_prev.min   () ; } 
      /// maximal value
      inline double    max   () const { return m_prev.max   () ; } 
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
      inline
      Moment_&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
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
      /// Reset the content
      void reset  () override { m_M = 0 ; m_prev.reset() ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      inline bool isfinite () const
      { return std::isfinite ( m_M ) && m_prev.isfinite () ; }
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
      inline data_type M ( const unsigned short k = N ) const
      { return 
        k <  N ? this->m_prev.M ( k ) :
        k == N ? this->m_M            : this->invalid_moment () ; }
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
      inline data_type M_ () const { return this->m_M ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(N>K),int>::type = 0 >
      inline data_type M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moment
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0)&&(K<=N),int>::type = 0 >
      inline data_type moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(N>=K),int>::type = 0 >
      inline data_type moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(K==N),int>::type = 0 >
      inline data_type moment_ () const
      { return !this->ok () ? this->invalid_moment () : this->m_M / this->size  () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError
      moment_ () const
      {
        if  ( !this -> ok () ) { return this->invalid_moment () ; }
        const size_type n = this->size() ;
        //
        const data_type muo  = this->template M_<K>   () / n ;
        //
        if ( this -> size() < 2 * K ) { return muo ; } // error estimate is impossible 
        //
        const data_type mu2o = this->template M_<2*K> () / n ;
        const data_type muop = this->template M_<K+1> () / n ;
        const data_type muom = this->template M_<K-1> () / n ;  
        const data_type mu2  = this->template M_<2>   () / n ;  
        //
        data_type cov2 = mu2o     ;
        cov2 -= 2 * K * muop * muom ;
        cov2 -=         muo  * muo  ;
        cov2 += K * K * mu2  * muom * muom ;
        cov2 /= n ;
        //
        return Ostap::Math::ValueWithError ( muo , cov2 ) ;
      }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(N<2*K)&&(N>K),int>::type = 0 >
      inline data_type moment_ () const
      { return this->m_prev.template moment_<K> () ; }      
      // =======================================================================
    public: // templated std_moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==2)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError
      std_moment_ () const
      {
        return !this->ok () ? Ostap::Math::ValueWithError ( this->invalid_moment () ) : 
          this->template moment_<K> () / std::pow ( double ( this->moment_<2>() ) , 0.5 * K ) ;
      }      
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(N<2*K)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const
      {
        return !this->ok () ? this->invalid_moment () :
          this->template moment_<K> () / std::pow ( double ( this->moment_<2>() ) , 0.5 * K ) ;
      }  
      // ======================================================================
    public:
      // ======================================================================      
      template <unsigned int K, typename std::enable_if<(K==0)&&(K<=N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double /* center */ ) const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(K<=N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center ) const { return this -> mu() - center ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K<N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center ) const 
      { return this->m_prev.template centralized_moment_<K> ( center ) ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center ) const 
      { 
        if ( !this->ok () ) { return this->invalid_moment() ; }
        //
        const data_type  delta = this -> mu () - center ;
        data_type result = 0 ;
        data_type deltai = 1 ;
        for ( unsigned short i = 0 ; i <= N ; ++ i )
        {
          result += s_Ck [ i ] * deltai * this->moment ( N - i ) ;
          deltai *= delta ;
        }
        return result ; 
      }
      // ======================================================================
    public: // cumulants 
      // ======================================================================      
      /// 1st cumulant 
      template <unsigned int K, typename std::enable_if<(1==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const { return this->mu() ; }
      /// 2nd cumulant 
      template <unsigned int K, typename std::enable_if<(2==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const { return this->template moment_<2> () ; }
      /// 3rd cumulant 
      template <unsigned int K, typename std::enable_if<(3==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const { return this->template moment_<3> ; }
      /// 4th cumulant 
      template <unsigned int K, typename std::enable_if<(4==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const 
      { 
        if ( !this->ok() ) { return this->invalid_moment () ; }
        //
        const data_type m4 = this->template M_<4> () / this->size() ;
        const data_type m2 = this->template M_<2> () / this->size() ;
        //
        return m4 - 3 * m2 * m2 ;
      }
      /// 5th cumulant 
      template <unsigned int K, typename std::enable_if<(5==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const 
      { 
        if ( !this->ok() ) { return this->invalid_moment () ; }
        //
        const data_type m5 = this->template M_<5> () / this->size() ;
        const data_type m3 = this->template M_<3> () / this->size() ;
        const data_type m2 = this->template M_<2> () / this->size() ;
        //
        return m5 - 10  *m3 * m2 ;
      }
      /// 6th cumulant 
      template <unsigned int K, typename std::enable_if<(6==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const 
      { 
        if ( !this->ok() ) { return this -> invalid_moment () ; }
        //
        const data_type m6 = this->template M_<6> () / this->size() ;       
        const data_type m4 = this->template M_<4> () / this->size() ;
        const data_type m3 = this->template M_<3> () / this->size() ;
        const data_type m2 = this->template M_<2> () / this->size() ;
        //
        return m6 - 15 * m4 * m2 - 10 * m3 * m3 + 30 * m2 * m2 * m2 ;
      }
      /// 7th cumulant 
      template <unsigned int K, typename std::enable_if<(7==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const 
      { 
        if ( !this->ok() ) { return this->invalid_moment () ; }
        //
        const data_type m7 = this->template M_<7> () / this->size() ;       
        const data_type m5 = this->template M_<5> () / this->size() ;       
        const data_type m4 = this->template M_<4> () / this->size() ;
        const data_type m3 = this->template M_<3> () / this->size() ;
        const data_type m2 = this->template M_<2> () / this->size() ;
        //
        return m7 - 21 * m5 * m2 - 35 * m4 * m3 + 210 *m3 * m2 * m2 ;
      }
      /// 8th cumulant 
      template <unsigned int K, typename std::enable_if<(8==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const 
      { 
        if ( !this->ok() ) { return this->invalid_moment () ; }
        //
        const data_type m8 = this->template M_<8> () / this->size() ;       
        const data_type m6 = this->template M_<6> () / this->size() ;       
        const data_type m5 = this->template M_<5> () / this->size() ;       
        const data_type m4 = this->template M_<4> () / this->size() ;
        const data_type m3 = this->template M_<3> () / this->size() ;
        const data_type m2 = this->template M_<2> () / this->size() ;
        //
        return m8 - 28 * m6 * m2 - 56 * m5 * m3 - 35 * m4 * m4 
        + 420 * m4 * m2 * m2 + 560 * m3 * m3 * m2 - 630 * m2 * m2 * m2 * m2 ;        
      }
      /// 9th cumulant 
      template <unsigned int K, typename std::enable_if<(9==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const 
      { 
        if ( !this->ok() ) { return this->invalid_moment () ; }
        //
        const data_type m9 = this->template M_<9> () / this->size() ;       
        const data_type m7 = this->template M_<7> () / this->size() ;       
        const data_type m6 = this->template M_<6> () / this->size() ;       
        const data_type m5 = this->template M_<5> () / this->size() ;       
        const data_type m4 = this->template M_<4> () / this->size() ;
        const data_type m3 = this->template M_<3> () / this->size() ;
        const data_type m2 = this->template M_<2> () / this->size() ;
        //
        return m9 - 36 * m7 * m2 - 84 * m6 * m3 - 126 * m5 * m4 
        + 756 * m5 * m2 * m2 + 2520 * m4 * m3 * m2 
        + 560 * m3 * m3 * m3 - 7560 * m3 * m2 * m2 * m2 ; // NB? typo here? 
        
      }
      /// 10th cumulant 
      template <unsigned int K, typename std::enable_if<(10==K)&&(K<=N),int>::type = 0 >
      inline data_type cumulant_ () const 
      { 
        if ( !this->ok() ) { return this->invalid_moment () ; }
        //
        const data_type m10 = this->template M_<10> () / this->size() ;       
        const data_type m8  = this->template M_<8>  () / this->size() ;       
        const data_type m7  = this->template M_<7>  () / this->size() ;       
        const data_type m6  = this->template M_<6>  () / this->size() ;       
        const data_type m5  = this->template M_<5>  () / this->size() ;       
        const data_type m4  = this->template M_<4>  () / this->size() ;
        const data_type m3  = this->template M_<3>  () / this->size() ;
        const data_type m2  = this->template M_<2>  () / this->size() ;
        //
        return m10 - 45 * m8 * m2 - 120 * m7 * m3 - 210 * m6 * m4 
        + 1260  * m6 * m2 * m2       - 126   * m5 * m5
        + 5040  * m5 * m3 * m2       + 3150  * m4 * m4 * m2 
        + 4200  * m4 * m3 * m3       - 18900 * m4 * m2 * m2 * m2  
        - 37800 * m3 * m3 * m2 * m2  + 22680 * m2 * m2 * m2 * m2 * m2 ;
      }
      // ======================================================================
    public:
      // ======================================================================
      ///  get the mean value with error estimate  
      inline Ostap::Math::ValueWithError mean () const 
      { return Ostap::Math::ValueWithError ( mu () , this->template moment_<2>() ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// counter of (N-1)th order
      Moment_<N-1> m_prev {   } ; // counter of (N-1)th order
      /// the current value of \f$ m_N \f$  
      data_type    m_M    { 0 } ; // the current value
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
      if ( !std::isfinite ( x ) ) { return *this ; } // NB: NO ACTION!!
      //
      const size_type nA = this->size() ;
      const size_type nB =      1 ;
      const size_type nN = nA + 1 ;
      //
      const data_type delta = x - this -> mu () ;
      const data_type b_n   =      -1.0L / nN ;
      const data_type a_n   =  nA * 1.0L / nN ;
      const data_type d_n   = - delta    / nN ;
      //
      m_M += ( nA * std::pow ( b_n , N ) + std::pow ( a_n , N ) ) * std::pow ( delta , N ) ;
      data_type d = 1 ;
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
      const size_type nA    = this->size () ;
      const size_type nB    =    x. size () ;
      const size_type nN    = nA + nB       ;
      //
      const data_type delta = x.mu() - this->mu() ;
      const data_type b_n   =  ( -1.0L * nB ) / nN ;
      const data_type a_n   =  (  1.0L * nA ) / nN ;
      //
      m_M += x.m_M ;
      m_M += nA * std::pow ( b_n * delta , N ) + nB * std::pow ( a_n * delta , N ) ;
      //
      data_type a = 1 ;
      data_type b = 1 ;
      data_type d = 1 ;
      //
      for ( unsigned short k = 1 ; k + 2 <= N ; ++k )
	    {
	      a   *= a_n   ;
	      b   *= b_n   ;
	      d   *= delta ;        
	      m_M += s_Ck [ k ] * d * ( this-> M ( N -k ) * b + x. M ( N - k ) * a ) ;
	    }
      /// update previous 
      this->m_prev += x.m_prev ; // update previous
      //
      return *this ;
    }
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
    protected:
      // ======================================================================
      /// return the value for invalid moment 
      double invalid_moment () const ;
      // ======================================================================
    } ;
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
    template <unsigned short N> class WMoment_    ;
    // ========================================================================    
    /// specialization for \f$ N=0 f\$
    template <>
    class WMoment_<0> : public WMoment 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// # of entries 
      typedef typename Moment_<0>::size_type size_type ;
      /// data type 
      typedef typename Moment_<0>::data_type data_type ; 
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 0 } ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor 
      WMoment_ () = default ;
      // ======================================================================
      /** full constructor
       *  @param size number of entries 
       *  @param sumw sum of weights 
       *  @param sumw sum of squared weights 
       *  @param wmin minimal weight 
       *  @param wmax maximal weight 
       */
      WMoment_
      ( const size_type size  ,
	      const double    sumw  ,
	      const double    sumw2 , 
	      const double    wmin  , 
	      const double    wmax  ) ;
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
       *  @param k the central moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const 
      { return 0 == k ? 1.0 : this->invalid_moment () ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 0 == k ? 1.0 : this->invalid_moment () ; }
      // ======================================================================
      /// get value of the kth centralized moment for \f$  k \le N \f$
      inline double centralized_moment
      ( const unsigned short k         , 
        const double      /* center */ ) const
      { return 0 == k ? 1.0 : this->invalid_moment () ; }
      // ======================================================================
    public :
      // ======================================================================
      /// get number of entries
      inline size_type size  () const { return m_size ; }
      /// get effective number of entries \f$  \frac{(\sum w_i)^2}{\sum w_i^2} \f$
      inline data_type nEff  () const { return m_w2 ? m_w * m_w / m_w2 : -1.0 ; }
      /// get sum of weights  \f$  \sum w_i \f$  
      inline data_type w     () const { return m_w    ; }
      /// get sum of weights squared  \f$  \sum w_i^2 \f$  
      inline data_type w2    () const { return m_w2   ; }
      /// empty ?
      inline bool      empty () const { return 0 == m_size ; }
      /// ok ?
      inline bool      ok    () const { return m_size && m_w && m_w2 ; }
      /// minimal weight 
      inline double    wmin  () const { return m_wmin ; } 
      /// maximal weight 
      inline double    wmax  () const { return m_wmax ; }       
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
      ( const double x     ,
        const double w = 1 )
      {
        /// INVALID values&weight and ZERO weight are ignred 
        if ( !w || !std::isfinite ( x ) || !std::isfinite ( w ) ) { return *this ;}
        //
        ++m_size ;
        //
        m_w   += w     ; 
        m_w2  += w * w ;
        //
        m_wmin = std::min ( m_wmin , w ) ;
        m_wmax = std::max ( m_wmax , w ) ;
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
        m_wmin = std::min ( m_wmin , x.m_wmin ) ;
        m_wmax = std::max ( m_wmax , x.m_wmax ) ;
        //
        return *this ;
      }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WMoment_&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
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
      void update
      ( const double x     ,
	      const double w = 1 ) override { add ( x , w ) ; }
      /// Reset the content
      void reset  () override
      {
	      m_size = 0 ;
	      m_w    = 0 ;
	      m_w2   = 0 ;
	      m_wmin =   std::numeric_limits<double>::max () ;
	      m_wmax = - std::numeric_limits<double>::max () ;
      }
      // ======================================================================      
    public :
      // ======================================================================      
      /// is finite ?
      inline bool isfinite () const
      { return
	        std::isfinite ( m_w    ) &&
	        std::isfinite ( m_w2   ) &&
	        std::isfinite ( m_wmin ) &&
	        std::isfinite ( m_wmax ) ; 
      }
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of \f$ M_k = \sum_i w_i \left( x_i - \bar{x} \right)^k \f$ 
       *  for \f$ k \le N \f$
       *  @param k the order   \f$  0 \le k \le N \f$
       *  @return the value of \$ M_k \f$  if \f$  0 \le k \le  \f$
       */   
      inline data_type M ( const unsigned short k = 0 ) const
      { return 0 == k ? m_w : this->invalid_moment () ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( WMoment_& right )
      {
        std::swap ( m_size , right.m_size ) ;
        std::swap ( m_w    , right.m_w    ) ;
        std::swap ( m_w2   , right.m_w2   ) ;
        std::swap ( m_wmin , right.m_wmin ) ;
        std::swap ( m_wmax , right.m_wmax ) ;
      }
      // ======================================================================
    public: //  templated M_ 
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type M_ () const { return this->m_w ; }
      // ======================================================================
    public:  // templated moment
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type moment_ () const { return 1 ; }
      // ======================================================================
    public: // templated std_moment 
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type std_moment_ () const { return 1 ; }
      // ======================================================================
    public: // templated centralized moment 
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double /* center */ ) const { return 1 ; }
      // ======================================================================
    private:
      // ======================================================================
      /// number of entries 
      size_type  m_size { 0 } ; // number of entries
      /// sum of weights \f$  \sum w_i \f$ 
      data_type  m_w    { 0 } ; // sum of weights
      /// sum of weights squared \f$  \sum w_i^2 \f$ 
      data_type  m_w2   { 0 } ; // sum of weights squared 
      // ======================================================================
      /// minimal weight 
      double     m_wmin {   std::numeric_limits<double>::max () } ;
      /// maximal weight 
      double     m_wmax { - std::numeric_limits<double>::max () } ;
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
      /// # of entries 
      typedef typename WMoment_<0>::size_type size_type ;
      /// data type 
      typedef typename WMoment_<0>::data_type data_type ; 
      // ======================================================================
    public:
      // ======================================================================
      enum _ { order = 1 } ;
      // ======================================================================
    public: 
      // ======================================================================
      /// default constructor 
      WMoment_ () = default ;
      // ======================================================================
      /// constructor from mu and previous moment 
      WMoment_
      ( const WMoment_<0>& prev ,
	      const double       mu   ,
	      const double       xmin ,
	      const double       xmax ) ;
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
      inline double moment ( const unsigned short k ) const
      { return 0 == k ? 1.0 : 1 == k ? 0.0 : this -> invalid_moment () ; }
      // ======================================================================
      /// get value of the kth standartized moment for \f$  k \le N \f$
      inline double std_moment ( const unsigned short k ) const
      { return 0 == k ? 1.0 : 1 == k ? 0.0 : this -> invalid_moment ()  ; }
      // ======================================================================
      /// get value of the kth centralized moment for \f$  k \le N \f$
      inline double centralized_moment 
      ( const unsigned short k      ,
        const double         center ) const
      { return 0 == k ? 1.0 : 1 == k ? mu () - center : this -> invalid_moment () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of entries
      inline size_type size  () const { return m_prev.size ()  ; }
      /// get effective number of entries \f$  \frac{(\sum w_i)^2}{\sum w_i^2} \f$
      inline data_type nEff  () const { return m_prev.nEff () ; }
      /// get sum of weights \f$ \sum w_i \f$
      inline data_type w     () const { return m_prev.w    () ; }
      /// get sum of weights squared 
      inline data_type w2    () const { return m_prev.w2   () ; }
      // get the mean value
      inline data_type mu    () const { return m_mu           ; } 
      /// empty ?
      inline bool      empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool      ok    () const { return m_prev.ok    () ; }
      /// minimal value
      inline double    min   () const { return m_min ; } 
      /// maximal value
      inline double    max   () const { return m_max ; }       
      // ======================================================================
      /// minimal weight 
      inline double    wmin  () const { return m_prev.wmin () ; } 
      /// maximal weight 
      inline double    wmax  () const { return m_prev.wmax () ; }
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
        /// INVALUID values and weights are ignored 
        if ( !std::isfinite ( x ) || !std::isfinite ( w ) ) { return *this ; }
        /// ZERO weights are ignored 
        if ( !w )                                           { return *this ; }  
        //
        const data_type wA = this -> w() ;
        const data_type wB = w ;
        //`
        m_mu = ( wA * m_mu + wB * x ) / ( wA + wB ) ; // update mean
        // minimal 
        m_min = std::min ( m_min , x ) ;
        // maximal
        m_max = std::max ( m_max , x ) ;
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
        const data_type wB = x.m_prev.w () ;
        if ( !wB ) { return *this ; }
        //
        const data_type wA = m_prev.w   () ;
        m_mu = ( wA * m_mu + wB * x.m_mu ) / ( wA + wB ) ; // update mean
        //
        // minimal 
        m_min = std::min ( m_min , x.m_min ) ;
        // maximal
        m_max = std::max ( m_max , x.m_max ) ;
        //
        m_prev += x.m_prev ;                                 // update previous
        //
        return *this ;
      }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WMoment_&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
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
      void update
      ( const double x     ,
        const double w = 1 ) override { add ( x , w ) ; }
      /// Reset the content
      void reset  () override
      {
        m_mu = 0 ;
        m_min =   std::numeric_limits<double>::max () ;
        m_max = - std::numeric_limits<double>::max () ;
        m_prev.reset () ;
      }
      // ======================================================================      
    public :
      // ======================================================================      
      /// is finite ?
      inline bool isfinite () const
      { return
          std::isfinite ( m_mu  ) &&
          std::isfinite ( m_min ) &&
          std::isfinite ( m_max ) && m_prev.isfinite () ;
      }      
      // ======================================================================
    public:      
      // ======================================================================
      /// get "previos" moment 
      const WMoment_<0>& previous () const { return this->m_prev ; }
      // ======================================================================
    public: 
      // ======================================================================
      inline data_type M ( const unsigned short k = 1 ) const
      { return 
          0 == k ? this->m_prev.M ( k ) :
          1 == k ? 0.0                  :
          this->invalid_moment () ; }
      // ======================================================================
    public:
      // ======================================================================
      void swap ( WMoment_& right )
      {
        m_prev.swap ( right.m_prev ) ;
        std::swap ( m_mu  , right.m_mu  ) ;
        std::swap ( m_min , right.m_min ) ;
        std::swap ( m_max , right.m_max ) ;
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
      inline data_type M_ () const { return 0 ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moments
      // ======================================================================
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline data_type moment_ () const { return 0 ; }
      // =======================================================================
    public: // templated std_moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline data_type std_moment_ () const { return 0 ; }
      // =======================================================================
    public: // templated centralized moment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double /* center */ ) const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center ) const { return this -> mu () - center ; }
      // ======================================================================
    public:
      // ======================================================================
      ///  get the mean value with error estimate  
      inline double mean () const{ return mu () ; }  
      // =======================================================================
    private:
      // ======================================================================
      WMoment_<0>  m_prev {   } ;
      /// mean value
      data_type m_mu      { 0 } ; // mean value
      // ======================================================================
      /// minimal value
      double    m_min  {   std::numeric_limits<double>::max () } ;
      /// maximal value
      double    m_max  { - std::numeric_limits<double>::max () } ;
      // ======================================================================
    private:
      // ======================================================================
      /// array of binomial coefficients 
      static const std::array<unsigned long long,2> s_Ck ;
      // ======================================================================
    } ;
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
      /// # of entries 
      typedef typename WMoment_<0>::size_type size_type ;
      /// data type 
      typedef typename WMoment_<0>::data_type data_type ; 
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
      ( const WMoment_<N-1>& prev , 
      	const double         mom  ) 
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
      inline double moment () const 
      { return this-> ok () ? this->M ( N ) / this-> w () : this->invalid_moment () ; }
      // ======================================================================
      /** get value of the kth moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{\sum w_i} \sum w_i \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double moment ( const unsigned short k ) const
      { return
          0 == k       ? 1.0 :
          1 == k       ? 0.0 :
          N >  k       ? this->m_prev.moment ( k ) :
          N <  k       ? this->invalid_moment   () :
          !this->ok () ? this->invalid_moment   () :
	      this->M ( k ) / this-> w () ; }
      // ======================================================================
      /** get value of the kth standartized moment for \f$  k \le N \f$
       *  \f[ \mu_k \equiv  \frac{1}{N} \sum \left( x_i - \bar{x} \right)^k \f]
       *  for \f$  k > N\f$ null is returned 
       *  @param k  the cenral moment order  \f$  0 \le k \le N \f$
       *  @return the value of the kth central moment if \f$  0 \le k \le N \f$, 0, otherwise 
       */
      inline double std_moment ( const unsigned short k ) const
      { return 
          0 == k       ? 1.0 :
          1 == k       ? 0.0 :
          2 == k       ? 1.0 : 
          N >  k       ? this->std_moment    ( k ) :
          N <  k       ? this->invalid_moment   () :
          !this->ok () ? this->invalid_moment   () :
          this->moment ( k ) / std::pow ( this->moment ( 2 ) , 0.5 * k ) ;
      }
      // ======================================================================
      /// get value of the kth centralized moment 
      inline double centralized_moment 
      ( const unsigned short k      , 
        const double         center ) const
      {
        if      ( 0 == k        ) { return 1 ; }
        else if ( 1 == k        ) { return mu () - center ; }
        else if ( N >  k        ) { return m_prev.centralized_moment ( k , center ) ; }
        else if ( N <  k        ) { return this->invalid_moment () ; }
        else if ( !this-> ok () ) { return this->invalid_moment () ; }
        //
        const data_type delta = mu () - center ;
        data_type result = 0 ;
        data_type deltai = 1 ; //  == delta**i 
        for ( unsigned short i = 0 ; i <= N ; ++i )
        {
          result += s_Ck [ i ] * deltai *  this->moment ( N - i ) ;
          deltai *= delta ; 
        }
        return result ;  
      }
      // ======================================================================
    public :
      // ======================================================================
      /// get number of entries
      inline size_type size  () const { return m_prev.size () ; }
      /// get effective number of entries \f$  \frac{(\sum w_i)^2}{\sum w_i^2} \f$
      inline data_type nEff  () const { return m_prev.nEff () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline data_type w     () const { return m_prev.w    () ; }
      /// get sum of weights squared 
      inline data_type w2    () const { return m_prev.w2   () ; }
      /// get the weighted mean value (if \f$ 1 \le N \$4)
      inline data_type mu    () const { return m_prev.mu   () ; } 
      /// empty ?
      inline bool      empty () const { return m_prev.empty () ; }
      /// ok ?
      inline bool      ok    () const { return m_prev.ok    () ; }
      // ======================================================================
      /// minimal value
      inline double    min  () const { return m_prev.min   () ; } 
      /// maximal value
      inline double    max  () const { return m_prev.max   () ; } 
      // ======================================================================
      /// minimal weight 
      inline double    wmin  () const { return m_prev.wmin () ; } 
      /// maximal weight 
      inline double    wmax  () const { return m_prev.wmax () ; }
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
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WMoment_&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
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
      void update
      ( const double x     ,
        const double w = 1 ) override { add ( x , w ) ; }
      /// Reset the content
      void reset  () override
      {
        m_M = 0 ;
        m_prev.reset() ;
      }      
      // ======================================================================      
    public :
      // ======================================================================      
      /// is finite ?
      inline bool isfinite () const
      { return std::isfinite ( m_M ) && m_prev.isfinite () ;}
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
       *  @return the value of \$ M_k \f$  if \f$  0 \le k \le N \f$ 
       */   
      inline data_type M ( const unsigned short k = N ) const
      { return 
          N > k ? this->m_prev.M    ( k ) :
          N < k ? this->invalid_moment () : this -> m_M ; }
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
      inline data_type M_ () const { return this->m_M ; }
      // =====================================================================
      template <unsigned  int K, typename std::enable_if<(N>K) ,int>::type = 0 >
      inline data_type M_ () const { return this->m_prev.template M_<K>() ; }
      // ======================================================================
    public:  // templated moment
      // ======================================================================      
      /// get the central moment of order K 
      template <unsigned int K, typename std::enable_if<(K==0)&&(K<=N),int>::type = 0 >
      inline data_type moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(K<=N),int>::type = 0 >
      inline data_type moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(K==N),int>::type = 0 >
      inline data_type moment_ () const
      { return !this->ok() ? this->invalid_moment () : this->m_M / this->w () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(N<2*K)&&(N>K),int>::type = 0 >
      inline data_type moment_ () const
      { return this->m_prev.template moment_<K> () ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<=K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError
      moment_ () const
      {
        //
        if ( !this->ok() ) { return this->invalid_moment() ; }
        //
		    const long double n = this->w() ; // ATENTION!
        //
        const data_type muo  = this->template M_ <K>  () / n ;
        const data_type mu2o = this->template M_<2*K> () / n ;
        const data_type muop = this->template M_<K+1> () / n ;
        const data_type muom = this->template M_<K-1> () / n ;  
        const data_type mu2  = this->         M_<2>   () / n ;  
        //
        data_type cov2 = mu2o     ;
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
      template <unsigned int K, typename std::enable_if<(K==0)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const { return 0 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==2)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(2<K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError
      std_moment_ () const
      {
        return !this->ok() ? Ostap::Math::ValueWithError ( this->invalid_moment () ) : 
          this->template moment_<K> () / std::pow ( double ( this->moment_<2>() ) , 0.5 * K ) ;
      }      
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(N<2*K)&&(K<=N),int>::type = 0 >
      inline data_type std_moment_ () const
      {
        return !this->ok() ? this->invalid_moment () :
          this->template moment_<K> () / std::pow ( double ( this->moment_<2>() ) , 0.5 * K ) ;
      }      
   // =======================================================================
    public: // templated centralzedmoment 
      // =======================================================================
      template <unsigned int K, typename std::enable_if<(K==0)&&(K<=N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double /* center */ ) const { return 1 ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==1)&&(K<=N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center ) const { return this->mu() - center  ; }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K<N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center  ) const { return this->m_prev.template centralized_moment_<K> ( center ); }
      // ======================================================================
      template <unsigned int K, typename std::enable_if<(K==N),int>::type = 0 >
      inline data_type centralized_moment_ 
      ( const double center ) const
      {
        if ( !this->ok () ) { return this->invalid_moment () ; }
        //
        const data_type delta = mu () - center ;
        data_type result = 0 ;
        data_type deltai = 1 ; //  == delta**i 
        for ( unsigned short i = 0 ; i <= N ; ++i )
          {
            result += s_Ck [ i ] * deltai *  this->moment ( N - i ) ;
            deltai *= delta ; 
          }
        return result ;  
      } 
       // ======================================================================
    public:
      // ======================================================================
      ///  get the mean value with error estimate  
      inline Ostap::Math::ValueWithError mean () const 
      { return Ostap::Math::ValueWithError ( mu () , this->template moment_<2>() ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// counter of (N-1)th order
      WMoment_<N-1> m_prev {}   ;  // counter of (N-1)th order
      /// the current value of \f$ M_N \f$  
      data_type     m_M   {0}  ; // the current value
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
      /// invalid value aand weights are ignored 
      if ( !std::isfinite ( x ) || !std::isfinite ( w ) ) { return *this ; }
      /// ZERO weights are also ignored 
      if ( !w )                                           { return *this ; }
      //
      const data_type wA    = this->w () ;
      const data_type wB    =       w    ;
      //
      const data_type wW    = wA + wB    ;
      const data_type delta = x - this->mu() ;
      //
      const data_type b_n =  -1.0L * wB / wW ;
      const data_type a_n =   1.0L * wA / wW ;
      const data_type d_n = - delta     / wW ;
      //
      m_M += ( wA * std::pow ( b_n , N ) + wB * std::pow ( a_n , N ) ) * std::pow ( delta , N ) ;
      data_type d = 1 ;
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
      const data_type wB    =  x. w () ;
      if ( !wB ) {               return *this ; }                    // RETURN! 
      //
      const data_type wA    = this->w () ;
      //
      const data_type wW    = wA + wB    ;
      const data_type delta = x.mu () - this->mu ()  ;
      const data_type b_n   =  ( -1.0L * wB ) / wW ;
      const data_type a_n   =  (  1.0L * wA ) / wW ;
      //
      m_M += x.m_M ;
      m_M += wA * std::pow ( b_n * delta , N ) + wB * std::pow ( a_n * delta , N ) ;
      //
      data_type a = 1 ;
      data_type b = 1 ;
      data_type d = 1 ;
      //
      for ( unsigned short k = 1 ; k + 2 <= N ; ++k )
      {
        a   *= a_n   ;
        b   *= b_n   ;
        d   *= delta ;        
        m_M += s_Ck [ k ] * d * ( this -> M ( N - k ) * b + x. M ( N - k ) * a ) ;
      }
      /// update previous 
      this->m_prev += x.m_prev ; // update previous
      //
      return *this ;
    }
    // =======================================================================

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
      typedef Moment_<2>         Counter   ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
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
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }
      inline Ostap::Math::ValueWithError value () const { return Ostap::Math::pow ( 2 , m_log.mean() ) ; }   
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
      inline
      GeometricMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      /// reset the cunter 
      void reset  () override { m_log.reset () ; }
      // ======================================================================      
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_log.size  () ; }
      /// empty ?
      inline bool      empty () const { return m_log.empty () ; } 
      /// ok ?
      inline bool      ok    () const { return m_log.ok    () ; } 
      // ======================================================================      
    private:
      // ======================================================================
      /// get the counter of log2(x) 
      Counter m_log {} ;
      // ======================================================================
    };
    // ========================================================================
    /** @class HarmonicMean 
     *  Calculate the harmonic mean 
     *  \f$ \frac{n}{ \frac{1}{x_1} + ... + \frac{1}{x_n}} \f$
     *  @see https://en.wikipedia.org/wiki/Harmonic_mean
     */
    class HarmonicMean : public Statistic
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef Moment_<2> Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
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
      inline Ostap::Math::ValueWithError value () const { return 1. / m_inv.mean() ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }
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
      /// add two counters together 
      inline HarmonicMean& add ( const HarmonicMean& x )
      {
        m_inv.add ( x.m_inv ) ;
        return *this ;
      }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline
      HarmonicMean& 
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of 1/x values 
      const Counter& counter () const { return m_inv ; }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      /// reset the cunter 
      void reset  () override { m_inv.reset () ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_inv.size  () ; }
      /// empty ?
      inline bool      empty () const { return m_inv.empty () ; } 
      /// ok ?
      inline bool      ok    () const { return m_inv.ok    () ; } 
      // ======================================================================      
    private:
      // ======================================================================
      /// get the counter of 1/x 
      Counter m_inv {} ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PowerMean 
     *  Calculate the power mean 
     *  \f$ \left(\frac{1}{n}\sum x_i^p \right)^{\frac{1}{p}}\f$
     *  @see https://en.wikipedia.org/wiki/Power_mean
     */
    class PowerMean : public Statistic 
    {
    public:
      // ======================================================================
      /// the actual underlying counter 
      typedef Moment_<2> Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;      
      // ======================================================================
    public:
      // ======================================================================
      /// default construictor
      PowerMean ( const double p = 1 ) ;
      /// constructor for the counters 
      PowerMean ( const double p , const Counter& cnt ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the power mean 
      inline Ostap::Math::ValueWithError value () const 
      { return Ostap::Math::pow ( m_pow.mean() , 1 / m_p ) ; }
      inline Ostap::Math::ValueWithError mean  () const 
      { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline PowerMean& operator+=( const double     x ) { return add ( x ) ; }
      inline PowerMean& operator+=( const PowerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      const Counter& counter () const { return m_pow ; }
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
      inline
      PowerMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      /// reset the cunter 
      void reset  () override { m_pow.reset () ; } 
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
   public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_pow.size  () ; }
      /// empty ?
      inline bool      empty () const { return m_pow.empty () ; } 
      /// ok ?
      inline bool      ok    () const { return m_pow.ok    () ; }
      /// power 
      inline double    p     () const { return m_p            ; }
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
      typedef Moment_<2> Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
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
      inline Ostap::Math::ValueWithError value () const 
      { return m_lp.mean() / m_lpm1.mean () ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline LehmerMean& operator+=( const double      x ) { return add ( x ) ; }
      inline LehmerMean& operator+=( const LehmerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      const Counter& counter1 () const { return m_lp   ; }
      /// the counter of x**(p-1)  values 
      const Counter& counter2 () const { return m_lpm1 ; }
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
      inline
      LehmerMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      /// reset the cunter 
      void reset  () override { m_lp.reset () ; m_lpm1.reset () ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_lp.size   () ; }
      /// empty ?
      inline bool      empty () const { return m_lp.empty  () ; } 
      /// ok ?
      inline bool      ok    () const { return m_lp.ok     () ; }
      /// power 
      inline double    p     () const { return m_p            ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the power
      double  m_p      {1} ;
      /// get the counter of x^p
      Counter m_lp     {}  ;
      /// get the counter of x^(p-1)
      Counter m_lpm1   {}  ;
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
      typedef WMoment_<2> Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
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
      inline Ostap::Math::ValueWithError value () const 
      { return Ostap::Math::pow ( 2 , m_log.mean() ) ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }
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
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WGeometricMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
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
      /// reset the cunter 
      void reset  () override { m_log.reset () ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_log.size  () ; }
      /// number of effective entries
      inline data_type nEff  () const { return m_log.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline data_type  w     () const { return m_log.w     () ; }
      /// get sum of weights squared 
      inline data_type  w2    () const { return m_log.w2    () ; }
      /// empty ?
      inline bool       empty () const { return m_log.empty () ; } 
      /// ok ?
      inline bool       ok    () const { return m_log.ok    () ; } 
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
      typedef WMoment_<2> Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
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
      inline Ostap::Math::ValueWithError value () const { return 1. / m_inv.mean() ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline WHarmonicMean& operator+=( const WHarmonicMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of 1/x values 
      const Counter& counter () const { return m_inv ; }
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
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WHarmonicMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a WMoment 
      void update
      ( const double x ,
        const double w ) override { add ( x , w ) ; }
      /// reset the cunter 
      void reset  () override { m_inv.reset () ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_inv.size  () ; }
      /// number of effective entries
      inline data_type nEff  () const { return m_inv.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline data_type w     () const { return m_inv.w     () ; }
      /// get sum of weights squared 
      inline data_type w2    () const { return m_inv.w2    () ; }
      /// empty ?
      inline bool      empty () const { return m_inv.empty () ; } 
      /// ok ?
      inline bool      ok    () const { return m_inv.ok    () ; } 
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
      typedef WMoment_<2> Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor
      WPowerMean ( const double = 1 ) ;
      /// constructor from the counter 
      WPowerMean ( const double p , const Counter& cnt ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the weighted power mean 
      inline Ostap::Math::ValueWithError value () const 
      { return Ostap::Math::pow ( m_pow.mean() , 1 / m_p ) ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline WPowerMean& operator+=( const WPowerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      const Counter& counter () const { return m_pow ; }
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
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WPowerMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update
      ( const double x ,
        const double w ) override { add ( x , w ) ; }
      /// reset the cunter 
      void reset  () override { m_pow.reset () ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_pow.size  () ; }
      /// number of effective entries
      inline data_type nEff  () const { return m_pow.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline data_type w     () const { return m_pow.w     () ; }
      /// get sum of weights squared 
      inline data_type w2    () const { return m_pow.w2    () ; }
      /// empty ?
      inline bool      empty () const { return m_pow.empty () ; } 
      /// ok ?
      inline bool      ok    () const { return m_pow.ok    () ; }
      /// power 
      inline double    p     () const { return m_p            ; }
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
      typedef WMoment_<2> Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
      // ======================================================================
    public:
      // ======================================================================
      /// default construictor
      WLehmerMean
      ( const double p = 1 ) ;
      /// constructor for the counters 
      WLehmerMean
      ( const double   p    ,
        const Counter& cnt1 , 
        const Counter& cnt2 ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the Lehmer  mean 
      inline Ostap::Math::ValueWithError value () const { return m_lp.mean () / m_lpm1.mean () ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; } 
      // ======================================================================
    public:
      // ======================================================================
      inline WLehmerMean& operator+=( const WLehmerMean& x ) { return add ( x ) ; }   
      // ======================================================================
    public:
      // ======================================================================
      /// the counter of x**p  values 
      inline const Counter& counter1 () const { return m_lp   ; }
      /// the counter of x**(p-1)  values 
      inline const Counter& counter2 () const { return m_lpm1 ; }
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
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WLehmerMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// to use it as a WMoment 
      void update
      ( const double x ,
        const double w ) override { add ( x , w ) ; }
      /// reset the cunter 
      void reset  () override { m_lp.reset() ; m_lpm1.reset() ; } 
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_lp.size   () ; }
      /// number of effective entries
      inline data_type nEff  () const { return m_lp.nEff   () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline data_type w     () const { return m_lp.w      () ; }
      /// get sum of weights squared 
      inline data_type w2    () const { return m_lp.w2     () ; }
      /// empty ?
      inline bool      empty () const { return m_lp.empty  () ; } 
      /// ok ?
      inline bool      ok    () const { return m_lp.ok     () ; }
      /// power 
      inline double    p     () const { return m_p            ; }
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
    class ArithmeticMean : public Statistic
    {
    public:
      // ======================================================================
      typedef Moment_<2>  Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
      // ======================================================================
    public:
      // ======================================================================
      ArithmeticMean () = default ;
      ArithmeticMean ( const Counter& cnt ) ; 
      // ======================================================================
    public :
      // ======================================================================
      /// accumulate only positive entries 
      inline ArithmeticMean& add ( const double          x )
      { m_cnt.add ( x       ) ; return *this ; } 
      /// add two counters togather
      inline ArithmeticMean& add ( const ArithmeticMean& x )
      { m_cnt.add ( x.m_cnt ) ; return *this ; }
      // ======================================================================
      /// add two counters together
      inline ArithmeticMean& operator+=( const double          x )
      { return add ( x ) ; }
      /// add two counters together
      inline ArithmeticMean& operator+=( const ArithmeticMean& x )
      { return add ( x ) ; }
      // ======================================================================
     /// add sequence of values  
      template <class ITERATOR>
      inline
      ArithmeticMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ====================================================================== 
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_cnt.size   () ; }
      /// empty ?
      inline bool      empty () const { return m_cnt.empty  () ; } 
      /// ok ?
      inline bool      ok    () const { return m_cnt.ok     () ; }
      // ======================================================================      
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update ( const double x ) override { add ( x ) ; }
      /// reset the cunter 
      void reset  () override { m_cnt.reset () ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const { return m_cnt.isfinite() ; } 
      // ======================================================================      
    public :
      // ======================================================================
      inline Ostap::Math::ValueWithError value () const { return m_cnt.mean () ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }      
      // ======================================================================
    public:
      // ======================================================================
      const Counter& counter() const { return m_cnt; }
      // ======================================================================
    private: 
      // ======================================================================
      Counter m_cnt {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class WArithmeticMean 
     *  Calculate the arithmetic mean 
     */
    class WArithmeticMean : public WStatistic
    {
    public:
      // ======================================================================
      typedef WMoment_<2>  Counter ;
      /// # entries
      typedef Counter::size_type size_type ;
      /// data type 
      typedef Counter::data_type data_type ;
      // ======================================================================
    public:
      // ======================================================================
      WArithmeticMean () = default ;
      WArithmeticMean ( const Counter& cnt ) ; 
      // ======================================================================
    public :
      // ======================================================================
      /// accumulate only positive entries 
      inline WArithmeticMean& add
      ( const double x          , 
        const double weight = 1 ) 
      { m_cnt.add ( x , weight ) ; return *this ; } 
      /// add two counters together
      inline WArithmeticMean& add ( const WArithmeticMean& x )
      { m_cnt.add ( x.m_cnt ) ; return *this ; }
      // ======================================================================
      /// add two counters together
      inline WArithmeticMean& operator+=( const WArithmeticMean& x )
      { return add ( x ) ; }
      // ======================================================================
      /// add sequence of values  
      template <class ITERATOR>
      inline
      WArithmeticMean&
      add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
        return *this ;
      }
      // ====================================================================== 
    public:
      // ======================================================================
      /// number of entries
      inline size_type size  () const { return m_cnt.size  () ; }
      /// number of effective entries
      inline data_type nEff  () const { return m_cnt.nEff  () ; }
      /// get sum of weighes \f$  \sum_i w_i \f$ 
      inline data_type  w     () const { return m_cnt.w     () ; }
      /// get sum of weights squared 
      inline data_type  w2    () const { return m_cnt.w2    () ; }
      /// empty ?
      inline bool       empty () const { return m_cnt.empty () ; } 
      /// ok ?
      inline bool       ok    () const { return m_cnt.ok    () ; } 
      // ======================================================================      
    public:
      // ======================================================================
      /// to use it as a Moment 
      void update
      ( const double x         , 
        const double weight = 1 ) override { add ( x , weight ) ; }
      /// reset the cunter 
      void reset  () override { m_cnt.reset () ; }
      // ======================================================================
    public :
      // ======================================================================      
      /// is finite ?
      bool isfinite () const { return m_cnt.isfinite() ; } 
      // ======================================================================      
    public :
      // ======================================================================
      inline Ostap::Math::ValueWithError value () const { return m_cnt.mean () ; }
      inline Ostap::Math::ValueWithError mean  () const { return value () ; }      
      // ======================================================================
    public:
      // ======================================================================
      const Counter& counter() const { return m_cnt; }
      // ======================================================================
    private: 
      // ======================================================================
      Counter m_cnt {} ;
      // ======================================================================
    } ;
    // =======================================================================

    // =======================================================================
    // More decorations 
    // ========================================================================
    
    // ========================================================================
    namespace Moments
    {
      // ======================================================================
      typedef Ostap::Math::ValueWithError VE ;
      // ======================================================================
      /// get the invalid moment 
      double invalid_moment () ;
      // ======================================================================
      
      /// get the central moment of order K=00
      template <unsigned short K ,
                unsigned short N ,
                typename std::enable_if<(K==0),int>::type = 0 >
      inline double
      moment ( const Moment_<N>& /* m */ ) { return 1 ; }
      /// get the central moment of order K=1
      template <unsigned short K ,
                unsigned short N ,
                typename std::enable_if<(K==1&&K<=N),int>::type = 0 >
      inline double
      moment ( const Moment_<N>& /* m */ ) { return 0 ; }
      /// get the central moment of order K=N 
      template <unsigned short K , 
                unsigned short N ,
                typename std::enable_if<(2<=K)&&(K==N),int>::type = 0 >
      inline Ostap::Math::ValueWithError      
      moment ( const Moment_<N>& m  ) { return m.template moment_<K> () ; }
      /// get the central moment of order 2<K && 2K<=N 
      template <unsigned short K , 
                unsigned short N ,
                typename std::enable_if<(2<=K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError      
      moment ( const Moment_<N>& m  ) { return m.template moment_<K> () ; }
      /// get the central moment of order K<N && N<=2K
      template <unsigned short K ,
                unsigned short N ,
                typename std::enable_if<(N<2*K)&&(N>K),int>::type = 0 >
      inline double
      moment ( const Moment_<N>& m  ) { return m.template moment_<K> () ; }
      
      // ======================================================================
      /** get the unbiased estimator for the 2nd order moment:
       *  \f[ \hat{\mu}_2 \equiv \frac{n}{n-1} \mu_2 \f] 
       *  @param  m input counter
       *  @return the unbiased estimate (if number of entries exceeds 2) and 
       *   <code>invalid_moment()</code> otherwise 
       */
      template <unsigned short N, typename std::enable_if<(N>=2),int>::type = 0 >
      inline double unbiased_2nd
      ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size() ;
        return !m.ok() || n < 2 ? invalid_moment () :
          m.template M_<2> () / ( n - 1 ) ;  
      }
      // ======================================================================
      /** get the unbiased estimator for the 3rd order moment:
       *  \f[ \hat{\mu}_3 \equiv \frac{n^2}{(n-1)(n-2)} \mu_3 \f] 
       *  @param  m input counter
       *  @return the unbiased estimate (if number of entries exceeds 3) and 
       *   <code>invalid_moment()</code> otherwise 
       */
      template <unsigned short N, typename std::enable_if<(N>=3),int>::type = 0 >
      inline double unbiased_3rd ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size() ;
        return !m.ok() || n < 3 ? invalid_moment() :
          m.template M_<3> * n / ( ( n - 1.0L ) * (  n - 2.0L  ) ) ;  
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
       *   <code>invalid_moment()</code> otherwise
       */
      template <unsigned short N, typename std::enable_if<(N>=4),int>::type = 0 > 
      inline double unbiased_4th ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( !m.ok() ||  (n < 4 ) ) { return invalid_moment()  ; }
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
       *   <code>invalid_moment()</code> otherwise
       */
      template <unsigned short N, typename std::enable_if<(N>=5),int>::type = 0 >
      inline double unbiased_5th ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( !m.ok() ||  ( n < 5 ) ) { return invalid_moment()  ; }
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
      
      // ======================================================================
      // MEAN
      // ======================================================================
      /// get the mean      
      template <unsigned short N, typename std::enable_if<(N>1),int>::type = 0 >
      inline VE     mean ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size () ;
        if ( !m.ok () || n  < 2  ) { return m.mu (); }  // ATTENTION! 
        const  double mu  = m.mu ()            ;
        const  double m2  = unbiased_2nd ( m ) ;        // ATTENTION! 
        return VE ( mu , m2 / n ) ;
      }
      // ======================================================================
      /// get the mean
      inline double mean ( const Moment_<1>& m ) { return m.mu () ; }
      // ======================================================================

      // ======================================================================
      // VARIANCE 
      // ======================================================================
      /// get the unbiased sample variance with uncertainty
      template <unsigned short N, typename std::enable_if<(N>3),int>::type = 0 >
      inline VE variance ( const Moment_<N>& m )
      {
        const unsigned long long n = m.size () ;
        if ( !m.ok() || n  < 2  ) { return invalid_moment() ; } // RETURN
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
      /// get the variance  
      template <unsigned short N, typename std::enable_if<(1<=N)&&(4>N),int>::type = 0 >
      inline double variance ( const Moment_<N>& m )
      {
        const unsigned long long n =  m.size() ;
        if ( !m.ok() || 2 > n ) { return invalid_moment() ; }        
        return unbiased_2nd ( m ) ;
      }
      // ======================================================================
      
      // ======================================================================
      // SKEWNESS 
      // ======================================================================
      /// get the estimate for the sample skewness  \f$ \frac{m_3}{\sigma^{3/2}}\f$
      template <unsigned short N, typename std::enable_if<(N>=3),int>::type = 0 >
      inline VE skewness ( const Moment_<N>& m ) 
      {
        const unsigned long long n = m.size() ;
        if ( !m.ok() || n < 3 ) { return invalid_moment() ; } // RETURN
        //
        const double m3   = unbiased_3rd ( m ) ;
        const double m2   = m.template M_<2>  () / n ;
        const double skew = m3 / std::pow ( m2 , 3.0/2 ) ;
        //
        const double cov2 = 6.0L * n * ( n - 1 ) / ( ( n - 2.0L ) * ( n + 1.0L ) * ( n + 3.0L ) ) ;
        return VE ( skew , cov2 ) ;
      }
      // ======================================================================
      
      // ======================================================================
      // KURTOSIS
      // ======================================================================
      /// get the estimate for the sample (excessive) kurtosis \f$ \frac{m_4}{\sigma^{4}}-3\f$
      template <unsigned short N, typename std::enable_if<(N>=4),int>::type = 0 >
      static inline VE kurtosis ( const Moment_<N>& m ) 
      {
        const unsigned long long n = m.size() ;
        if ( !m.ok() || n < 4 ) { return invalid_moment() ; }
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
      inline double central_moment ( const Moment_<N>& /* m */ ) { return 1 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      inline double central_moment ( const Moment_<N>& /* m */ ) { return 0 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K>=2)&&(2*K>N)&&(K<=N),int>::type = 1 >
      inline double central_moment ( const Moment_<N>& m )
      { return m.template moment_<K> () ; }
      /// get the central moment of order \f$ N \f$      
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(2<=K) && (N>=2*K),int>::type = 0 >
      inline VE  central_moment ( const Moment_<N>& m )
      { return m.template moment_<K> () ; }
      // ======================================================================
      
      // ======================================================================
      /// get the standartized central moment of order 1
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==0),int>::type = 1 >
      inline double std_moment ( const Moment_<N>& /* m */ ) { return 1 ; }
      /// get the standartised central moment of order \f$ K \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      inline double std_moment ( const Moment_<N>& /* m */ ) { return 0 ; }
      /// get the standartised central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K>=2)&&(K<=N),int>::type = 1 >
      inline double std_moment ( const Moment_<N>& m )
      { return m.template std_moment_<K> () ; }
      // ======================================================================
      
      // ======================================================================
      /// get the centralized moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K<=N),int>::type = 1 >
      inline double centralized_moment
      ( const Moment_<N>& m      ,
        const double      center )
      { return m.template centralized_moment_<K> ( center ) ; }
      
      // ======================================================================
      /// get the first cumulant, well, it is actually the first central moment 
      template <unsigned short N,
                typename std::enable_if<(N==1),int>::type = 1 >
      inline double cumulant_1st ( const Moment_<N>& m ) 
      { return m.mean () ; }
      /// get the first cumulant, well, it is actually the first central moment 
      template <unsigned short N,
                typename std::enable_if<(N>1),int>::type = 1 >
      inline VE     cumulant_1st ( const Moment_<N>& m ) 
      { return m.mean () ; }
      /// get the second unbiased cumulant, well it is actually the second unbiased central moment 
      template <unsigned short N,
                typename std::enable_if<(N>=2)&&(N<=3),int>::type = 1 >
      inline double cumulant_2nd ( const Moment_<N>& m ) { return unbiased_2nd ( m ) ; }
      /// get the second unbiased cumulant, well it is actually the second unbiased central moment 
      template <unsigned short N,
                typename std::enable_if<(N>=4),int>::type = 1 >
      inline VE     cumulant_2nd ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 2 ) ) { return invalid_moment() ;} 
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
      inline double cumulant_3rd  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 3 ) ) { return invalid_moment() ; }
        //
        const auto n = m.size() ;
        const double m3 = m.template M_<3> () / n ; 
        //
        return m3 * n * n  / ( ( n - 1 ) * ( n - 2 ) ) ; 
      }
      /// get the third unbiased cumulant, well it is actually the third unbiased central moment 
      template <unsigned short N,
                typename std::enable_if<(N>=6),int>::type = 1 >
      inline VE cumulant_3rd  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 3 ) ) { return invalid_moment() ; }
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
      inline double cumulant_4th  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 4 ) ) { return invalid_moment() ; }
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
      inline VE cumulant_4th  ( const Moment_<N>& m ) 
      { 
        if ( !m.ok() || ( m.size() < 4 ) ) { return invalid_moment() ; }
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

      /// get the central moment of order K=00
      template <unsigned short K ,
                unsigned short N ,
                typename std::enable_if<(K==0),int>::type = 0 >
      inline double
      moment ( const WMoment_<N>& /* m */ ) { return 1 ; }
      /// get the central moment of order K=1
      template <unsigned short K ,
                unsigned short N ,
                typename std::enable_if<(K==1&&K<=N),int>::type = 0 >
      inline double
      moment ( const WMoment_<N>& /* m */ ) { return 0 ; }
      /// get the central moment of order K=N 
      template <unsigned short K , 
                unsigned short N ,
                typename std::enable_if<(2<=K)&&(K==N),int>::type = 0 >
      inline Ostap::Math::ValueWithError      
      moment ( const WMoment_<N>& m  ) { return m.template moment_<K> () ; }
      /// get the central moment of order 2<K && 2K<=N 
      template <unsigned short K , 
                unsigned short N ,
                typename std::enable_if<(2<=K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError      
      moment ( const WMoment_<N>& m  ) { return m.template moment_<K> () ; }
      /// get the central moment of order K<N && N<=2K
      template <unsigned short K ,
                unsigned short N ,
                typename std::enable_if<(N<2*K)&&(N>K),int>::type = 0 >
      inline double
      moment ( const WMoment_<N>& m  ) { return m.template moment_<K> () ; }      
      
      // ======================================================================
      // MEAN
      // ======================================================================
      /// get the mean      
      template <unsigned short N, typename std::enable_if<(N>1),int>::type = 0 >
      inline VE     mean ( const WMoment_<N>& m )
      {
        if ( !m.ok() || m.size() < 2 ) { return m.mu() ; }
        //
        const auto n = m.nEff () ;
        const  double mu  = m.mu     (   ) ;
        const  double m2  = m.template moment_<2> () ;
        return VE ( mu , m2 / n ) ;
      }
      /// get the mean
      inline double mean ( const WMoment_<1>& m ) { return m.mu () ; }
      // ======================================================================
      
      // ======================================================================
      // VARIANCE 
      // ======================================================================
      /// get the  sample variance with uncertainty
      template <unsigned short N, typename std::enable_if<(N>3),int>::type = 0 >
      inline VE variance ( const WMoment_<N>& m )
      {
        if ( !m.ok() || m.size() < 2  ) { return invalid_moment() ; } // RETURN
        //
        const double m2 = m.template M_<2> () / m.w()  ;
	    // 
        if ( 0 >  m2 ) { return invalid_moment() ; } // RETURN
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
      /// get the variance  
      template <unsigned short N, typename std::enable_if<(1<=N)&&(4>N),int>::type = 0 >
      inline double variance ( const WMoment_<N>& m )
      { return !m.ok() || ( m.size() < 2 ) ? invalid_moment() : m.template moment_<2>() ; }      
      // ======================================================================
      
      // ======================================================================
      // SKEWNESS 
      // ======================================================================
      /// get the estimate for the sample skewness  \f$ \frac{m_3}{\sigma^{3/2}}\f$
      template <unsigned short N, typename std::enable_if<(N>=3),int>::type = 0 >
      inline VE skewness ( const WMoment_<N>& m ) 
      {
        if ( !m.ok() || ( m.size() < 3 ) ) { return invalid_moment() ; }
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

      // ======================================================================
      // KURTOSIS 
      // ======================================================================
      /// get the estimate for the sample (excessive) kurtosis \f$ \frac{m_4}{\sigma^{4}}-3\f$
      template <unsigned short N, typename std::enable_if<(N>3),int>::type = 0 >
      inline VE kurtosis ( const WMoment_<N>& m ) 
      {
        if ( !m.ok() || ( m.size() < 4 ) ) { return invalid_moment() ; }
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
      inline double central_moment ( const WMoment_<N>& /* m */ ) { return 1 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      inline double central_moment ( const WMoment_<N>& /* m */ ) { return 0 ; }
      /// get the central moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K>=2)&&(2*K>N)&&(K<=N),int>::type = 1 >
      inline double central_moment ( const WMoment_<N>& m )
      { return m.template moment_<K> () ; }
      /// get the central moment of order \f$ N \f$      
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(2<=K) && (N>=2*K),int>::type = 0 >
      inline VE     central_moment ( const WMoment_<N>& m )
      { return m.template moment_<K> () ; }

      
      // ======================================================================
      /// get the standartized moment of order 1
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==0),int>::type = 1 >
      inline double std_moment ( const WMoment_<N>& /* m */ ) { return 1 ; }
      /// get the standartized moment of order 1
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==1),int>::type = 1 >
      inline double std_moment ( const WMoment_<N>& /* m */ ) { return 0 ; }
      /// get the standartized moment of order 2
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K==2),int>::type = 1 >
      inline double std_moment ( const WMoment_<N>& /* m */ ) { return 1 ; }
      /// get the standartized moment of order 1
      template <unsigned short K, unsigned short N, 
                typename std::enable_if<(2<K)&&(2*K<=N),int>::type = 0 >
      inline Ostap::Math::ValueWithError
      std_moment ( const WMoment_<N>& m ) { return m.template std_moment_<K> () ; }
      /// get the standartized moment of order 1
      template <unsigned short K, unsigned short N, 
                typename std::enable_if<(N<2*K)&&(K<=N),int>::type = 0 >
      inline double
      std_moment ( const WMoment_<N>& m ) { return m.template std_moment_<K> () ; }

      // ======================================================================
      /// get the centralized moment of order \f$ N \f$  
      template <unsigned short K, unsigned short N,
                typename std::enable_if<(K<=N),int>::type = 1 >
      inline double centralized_moment
      ( const WMoment_<N>& m      ,
        const double       center )
      { return m.template centralized_moment_<K> ( center ) ; }
      
     
      // ======================================================================
    } ; //                            The end of namespace Ostap::Math::Moments
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END  
// ============================================================================
#endif // OSTAP_MOMENTS_H
// ============================================================================

