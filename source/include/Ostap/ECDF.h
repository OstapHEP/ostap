// ============================================================================
#ifndef OSTAP_ECDF_H 
#define OSTAP_ECDF_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STD
// ============================================================================
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/Statistic.h"
#include "Ostap/Moments.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class ECDF Ostap/ECDF.h
     *  Empirical cumulative distribution function 
     *  @author Vanya Belyaev
     *  @date   2023-09-17
     */
    class ECDF : public Statistic 
    {
    public: 
      // ======================================================================
      /// the actual type of data 
      typedef std::vector<double>           Data     ;
      /// iterator type
      typedef Data::const_iterator          iterator ;
      /// the actual type of indices 
      typedef std::vector<Data::size_type>  Indices  ;
      // ======================================================================
      /** quantile Hyndman Fan taxonomy
       *  see Table 1 in 
       *  @see https://arxiv.org/abs/2304.07265
       *  @see Andrey Akinshin, "Weighted quantile estimators", arXiv:2304.07265
       */
      enum QType {
	One = 1 ,
	Two     ,
	Three   ,
	Four    ,
	Five    ,
	Six     ,
	Seven   ,
	Eight   ,
	Nine    
      } ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Constructor from  data
       *  data must be non-empty!
       */ 
      ECDF
      ( const Data&  data                  ,
        const bool   complementary = false ) ;
      // =======================================================================
      /** Constructor from data
       *  data must be non-empty!
       */ 
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename            = std::enable_if<std::is_convertible<value_type,double>::value> >
      ECDF
      ( ITERATOR   begin                 ,
        ITERATOR   end                   ,
        const bool complementary = false ) 
        : m_data          ( begin , end   )
        , m_complementary ( complementary )
	      , m_counter       () 
      {
	/// (1) adjust the content 
	this -> cleanup () ;
	/// (2) sort it if needed 
        if ( !std::is_sorted ( m_data.begin () , m_data.end () ) )
          { std::sort ( m_data.begin () , m_data.end () ) ; }
	/// (3) update counters 
	for ( auto v : m_data ) { m_counter.add ( v ) ; }
      }
      /// constructor to create complementary/ordinary ECDF
      ECDF
      ( const ECDF&  right         ,
        const bool   complementary ) ;
      /// copy constructor
      ECDF ( const ECDF&  right ) = default ;
      /// move constructor 
      ECDF (       ECDF&& right ) = default ;
      /// default constructor
      ECDF () = default ;
      // ======================================================================
    public: 
      // ======================================================================
      /// copy assignement 
      ECDF& operator=( const ECDF&  right ) ;
      /// move assignement 
      ECDF& operator=(       ECDF&& right ) ;
      // ======================================================================
    public: // the main method 
      // ======================================================================
      /// the main method 
      double        evaluate   ( const double x ) const ;
      /// the main method 
      inline double operator() ( const double x ) const
      { return evaluate ( x ) ; }
      // ======================================================================
      /** get the value with binomial estimate for uncertainty
       *  @see Ostap::Math::binomEff 
       */ 
      Ostap::Math::ValueWithError estimate ( const double x ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// add a value to data container  
      ECDF& add ( const double value  ) ;
      /// add more values to data constainer 
      ECDF& add ( const Data&  values ) ; 
      /// add more values to data container 
      ECDF& add ( const ECDF&  values ) ;
      /// add a bunch of values 
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename            = std::enable_if<std::is_convertible<value_type,double>::value> >
      ECDF& add
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        if ( std::is_sorted ( begin , end ) ) { return this -> add_sorted ( begin , end ) ; }
        /// copy input data (only valid/finite entries)  
        Data tmp ( begin, end ) ;	
        /// sort input data 
        std::sort ( tmp.begin () , tmp.end () ) ;
	/// add sorted input 
        return this->add_sorted ( tmp.begin() , tmp.end() ) ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// Ostap::Math::Statistic::update
      void update ( const double x ) override { this -> add  ( x ) ; }
      /// Reset
      void reset  () override { m_data.clear() ; m_counter.reset() ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// add a bunch of sorted values 
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename            = std::enable_if<std::is_convertible<value_type,double>::value> >
      ECDF& add_sorted 
      ( ITERATOR begin ,
        ITERATOR end   )
      {
        const std::size_t d = std::distance ( begin , end ) ;
        /// nothing to add ?
        if ( !d  ) { return  *this  ; } // nothing to add 
        /// prepare the output 
        Data tmp  ( d + m_data.size () )  ;
        /// merge two sorted containers 
        std::merge ( begin           ,                     
                     end             ,
                     m_data.begin () ,
                     m_data.end   () ,                     
                     tmp.begin    () ) ;
        /// swap the merged result  with own data 
        std::swap ( m_data , tmp ) ;
	/// (1) remove bad entries (if needed) 
	this -> cleanup () ;
	/// (2) update counter
	for ( ITERATOR v = begin ; end != v ; ++v) { m_counter.add ( *v ) ; }
	return *this ; 
      }
      // ======================================================================
    public :
      // ======================================================================
      inline ECDF& __iadd__ ( const double      x ) { return add ( x ) ; }
      inline ECDF& __iadd__ ( const ECDF&       x ) { return add ( x ) ; }
      inline ECDF& __iadd__ ( const ECDF::Data& x ) { return add ( x ) ; }
      ECDF         __add__  ( const double      x ) const ;      
      ECDF         __add__  ( const ECDF&       x ) const ;      
      ECDF         __add__  ( const ECDF::Data& x ) const ;      
      // ======================================================================
    public :
      // ======================================================================
      inline ECDF& operator+= ( const double      x ) { return add ( x ) ; }
      inline ECDF& operator+= ( const ECDF&       x ) { return add ( x ) ; }      
      inline ECDF& operator+= ( const ECDF::Data& x ) { return add ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// ok ?
      inline bool            ok            () const
      { return !m_data.empty() && m_data.size() == m_counter.nEntries () ; }
      /// empty?
      inline bool            empty         () const { return m_data.empty () ; } 
      /// data size 
      inline Data::size_type N             () const { return m_data.size  () ; } 
      /// data size 
      inline Data::size_type size          () const { return m_data.size  () ; }
      /// number of effective entries
      inline Data::size_type nEff          () const { return m_data.size  () ; }      
      // ======================================================================
      /// access to data
      inline const Data&     data          () const { return m_data         ; }
      // access to data 
      inline double          data   ( const unsigned int index ) const
      { return index < m_data.size() ?  m_data [ index ] : xmax () ; }
      // ======================================================================
      /// complementary?
      inline bool            complementary () const { return m_complementary ; }
      /// minimal x-value
      inline double          xmin          () const { return m_counter.min () ; } 
      /// maximal x-value
      inline double          xmax          () const { return m_counter.max () ; } 
      // ======================================================================
      /// get the abscissa value with the given index 
      inline double operator[] ( const unsigned int index ) const
      { return index < m_data.size() ? m_data[index] : xmax ()  ; }
      // ======================================================================
      // get the value of F_k
      inline double Fk ( const unsigned int k ) const
      {
        const Data::size_type n = m_data.size() ;
        return 0 == k ? 0.0 : n <= k ? 1.0 : k * 1.0 / n  ;
      }
      // ======================================================================
      /// expose the data: begin iterator 
      inline iterator begin () const { return m_data.begin () ; }
      /// expose the data: end-iterator 
      inline iterator end   () const { return m_data.end   () ; }      
      // ======================================================================
    public:
      // ======================================================================
      /** assuming that x comes from the same distribution 
       *  return transformed value  \f$ g = f(x) \f$, such that 
       *  \f$ g \f$  has Gaussian distribution 
       */
      double gauss   ( const double x ) const ;
      // =======================================================================
      /** assuming that x comes from the same distribution 
       *  return transformed value  \f$ u = f(x) \f$, such that 
       *  \f$ u \f$  has uniform distribution for \f$ 0 \le  u \le 1 \f$ 
       */
      double uniform ( const double x ) const ;
      // ======================================================================
      /** get p-quantile of distribution: \f$ 1 \le p \le1  \f$
       *  @see scipy.stats.mstats
       *
       * Typical avleus for alphap, betap are:
       * 
       * - (0,1) : p(k) = k/n : linear interpolation of cdf (R type 4)
       * - (.5,.5) : p(k) = (k - 1/2.)/n : piecewise linear function (R type 5)
       * - (0,0) : p(k) = k/(n+1) : (R type 6)
       * - (1,1) : p(k) = (k-1)/(n-1): p(k) = mode[F(x[k])]. (R type 7, R default)
       * - (1/3,1/3): p(k) = (k-1/3)/(n+1/3): Then p(k) ~ median[F(x[k])].
       *   The resulting quantile estimates are approximately median-unbiased 
       *   regardless of the distribution of x. (R type 8)
       * - (3/8,3/8): p(k) = (k-3/8)/(n+1/4): Blom.
       *   The resulting quantile estimates are approximately 
       *   unbiased if x is normally distributed (R type 9)
       * - 0(.4,.4) : approximately quantile unbiased (Cunnane)
       * - (.35,.35): APL, used with PWM
       *
       *  @param p      (INPUT) quantile
       *  @param alphap (INPUT) parameter alphap \f$ 0 \le\alpha_p \le 1 \f$ 
       *  @param abetap (INPUT) parameter betap \f$ 0 \le\beta_p \le 1 \f$ 
       */    
      double quantile
      ( const double p            ,  
        const double alphap = 0.4 , 
        const double betap  = 0.4 ) const ;
      // ======================================================================
      /** get p-quantile
       *  @see https://arxiv.org/abs/2304.07265
       *  @see Andrey Akinshin, "Weighted quantile estimators", arXiv:2304.07265     
       */
      double quantile_HF
      ( const double p         ,
	const QType  t = Seven ) ;
      // ======================================================================
      /** Get Harrel-Davis estimator for quantile function
       *  @param p  (INPUT) quantile 
       *  @retiurn Harrel-Davis quantile estimator
       *  @see https://doi.org/10.1093/biomet/69.3.635
       *  @see F.E. Harrel and C.E.Davis,  "A new distribution-free quantile estimator",
       *       Biometrika 63.9 (Dec 1982), pp. 635-640
       */
      double quantile_HD ( const double p ) const ;
      // ======================================================================
    public:
      //=======================================================================
      /// get simple statistics 
      const Ostap::StatEntity& counter() const { return m_counter ; }
      // ======================================================================
      /** statistics (as statistics)
       * @param stat (UPDATE) input statistic object
       * @return updated statistics object
       */
      Ostap::Math::Statistic& 
      statistics
      ( Ostap::Math::Statistic& stat ) ;
      // ======================================================================
      /// get the moments 
      template <unsigned short K>
      Ostap::Math::Moment_<K> moment_ () const
     {
       Ostap::Math::Moment_<K> m{} ;
       for ( auto value : m_data ) { m.add ( value ) ; }
       return m ;
     } 
      // ======================================================================
    public:
      //=======================================================================
      /** number of elements that are less or equal to x
       *  "rank of x" in the ordered sample
       */
      inline Data::size_type rank ( const double x ) const
      { return std::upper_bound ( m_data.begin () , m_data.end () , x ) - m_data.begin() ; }
      // ======================================================================
      /// get ranks for all elements from another sample 
      Indices ranks ( const ECDF& sample ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( ECDF& right ) ;
      // ======================================================================
    private :
      // ======================================================================
      /// Remove elments <code>!std::isfinite</code>
      ECDF& cleanup () ; // Remove elments <code>!std::isfinite</code>
      // ======================================================================
    private:
      // ======================================================================
      /// container of sorted data 
      Data              m_data          {       } ; // container of sorted data
      /// complementary CDF?
      bool              m_complementary { false } ; // complementary CDF ?
      /// counter
      Ostap::StatEntity m_counter       {} ;
      // ======================================================================
    }; //                                The end of the class Ostap::Math::ECDF
    // ========================================================================
    /// swap two objects 
    inline void swap
    ( ECDF& left ,
      ECDF& right )
    { left.swap ( right ) ; }
    /// add two ECDFs
    inline ECDF
    operator+
    ( const ECDF& a ,
      const ECDF& b )
    { ECDF c { a } ; c += b ; return c ; }
    // ========================================================================
    /** @class WECDF Ostap/ECDF.h
     *  Empirical cumulative distribution function for weighted data 
     *  @author Vanya Belyaev
     *  @date   2024-11-28
     */
    class WECDF : public WStatistic
    {
    public: 
      // ======================================================================
      /// The actual type of entry: (y,w)
      typedef std::pair<double,double>    Entry    ;
      /// the actual type of data
      typedef std::vector<Entry>          Data     ;
      /// iterator type
      typedef Data::const_iterator        iterator ;
      /// the actual type of indices 
      typedef ECDF::Indices               Indices  ;
      // ======================================================================
      /// ordering criteria - order by the first component/abscissas  
      struct COMPARE
      {
	/// comparison criteria: compare abscissas 
	inline bool operator () ( const Entry& a , const Entry& b ) const
	{ return a.first < b.first ; }
	/// comparison criteria: compare abscissas 
	inline bool operator () ( const Entry& a , const double b ) const
	{ return a.first < b ; }
	/// comparison criteria: compare abscissas
	inline bool operator () ( const double a , const Entry& b ) const
	{ return a < b.first ; }
      } ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Constructor from  data
       *  data must be non-empty!
       */ 
      WECDF
      ( const Data&  data                  ,
        const bool   complementary = false ) ;
      // ======================================================================
      /** Constructor from data
       *  data must be non-empty!
       */ 
      WECDF
      ( const ECDF::Data&  data                  ,
        const ECDF::Data&  weights               ,
        const bool         complementary = false ) ;
      // =======================================================================
      /** Constructor from  data
       *  data must be non-empty!
       */ 
      WECDF
      ( const ECDF::Data&  data                  ,
        const bool         complementary = false ) ;      
      // =======================================================================
      WECDF
      ( const WECDF&  right         ,
        const bool    complementary ) ;
      // ======================================================================
      WECDF
      ( const ECDF&   right         ,
        const bool    complementary ) ;
      // ======================================================================
      WECDF
      ( const ECDF&   right         ) ;
      // ======================================================================
      /// copy constructor
      WECDF ( const WECDF&  right ) = default ;
      /// move constructor 
      WECDF (       WECDF&& right ) = default ;
      // ======================================================================
      /// default constructor
      WECDF () = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// copy assignement 
      WECDF& operator=( const WECDF&  right ) ;
      /// move assignement 
      WECDF& operator=(       WECDF&& right ) ;
      // ======================================================================
    public: // the main method 
      // ======================================================================
      /// the main method 
      double        evaluate   ( const double x ) const ;
      /// the main method 
      inline double operator() ( const double x ) const
      { return evaluate ( x ) ; }
      // ======================================================================
      /** get the value with estimate for uncertainty
       *  @see Ostap::Math::binomEff 
       */ 
      Ostap::Math::ValueWithError estimate ( const double x ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// add a value to data container  
      inline WECDF&
      add
      ( const double value        ,
        const double weight = 1.0 ) { return add ( Entry ( value , weight ) ) ; }
      /// add a value to data container  
      WECDF& add ( const Entry&       entry  ) ;      
      /// add more values to data container 
      WECDF& add ( const WECDF&       values ) ;
      /// add more values to data container       
      WECDF& add ( const WECDF::Data& values ) ;
      /// add more values to data container 
      WECDF& add ( const ECDF&        values ) ;
      /// add more values to data container 
      WECDF& add ( const ECDF::Data&  values ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// Ostap::Math::WStatistic::update
      void update
      ( const double x     , 
        const double w = 1 ) override { this->add  ( x , w ) ; }
      /// Reset
      void reset  () override { m_data.clear() ; m_counter.reset () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline WECDF& operator+= ( const double       x ) { return add ( x ) ; }
      inline WECDF& operator+= ( const Entry&       x ) { return add ( x ) ; }
      inline WECDF& operator+= ( const WECDF&       x ) { return add ( x ) ; }      
      inline WECDF& operator+= ( const WECDF::Data& x ) { return add ( x ) ; }      
      inline WECDF& operator+= ( const ECDF&        x ) { return add ( x ) ; }
      inline WECDF& operator+= ( const ECDF::Data&  x ) { return add ( x ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      inline WECDF& __iadd__ ( const double       x ) { return add ( x ) ; }
      inline WECDF& __iadd__ ( const Entry&       x ) { return add ( x ) ; }
      inline WECDF& __iadd__ ( const WECDF&       x ) { return add ( x ) ; }
      inline WECDF& __iadd__ ( const WECDF::Data& x ) { return add ( x ) ; }
      inline WECDF& __iadd__ ( const ECDF&        x ) { return add ( x ) ; }
      inline WECDF& __iadd__ ( const ECDF::Data&  x ) { return add ( x ) ; }
      // ======================================================================      
      WECDF  __add__  ( const WECDF&  x ) const ;      
      WECDF  __add__  ( const  ECDF&  x ) const ;      
      // ======================================================================
    public:
      // ======================================================================
      /// ok ?
      inline bool            ok            () const
      { return !m_data.empty()
	  && m_data.size() == m_counter.nEntries ()
	  && 0 < m_counter.sumw () ; }      
      /// empty?
      inline bool            empty         () const { return m_data.empty () ; } 
      /// data size 
      inline Data::size_type N             () const { return m_data.size  () ; } 
      /// data size 
      inline Data::size_type size          () const { return m_data.size  () ; }
      /// number of effective entries
      inline double          nEff          () const { return m_counter.nEff () ;  }
      /// sum of all weights
      inline double          sumw          () const { return m_counter.sumw  () ; }      
      /// sum of all squared weights
      inline double          sumw2         () const { return m_counter.sumw2 () ; }      
      // ======================================================================      
      /// access to data
      inline const Data&     data          () const { return m_data         ; }
      // access to data 
      inline double data   ( const unsigned int index ) const
      { return index < m_data.size() ?  m_data [ index ].first  : xmax ()   ; }
      // access to weight 
      inline double weight ( const unsigned int index ) const
      { return index < m_data.size() ?  m_data [ index ].second : 0.0       ; }      
      // ======================================================================
      /// complementary?
      inline bool            complementary () const { return m_complementary  ; }
      /// minimal x-value
      inline double          xmin          () const { return m_counter.min () ; } 
      /// maximal x-value
      inline double          xmax          () const { return m_counter.max () ; } 
      // ======================================================================
      /// get the abscissa value with the given index 
      inline const Entry& operator[] ( const unsigned int index ) const
      { return index < m_data.size() ? m_data [ index ] : m_data.back() ; }
      // ======================================================================      
      /// expose the data: begin iterator 
      inline iterator begin () const { return m_data.begin () ; }
      /// expose the data: end-iterator 
      inline iterator end   () const { return m_data.end   () ; }      
      // ======================================================================
   public:
      //=======================================================================
      /// statistics (as Counter)
      const Ostap::WStatEntity& counter() const { return m_counter ; } 
      /** statistics (as statistics)
       * @param stat (UPDATE) input statistic object
       * @return updated statistics object
       */
      Ostap::Math::WStatistic& 
      statistics
      ( Ostap::Math::WStatistic& stat ) ;
      /// get the moments 
      template <unsigned short K>
      Ostap::Math::WMoment_<K> moment_ () const
      {
        Ostap::Math::WMoment_<K> m{} ;
	for ( const auto& value : m_data ) { m.add ( value.first, value.second ) ; }
	return m ;
      } 
      // =======================================================================
    public: 
      // =======================================================================
      /** Get Harrel-Davis estimator for quantile function
       *  @param p  (INPUT) quantile 
       *  @retiurn Harrel-Davis quantile estimator
       *
       *  @see https://arxiv.org/abs/2304.07265
       *  @see Andrey Akinshin, "Weighted quantile estimators", arXiv:2304.07265
       *   
       *  @see https://doi.org/10.1093/biomet/69.3.635
       *  @see F.E. Harrel and C.E.Davis,  "A new distribution-free quantile estimator",
       *       Biometrika 63.9 (Dec 1982), pp. 635-640
       */
      double quantile_HD ( const double p ) const ;
      // ======================================================================
    public: 
      // ======================================================================
      /// swap two objects 
      void swap ( WECDF& right ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// calculate \f$ \sum_i^{n} w_i \f$
      inline double sumw  ( const unsigned long n ) const
      { return std::accumulate ( m_data.begin () ,
				 m_data.begin () + std::min ( n , m_data.size() ) , 
				 0.0             ,
				 [] ( const double s , const Entry& entry ) -> double
				 { return s + entry.second ; } ) ;	}
      /// calculate \f$ \sum_i^{n} w_i \f$
      inline double sumw2 ( const unsigned long n ) const
      { return std::accumulate ( m_data.begin () ,
                                 m_data.begin () + std::min ( n , m_data.size() ) , 
                                 0.0             ,
                                 [] ( const double s , const Entry& entry ) -> double
                                 { return s + entry.second * entry.second ; } ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      /** number of elements that are less or equal to x
       *  "rank of x" in the ordered sample
       */
      inline Data::size_type rank ( const double x ) const
      { return std::upper_bound ( m_data.begin ()   ,
                                  m_data.end   ()   ,
                                  Entry ( x , 1.0 ) ,
				  COMPARE      ()   ) - m_data.begin() ; }
      // ======================================================================
      /// get ranks for all elements from another sample 
      ECDF::Indices ranks ( const  ECDF& sample ) const ;
      /// get ranks for all elements from another sample 
      ECDF::Indices ranks ( const WECDF& sample ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// remove elements <code>!std::isfinite</code>
      WECDF& cleanup () ; // check that ECDF is OK: there are some entries 
      // ======================================================================
    private: 
      // ======================================================================
      /// container of sorted data 
      Data               m_data    {} ; // container of sorted data
      /// complementary CDF?
      bool               m_complementary { false } ; // complementary CDF ? 
      /// Basic Statistics
      Ostap::WStatEntity m_counter {} ;        // Basic statistics
      // ======================================================================
    } ; //                              The end of the class Ostap::Math::WECDF 
    // ========================================================================
    /// swap two objects 
    inline void swap ( WECDF& left , WECDF& right ) { left.swap ( right ) ; }
    /// add two WECDFs
    inline WECDF
    operator+
    ( const WECDF& a ,
      const WECDF& b )
    { WECDF c { a } ; c += b ; return c ; }
    /// add WECDF and ECDF 
    inline WECDF
    operator+
    ( const WECDF& a ,
      const  ECDF& b )
    { WECDF c { a } ; c += b ; return c ; }
    /// add WECDF and ECDF 
    inline WECDF
    operator+
    ( const  ECDF& a ,
      const WECDF& b )
    { WECDF c { b } ; c += a ; return c ; }
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ECDF_H
// ============================================================================
