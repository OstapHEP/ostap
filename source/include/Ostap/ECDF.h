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
#include "Ostap/Power.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    class  ECDF ;
    class WECDF ;
    // ========================================================================
    // Kernels
    // https://en.wikipedia.org/wiki/Kernel_(statistics)
    // ========================================================================
    /// uniform Kernel 
    inline double k_uniform      ( const double u ) { return 0.5 * ( std::abs ( u ) <= 1 ) ; }
    /// Triangular kernel
    inline double k_triangular   ( const double u ) { return std::abs ( u ) <= 1 ? ( 1 - std::abs ( u ) ) : 0.0 ; }
    /// Epanechnnikov/parabolic  kernel 
    inline double k_epanechnikov ( const double u ) { return std::abs ( u ) <= 1 ? ( 1 - u * u )          : 0.0 ; }
    /// Epanechnikov Parabolic kernel
    inline double k_parabolic    ( const double u ) { return k_epanechnikov ( u ) ; }
    /// Quartic/biweight kernel  
    inline double k_quartic      ( const double u ) { return std::abs ( u ) <= 1 ? 15 * Ostap::Math::POW ( 1.0 - u * u ,  2  ) / 16 : 0.0 ; } 
    /// Quartic/biweight kernel
    inline double k_biweight     ( const double u ) { return k_quartic ( u ) ; }
    /// Triweight kernel  
    inline double k_triweight    ( const double u ) { return std::abs ( u ) <= 1 ? 35 * Ostap::Math::POW ( 1.0 - u * u ,  3  ) / 32 : 0.0 ; } 
    /// Tricube kernel
    inline double k_tricube      ( const double u ) { return std::abs ( u ) <= 1 ? 70 * Ostap::Math::POW ( 1.0 - std::abs ( u * u * u ) , 3 ) / 81 : 0.0 ; } 
    /// Gaussian kernel 
    double        k_gaussian     ( const double u ) ; 
    /// Cosine kernel 
    double        k_cosine       ( const double u ) ; 
    /// Logistic Kernel
    double        k_logistic     ( const double u ) ;
    /// sigmoid kernel 
    double        k_sigmoid      ( const double u ) ; 
    // ==========================================================================
    /** @class DensityEstimator
     *  Helper class for non-parametetric density estimators
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     */ 
    class DensityEstimator
    {
      //====================================================================
    public: 
      //====================================================================
      enum Kernel
	{
	  Uniform      , 
	  Rectangular  = Uniform      , 
	  Boxcar       = Uniform      , 
	  Triangular   ,
	  Epanechnikov , 
	  Parabolic    = Epanechnikov , 
	  Quartic      , 
	  Biweight     = Quartic      , 
	  Triweight    , 
	  Tricube      ,
	  Gaussian     , 
	  Cosine       , 
	  Logistic     , 
	  Sigmoid      , 
	  Last         = Sigmoid    
	};
      // ======================================================================
    public:
      // ======================================================================
      /// get the kernel estimate
      static double kernel ( const double u , const Kernel k ) ;
      // ======================================================================
      /// get the "optimal" value for smoothing parameter 
      static double hopt   ( const  ECDF& data ) ;
      /// get the "optimal" value for smoothing parameter 
      static double hopt   ( const WECDF& data ) ;
      // ======================================================================
    };
    // ========================================================================
    /** @class ECDF Ostap/ECDF.h
     *  Empirical cumulative distribution function 
     *  @author Vanya Belyaev
     *  @date   2023-09-17
     */
    class ECDF
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
      {
        if ( !std::is_sorted ( m_data.begin() , m_data.end() ) )
          { std::sort ( m_data.begin() , m_data.end() ) ; }
        this->check_me() ;
      }
      /// constructor to create complementary/ordinary ECDF
      ECDF
      ( const ECDF&  right         ,
        const bool   complementary ) ;
      /// copy constructor
      ECDF ( const ECDF&  right ) = default ;
      /// move constructor 
      ECDF (       ECDF&& right ) = default ;
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
        return this->add_sorted ( tmp.begin() , tmp.end() ) ;
      }
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
        return this->check_me() ;
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
      /// data size 
      inline Data::size_type N             () const { return m_data.size () ; } 
      /// data size 
      inline Data::size_type size          () const { return m_data.size () ; }
      /// number of effective entries
      inline Data::size_type nEff          () const { return m_data.size () ; }      
      // ======================================================================
      /// access to data
      inline const Data&     data          () const { return m_data         ; }
      // access to data 
      inline double          data   ( const unsigned int index ) const
      { return index < m_data.size() ?  m_data[index] : m_data.back() ; }
      // ======================================================================
      /// complementary?
      inline bool            complementary () const { return m_complementary ; }
      /// minimal x-value
      inline double          xmin          () const { return m_data.front () ; } 
      /// maximal x-value
      inline double          xmax          () const { return m_data.back  () ; }
      // ======================================================================
      /// get the abscissa value with the given index 
      inline double operator[] ( const unsigned int index ) const
      { return index < m_data.size() ? m_data[index] : m_data.back() ; }
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
    public:
      // ======================================================================
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
      /** check that ECDF is OK: there are some entries
       *  - remove elments <code>!std::isfinite</code>
       */
      ECDF& check_me () ; // check that ECDF is OK: there are some entries 
      // ======================================================================
    private:
      // ======================================================================
      /// container of sorted data 
      Data m_data          {       } ; // container of sorted data
      /// complementary CDF?
      bool m_complementary { false } ; // complementary CDF ? 
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
    class WECDF
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
      /// data size 
      inline Data::size_type N             () const { return m_data.size () ; } 
      /// data size 
      inline Data::size_type size          () const { return m_data.size () ; }
      /// number of effective entries
      inline double          nEff          () const { return m_sumw * m_sumw / m_sumw2 ;  }
      /// sum of all weights
      inline double          sumw          () const { return m_sumw  ; }      
      /// sum of all squared weights
      inline double          sumw2         () const { return m_sumw2 ; }      
      // ======================================================================      
      /// access to data
      inline const Data&     data          () const { return m_data         ; }
      // access to data 
      inline double data   ( const unsigned int index ) const
      { return index < m_data.size() ?  m_data[index].first  : m_data.back().first  ; }
      // access to weight 
      inline double weight ( const unsigned int index ) const
      { return index < m_data.size() ?  m_data[index].second : m_data.back().second ; }      
      // ======================================================================
      /// complementary?
      inline bool            complementary () const { return m_complementary ; }
      /// minimal x-value
      inline double          xmin          () const { return m_data.front ().first ; } 
      /// maximal x-value
      inline double          xmax          () const { return m_data.back  ().first ; }
      // ======================================================================
      /// get the abscissa value with the given index 
      inline const Entry& operator[] ( const unsigned int index ) const
      { return index < m_data.size() ? m_data[index] : m_data.back() ; }
      // ======================================================================      
      /// expose the data: begin iterator 
      inline iterator begin () const { return m_data.begin () ; }
      /// expose the data: end-iterator 
      inline iterator end   () const { return m_data.end   () ; }      
      // ======================================================================
    public: 
      // ======================================================================
      /// swap two objects 
      void swap ( WECDF& right ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// calculate \f$ \sum_i^{n} w_i \f$
      inline double calc_sumw ( const unsigned long n ) const
      { return std::accumulate  ( m_data.begin () ,
                                  m_data.begin () + std::min ( n , m_data.size() ) , 
                                  0.0             ,
                                  [] ( const double s , const Entry& entry ) -> double
                                  { return s + entry.second ; } ) ;	}
      /// calculate \f$ \sum_i^{n} w_i \f$
      inline double calc_sumw2 ( const unsigned long n ) const
      { return std::accumulate ( m_data.begin () ,
                                 m_data.begin () + std::min ( n , m_data.size() ) , 
                                 0.0             ,
                                 [] ( const double s , const Entry& entry ) -> double
                                 { return s + entry.second * entry.second ; } ) ; }
      // ======================================================================
      /// calculate \f$ \sum_i w_i   \f$
      inline double calc_sumw  () const { return calc_sumw  ( m_data.size() ) ; }
      /// calculate \f$ \sum_i w^2_i \f$
      inline double calc_sumw2 () const { return calc_sumw2 ( m_data.size() ) ; }
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
      /** check that WECDF is OK: 
       *  - there are some entries
       *  - sum of weigths is positive 
       *  - sum of squared weigths is positive 
       *  - remove elements <code>!std::isfinite</code>
       */
      WECDF& check_me () ; // check that ECDF is OK: there are some entries 
      // ======================================================================
    private: 
      // ======================================================================
      /// container of sorted data 
      Data m_data          {       } ; // container of sorted data
      /// container of weights 
      Data m_weights       {       } ; // container of weights
      /// sum of all weights
      double m_sumw        { 0     } ; // sum opf all weights 
      /// sum of all squared weights
      double m_sumw2       { 0     } ; // sum opf all squared weights 
      /// complementary CDF?
      bool m_complementary { false } ; // complementary CDF ? 
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
    // =======================================================================
    /** @class EPDF
     *  Helper utility to eatimate the PDF for emprical data using 
     *  Kernel estimators
     */
    class EPDF 
    {
    public:
      // =====================================================================
      /// create the emppirical PDF from empirical CDF 
      EPDF
      ( const ECDF&                                 ecdf  ,
	const Ostap::Math::DensityEstimator::Kernel k     ,
	const double                                h = 0 ) ;
      // =====================================================================
    public:
      // =====================================================================
      /// get the PDF 
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      double evaluate ( const double x ) const ;
      /// get the CDF 
      inline double cdf         ( const double x ) const { return m_cdf ( x ) ; }
      // =====================================================================
    public:
      // =====================================================================
      /// get ECDF 
      const  ECDF&   cdf () const { return m_cdf  ; }
      // =====================================================================
      /// get the kernel 
      inline Ostap::Math::DensityEstimator::Kernel kernel () const { return m_k ; }
      // =====================================================================
      /// get the smoothing parameter
      inline double h    () const { return m_h ; }
      // =====================================================================
    public:
      // =====================================================================
      /// update smoothing parameters  (
      bool setH      ( const double h ) ;
      /// set hew kernel
      bool setKernel ( const Ostap::Math::DensityEstimator::Kernel k ) ;
      // =====================================================================
    private:
      // =====================================================================
      ECDF                                  m_cdf ;
      Ostap::Math::DensityEstimator::Kernel m_k    { Ostap::Math::DensityEstimator::Epanechnikov } ;
      double                                m_h    { 0       } ; 
      // =====================================================================
    } ;
    // ========================================================================
    /** @class EPDF
     *  Helper utility to eatimate the PDF for emprical data using 
     *  Kernel estimators
     */
    class WEPDF 
    {
    public:
      // =====================================================================
      /// create the emppirical PDF from empirical CDF 
      WEPDF
      ( const WECDF&                                ecdf  ,
	const Ostap::Math::DensityEstimator::Kernel k     ,
	const double                                h = 0 ) ;
      // =====================================================================
    public:
      // =====================================================================
      /// get the PDF 
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      double evaluate ( const double x ) const ;
      /// get the CDF 
      inline double cdf         ( const double x ) const { return m_cdf ( x ) ; }
      // =====================================================================
    public:
      // =====================================================================
      /// get ECDF 
      const WECDF& cdf () const { return m_cdf  ; }
      // =====================================================================
      /// get the kernel 
      inline Ostap::Math::DensityEstimator::Kernel kernel () const { return m_k ; }
      // =====================================================================
      /// get the smoothing parameter
      inline double h  () const { return m_h ; }
      // =====================================================================
    public:
      // =====================================================================
      /// update smoothing parameters  (
      bool setH      ( const double h ) ;
      /// set hew kernel
      bool setKernel ( const Ostap::Math::DensityEstimator::Kernel k ) ;
      // =====================================================================
    private: 
      // =====================================================================
      WECDF                                 m_cdf ;
      Ostap::Math::DensityEstimator::Kernel m_k    { Ostap::Math::DensityEstimator::Epanechnikov } ;
      double                                m_h    { 0       } ; 
      // =====================================================================
    } ;    
    // =======================================================================    
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ECDF_H
// ============================================================================
