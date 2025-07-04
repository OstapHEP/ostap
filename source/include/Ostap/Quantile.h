// ============================================================================
#ifndef OSTAP_QUANTILE_H 
#define OSTAP_QUANTILE_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <array>
#include <vector>
// ============================================================================
// Ostap
// =============================================================================
#include "Ostap/Statistic.h"
#include "Ostap/StatEntity.h"
// =============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace  Math
  {
    // ========================================================================
    /** @class Quantile 
     *  P2 algorithm for (approximate) quantile estimamtion
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
     * 
     *  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
     */
    class Quantile : public Ostap::Math::Statistic
    {
    public :
      // =======================================================================
      /// the type of helper counter 
      typedef typename Ostap::StatEntity  Counter   ;
      /// #of entries 
      typedef typename Counter::size_type size_type ; 
      // =======================================================================      
    public :
      // =======================================================================
      /// initialization strategy 
      enum Initialization
    	{
          Classic = 0 ,
          Adaptive 
        } ;
      // =======================================================================      
    public :
      // =======================================================================
      /** constructor from  p
       *  @parameter p (INPUT) p: \f$ 0 \le p \le 1 \f$  
       *  - median : p = 0.5
       */  
      Quantile
      ( const double         p = 0.5      ,
        const Initialization s = Adaptive ) ;
      // ======================================================================
    public: // Ostap::Math::Statistic
      // =====================================================================
      // Generic counter interface
      void update ( const double value ) override { add  ( value ) ; }
      // reset quantile 
      void reset  ()                     override ; 
      // ======================================================================
    public :
      // ======================================================================
      /// add one more measurable, update quantiles    
      Quantile& add ( const double v ) ;
      /// add a range of values
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      Quantile& add
      ( ITERATOR first ,
        ITERATOR last  )
      { for ( ; first != last ; ++first ) { add ( *first ) ; } ; return *this ; }
      // =====================================================================
    public:
      // ======================================================================
      /// sample size 
      inline size_type N          () const { return m_N  ; }
      /// sample size 
      inline size_type size       () const { return N () ; }
      /// valid counter?
      inline bool      valid      () const { return m_N  ; } 
      /// valid counter?
      inline bool      ok         () const { return m_N  ; }
      /// probability
      inline double    p          () const { return m_p  ; } 
      /// probability
      inline double    probabilty () const { return m_p  ; } 
      // ======================================================================
    public:
      // ======================================================================
      /// minimal value ( == quantile for  p = 0 )
      inline double min () const { return m_counter.min () ; } 
      /// maximal value ( == quantile for  p = 1 )
      inline double max () const { return m_counter.max () ; } 
      // ======================================================================
    public:
      // ======================================================================
      /// get the triplet: { min , quantile , max }
      std::array<double,3> quantiles () const ; 
      /// get the quantile value
      double quantile () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects
      void swap ( Quantile& right ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// helper counter 
      inline const Counter& counter() const { return m_counter ; }
      // ======================================================================
    private:
      // ======================================================================
      /// initialiation strategy
      Initialization          m_init { Adaptive } ; // initialiation strategy
      /// quantile
      double                  m_p       { 0.5 } ; // quantile 
      /// sample size 
      size_type               m_N       { 0   } ; // sample size
      //
      std::array<double,5>    m_q       {     } ;
      std::array<double,5>    m_ns      {     } ;
      std::array<size_type,5> m_n       { 0 , 1 , 2 , 3 , 4 } ;
      // ======================================================================= 
    public:
      // =======================================================================
      /// helper counter (non needed in the original algorithm)
      Counter                 m_counter {} ;
      // =======================================================================
    };
    // ========================================================================
    // swap two objects 
    inline void swap ( Quantile& a , Quantile& b ) { a.swap ( b ) ; } 
    // ========================================================================
    /** @class Quantiles 
     *  Extended P2 algorithm for (approximate) quantile estimation
     *  @see https://aakinshin.net/posts/ex-p2-quantile-estimator
     */
    class Quantiles : public Ostap::Math::Statistic
    {
    public:
      // ======================================================================
      /// helper countter 
      typedef typename Quantile::Counter  Counter   ;
      /// #of entri4s
      typedef typename Counter::size_type size_type ;
      // ======================================================================
    public:
      // ======================================================================
      /// from vector of probabiliies 
      Quantiles
      ( const std::vector<double>& ps ) ;
      /// get N-quantiles 
      Quantiles ( const std::size_t N ) ; 
      /// from soem range
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >      
      Quantiles
      ( ITERATOR first ,
        ITERATOR last  )
        : Quantiles ( std::vector<double> ( first , last ) )
      {}
      // ======================================================================
    public: // 
      // ======================================================================
      Quantiles& add ( const double value ) ;
      /// add a rageg of values
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >            
      Quantiles& add
      ( ITERATOR first ,
        ITERATOR last  )
      { for ( ; first != last ; ++first ) { add ( *first ) ; } ; return *this ; }
      // ======================================================================
    public: // Ostap::Math::Statistic
      // =====================================================================
      // Generic counter interface
      void update ( const double value ) override { add  ( value ) ; }
      // reset quantile 
      void reset  ()                     override ; 
      // ======================================================================
    public:
      // ======================================================================
      /// sample size 
      inline size_type   N     () const { return m_N  ; }
      /// sample size 
      inline size_type   size  () const { return N () ; }
      /// valid counter?
      inline bool        valid () const { return m_N  ; } 
      /// valid counter?
      inline bool        ok    () const { return m_N  ; }
      /// # of probabilities 
      inline std::size_t NP    () const { return m_p.size() ; } 
      /// # of quantiles
      inline std::size_t NQ    () const { return m_p.size() ; } 
      // ======================================================================
    public:
      // ======================================================================
      /// minimal  value ( quantile for p = 0 )
      inline double min () const { return m_counter.min () ; }
      /// maximal value  ( quantile for p = 1 )
      inline double max () const { return m_counter.max () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the probabilty value
      inline double probability ( const std::size_t index ) const
      { return index < m_p.size() ? m_p [ index ] : 1.0 ; }
      /// get the probabilty value
      inline double p           ( const std::size_t index ) const
      { return probability ( index ) ; }      
      /// get the quantile value
      double quantile ( const std::size_t index ) const ;
      /// get all quantiles 
      std::vector<double> quantiles () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects
      void swap ( Quantiles& right ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// helper counter 
      inline const Counter& counter() const { return m_counter ; }
      // ======================================================================
    private:
      // ======================================================================
      /// vector of probabilities 
      std::vector<double>    m_p  {   } ; // vector of probabillities 
      /// sample size 
      size_type              m_N  { 0 } ; // sample size
      //
      std::vector<double>    m_q  {   } ;
      std::vector<double>    m_ns {   } ;
      std::vector<size_type> m_n  {   } ;
      // ======================================================================
    private: 
      // ======================================================================
      /// helper counter (non needed in the original algorithm)
      Counter                m_counter {} ; // helper counter 
      // =======================================================================
    } ;
    // ========================================================================
    // swap two objects 
    inline void swap ( Quantiles& a , Quantiles& b ) { a.swap ( b ) ; } 
    // ========================================================================
    /** @class Quantiles_
     *  Extended P2 algorithm for (approximate) equidistant quantiles 
     *  @see https://aakinshin.net/posts/ex-p2-quantile-estimator
     */
    template <unsigned int N,
              typename std::enable_if<(1<=N),bool>::type = true >
    class Quantiles_ : public Ostap::Math::Statistic
    {
      // ====================================================================
    public :
      // ====================================================================
      /// update the counter 
      Quantiles_& add  ( const  double v ) { m_qs.add ( v ) ; return *this ; }
      /// add a range of of values
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      Quantiles_& add
      ( ITERATOR first ,
        ITERATOR last  )
      { for ( ; first != last ; ++first ) { this->add ( *first ) ; } ; return *this ; }
      // ====================================================================
    public: // Ostap::Math::Statistic
      // =====================================================================
      // Generic counter interface
      void update ( const double value ) override { add  ( value ) ; }
      // reset quantile 
      void reset  ()                     override { m_qs.reset ()  ; } 
      // ======================================================================
    public :
      // ======================================================================
      std::array<double,N+2> quantiles ()  const 
      {
        const std::vector<double> _output { m_qs.quantiles () } ;
        std::array<double,N+2> result ;
        /// 
        const std::size_t  o = _output.size() ; 
        const std::size_t  L = o < N + 2 ? o : N + 2 ;
        for  ( std::size_t i = 0 ; i < L ; ++i )
          { result [ i ] = _output [ i ] ; }  
        return result ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// minimal  value ( quantile for p = 0 )
      inline double min () const { return m_qs.min () ; }
      /// maximal value  ( quantile for p = 1 )
      inline double max () const { return m_qs.max () ; }
      // ======================================================================
    private:
      // ======================================================================
      /// actual quantile counter  
      Quantiles m_qs { N } ; // actual quantile counter 
      // ======================================================================
    }; //                                             The end of class Qantile_ 
    // ========================================================================
    /// get minmax:  trivial "quantiles" for p=0 & p=1 
    typedef Quantiles_<1>   QMinMax      ;
    ///              get ( min, median, max)    
    typedef Quantiles_<2>   QMedian      ;
    /// Terciles:    get ( min, t1 , t2 , max )    
    typedef Quantiles_<3>   QTerciles    ;
    /// Quartiles:   get ( min, q1 , median  , q3   , max )        
    typedef Quantiles_<4>   QQuartiles   ;
    /// Quintiles:   get ( min, q1 , q2 , a3 , q4   , max )            
    typedef Quantiles_<5>   QQuintiles   ;
    /// Sextiles:    get ( min, q1 , q2 , ... , q5  , max )                
    typedef Quantiles_<6>   QSextiles    ;
    /// Septiles:    get ( min, q1 , q2 , ... , q6  , max )            
    typedef Quantiles_<7>   QSeptiles    ;
    /// Octiles:     get ( min, q1 , q2 , ... , q7  , max )                
    typedef Quantiles_<8>   QOctiles     ;
    /// Deciles:     get ( min, q1 , q2 , ... , q9  , max )                
    typedef Quantiles_<10>  QDeciles     ;
    /// Ventiles:    get ( min, q1 , q2 , ... , q19 , max )                
    typedef Quantiles_<20>  QVentiles    ;
    /// Percentiles: get ( min, q1 , q2 , ... , q99 , max )                
    typedef Quantiles_<100> QPercentiles ;
    // ========================================================================
  } //                                         The end od namesapce Ostap::Math
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
#endif // OSTAP_QUANTILE_H
// ============================================================================
//                                                                      The END 
// ============================================================================
