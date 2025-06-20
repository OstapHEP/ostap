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
// =============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace  Math
  {
    // ========================================================================
    /** @class Qantile 
     *  P2 algorithm for (approximate) quantile estimamtion
     *
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
     * 
     *  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
     */
    class  Quantile : public Ostap::Math::Statistic
    {
    public :
      // =======================================================================
      /// # of events 
      typedef unsigned long long  size_type ; 
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
      void reset  ()                     override { m_N = 0 ; } ; 
      // ======================================================================
    public :
      // ======================================================================
      /// add one more measurable, update quantiles    
      Quantile& add ( const double v ) ;
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
      /// minimal  value ( quantile for p=0)
      double min () const ;
      /// maximal value ( quantile for p=1)
      double max () const ; 	
      // ======================================================================
    public:
      // ======================================================================
      /// get the quantile value
      double quantile () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects
      void swap ( Quantile& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// initialiation strategy
      Initialization          m_init { Adaptive } ; // initialiation strategy
      /// quantile
      double                  m_p    { 0.5 } ; // quantile 
      /// sample size 
      size_type               m_N    { 0   } ; // sample size
      //
      std::array<double,5>    m_q    {     } ;
      std::array<double,5>    m_ns   {     } ;
      std::array<size_type,5> m_n    { 0 , 1 , 2 , 3 , 4 } ;
      // =======================================================================
    };
    // ========================================================================
    // swap two objects 
    inline void swap ( Quantile& a , Quantile& b ) { a.swap ( b ) ; } 
    // ========================================================================
    /** @class Quantiles 
     *  Extended P2 algorithm for (approximate) quantile estimamtion
     *  @see https://aakinshin.net/posts/ex-p2-quantile-estimator/
     */
    class Quantiles : public Ostap::Math::Statistic
    {
    public:
      // ======================================================================
      /// #of entri4s
      typedef typename Quantile::size_type  size_type ;
      // ======================================================================
    public:
      // ======================================================================
      /// from vector oof probabiliies 
      Quantiles
      ( const std::vector<double>& ps ) ;
      /// from soem range
      template <class ITERATOR,
		typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
		typename = std::enable_if<std::is_convertible<value_type,double>::value> >
      Quantiles
      ( ITERATOR first ,
	ITERATOR last  )
	: Quantiles ( std::vector<double> ( first , last ) )
      {}
      // ======================================================================
    public: // 
      // ======================================================================
      Quantiles& add ( const double value ) ;
      // ======================================================================
    public: // Ostap::Math::Statistic
      // =====================================================================
      // Generic counter interface
      void update ( const double value ) override { add  ( value ) ; }
      // reset quantile 
      void reset  ()                     override { m_N = 0 ; } ; 
      // ======================================================================
    public:
      // ======================================================================
      /// sample size 
      inline size_type    N     () const { return m_N  ; }
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
      double min () const ;
      /// maximal value  ( quantile for p = 1 )
      double max () const ; 	
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
    private:
      // ======================================================================
      /// vector of probabilities 
      std::vector<double>    m_p  {   } ; // vector of probabilities
      /// sample size 
      size_type              m_N  { 0 } ; // sample size
      //
      std::vector<double>    m_q  {} ;
      std::vector<double>    m_ns {} ;
      std::vector<size_type> m_n  {} ;
      // ======================================================================
    } ;
    // ========================================================================
    // swap two objects 
    inline void swap ( Quantiles& a , Quantiles& b ) { a.swap ( b ) ; } 
    // ========================================================================
  } //                                         The end od namesapce Ostap::Math
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_P2QUANTILE_H
// ============================================================================
