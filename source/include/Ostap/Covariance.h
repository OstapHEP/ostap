// ============================================================================
#ifndef OSTAP_COVARIANCE_H
#define OSTAP_COVARIANCE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/SymmetricMatrixTypes.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Covariance
     *  Counter for two variables which also counts the covariance 
     *  @see Ostap::StatEntity
     */
    class Covariance
    {
      // ======================================================================
    public: 
      // ======================================================================
      /// the actual type of counter 
      typedef Ostap::StatEntity                       Counter  ;
      /// covarinace/correlation matrices 
      typedef Ostap::SymMatrix2x2                     Matrix   ;
      /// #of entries 
      typedef Counter::size_type                      size_type ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor 
      Covariance () = default ;
      // ======================================================================
      // constructor from two counters and the correlation coefficient 
      Covariance
      ( const Counter& c1       , 
      	const Counter& c2       ,
	      const double   corr = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// the first counter 
      const Counter& counter1   () const { return m_cnt1   ; }
      /// the second counter 
      const Counter& counter2   () const { return m_cnt2   ; }
      // ======================================================================
      /// get the moment \f$ \sum_i (x_i - \bar{x} ) ( y_i - \bar{y}) \f$
      inline double cov2m       () const  { return m_cov2m ; }
      /// get the true covariance 
      inline double covariance  () const  { return empty() ? 0.0 : m_cov2m / n() ; }
      /// get the correlation coefficient 
      double        correlation () const ;
      // ======================================================================
      /// number of entries 
      inline size_type  n     () const { return m_cnt1.n     () ; }
      /// effective number of entries 
      inline size_type  nEff  () const { return m_cnt1.nEff  () ; }
      ///  number of "good" (non-zero) entries
      inline size_type  nGood () const { return m_cnt1.nGood () ; }
      /// empty ?
      inline bool       empty () const { return m_cnt1.empty () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add two values to the counters 
      Covariance& add	
      ( const double x ,
      	const double y ) ;
      /// add another counter 
      Covariance& add ( const Covariance& right ) ;
      // ======================================================================
      /// add another counter 
      inline Covariance& operator+= ( const Covariance& right )
      { return add ( right ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add x,y
      inline void update
      ( const double x ,
	      const double y ) { add ( x , y ) ; }
      // ======================================================================
      /// reset counters 
      void reset () ;
      // ======================================================================
      // everything is finite?
      inline bool isfinite () const
      { return std::isfinite ( m_cov2m ) && m_cnt1.isfinite() && m_cnt2.isfinite() ; }
      // ======================================================================
    private:
      // ======================================================================
      Ostap::StatEntity m_cnt1  {   } ;
      Ostap::StatEntity m_cnt2  {   } ;
      double            m_cov2m { 0 } ;
      // ======================================================================      
    };
    // ========================================================================      
    /// external operator for addition of two covariance objects
    inline Covariance operator+ ( Covariance a , const Covariance& b ) 
    { a += b ; return a ; }
    // ========================================================================
    /// get the covariance matrix
    Covariance::Matrix   covariance   ( const Covariance& ) ;
    // ========================================================================
    /// get the correlation matrix
    Covariance::Matrix   correlation  ( const Covariance& ) ;
    // ========================================================================          
    /** @class WCovariance
     *  Counter for two variables which also counts the covariance 
     *  @see Ostap::WStatEntity
     */
    class WCovariance
    {
      // ======================================================================
    public: 
      // ======================================================================
      /// the actual type of counter 
      typedef Ostap::WStatEntity                Counter  ;
      /// covarinace/correlation matrices 
      typedef Ostap::Math::Covariance::Matrix   Matrix   ;
      /// #of entries 
      typedef Counter::size_type                size_type ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor
      WCovariance () = default ;
      // ======================================================================
      // constructor from two counters and the correlation coefficient
      WCovariance
      ( const Counter& c1       , 
	      const Counter& c2       ,
	      const double   corr = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// the first counter 
      const Counter& counter1   () const { return m_cnt1   ; }
      /// the second counter 
      const Counter& counter2   () const { return m_cnt2   ; }
      // ======================================================================
      /// get the moment \f$ \sum_i (x_i - \bar{x} ) ( y_i - \bar{y}) \f$
      inline double cov2m       () const  { return m_cov2m ; }
      /// get the true covariance 
      inline double covariance  () const  { return empty() ? 0.0 : m_cov2m / w () ; }
      /// get the correlation coefficient 
      double        correlation () const ;
      // ======================================================================
      /// number of entries 
      inline size_type  n     () const { return m_cnt1.n     () ; }
      /// effective number of entries 
      inline double     nEff  () const { return m_cnt1.nEff  () ; }
      ///  number of "good" (non-zero) entries
      inline size_t     nGood () const { return m_cnt1.nGood () ; }
      /// empty ?
      inline bool       empty () const { return m_cnt1.empty () ; }
      /// sum of weights 
      inline double     w     () const { return m_cnt1.sumw  () ; }
      /// sum of weights 
      inline double     sumw  () const { return m_cnt1.sumw  () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add two values to the counters 
      WCovariance& add	
      ( const double x     ,
        const double y     ,
        const double w = 1 ) ; 
      /// add another counter 
      WCovariance& add ( const WCovariance& right ) ;
      // ======================================================================
      /// add another counter 
      inline WCovariance& operator+= ( const WCovariance& right )
      { return add ( right ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// add x,y
      inline void update
      ( const double x     ,
        const double y     ,
        const double w = 1 ) { add ( x , y , w ) ; }
      // ======================================================================
      /// reset counters 
      void reset () ;
      // ======================================================================
      // everything is finite?
      inline bool isfinite () const
      { return std::isfinite ( m_cov2m ) && m_cnt1.isfinite() && m_cnt2.isfinite() ; }
      // ======================================================================
    private:
      // ======================================================================
      Ostap::WStatEntity m_cnt1  {   } ;
      Ostap::WStatEntity m_cnt2  {   } ;
      double             m_cov2m { 0 } ;
      // ======================================================================      
    };
    // ========================================================================      
    /// external operator for addition of two covariance objects
    inline WCovariance operator+ ( WCovariance a , const WCovariance& b ) 
    { a += b ; return a ; }
    // ========================================================================
    /// get the covariance matrix
    WCovariance::Matrix    covariance   ( const WCovariance& ) ;
    /// get the correlation matrix
    WCovariance::Matrix    correlation  ( const WCovariance& ) ;
    // ========================================================================
  } //                                         The end of nameapace Ostap::Math
  // ==========================================================================
} //                                                    The end namespace Ostap
// ============================================================================
#endif // OSTAP_COVARIANCE_H
// ============================================================================
