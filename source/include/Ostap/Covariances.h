// ============================================================================
#ifndef OSTAP_COVARIANCES_H
#define OSTAP_COVARIANCES_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "TMatrixTSym.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Covariances 
     *  Covariance N=-dimensional data 
     *  @author Vanya BELYAEV Ivan.Belyae@cern.ch
     *  @date 2025-06-17
     */
    class Covariances
    {
    public:
      // ======================================================================
      typedef Ostap::StatEntity    Counter    ;
      typedef std::vector<Counter> Counters   ;
      typedef TMatrixTSym<double>  CovMatrix  ;
      typedef Counter::size_type   size_type  ;      
      // ======================================================================
    public:
      // ======================================================================
      /// Constructor from the dimension (2<=N)
      Covariances ( const unsigned short N = 2 ) ;
      // ======================================================================
      /// constructor from  the content
      Covariances
      ( const Counters&  counters , 
        const CovMatrix& cov2     ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// number of variables 
      inline std::size_t       N           () const { return m_counters.size() ; }
      /// counters
      inline const Counters&   counters    () const { return m_counters ; } 
      /// covariance
      inline const CovMatrix&  covariance  () const { return m_cov2     ; } 
      /// covariance
      inline const CovMatrix&  cov2        () const { return m_cov2     ; } 
      // ======================================================================
    public :   // derived statistics 
      // =======================================================================
      /// number of entries
      inline size_type n    () const { return m_counters.front().n    ()  ; }
      /// effective number of entries
      inline size_type nEff () const { return m_counters.front().nEff () ; }
      // ======================================================================
    public:  // the main operation 
      // ======================================================================
      /// update the correlation counter 
      Covariances& add ( const std::vector<double>& input ) ; 
      // ======================================================================
    public: // add to conuters of the sam size 
      // ======================================================================
      Covariances& add
      ( const Covariances& right ) ;
      // ======================================================================
      Covariances& operator+=( const Covariances& right )
      { return add ( right ) ; }
      // ======================================================================      
    public:
      // ======================================================================
      Covariances& __iadd__ ( const Covariances& right )
      { return add ( right ) ; }
      // ======================================================================
      Covariances  __add__  ( const Covariances& right ) const
      { Covariances result { *this } ; result.add ( right ) ; return result ; }      
      // ======================================================================
    public: 
      // ======================================================================
      /// swap two objects
      void swap ( Covariances& c2 ) ;
      // ======================================================================
    private :
      // ======================================================================
      /// Counters fom variables  
      Counters   m_counters   { 2 } ; // Counters fom variables  
      /// Covariance matrix 
      CovMatrix  m_cov2       { 2  } ; // Covariance matrix  
      // ======================================================================
    private:
      // ======================================================================
      mutable std::vector<double>  m_delta {} ;
      // ======================================================================
    } ;
    // ========================================================================      
    /// external operator for addition of two covariance objects
    inline Covariances operator+( Covariances a , const Covariances& b ) 
    { a += b ; return a ; }
    // ========================================================================
    /// swap two objects 
    inline void swap
    ( Covariances& a ,
      Covariances& b ) { a.swap ( b ) ; } 
    // ========================================================================
    /** @class WCovariances 
     *  Covariance for  N-dimensional data (with weight) 
     *  @author Vanya BELYAEV Ivan.Belyae@cern.ch
     *  @date 2025-06-17
     */
    class WCovariances
    {
    public:
      // ======================================================================
      typedef Ostap::WStatEntity   Counter    ;
      typedef std::vector<Counter> Counters   ;
      typedef TMatrixTSym<double>  CovMatrix  ;
      typedef Counter::size_type   size_type  ;
      // ======================================================================
    public:
      // ======================================================================
      /// Constructor from the dimension (2<=N)
      WCovariances ( const unsigned short N = 2 ) ;
      // ======================================================================
      /// constructor from  the content
      WCovariances
      ( const Counters&  counters , 
        const CovMatrix& cov2     ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// number of variables 
      inline std::size_t       N           () const { return m_counters.size() ; }
      /// counters
      inline const Counters&   counters    () const { return m_counters ; } 
      /// covariance
      inline const CovMatrix&  covariance  () const { return m_cov2     ; } 
      /// covariance
      inline const CovMatrix&  cov2        () const { return m_cov2     ; } 
      // ======================================================================
    public :   // derived statistics 
      // ======================================================================
      /// number of entries
      inline size_type n    () const { return m_counters.front().n    ()  ; }
      /// effective number of entries
      inline double    nEff () const { return m_counters.front().nEff () ; }
      /// sum of weights 
      inline double    sumw () const { return m_counters.front().sumw () ; }
      // ======================================================================
    public:  // the main operation 
      // ======================================================================
      /// update the correlation counter 
      WCovariances& add 
      ( const std::vector<double>& input      , 
        const  double              weight = 1) ; 
      // ======================================================================
    public: // add to conuters of the sam size 
      // ======================================================================
      WCovariances& add
      ( const WCovariances& right ) ;
      // ======================================================================
      WCovariances& operator+=( const WCovariances& right )
      { return add ( right ) ; }
      // ======================================================================      
    public:
      // ======================================================================
      WCovariances& __iadd__ ( const WCovariances& right )
      { return add ( right ) ; }
      // ======================================================================
      WCovariances  __add__  ( const WCovariances& right ) const
      { WCovariances result { *this } ; result.add ( right ) ; return result ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects
      void swap ( WCovariances& c2 ) ;
      // ======================================================================
    private :
      // ======================================================================
      /// Counters fom variables  
      Counters   m_counters   { 2 } ; // Counters fom variables  
      /// Covariance matrix 
      CovMatrix  m_cov2       { 2 } ; // Covariance matrix  
      // ======================================================================
    private:
      // ======================================================================
      mutable std::vector<double>  m_delta {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /// external operator for addition of two covariance objects
    inline WCovariances operator+( WCovariances a , const WCovariances& b ) 
    { a += b ; return a ; }
    // ========================================================================
    /// swap two objects 
    inline void swap
    ( WCovariances& a ,
      WCovariances& b ) { a.swap ( b ) ; } 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
#endif // OSTAP_NCOVARIANCE_H
// ============================================================================
//                                                                      The END 
// ============================================================================
