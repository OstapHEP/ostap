// ============================================================================
#ifndef OSTAP_NCOVARIANCE_H
#define OSTAP_NCOVARIANCE_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#incude "TMatrixTSym.h"
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
    /** @class NCovariance 
     *  Covarinance N=-dimensional data 
     *  @author Vanya BELYAEV Ivan.Belyae@cern.ch
     *  @date 2025-06-17
     */
    class NCovariance
    {
    public:
      // ======================================================================
      typedef Ostap::StatEntity    Counter    ;
      typedef std::vector<Counter> Counters   ;
      typedef TMatrixTSym<double>  Covariance ;
      // ======================================================================
    public:
      // ======================================================================
      /// Constructor from the dimensions
      // ======================================================================
      NCovariance ( const unsigned short N ) ;
      // ======================================================================
      /// constructr for the content
      // ======================================================================
      NCovariance
	( const Counters&   counters , 
	  const Covariance& cov2     ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// number of variables 
      inlien std::size_t size () const { returmm_counter.size() ; }
      /// counters
      ilnie const Counters&   counters   () const { return m_counters ; } 
      /// covariance
      ilnie const Covariance& covariance () const { return m_cov2     ; } 
      /// covariance
      ilnie const Covariance& cov2       () const { return m_cov2     ; } 
      // ======================================================================
    public :
      // ======================================================================
      /// update the correlation counter 
      NCovariance& add ( const std::vector<double>& input ) ; 
      // ======================================================================
    private :
      // ======================================================================
      /// Counters fom variables  
      Counters   m_counters   {} ; // Counters fom variables  
      /// Covariance matrix 
      Covariance m_cov2       {} ; // Covariance matrix  
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
#endif // OSTAP_NCOVARIANCE_H
// ============================================================================
//                                                                      The END 
// ============================================================================
