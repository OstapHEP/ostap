// ============================================================================
#ifndef OSTAP_WSTATENTITY_H 
#define OSTAP_WSTATENTITY_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
// ============================================================================
namespace Ostap
{
  // ========================================================================
  /** @class WStatEntity  WStatEntity.h
   *  Statistic with weight
   *  @see StatEntity 
   *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
   *  @date   2014-04-07
   */
  class WStatEntity 
  {
  public:
    // ======================================================================
    /// empty constructor 
    WStatEntity () ;
    /// copy constructor 
    WStatEntity ( const WStatEntity& ) ;
    /// constructor from StatEntity of values 
    WStatEntity ( const StatEntity& values ) ;
    // ======================================================================
  public:
    // ======================================================================
    /// the only one important method: add value with the given weight 
    WStatEntity& add   
    ( const double value      ,  
      const double weight = 1 ) ;
    /// ditto
    WStatEntity& update 
    ( const double value      ,  
      const double weight = 1 ) { return add ( value , weight ) ; }
    // ======================================================================
    WStatEntity& operator+= ( const double value ) { return add ( value ) ; }
    // ======================================================================
    // reset statistic 
    void reset() ;
    // ======================================================================      
  public: // the basic getters 
    // ======================================================================
    /// get the mean-value 
    double             mean       () const ;  
    double             meanErr    () const ;  
    double             dispersion () const ;
    double             rms        () const ;
    /// get the effective number of entries 
    double             nEff       () const ;
    /// get the actual number of entries 
    unsigned long long nEntries   () const { return m_weights.nEntries() ; }
    // ======================================================================
    /// sum_i weight_i*value_i
    double        sum  () const { return m_sum  ; } // sum_i weight_i*value_i
    /// sum_i weight_i*value_i**2
    double        sum2 () const { return m_sum2 ; } // sum_i weight_i*value_i**2
    // ======================================================================
  public:
    // ======================================================================
    /// get the statistic of weights 
    const Ostap::StatEntity& weights () const { return m_weights ; }
    /// get the statistic of values with non-zero weight 
    const Ostap::StatEntity& values  () const { return m_values  ; }
    // ======================================================================
  public: // printout 
    // ======================================================================
    std::ostream& fillStream ( std::ostream& s ) const ;      
    std::string   toString   () const ;
    // ======================================================================
  private: /// the basic statistic 
    // ======================================================================
    /// sum_i weight_i*value_i
    double m_sum    ;  // sum_i weight_i*value_i
    /// sum_i weight_i*value_i**2 
    double m_sum2   ;  // sum_i weight_i*value_i**2 
    // ======================================================================     
  private: /// helper statistics 
    // ======================================================================
    /// statistic of values with non-zero weight 
    Ostap::StatEntity m_values  ;
    /// statistic of weights 
    Ostap::StatEntity m_weights ;
    // ======================================================================
  };
  // ========================================================================
  inline std::ostream& operator<<( std::ostream& s, const WStatEntity& e )
  { return e.fillStream ( s ) ; }
  // ==========================================================================
  /// conversion to string 
  inline std::string to_string ( const WStatEntity& e ) { return e.toString() ;}
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_WSTATENTITY_H
// ============================================================================
