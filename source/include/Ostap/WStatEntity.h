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
  /** @class WStatEntity  Ostap/WStatEntity.h
   *  Statistics with "weight"
   *  @see Ostap::StatEntity 
   *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
   *  @date   2014-04-07
   */
  class WStatEntity 
  {
  public:
    // ======================================================================
    /// empty constructor 
    WStatEntity () = default ;
    /// constructor from StatEntity of values 
    WStatEntity					
    ( const StatEntity& values ) ;
    /// full consructor
    WStatEntity
    ( const double      mu      ,
      const double      mu2     ,
      const StatEntity& values  ,
      const StatEntity& weights ) ;
    // ======================================================================
  public: // the basic getters 
    // ======================================================================
    /// total number of entries 
    unsigned long long       n   () const { return m_weights.n () ; }
    /// the first weighted moment/mean-value 
    double                   mu  () const { return m_mu  ; }
    /// the second central weighted moment/dispersion/variance  
    double                   mu2 () const { return m_mu2 ; }
    // ======================================================================
  public: // derived getters and aliases 
    // ======================================================================
    /// empty ?
    bool                     empty      () const { return m_weights.empty()  ; }
    /// get the actual number of entries 
    unsigned long long       nEntries   () const { return n () ; }
    /// mean-value 
    double                   mean       () const { return m_mu  ; }
    /// error im mean-value 
    double                   meanErr    () const ;  
    /// dispersion 
    double                   dispersion () const { return m_mu2 ; }
    /// variance 
    double                   variance   () const { return m_mu2 ; }
    /// RMS 
    double                   rms        () const ;
    /// get the effective number of entries 
    double                   nEff       () const ;
    /// get number of "good" (non-zero) entries
    unsigned long long       nGood      () const { return m_values. n   () ; }
    /// minimal value (for non-zero weights) 
    double                   min        () const { return m_values. min () ; }
    /// maximal value (for non-zero weights) 
    double                   max        () const { return m_values. max () ; }       
    // ======================================================================
  public: // helper sums 
    // ======================================================================
    /// sum_i weight_i*value_i
    double                   sum   () const ; // sum_i weight_i * value_i
    /// sum_i weight_i*value_i**2
    double                   sum2  () const ; // sum_i weight_i * value_i**2
    /// sum_i weight_i  
    double                   sumw  () const {  return m_weights.sum  () ; }
    /// sum_i weight_i*wight_i   
    double                   sumw2 () const {  return m_weights.sum2 () ; }
    // ======================================================================
  public:  // statistics &  weights and values 
    // ======================================================================
    /// get the statistic of weights 
    const Ostap::StatEntity& weights () const { return m_weights ; }
    /// get the statistic of values with non-zero weight 
    const Ostap::StatEntity& values  () const { return m_values  ; }
    // ======================================================================
  public:
    // ======================================================================
    /// add   the  value      (with weight=1) to the counter 
    WStatEntity& operator+= ( const double value ) { return add (  value ) ; }
    /// add   the  value +1   (with weight=1) to the counter 
    WStatEntity& operator++ ()      { return   (*this) += 1 ; }
    /// add   the  value +1   (with weight=1) to the counter 
    WStatEntity& operator++ ( int ) { return ++(*this)      ; }    
    /// add   the -value      (with weight=1) to the counter 
    WStatEntity& operator-= ( const double value ) { return add ( -value ) ; }
    /// add   the  value -1   (with weight=1) to the counter 
    WStatEntity& operator-- ()      { return   (*this) -= 1 ; }
    /// add   the  value -1   (with weight=1) to the counter 
    WStatEntity& operator-- ( int ) { return --(*this)      ; }    
    // ======================================================================
  public: // ordering operators 
    // ======================================================================
    /// basic comparison
    bool operator< ( const WStatEntity& s ) const ;
    bool operator==( const WStatEntity& s ) const ;
    /// derived comparisons:
    bool operator> ( const WStatEntity& s ) const { return s < *this ; }
    bool operator<=( const WStatEntity& s ) const { return    (*this) == s || (*this) < s ; }
    bool operator>=( const WStatEntity& s ) const { return    (*this) == s || (*this) > s ; }
    bool operator!=( const WStatEntity& s ) const { return   !(*this  == s ) ; }
    // ======================================================================
  public: // various technical helper methods  
    // ========================================================================
    /// reset the counters
    void reset () ;
    /// swap two cunters
    void swap ( WStatEntity& right ) ;
    /// representation as string
    std::string   toString   () const;
    /// printout  to std::ostream
    std::ostream& fillStream ( std::ostream& o ) const ;
    // ======================================================================
    // all finite values ?
    inline bool isfinite () const
    {
      return
	std::isfinite ( m_mu  ) &&
	std::isfinite ( m_mu2 ) &&
	m_values .isfinite ()   &&
	m_weights.isfinite ()   ;
    }
    // ======================================================================
  public: // the main methods 
    // ======================================================================
    /** the only one important method: add value with the given weight 
     *  @param value  value to be added 
     *  @param weight associated weight 
     *  @attention non-finite values  are ignored! 
     *  @attention non-finite weights are ignored! 
     */
    WStatEntity& add   
    ( const double value      ,  
      const double weight = 1 ) ;
    /// ditto
    WStatEntity& update 
    ( const double value      ,  
      const double weight = 1 ) { return add ( value , weight ) ; }
    /** add another counter 
     *  @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    WStatEntity& add ( const WStatEntity& stat ) ;
    /// add another counter 
    WStatEntity& add ( const  StatEntity& stat ) 
    { return add ( WStatEntity (   stat ) ) ; }
    /** add another counter 
     *  @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    WStatEntity& update ( const WStatEntity& stat ) { return add ( stat ) ; }
    /// add another counter 
    WStatEntity& update ( const  StatEntity& stat ) { return add ( stat ) ; }
    // ======================================================================
  public:
    // ======================================================================
    /** add another counter 
     *  @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    WStatEntity& operator+= ( const WStatEntity& other ) { return add ( other ) ; }
    /// add another counter 
    WStatEntity& operator+= ( const  StatEntity& other ) { return add ( other ) ; }
    // ======================================================================
  private: /// the basic quantities 
    // ======================================================================
    /// The  first weighted moment/mean value
    double m_mu  { 0 } ; // The  first weighted moment/mean value
    double m_mu2 { 0 } ; // The second central weighted moment/variance/dispersion
    // ======================================================================     
  private: /// helper statistics for the values and  weights 
    // ======================================================================
    /// statistic of values with non-zero weight 
    Ostap::StatEntity m_values  {} ;
    /// statistic of weights 
    Ostap::StatEntity m_weights {} ;
    // ======================================================================
  };
  // ========================================================================
  /// add two counters 
  inline WStatEntity operator+ ( WStatEntity a , const WStatEntity& b ) 
  { a += b ; return a ; }
  // ========================================================================
  inline std::ostream& operator<<( std::ostream& s, const WStatEntity& e )
  { return e.fillStream ( s ) ; }
  // ==========================================================================
  /// conversion to string 
  inline std::string to_string ( const WStatEntity& e ) { return e.toString() ;}
  // ==========================================================================
  /// swap two counters 
  inline void swap ( WStatEntity&  a , WStatEntity&  b  ) { a.swap( b ) ; }
  // =========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_WSTATENTITY_H
// ============================================================================
