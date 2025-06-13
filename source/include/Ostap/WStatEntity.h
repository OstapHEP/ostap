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
  class WStatEntity : public Ostap::Math::WStatistic
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
    inline unsigned long long       n   () const { return m_weights.n () ; }
    /// the first weighted moment/mean-value 
    inline double                   mu  () const { return m_mu  ; }
    /// the second central weighted moment/dispersion/variance  
    inline double                   mu2 () const { return m_mu2 ; }
    // ======================================================================
  public: // derived getters and aliases 
    // ======================================================================
    /// empty ?
    inline bool                     empty      () const { return m_weights.empty()  ; }
    /// get the actual number of entries 
    inline unsigned long long       nEntries   () const { return n () ; }
    /// mean-value 
    inline double                   mean       () const { return m_mu  ; }
    /// error im mean-value 
    double                          meanErr    () const ;  
    /// dispersion 
    inline double                   dispersion () const { return m_mu2 ; }
    /// variance 
    inline double                   variance   () const { return m_mu2 ; }
    /// RMS 
    double                          rms        () const ;
    /// get the effective number of entries 
    double                          nEff       () const ;
    /// get number of "good" (non-zero) entries
    inline unsigned long long       nGood      () const { return m_values. n   () ; }
    /// minimal value (for non-zero weights) 
    inline double                   min        () const { return m_values. min () ; }
    /// maximal value (for non-zero weights) 
    inline double                   max        () const { return m_values. max () ; }       
    // ======================================================================
  public: // helper sums 
    // ======================================================================
    /// sum_i weight_i*value_i
    double                   sum   () const ; // sum_i weight_i * value_i
    /// sum_i weight_i*value_i**2
    double                   sum2  () const ; // sum_i weight_i * value_i**2
    /// sum_i weight_i  
    inline double            sumw  () const {  return m_weights.sum  () ; }
    /// sum_i weight_i*wight_i   
    inline double            sumw2 () const {  return m_weights.sum2 () ; }
    // ======================================================================
  public:  // statistics &  weights and values 
    // ======================================================================
    /// get the statistic of weights 
    inline const Ostap::StatEntity& weights () const { return m_weights ; }
    /// get the statistic of values with non-zero weight 
    inline const Ostap::StatEntity& values  () const { return m_values  ; }
    // ======================================================================
  public:
    // ======================================================================
    /// add   the  value      (with weight=1) to the counter 
    inline WStatEntity& operator+= ( const double value ) { return add (  value ) ; }
    /// add   the  value +1   (with weight=1) to the counter 
    inline WStatEntity& operator++ ()      { return   (*this) += 1 ; }
    /// add   the  value +1   (with weight=1) to the counter 
    inline WStatEntity& operator++ ( int ) { return ++(*this)      ; }    
    /// add   the -value      (with weight=1) to the counter 
    inline WStatEntity& operator-= ( const double value ) { return add ( -value ) ; }
    /// add   the  value -1   (with weight=1) to the counter 
    inline WStatEntity& operator-- ()      { return   (*this) -= 1 ; }
    /// add   the  value -1   (with weight=1) to the counter 
    inline WStatEntity& operator-- ( int ) { return --(*this)      ; }    
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
    // ===============================================================================
  public:
    // ===============================================================================
    /// Ostap::Math::WStatistics
    void update 
      ( const double value      ,  
	const double weight = 1 ) override { this -> add ( value , weight ) ; }
    // ===============================================================================
  public:
    // =============================================================================== 
    /** add another counter 
     *  @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    WStatEntity& add ( const WStatEntity& stat ) ;
    /// add another counter 
    WStatEntity& add ( const  StatEntity& stat ) 
    { return add ( WStatEntity ( stat ) ) ; }
    // ======================================================================
  public:
    // ======================================================================
    /** add another counter 
     *  @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    inline WStatEntity& operator+= ( const WStatEntity& other ) { return add ( other ) ; }
    /// add another counter 
    inline WStatEntity& operator+= ( const  StatEntity& other ) { return add ( other ) ; }
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
