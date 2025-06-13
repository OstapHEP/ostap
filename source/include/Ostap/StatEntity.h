// ============================================================================
#ifndef OSTAP_STATENTITY_H
#define OSTAP_STATENTITY_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <limits>
#include <string>
#include <ostream>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Statistic.h"
// ============================================================================
namespace Ostap
{
  // ========================================================================
  /** @class StatEntity Ostap/StatEntity.h
   *  A slghtly modified version of StatEntity from Gaudi project.
   *  Essentially the generic counter could be considered as
   *  the trivial 1-bin histogram 
   *  @code
   *  // get all tracks
   *  const Tracks* tracks = ... ;
   *  // create the counter
   *  StatEntity chi2 ;
   *  // make a loop over all tracks:
   *  for ( auto& track : *tracks ) 
   *  {
   *     const Track* track = *itrack ;
   *     // use the counters to accumulate information:
   *    chi2 += track -> chi2  () ;
   *  } ;
   *  // Extract the information from the counter:
   *  // get number of entries (== number of tracks)
   *  int nEntries = chi2.n       () ;
   *  // get the minimal value of chi2
   *  double chi2Min = chi2.min   () ;
   *  // get the average value of chi2:
   *  double meanChi2 = chi2.mean () ;
   *  // get R.M.S. for chi2-distribution:
   *  double rmsChi2   = chi2.rms () ;
   *  .. etc...
   *  @endcode
   *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date    26/11/1999
   *  @date    2005-08-02
   */
  class StatEntity : public Ostap::Math::Statistic 
  {
  public:
    // ======================================================================
    /// the default constructor
    StatEntity  () = default ;
    // ======================================================================
    /* The full constructor from all important values:
     *  @param entries number of entries
     *  @param mu   the mean value 
     *  @param mu2  the second central moment/variance/dispersion 
     *  @param minv the minimum value 
     *  @param maxv the maximum value 
     */
    StatEntity 
    ( const unsigned long entries ,
      const double        mu      ,
      const double        mu2     ,
      const double        minv    ,
      const double        maxv    ) ;
    // ======================================================================
  public: // the basic accessors 
    // ======================================================================
    /// number of entries 
    unsigned long long   n          () const { return m_n    ; }
    /// effective number of entries 
    unsigned long long   nEff       () const { return m_n    ; }
    /// mean value 
    double               mu         () const { return m_mu   ; }
    /// the second central moment/dispersion/variance  
    double               mu2        () const { return m_mu2  ; }    
    /// minimal value
    double               min        () const { return m_min  ; }
    /// maximal valu e
    double               max        () const { return m_max  ; }
    /// get number of "good" (mnon-zero) entries
    unsigned long long   nGood      () const { return n   () ; }
    // ======================================================================
  public: // derived quantities & aliases  
    // ======================================================================
    /// empty ?
    bool                 empty      () const { return 0 == m_n ; }
    /// number of entries 
    unsigned long long   nEntries   () const { return m_n    ; }
    /// variance
    double               variance   () const { return m_mu2  ; }
    /// dispersion 
    double               dispersion () const { return m_mu2  ; }
    /// r.m.s of value
    double               rms        () const ;
    /// mean value of counter
    double               mean       () const { return m_mu   ; }
    /// error in mean value of counter
    double               meanErr    () const ;
    // ======================================================================
  public: // various helper sums 
    // ======================================================================
    /// sum of all  values 
    double               sum        () const ; 
    /// accumulated value**2
    double               sum2       () const ;
    // ======================================================================
  public:
    // ======================================================================
    /** Interpret the counter as some measure of efficiency
     *  The efficiency is calculated as the ratio of the weight
     *  over the number of entries
     *  One gets the correct interpretation in the case of
     *  filling the counters only with 0 and 1.
     *  Some checks are performed:
     *   - number of counts must be positive
     *   - "value" must be non-negative
     *   - "value" does not exceed the overall number of counts
     *
     *  If these conditions are not satisfied the method returns -1,
     *  otherwise it returns the value of "mean"
     *
     *  @code
     *  StatEntity& stat = ... ;
     *  for ( auto& track : tracks )  
     *  {
     *    const bool good = PT( track ) > 1 * GeV ;
     *    stat += good ;
     *  }
     *  std::cout << " Efficiency of PT-cut is "
     *            << stat.efficiency() << std::endl ;
     *  @endcode
     *  @see StatEntity::efficiency 
     */
    double efficiency () const ;
    // ======================================================================
    /** Interpret the counter as some measure of efficiency and evaluate the
     *  uncertainty of this efficiency using <b>binomial</b> estimate.
     *  The efficiency is calculated as the ratio of the weight
     *  over the number of entries
     *  One gets the correct interpretation in the case of
     *  filling the counters only with 0 and 1.
     *  Some checks are performed:
     *   - number of counts must be positive
     *   - "value" must be non-negative
     *   - "value" does not exceed the overall number of counts
     *
     *  If these conditions are not satisfied the method returns -1.
     *
     *  @attention The action of this method is <b>DIFFERENT</b> from the action of
     *             the method StatEntity::meanErr
     *
     *  @code
     *  StatEntity& stat = ... ;
     *  for ( auto& track : tracks )
     *  {
     *    const bool good = PT( track ) > 1 * GeV ;
     *    stat += good ;
     *  }
     *  std::cout << " Efficiency of PT-cut is "
     *            << stat.efficiency    () << "+-"
     *            << stat.efficiencyErr () << std::endl ;
     *  @endcode
     *  @see StatEntity::efficiency
     *  @see StatEntity::efficiencyErr
     */
    double efficiencyErr () const ;
    /// shortcut, @see StatEntity::efficiency
    double eff           () const { return efficiency    () ; }
    /// shortcut, @see StatEntity::efficiencyErr
    double effErr        () const { return efficiencyErr () ; }
    // ======================================================================
  public: // operators
    // ======================================================================
    /** General increment operator for the flag
     *  Could be used for easy manipulation with StatEntity object:
     *  @code
     *  StatEntity stat ;
     *  for ( int i =    ... )
     *   {
     *    double eTotal = .... ;
     *    // increment the counter
     *    stat += eTotal ;
     *  }
     *  @endcode
     *  @param f counter increment
     */
    StatEntity& operator+= ( const double f ) { return add ( f )  ; }
    // ======================================================================
    /** Pre-increment operator for the flag
     *  Could be used for easy manipulation with StatEntity object:
     *  @code
     *  StatEntity stat ;
     *  for ( int i =    ... )
     *   {
     *     double eTotal = .... ;
     *     // increment the counter
     *     if ( eTotal > 1 * TeV ) { ++stat ; } ;
     *  }
     *  @endcode
     */
    StatEntity& operator++ ()    { return   (*this)+= 1  ; }
    // ======================================================================
    /** Post-increment operator for the flag.
     *  Actually it is the same as pre-increment.
     *  Could be used for easy manipulation with StatEntity object:
     *  @code
     *  StatEntity stat ;
     *  for ( int i =    ... )
     *   {
     *    double eTotal = .... ;
     *    // increment the counter
     *    if ( eTotal > 1 * TeV ) { stat++ ; } ;
     *  }
     *  @endcode
     */
    StatEntity& operator++ (int) { return ++(*this) ; }
    // ======================================================================
    /** General decrement operator for the flag
     *  Could be used for easy manipulation with StatEntity object:
     *  @code
     *  StatEntity stat ;
     *  for ( int i =    ... )
     *   {
     *    double eTotal = .... ;
     *    // decrement the counter
     *    stat -= eTotal ;
     *  }
     *  @endcode
     *  @param f counter increment
     */
    StatEntity& operator-= ( const double   f ) { return add ( -f ) ; }
    // ======================================================================
    /** Pre-decrement operator for the flag
     *  Could be used for easy manipulation with StatEntity object:
     *  @code
     *  StatEntity stat ;
     *  for ( int i =    ... )
     *  {
     *    double eTotal = .... ;
     *    // increment the counter
     *    if ( eTotal > 1 * TeV ) { --stat ; } ;
     *  }
     *  @endcode
     */
    StatEntity& operator-- () { return (*this)-=1  ; }
    // ======================================================================
    /** Post-decrement operator for the flag
     *  Could be used for easy manipulation with StatEntity object:
     *  @code
     *  StatEntity stat ;
     *  for ( int i = ... )
     *  {
     *    double eTotal = .... ;
     *    // increment the counter
     *    if ( eTotal > 1 * TeV ) { stat-- ; } ;
     *  }
     *  @endcode
     */
    StatEntity& operator-- (int) { return --(*this) ; }
    // ======================================================================
    /** increment with other counter
     *  @code
     *  const StatEntity second = ... ;
     *  StatEntity first = ... ;
     *  first += second ;
     *  @endcode
     *  @param other counter to be added
     *  @return self-reference
     *  @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
     *  @see https://doi.org/10.1007/s00180-015-0637-z
     */
    StatEntity& operator+= ( const StatEntity& other ) { return add  ( other ) ; }
    // ======================================================================
  public:
    // ======================================================================
    /// basic comparison
    bool operator< ( const StatEntity& s ) const ;
    bool operator==( const StatEntity& s ) const ;
    /// derived comparisons:
    bool operator> ( const StatEntity& s ) const { return s < *this ; }
    bool operator<=( const StatEntity& s ) const { return    (*this) == s || (*this) < s ; }
    bool operator>=( const StatEntity& s ) const { return    (*this) == s || (*this) > s ; }
    bool operator!=( const StatEntity& s ) const { return   !(*this  == s ) ; }
    // ======================================================================
  public:
    // ======================================================================
    /** add a value : the main method 
     *  @see https://arxiv.org/abs/1510.04923v1
     *  @param value (INPUT) value to be added
     *  @return self-reference 
     *  @attention non-finite values ar eignored! 
     */
    StatEntity& add    ( const double      value ) ;
    // ========================================================================
    /** increment with other counter
     *  @code
     *  const StatEntity second = ... ;
     *  StatEntity first = ... ;
     *  first += second ;
     *  @endcode
     *  @param other counter to be added
     *  @return self-reference
     *  @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
     *  @see https://doi.org/10.1007/s00180-015-0637-zu
     */
    StatEntity& add    ( const StatEntity& value ) ;
    // ========================================================================
  public:
    // ========================================================================
    /// update counter Ostap::Math::Statistic
    void update ( const double value ) override{ add ( value ) ; } ;
    // ========================================================================
  public: // various technical helper methods  
    // ========================================================================
    /// reset the counters
    void reset () ;
    /// swap two counters
    void swap ( StatEntity& right ) ;
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
	std::isfinite ( m_min ) &&
	std::isfinite ( m_max ) ;	
    }
    // ======================================================================
  private: // data members 
    // ======================================================================
    /// number of calls
    unsigned long long m_n   { 0 } ;
    /// accumulated flag
    double             m_mu  { 0 } ;
    double             m_mu2 { 0 } ;
    double             m_min {   std::numeric_limits<double>::max() } ;
    double             m_max { - std::numeric_limits<double>::max() } ;
    // ======================================================================
  };
  // ========================================================================
  /// external operator for addition of StatEntity and a number
  inline StatEntity    operator+ ( StatEntity   e , const double      v ) 
  { e += v ; return e ; }
  /// external operator for addition of StatEntity and a number
  inline StatEntity    operator+ ( const double v , const StatEntity& e )
  { return e + v ; }
  /// external operator for addition of StatEntity and StatEntity
  inline StatEntity    operator+ ( StatEntity   a , const StatEntity& b )
  {  a += b ; return a ; }
  /// external operator for subtraction of StatEntity and a number
  inline StatEntity    operator- ( StatEntity   e , const double      v ) 
  { e -= v ; return e ; }
  /// external printout operator to std::ostream
  inline std::ostream& operator<<( std::ostream& o , const StatEntity& e ) 
  { return e.fillStream ( o ) ; }
  // ==========================================================================
  /// conversion to string 
  inline std::string to_string ( const StatEntity& s ) { return s.toString() ; }
  // ==========================================================================
  /// swap two counters 
  inline void swap ( StatEntity&  a , StatEntity&  b  ) { a.swap( b ) ; }
  // =========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_STATENTITY_H
// ============================================================================
