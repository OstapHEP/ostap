// ============================================================================
#ifndef OSTAP_STATENTITY_H
#define OSTAP_STATENTITY_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <string>
#include <iostream>
// ============================================================================
namespace Ostap
{
  // ========================================================================
  /** @class StatEntity StatEntity.h
   *
   *  A slghtly modified version of StatEntity from Gaudi project
   *
   *  The basic counter used for Monitoring purposes.
   *
   *  Essentially the generic counter could be considered as
   *  the trivial 1-bin
   *
   *  @code
   *
   *   // get all tracks
   *   const Tracks* tracks = ... ;
   *   // create the counter
   *   StatEntity chi2 ;
   *   // make a loop over all tracks:
   *   for ( Tracks::const_iterator itrack = tracks->begin() ;
   *         tracks->end() != itrack ; ++itrack )
   *    {
   *        const Track* track = *itrack ;
   *
   *        // use the counters to accumulate information:
   *        chi2 += track -> chi2  () ;
   *    } ;
   *
   *   // Extract the information from the counter:
   *
   *   // get number of entries (== number of tracks)
   *   int nEntries = chi2.nEntries() ;
   *
   *   // get the minimal value of chi2
   *   double chi2Min = chi2.flagMin() ;
   *
   *   // get the average value of chi2:
   *   double meanChi2 = chi2.flagMean() ;
   *
   *   // get R.M.S. for chi2-distribution:
   *   double rmsChi2   = chi2.flagRMS() ;
   *
   *   .. etc...
   *
   *  @endcode
   *
   *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date    26/11/1999
   *  @date    2005-08-02
   */
  class StatEntity
  {
  public:
    // ======================================================================
    /// the default constructor
    StatEntity  () { reset() ; }
    // ======================================================================
    /* The full constructor from all important values:
     * @attention it need to be coherent with
     *            the actual structure of the class
     *            and the format description
     *  @see StatEntity::format
     *  @param entries number of entries
     *  @param flag  accumulated flag
     *  @param flag2 accumulated statistics: flag**2
     *  @param minFlag the minimum value for flag
     *  @param maxFlag the maximum value for flag
     */
    StatEntity ( const unsigned long entries ,
                 const double        flag    ,
                 const double        flag2   ,
                 const double        minFlag ,
                 const double        maxFlag ) ;
    /// copy 
    StatEntity  ( const StatEntity& ) = default ;    
    /// destructor
    ~StatEntity () = default;
    // ======================================================================
  public: // the basic quantities
    // ======================================================================
    /// getters
    unsigned long long   nEntries () const { return m_se_nEntries    ; }
    /// accumulated value
    double               sum      () const { return m_se_accumulatedFlag  ; }
    /// accumulated value**2
    double               sum2     () const { return m_se_accumulatedFlag2 ; }
    /// mean value of counter
    double               mean     () const ;
    /// r.m.s of value
    double               rms      () const ;
    /// error in mean value of counter
    double               meanErr  () const ;
    /// minimal value
    double               min      () const { return m_se_minimalFlag ; }
      /// maximal value
    double               max      () const { return m_se_maximalFlag ; }
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
     *   - "flag" must be non-negative
     *   - "flag" does not exceed the overall number of counts
     *
     *  If these conditions are not satisfied the method returns -1,
     *  otherwise it returns the value of "flagMean"
     *
     *  @code
     *
     *  StatEntity& stat = ... ;
     *
     *  for ( TRACKS::iterator itrack = tracks.begin() ;
     *        tracks.end() != itrack ; itrack )
     *   {
     *     const bool good = PT( *itrack ) > 1 * Gaudi::Units::GeV ;
     *     stat += good ;
     *   }
     *
     *  std::cout << " Efficiency of PT-cut is "
     *            << stat.efficiency() << std::endl ;
     *
     *  @endcode
     *  @see StatEntity::flagMean
     *  @see StatEntity::eff
     *  @see StatEntity::efficiencyErr
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
     *   - "flag" must be non-negative
     *   - "flag" does not exceed the overall number of counts
     *
     *  If these conditions are not satisfied the method returns -1.
     *
     *  @attention The action of this method is <b>DIFFERENT</b> from the action of
     *             the method StatEntity::flagMeanErr
     *
     *  @code
     *
     *  StatEntity& stat = ... ;
     *
     *  for ( TRACKS::iterator itrack = tracks.begin() ;
     *        tracks.end() != itrack ; itrack )
     *   {
     *     const bool good = PT( *itrack ) > 1 * Gaudi::Units::GeV ;
     *     stat += good ;
     *   }
     *
     *  std::cout << " Efficiency of PT-cut is "
     *            << stat.efficiency    () << "+-"
     *            << stat.efficiencyErr () << std::endl ;
     *
     *  @endcode
     *  @see StatEntity::efficiency
     *  @see StatEntity::effErr
     *  @see StatEntity::flagMeanErr
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
     *
     *  @code
     *
     *   StatEntity stat ;
     *
     *   for ( int i =    ... )
     *    {
     *       double eTotal = .... ;
     *
     *       // increment the counter
     *       stat += eTotal ;
     *    }
     *
     *  @endcode
     *
     *  @param f counter increment
     */
    StatEntity& operator+= ( const double f ) { add ( f ) ; return *this ; }
    // ======================================================================
    /**  Pre-increment operator for the flag
     *   Could be used for easy manipulation with StatEntity object:
     *
     *  @code
     *
     *   StatEntity stat ;
     *
     *   for ( int i =    ... )
     *    {
     *       double eTotal = .... ;
     *
     *       // increment the counter
     *       if ( eTotal > 1 * TeV ) { ++stat ; } ;
     *    }
     *
     *  @endcode
     */
    StatEntity& operator++ ()    { return   (*this)+= 1  ; }
    // ======================================================================
    /** Post-increment operator for the flag.
     *  Actually it is the same as pre-increment.
     *  Could be used for easy manipulation with StatEntity object:
     *
     *  @code
     *
     *   StatEntity stat ;
     *
     *   for ( int i =    ... )
     *    {
     *       double eTotal = .... ;
     *
     *       // increment the counter
     *       if ( eTotal > 1 * TeV ) { stat++ ; } ;
     *    }
     *
     *  @endcode
     */
    StatEntity& operator++ (int) { return ++(*this) ; }
    // ======================================================================
    /** General decrement operator for the flag
     *  Could be used for easy manipulation with StatEntity object:
     *
     *  @code
     *
     *   StatEntity stat ;
     *
     *   for ( int i =    ... )
     *    {
     *       double eTotal = .... ;
     *
     *       // decrement the counter
     *       stat -= eTotal ;
     *    }
     *
     *  @endcode
     *
     *  @param f counter increment
     */
    StatEntity& operator-= ( const double   f ) { add ( -f ) ; return *this ; }
    // ======================================================================
    /** Pre-decrement operator for the flag
     *   Could be used for easy manipulation with StatEntity object:
     *
     *  @code
     *
     *   StatEntity stat ;
     *
     *   for ( int i =    ... )
     *    {
     *       double eTotal = .... ;
     *
     *       // increment the counter
     *       if ( eTotal > 1 * TeV ) { --stat ; } ;
     *    }
     *
     *  @endcode
     */
    StatEntity& operator-- () { return (*this)-=1  ; }
    // ======================================================================
    /** Post-decrement operator for the flag
     *  Could be used for easy manipulation with StatEntity object:
     *
     *  @code
     *
     *   StatEntity stat ;
     *
     *   for ( int i =    ... )
     *    {
     *       double eTotal = .... ;
     *
     *       // increment the counter
     *       if ( eTotal > 1 * TeV ) { stat-- ; } ;
     *    }
     *
     *  @endcode
     */
    StatEntity& operator-- (int) { return --(*this) ; }
    // ======================================================================
    /** increment with other counter (useful for Monitoring/Adder )
     *
     *  @code
     *
     *   const StatEntity second = ... ;
     *
     *   StatEntity first = ... ;
     *
     *   first += second ;
     *
     *  @endcode
     *
     *  @param other counter to be added
     *  @return self-reference
     */
    StatEntity& operator+= ( const StatEntity& other ) ;
    // ======================================================================
  public:
    // ======================================================================
    /// basic comparison
    bool operator< ( const StatEntity& s ) const ;
    bool operator==( const StatEntity& s ) const ;
    /// derived comparisons:
    bool operator> ( const StatEntity& s ) const { return s < *this ; }
    bool operator<=( const StatEntity& s ) const { return   (*this) == s || (*this) < s ; }
    bool operator>=( const StatEntity& s ) const { return   (*this) == s || (*this) > s ; }
    bool operator!=( const StatEntity& s ) const { return !( *this  == s ) ; }
    // ======================================================================
  public:
    // ======================================================================
    /** add a value
     *  @param Flag value to be added
     *  @return number of entries
     */
    unsigned long long add ( const double Flag ) ;
    /// reset the counters
    void reset () ;
    /// representation as string
    std::string   toString () const;
    /// printout  to std::ostream
    std::ostream& fillStream ( std::ostream& o ) const ;
    // ======================================================================
  private:
    // ======================================================================
    /// number of calls
    unsigned long long           m_se_nEntries          ;
    /// accumulated flag
    double                       m_se_accumulatedFlag   ;
    double                       m_se_accumulatedFlag2  ;
    double                       m_se_minimalFlag       ;
    double                       m_se_maximalFlag       ;
    // ======================================================================
  };
  // ========================================================================
  /// external operator for addition of StatEntity and a number
  inline StatEntity    operator+ ( StatEntity e , const double v ) { e += v ; return e ; }
  /// external operator for addition of StatEntity and a number
  inline StatEntity    operator+ ( const double v , const StatEntity& e ) { return e + v ; }
  /// external operator for addition of StatEntity and StatEntity
  inline StatEntity    operator+ ( StatEntity   a , const StatEntity& b ) {  a += b ; return a ; }
  /// external operator for subtraction of StatEntity and a number
  inline StatEntity    operator- ( StatEntity e , const double v ) { e -= v ; return e ; }
  /// external printout operator to std::ostream
  inline std::ostream& operator<<( std::ostream& o , const StatEntity& e ) 
  { return e.fillStream ( o ) ; }
  // ==========================================================================
  /// conversion to string 
  inline std::string to_string ( const StatEntity& s ) { return s.toString() ; }
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
// The END
// ============================================================================
#endif // OSTAP_STATENTITY_H
// ============================================================================
