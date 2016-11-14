// $Id$ 
// ============================================================================
#ifndef OSTAP_NSTATENTITY_H 
#define OSTAP_NSTATENTITY_H 1
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
  /** @class NStatEntity LHCbMath/NStatEntity.h
   *  A good and efficient approximation to the running
   *  statistic for last N-events.Actually statiustic is calcualetd for 
   *  last n-events, where   N<= n < 2*N.
   *  Implemeneted as two sliding counters, each of them is reset 
   *  for every 2*N events.
   *  Something similar I've written long time ago for Velo.
   *  It can be useful for e.g. "current rate" counters 
   *  @see Ostap::StatEntity 
   *  Coplied from LHCbMath 
   *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
   *  @date   2015-04-03
   */
  class NStatEntity 
  {
    // ======================================================================
  public:
    // ======================================================================
    /// constructor with N-parameter  
    NStatEntity ( const unsigned long N = 1000 ) ;
    // ======================================================================
  public:
    // ======================================================================
    /// get N, such that statistic is calculated for n-events  N<= n < 2*N 
    unsigned long N () const { return m_N ; }
    // ======================================================================
  public: // conversion to the counter with longer history: N<= n < 2*N 
    // ======================================================================
    /// get the counter with longer history 
    inline const StatEntity& counter () const 
    { return m_cnt1.nEntries() > m_cnt2.nEntries() ? m_cnt1 : m_cnt2 ; }
    /// implicit cast to counters 
    inline operator const StatEntity&() const { return counter () ; }
    // ======================================================================
  public: // mimic regular counter interface, the basic quantities
    // ======================================================================
    /// numbner of entries (  N <= n < 2*N ) 
    unsigned long nEntries      () const { return counter().nEntries     () ; }
    /// accumulated value
    double        sum           () const { return counter().sum          () ; }
    /// accumulated value**2
    double        sum2          () const { return counter().sum2         () ; }
    /// mean value of counter
    double        mean          () const { return counter().mean         () ; }
    /// r.m.s of value
    double        rms           () const { return counter().rms          () ; }
    /// error in mean value of counter
    double        meanErr       () const { return counter().meanErr      () ; }
    /// minimal value
    double        min           () const { return counter().min          () ; }
    /// maximal value
    double        max           () const { return counter().max          () ; }
    // efficiency 
    double        efficiency    () const { return counter().efficiency   () ; }
    // efficiency error
    double        efficiencyErr () const { return counter().efficiencyErr() ; }
    // ======================================================================
  public: // increment and decrement 
    // ======================================================================
    /// pre-increment
    NStatEntity& operator++ ()    { return add (  1 ) ; }
    /// post-increment (actually udentical to pre-increment 
    NStatEntity& operator++ (int) { return add (  1 ) ; }
    /// pre-decrement
    NStatEntity& operator-- ()    { return add ( -1 ) ; }
    /// post-decrement (actually udentical to pre-decrement 
    NStatEntity& operator-- (int) { return add ( -1 ) ; }
    // ======================================================================
  public: // another form of increment/decrement 
    // ======================================================================
    /// another form of increment 
    NStatEntity& operator+= ( const double   f ) { return add (  f ) ; }
    /// another form of decrement 
    NStatEntity& operator-= ( const double   f ) { return add ( -f ) ; }
    // ======================================================================
  public:
    // ======================================================================
    /** printout  to std::ostream
     *  @param s the reference to the output stream
     */
    std::ostream& fillStream ( std::ostream& o ) const 
    { return counter().fillStream ( o ) ; }
    // =====================================================================
    /// conversion to string
    std::string toString() const { return counter().toString () ; }
    // =====================================================================
    /// reset method (likely not needed at all) 
    void reset() 
    {
      m_cnt1.reset() ;
      m_cnt2.reset() ;        
    }
    // =====================================================================
  public: // the main method without decorations  
    // ======================================================================
    /// the main method without decorations  
    NStatEntity& add ( const double f ) 
    {
      /// inncrement both counters 
      m_cnt1.add ( f ) ;
      m_cnt2.add ( f ) ;
      // reset them when needed 
      if ( 0    == m_cnt1.nEntries() % ( 2 * m_N ) ) { m_cnt1.reset() ; }
      if ( m_N  == m_cnt1.nEntries() % ( 2 * m_N ) ) { m_cnt2.reset() ; }
      //
      return *this ;
    }
    // ======================================================================
  private:
    // ======================================================================
    /// the first  counter 
    StatEntity    m_cnt1 ;                            // the first  counter 
    /// the second counter 
    StatEntity    m_cnt2 ;                            // the second counter 
    /// the sliding window 
    unsigned long m_N    ;                            // the sliding window 
    // ======================================================================
  };
  // ========================================================================
  /// printout operator to std::ostream
  inline std::ostream& operator<<
  ( std::ostream&      stream ,
    const NStatEntity& entity ) { return entity.fillStream ( stream ) ; }
  // ========================================================================
  /// operator for addition of NStatEntity and a number
  inline NStatEntity operator+( NStatEntity        e , 
                                const double       v ) { return e.add ( v ) ; }
  /// operator for addition of NStatEntity and a number
  inline NStatEntity operator+( const double       v , 
                                const NStatEntity& e ) { return e + v       ; }
  /// operator for addition of NStatEntity and a number
  inline NStatEntity operator-( const NStatEntity& e , 
                                const double       v ) { return e + (-1*v)  ; }
  // ==========================================================================
  /// conversion to string 
  inline std::string to_string ( const NStatEntity& e ) { return e.toString() ;}
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_NSTATENTITY_H
// ============================================================================
