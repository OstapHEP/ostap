// ============================================================================
#ifndef OSTAP_STATISTIC_H 
#define OSTAP_STATISTIC_H 1
// =============================================================================
namespace  Ostap
{
  // ===========================================================================
  /// The actual  type  ofevent undex(ROOT::Int_64)
  typedef unsigned long EventIndex; 
  // ===========================================================================
  /** @var FirstEvent  
   *  Index for the first event in the loop
   */
  extern const EventIndex FirstEvent ;
   // ===========================================================================
  /** @var LastEvent  
   *  Index for the last (exclusive) event in the loop
   */
  extern const EventIndex LastEvent ;
  // ===========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Statistic
     *  Helper abstract base class for statistic counters
     */
    class Statistic
    {
    public :
      // ======================================================================
      virtual ~Statistic() ;
      /// add new value to the counter
      virtual void update ( const double x ) = 0 ;
      // ======================================================================
    } ; //                          The end of the class Ostap::Math::Statistic
    // ========================================================================
    /** @class WStatistic
     *  Helper (empty) base class for weighted statistics
     */
    class WStatistic
    {
    public :
      // ======================================================================
      virtual ~WStatistic () ;
      /// add new value to the counter
      virtual void update ( const double x , const double w = 1 ) = 0 ;
      // ======================================================================
    } ; //                         The end of the class Ostap::Math::WStatistic 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
#endif // OSTAP_STATISTIC_H
// ============================================================================
