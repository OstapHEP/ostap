// ============================================================================
#ifndef OSTAP_STATISTIC_H 
#define OSTAP_STATISTIC_H 1
// =============================================================================
namespace  Ostap
{
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
      /// reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================
    } ; //                          The end of the class Ostap::Math::Statistic
    // ========================================================================
   /** @class Statistic2
     *  Helper abstract base class for statistic counters
     */
    class Statistic2
    {
    public :
      // ======================================================================
      virtual ~Statistic2() ;
      /// add new value to the counter
      virtual void update 
      ( const double x , 
        const double y ) = 0 ; 
      // reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================
    } ; //                         The end of the class Ostap::Math::Statistic2
    // ========================================================================
   /** @class Statistic3
     *  Helper abstract base class for statistic counters
     */
    class Statistic3
    {
    public :
      // ======================================================================
      virtual ~Statistic3() ;
      /// add new value to the counter
      virtual void update 
      ( const double x , 
        const double y , 
        const double z ) = 0 ;
      /// reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================
    } ; //                         The end of the class Ostap::Math::Statistic3
    // ========================================================================
   /** @class Statistic4
     *  Helper abstract base class for statistic counters
     */
    class Statistic4
    {
    public :
      // ======================================================================
      virtual ~Statistic4() ;
      /// add new value to the counter
      virtual void update 
      ( const double x , 
        const double y , 
        const double z , 
        const double t ) = 0 ;
      /// reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================
    } ; //                         The end of the class Ostap::Math::Statistic4
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
      /// reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================
    } ; //                         The end of the class Ostap::Math::WStatistic 
    // ========================================================================
   /** @class WStatistic2
     *  Helper (empty) base class for weighted statistics
     */
    class WStatistic2
    {
    public :
      // ======================================================================
      virtual ~WStatistic2 () ;
      /// add new value to the counter
      virtual void update 
      ( const double x , 
        const double y , 
        const double w = 1 ) = 0 ;
      /// reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================
    } ; //                        The end of the class Ostap::Math::WStatistic2
    // ========================================================================
   /** @class WStatistic3
     *  Helper (empty) base class for weighted statistics
     */
    class WStatistic3
    {
    public :
      // ======================================================================
      virtual ~WStatistic3 () ;
      /// add new value to the counter
      virtual void update 
      ( const double x , 
        const double y , 
        const double z , 
        const double w = 1 ) = 0 ;
      /// reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================
    } ; //                        The end of the class Ostap::Math::WStatistic3
    // ========================================================================
    /** @class WStatistic4
     *  Helper (empty) base class for weighted statistics
     */
    class WStatistic4
    {
    public :
      // ======================================================================
      virtual ~WStatistic4 () ;
      /// add new value to the counter
      virtual void update 
      ( const double x , 
        const double y , 
        const double z , 
        const double t , 
        const double w = 1 ) = 0 ;
      /// reset  the context
      virtual void reset () = 0 ; 
      // ======================================================================      
    } ; //                        The end of the class Ostap::Math::WStatistic4
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
#endif // OSTAP_STATISTIC_H
// ============================================================================
