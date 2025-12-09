// ============================================================================
#ifndef OSTAP_EFFENTITY_H
#define OSTAP_EFFENTITY_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <utility>
#include <iosfwd>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
// ============================================================================
namespace Ostap
{
  // ========================================================================
  /** @class EffEntity Ostap/EffEntity.h
   *  Trivial efficency counter 
   */
  class EffEntity 
  {
  public:
    // ======================================================================
    /// the  content type: Number of entries  
    typedef unsigned long long       size_type ;
    /// interval type
    typedef std::pair<double,double> Interval  ; 
    // ======================================================================
  public:
    // ======================================================================
    /** constructor
     *  @param accepted number of accepted events 
     *  @param rejected number of rejected events 
     */
    EffEntity
      ( const size_type accepted = 0 ,
	const size_type rejected = 0 ) ;
    // ======================================================================
  public: // the basic accessors 
    // ======================================================================
    /// number of accepted events 
    inline size_type accepted () const { return m_accepted ; }
    /// number of rejected events 
    inline size_type rejected () const { return m_rejected ; }
    /// total numebr of events: accepted + rejected 
    inline size_type total    () const { return m_accepted + m_rejected ; }    
    // ======================================================================
  public: // the main method 
    // ======================================================================
    /** convert counter to binomial efficiency
     *  - evaluate the binomial efficiency for Bernulli scheme    
     *  @see Ostap::Math::ValueWithError 
     *  @see Ostap::Math::binomEff 
     *  @see Ostap::EffEntity::binomEff 
     */
    inline Ostap::Math::ValueWithError efficiency () const
    { return binomEff () ; }
    // ======================================================================
    /** convert counter to binomial efficiency
     *  - evaluate the binomial efficiency for Bernulli scheme    
     *  @see Ostap::Math::ValueWithError 
     *  @see Ostap::Math::binomEff 
     */
    Ostap::Math::ValueWithError binomEff         () const ;
    // ======================================================================
    /** evaluate the binomial efficiency interval using Wilson's prescription
     *  @see Ostap::Math::ValueWithError 
     *  @see Ostap::Math::wilsonEff 
     *  @return the binomial efficiency 
     */
    Ostap::Math::ValueWithError wilsonEff        () const ;
    // ======================================================================
    /** evaluate the binomial efficiency interval 
     *  using Agresti-Coull's prescription
     *  @return the binomial efficiency 
     */
    Ostap::Math::ValueWithError agrestiCoullEff  () const ;
    // ======================================================================    
  public: // get the confidence intervals 
    // ======================================================================
    /** normal approximation interval for binomial proportion/efficiency 
     *  ( "Wald test")
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    Interval  wald_interval (  const double        conflevel ) const ;
    // ========================================================================
    /** Wilson score interval for binomial proportion/efficiency 
     *  @param  accepted  number of accepted events
     *  @param  rejected  number of rejected events
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    Interval wilson_score_interval (  const double        conflevel ) const ;
    // ========================================================================
    /** Wilson score interval with continuity correction for binomial proportion/efficiency 
     *  @param  accepted  number of accepted events
     *  @param  rejected  number of rejected events
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    Interval wilson_score_continuity_interval ( const double        conflevel ) const ;
    // ========================================================================
    /** ArcSin interval with continuity correction for binomial proportion/efficiency 
     *  @param  accepted  number of accepted events
     *  @param  rejected  number of rejected events
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    Interval arcsin_interval ( const double        conflevel ) const ;
    // ========================================================================
    /** Agresti-Coull interval with continuity correction for binomial proportion/efficiency 
     *  @param  accepted  number of accepted events
     *  @param  rejected  number of rejected events
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    Interval agresti_coull_interval ( const double        conflevel ) const ;
    // ========================================================================
    /** Jeffreys interval for binomial proportion/efficiency 
     *  @param  accepted  number of accepted events
     *  @param  rejected  number of rejected events
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    Interval jeffreys_interval ( const double        conflevel ) const ;
    // ========================================================================
    /** Clopper-Pearson interval for binomial proportion/efficiency 
     *  @param  accepted  number of accepted events
     *  @param  rejected  number of rejected events
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    Interval clopper_pearson_interval ( const double         conflevel ) const ;
    // ========================================================================
    /** Bayes' theorem based interval
     *  @see M.Paterno, "Calculationg efficiencies and their uncertainties", 
     *                   FERMILAB-TM-2286-CD
     *  @see https://inspirehep.net/literature/669498
     *  @see DOI: 10.2172/15017262
     */
    Interval bayes_interval ( const double conflevel ) const ;
    // ========================================================================
  public:
    // ======================================================================
    /// update the counter 
    EffEntity& add ( const bool value )
    {
      if ( value ) { ++m_accepted ; }
      else         { ++m_rejected ; }
      return *this ;
    }
    /// update the counter 
    EffEntity& add ( const EffEntity& right  )
    {
      m_accepted += right.m_accepted ;
      m_rejected += right.m_rejected ;
      return *this ;
    }
    // ======================================================================
    /// update the counter 
    EffEntity& operator+= ( const bool        value ) { return add ( value ) ; }
    /// update the counter 
    EffEntity& operator+= ( const EffEntity&  value ) { return add ( value ) ; }
    // ======================================================================    
  private :
    // ======================================================================
    /// number of accepted events 
    size_type m_accepted { 0 } ; // number of accepted events 
    /// number of rejected events 
    size_type m_rejected { 0 } ; // number of rejected events 
    // ======================================================================
  }; //                                     The end of class Ostap::EffEntity
  // ========================================================================
  /// add two counters 
  inline EffEntity operator+( const EffEntity& a , const EffEntity& b )
  { EffEntity c{ a } ; c+= b ; return c ; }
  // ========================================================================
  /// add two counters 
  inline EffEntity operator+( const EffEntity& a , const bool       b )
  { EffEntity c { a } ; c+= b ; return c ; }
  // ========================================================================
  /// add two counters 
  inline EffEntity operator+( const bool       b , const EffEntity& a ) 
  { EffEntity c { a } ; c+= b ; return c ; }
  // ========================================================================
} //                                               The end of namespace Ostap 
// ==========================================================================
//                                                                    The END 
// ==========================================================================
#endif // OSTAP_EFFENTITY_H
// ==========================================================================
