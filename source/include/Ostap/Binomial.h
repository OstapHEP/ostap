// ============================================================================
#ifndef OSTAP_BINOMIAL_H 
#define OSTAP_BINOMIAL_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <utility>
// ============================================================================
/** @file Ostap/Binomial.h
 *  Collection of functions to estimate the confidence intervals for 
 *  binomial proportion/efficiency
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** normal approximation interval for binomial proportion/efficiency 
     *  ( "Wald test")
     *  @param  accepted  number of accepted events
     *  @param  rejected  number of rejected events
     *  @param  conflevel the confidence level:    0<=CL<=1 
     *  @return the confidence interval 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-09-17
     *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
     */
    std::pair<double,double>
    wald_interval
    ( const unsigned long accepted  ,
      const unsigned long rejected  ,
      const double        conflevel ) ;
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
    std::pair<double,double>
    wilson_score_interval
    ( const unsigned long accepted  ,
      const unsigned long rejected  ,
      const double        conflevel ) ;
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
    std::pair<double,double>
    wilson_score_continuity_interval
    ( const unsigned long accepted  ,
      const unsigned long rejected  ,
      const double        conflevel ) ;
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
    std::pair<double,double>
    arcsin_interval
    ( const unsigned long accepted  ,
      const unsigned long rejected  ,
      const double        conflevel ) ;
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
    std::pair<double,double>
    agresti_coull_interval
    ( const unsigned long accepted  ,
      const unsigned long rejected  ,
      const double        conflevel ) ;
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
    std::pair<double,double>
    jeffreys_interval
    ( const unsigned long accepted  ,
      const unsigned long rejected  ,
      const double        conflevel ) ;
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
    std::pair<double,double>
    clopper_pearson_interval
    ( const unsigned long  accepted  ,
      const unsigned long  rejected  ,
      const double         conflevel ) ;
    // ========================================================================
    /** Bayes' theorem based interval
     *  @see M.Paterno, "Calculationg efficiencies and their uncertainties", 
     *                   FERMILAB-TM-2286-CD
     *  @see https://inspirehep.net/literature/669498
     *  @see DOI: 10.2172/15017262
     */
    std::pair<double,double>
    bayes_interval
    ( const unsigned long accepted  ,
      const unsigned long rejected  ,
      const double        conflevel ) ;
    // ========================================================================
  } //                                            end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_BINOMIAL_H
// ============================================================================
