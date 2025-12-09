// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
#include "Ostap/EffEntity.h"
#include "Ostap/Binomial.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local 
// ============================================================================
#include "format.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::EffEntity
 *  @date 26/06/2001
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 */
// ============================================================================
/** constructor
 *  @param accepted number of accepted events 
 *  @param rejected number of rejected events 
 */
// ============================================================================
Ostap::EffEntity::EffEntity
( const size_type accepted ,
  const size_type rejected )
  : m_accepted ( accepted )
  , m_rejected ( rejected )
{}
// ============================================================================
/* convert counter to binomial efficiency
 *  - evaluate the binomial efficiency for Bernulli scheme    
 *  @see Ostap::Math::ValueWithError 
 *  @see Ostap::Math::binomEff 
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::EffEntity::binomEff () const
{ return Ostap::Math::binomEff ( accepted () , total () ) ; }
// ============================================================================
/*  evaluate the binomial efficiency interval using Wilson's prescription
 *  @return the binomial efficiency 
 *  @see Ostap::Math::ValueWithError 
 *  @see Ostap::Math::wilsonEff 
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::EffEntity::wilsonEff  () const
{ return Ostap::Math::wilsonEff ( accepted () , total () ) ; }
// ============================================================================
/* evaluate the binomial efficiency interval 
 *  using Agresti-Coull's prescription
 *  @return the binomial efficiency 
 *  @see Ostap::Math::agrestiCoullEff 
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::EffEntity::agrestiCoullEff  () const
{ return Ostap::Math::agrestiCoullEff ( accepted () , total () ) ; }
// ============================================================================
// Intervals 
// ============================================================================
/*  normal approximation interval for binomial proportion/efficiency 
 *  ( "Wald test")
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::wald_interval
( const double  conflevel ) const
{ return Ostap::Math::wald_interval ( accepted () ,
				      rejected () ,
				      conflevel   ) ; }
// ============================================================================
/* Wilson score interval for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::wilson_score_interval
( const double  conflevel ) const
{ return Ostap::Math::wilson_score_interval ( accepted () ,
					      rejected () ,
					      conflevel   ) ; }
// ============================================================================
/* Wilson score interval with continuity correction for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::wilson_score_continuity_interval
( const double  conflevel ) const
{ return Ostap::Math::wilson_score_continuity_interval ( accepted () ,
							 rejected () ,
							 conflevel   ) ; }
// ============================================================================
/*  ArcSin interval with continuity correction for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::arcsin_interval
( const double  conflevel ) const
{ return Ostap::Math::arcsin_interval ( accepted () ,
					rejected () ,
					conflevel   ) ; }
// ============================================================================
/*  Agresti-Coull interval with continuity correction for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::agresti_coull_interval
( const double  conflevel ) const
{ return Ostap::Math::agresti_coull_interval ( accepted () ,
					       rejected () ,
					       conflevel   ) ; }
// ============================================================================
/* Jeffreys interval for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::jeffreys_interval
( const double  conflevel ) const
{ return Ostap::Math::jeffreys_interval ( accepted () ,
					  rejected () ,
					  conflevel   ) ; }
// ============================================================================
/* Clopper-Pearson interval for binomial proportion/efficiency 
 *  @param  accepted  number of accepted events
 *  @param  rejected  number of rejected events
 *  @param  conflevel the confidence level:    0<=CL<=1 
 *  @return the confidence interval 
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 *  @see http://en.wikipedia.org.wiki/Binomial_proportion_connfidence_interval
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::clopper_pearson_interval
( const double  conflevel ) const
{ return Ostap::Math::clopper_pearson_interval ( accepted () ,
						 rejected () ,
						 conflevel   ) ; }
// ============================================================================
/*  Bayes' theorem based interval
 *  @see M.Paterno, "Calculationg efficiencies and their uncertainties", 
 *                   FERMILAB-TM-2286-CD
 *  @see https://inspirehep.net/literature/669498
 *  @see DOI: 10.2172/15017262
 */
// ============================================================================
Ostap::EffEntity::Interval
Ostap::EffEntity::bayes_interval
( const double  conflevel ) const
{ return Ostap::Math::bayes_interval ( accepted () ,
				       rejected () ,
				       conflevel   ) ; }
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
