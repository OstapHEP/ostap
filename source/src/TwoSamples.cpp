
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <vector>
#include <algorithm>
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ECDF.h"
#include "Ostap/TwoSamples.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// local 
// ============================================================================
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implement Two-Samples (weighted) Test statistics and distances.
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2026-02-22
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  using namespace Ostap::Math ;
  // ==========================================================================
  /// Internal structure storing cumulative probability state at pooled points
  struct PooledPoint 
  {
    double x    { 0.0 } ;
    double cdf1 { 0.0 } ;
    double cdf2 { 0.0 } ;
    double H    { 0.0 } ;
    double dH   { 0.0 } ;
    double dx   { 0.0 } ;
  } ;

  /// Internal helper to calculate interval widths dx for integration
  inline void _fill_dx ( std::vector<PooledPoint>& pts ) 
  {
    if ( pts.size() < 2 ) { return ; }
    for ( std::size_t i = 0 ; i + 1 < pts.size() ; ++i ) 
    {
      pts[i].dx = pts[i+1].x - pts[i].x ;
    }
    pts.back().dx = 0.0 ;
  }

  // ==========================================================================
  /// Internal helper for clipped Kullback-Leibler divergence
  inline double kl_div
  ( const double p   ,
    const double q   ,
    const double eps ) 
  {
    const long double p_c = std::min ( std::max ( p , eps ) , 1.0 - eps ) ;
    const long double q_c = std::min ( std::max ( q , eps ) , 1.0 - eps ) ;
    return p_c * std::log ( p_c / q_c ) + ( 1.0 - p_c ) * std::log ( ( 1.0 - p_c ) / ( 1.0 - q_c ) ) ;
  }
  // ==========================================================================
  /// Two-pointer sweep for unweighted ECDF O(N1 + N2)
  std::vector<PooledPoint> _prepare_ecdfs
  ( const ECDF& ecdf1 ,
    const ECDF& ecdf2 ) 
  {
    //
    Ostap::Assert ( ecdf1.ok () && ecdf2.ok () ,
                    "Invalid ECDF(s)"          ,
                    "Ostap::Math::TwoSamples::prepare_ecdfs" ,
                    INVALID_ECDF , __FILE__ , __LINE__       ) ;
    
    std::vector<PooledPoint> result ;
    if ( !ecdf1.ok() || !ecdf2.ok() ) { return result ; }

    const double N1    = static_cast<double>( ecdf1.size() ) ;
    const double N2    = static_cast<double>( ecdf2.size() ) ;
    const double W_tot = N1 + N2 ;
    
    result.reserve ( ecdf1.size() + ecdf2.size() ) ;

    auto it1 = ecdf1.begin () , end1 = ecdf1.end () ;
    auto it2 = ecdf2.begin () , end2 = ecdf2.end () ;

    std::size_t cnt1 = 0 , cnt2 = 0 ;

    while ( it1 != end1 || it2 != end2 ) 
    {
      double x = 0.0 ;
      if      ( it1 == end1 ) { x = *it2 ; }
      else if ( it2 == end2 ) { x = *it1 ; }
      else                    { x = std::min ( *it1 , *it2 ) ; }

      std::size_t c1 = 0 , c2 = 0 ;
      while ( it1 != end1 && *it1 == x ) { ++c1 ; ++it1 ; }
      while ( it2 != end2 && *it2 == x ) { ++c2 ; ++it2 ; }

      cnt1 += c1 ;
      cnt2 += c2 ;

      PooledPoint pt ;
      pt.x    = x ;
      pt.cdf1 = cnt1 / N1 ;
      pt.cdf2 = cnt2 / N2 ;
      pt.dH   = static_cast<double>( c1 + c2 ) / W_tot ;
      pt.H    = static_cast<double>( cnt1 + cnt2 ) / W_tot ;
      
      result.push_back ( pt ) ;
    }
    _fill_dx ( result ) ;
    return result ;
  }
  // ==========================================================================
  /// Two-pointer sweep for weighted WECDF O(N1 + N2)
  std::vector<PooledPoint> _prepare_ecdfs
  ( const WECDF& ecdf1 ,
    const WECDF& ecdf2 ) 
  {
    //
    Ostap::Assert ( ecdf1.ok () && ecdf2.ok () ,
                    "Invalid WECDF(s)"         ,
                    "Ostap::Math::TwoSamples::prepare_ecdfs" ,
                    INVALID_WECDF , __FILE__ , __LINE__      ) ;
    
    std::vector<PooledPoint> result ;
    if ( !ecdf1.ok() || !ecdf2.ok() ) { return result ; }

    const double sum_w1 = ecdf1.sumw() ;
    const double sum_w2 = ecdf2.sumw() ;
    const double W_tot  = sum_w1 + sum_w2 ;
    if ( W_tot <= 0.0 || sum_w1 <= 0.0 || sum_w2 <= 0.0 ) { return result ; }
    
    result.reserve ( ecdf1.size() + ecdf2.size() ) ;

    auto it1 = ecdf1.begin () , end1 = ecdf1.end () ;
    auto it2 = ecdf2.begin () , end2 = ecdf2.end () ;

    double acc_w1 = 0.0 , acc_w2 = 0.0 ;

    while ( it1 != end1 || it2 != end2 ) 
    {
      double x = 0.0 ;
      if      ( it1 == end1 ) { x = it2->first ; }
      else if ( it2 == end2 ) { x = it1->first ; }
      else                    { x = std::min ( it1->first , it2->first ) ; }

      double step_w1 = 0.0 , step_w2 = 0.0 ;
      while ( it1 != end1 && it1->first == x ) { step_w1 += it1->second ; ++it1 ; }
      while ( it2 != end2 && it2->first == x ) { step_w2 += it2->second ; ++it2 ; }

      acc_w1 += step_w1 ;
      acc_w2 += step_w2 ;

      PooledPoint pt ;
      pt.x    = x ;
      pt.cdf1 = acc_w1 / sum_w1 ;
      pt.cdf2 = acc_w2 / sum_w2 ;
      pt.dH   = ( step_w1 + step_w2 ) / W_tot ;
      pt.H    = ( acc_w1 + acc_w2 ) / W_tot ;

      result.push_back ( pt ) ;
    }
    _fill_dx ( result ) ;
    return result ;
  }

  // ==========================================================================
  // Computation engines
  // ==========================================================================

  double compute_ks
  ( const std::vector<PooledPoint>& pts ) 
  {
    double max_diff = 0.0 ;
    for ( const auto& pt : pts ) 
    { max_diff = std::max ( max_diff , std::abs ( pt.cdf1 - pt.cdf2 ) ) ; }
    return max_diff ;
  }

  double compute_kuiper
  ( const std::vector<PooledPoint>& pts ) 
  {
    double d_plus  = 0.0 ;
    double d_minus = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double diff = pt.cdf1 - pt.cdf2 ;
      d_plus  = std::max ( d_plus  ,  diff ) ;
      d_minus = std::max ( d_minus , -diff ) ;
    }
    return std::max ( 0.0 , d_plus ) + std::max ( 0.0 , d_minus ) ;
  }

  double compute_cvm
  ( const std::vector<PooledPoint>& pts ) 
  {
    double cvm = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double diff = pt.cdf1 - pt.cdf2 ;
      cvm += diff * diff * std::abs ( pt.dH ) ;
    }
    return cvm ;
  }

  double compute_ad
  ( const std::vector<PooledPoint>& pts , const double eps ) 
  {
    double ad = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double H_c  = std::min ( std::max ( pt.H , eps ) , 1.0 - eps ) ;
      const double var  = H_c * ( 1.0 - H_c ) ;
      const double diff = pt.cdf1 - pt.cdf2 ;
      ad += ( diff * diff / var ) * std::abs ( pt.dH ) ;
    }
    return ad ;
  }

  double compute_bj
  ( const std::vector<PooledPoint>& pts , const double eps ) 
  {
    double bj = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double kl1 = kl_div ( pt.cdf1 , pt.H , eps ) ;
      const double kl2 = kl_div ( pt.cdf2 , pt.H , eps ) ;
      bj = std::max ( bj , std::max ( kl1 , kl2 ) ) ;
    }
    return bj ;
  }

  double compute_za
  ( const std::vector<PooledPoint>& pts , const double eps ) 
  {
    double za = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double kl1 = kl_div ( pt.cdf1 , pt.H , eps ) ;
      const double kl2 = kl_div ( pt.cdf2 , pt.H , eps ) ;
      za += ( kl1 + kl2 ) * std::abs ( pt.dH ) ;
    }
    return za ;
  }

  double compute_zk
  ( const std::vector<PooledPoint>& pts , const double eps ) 
  {
    double zk = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double kl1 = kl_div ( pt.cdf1 , pt.H , eps ) ;
      const double kl2 = kl_div ( pt.cdf2 , pt.H , eps ) ;
      zk = std::max ( zk , kl1 + kl2 ) ;
    }
    return zk ;
  }

  double compute_w1
  ( const std::vector<PooledPoint>& pts ) 
  {
    double w1 = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      if ( pt.dx <= 0.0 ) { continue ; }
      w1 += std::abs ( pt.cdf1 - pt.cdf2 ) * pt.dx ;
    }
    return w1 ;
  }

  double compute_w2
  ( const std::vector<PooledPoint>& pts ) 
  {
    double w2_sq = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      if ( pt.dx <= 0.0 ) { continue ; }
      const double diff = pt.cdf1 - pt.cdf2 ;
      w2_sq += ( diff * diff ) * pt.dx ;
    }
    return std::sqrt ( w2_sq ) ;
  }
  
  // ==========================================================================
  /// Internal computation engine for Rosenblatt-type two-sample test statistic.
  /// Accumulates squared pointwise differences between two empirical CDFs.
  double compute_rosenblatt
  ( const std::vector<PooledPoint>& pts ) 
  {
    double result = 0.0 ;
    // Iterate over all pooled points generated by the two-pointer sweep
    for ( const auto& pt : pts ) 
    {
      // Compute the difference between the first and second empirical CDFs
      const double diff = pt.cdf1 - pt.cdf2 ;
      
      // Accumulate the squared difference as an integral analogue component
      result += diff * diff ; 
    }
    return result ;
  }

  // ==========================================================================
  /// Internal computation engine for Watson's U^2 two-sample test statistic.
  /// Computes Cramér-von Mises type discrepancy adjusted for circular/cyclic invariance.
  double compute_watson_u2
  ( const std::vector<PooledPoint>& pts ) 
  {
    if ( pts.empty() ) { return 0.0 ; }

    // First pass: compute the weighted mean difference (D_bar) over the pooled distribution
    double d_bar = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double diff = pt.cdf1 - pt.cdf2 ;
      d_bar += diff * std::abs ( pt.dH ) ;
    }

    // Second pass: accumulate squared deviations from the mean difference
    double u2 = 0.0 ;
    for ( const auto& pt : pts ) 
    {
      const double diff = ( pt.cdf1 - pt.cdf2 ) - d_bar ;
      u2 += diff * diff * std::abs ( pt.dH ) ;
    }
    
    return u2 ;
  }
  // ==========================================================================
} // The end of anonymous namespace
// ============================================================================
// Public API Implementation
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
#define IMPLEMENT_TWOSAMPLE_NO_EPS( NAME , ENGINE )                           \
    double NAME ( const ECDF&  e1 , const ECDF&  e2 )                         \
    { return ENGINE ( _prepare_ecdfs ( e1 , e2 ) ) ; }                        \
    double NAME ( const WECDF& e1 , const WECDF& e2 )                         \
    { return ENGINE ( _prepare_ecdfs ( e1 , e2 ) ) ; }                        \
    double NAME ( const ECDF&  e1 , const WECDF& e2 )                         \
    { return NAME ( WECDF ( e1 ) , e2 ) ; }                                   \
    double NAME ( const WECDF& e1 , const ECDF&  e2 )                         \
    { return NAME ( e1 , WECDF ( e2 ) ) ; }

#define IMPLEMENT_TWOSAMPLE_EPS( NAME , ENGINE )                              \
    double NAME ( const ECDF&  e1 , const ECDF&  e2 , const double eps )      \
    { return ENGINE ( _prepare_ecdfs ( e1 , e2 ) , eps ) ; }                  \
    double NAME ( const WECDF& e1 , const WECDF& e2 , const double eps )      \
    { return ENGINE ( _prepare_ecdfs ( e1 , e2 ) , eps ) ; }                  \
    double NAME ( const ECDF&  e1 , const WECDF& e2 , const double eps )      \
    { return NAME ( WECDF ( e1 ) , e2 , eps ) ; }                             \
    double NAME ( const WECDF& e1 , const ECDF&  e2 , const double eps )      \
    { return NAME ( e1 , WECDF ( e2 ) , eps ) ; }

    // 1. Kolmogorov-Smirnov
    IMPLEMENT_TWOSAMPLE_NO_EPS ( kolmogorov_smirnov , compute_ks )

    // 2. Kuiper
    IMPLEMENT_TWOSAMPLE_NO_EPS ( kuiper , compute_kuiper )

    // 3. Cramér-von Mises
    IMPLEMENT_TWOSAMPLE_NO_EPS ( cramer_von_mises , compute_cvm )

    // 4. Anderson-Darling
    IMPLEMENT_TWOSAMPLE_EPS    ( anderson_darling , compute_ad )

    // 5. Berk-Jones
    IMPLEMENT_TWOSAMPLE_EPS    ( berk_jones , compute_bj )

    // 6. Zhang's ZA
    IMPLEMENT_TWOSAMPLE_EPS    ( ZA , compute_za )

    // 7. Zhang's ZC
    IMPLEMENT_TWOSAMPLE_EPS    ( ZC , compute_ad )

    // 8. Zhang's ZK
    IMPLEMENT_TWOSAMPLE_EPS    ( ZK , compute_zk )

    // 9. Wasserstein-1
    IMPLEMENT_TWOSAMPLE_NO_EPS ( wasserstein1 , compute_w1 )

    // 10. Wasserstein-2
    IMPLEMENT_TWOSAMPLE_NO_EPS ( wasserstein2 , compute_w2 )
    
    // 11. Rosenblatt-type two-sample test statistic
    IMPLEMENT_TWOSAMPLE_NO_EPS ( rosenblatt , compute_rosenblatt )
    
    // 12. Watson's U^2 two-sample test statistic
    IMPLEMENT_TWOSAMPLE_NO_EPS ( watson_u2 , compute_watson_u2 )
    
#undef IMPLEMENT_TWOSAMPLE_NO_EPS
#undef IMPLEMENT_TWOSAMPLE_EPS

    // ========================================================================    
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================

// ============================================================================
// Rank-based estimators (for unweighted data only)
// ============================================================================

// ===========================================================================
/*  Mann-Whitney U two-sample test statistic.
 *
 *  Calculates the Mann-Whitney U statistic (also known as Wilcoxon rank-sum statistic)
 *  for two unweighted empirical cumulative distribution functions (ECDF).
 *  To optimize performance, the algorithm iterates over the shorter sample and
 *  queries ranks in the longer sample using binary search.
 *
 *  @param e1 First ECDF sample
 *  @param e2 Second ECDF sample
 *  @return Minimum Mann-Whitney U statistic value (\f$ \min(U_1, U_2) \f$)
 *
 *  @see https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
 *  @note Mann, H.B., Whitney, D.R. (1947). "On a test of whether one of two random variables is stochastically larger than the other".
 *        The Annals of Mathematical Statistics, 18(1): 50–60.
 *  @note Wilcoxon, F. (1945). "Individual comparisons by ranking methods".
 *        Biometrics Bulletin, 1(6): 80–83.
 */
// ===========================================================================
double Ostap::Math::mann_whitney
( const Ostap::Math::ECDF& e1 ,
  const Ostap::Math::ECDF& e2 ) 
{
  if ( !e1.ok() || !e2.ok() ) { return 0.0 ; }
  
  // Determine the shorter and longer samples to optimize binary search complexity
  const bool               e1_shorter   = ( e1.size() <= e2.size() ) ;
  const Ostap::Math::ECDF& short_sample = e1_shorter ? e1 : e2 ;
  const Ostap::Math::ECDF& long_sample  = e1_shorter ? e2 : e1 ;

  double u_short = 0.0 ;
  // Sum the ranks of elements from the short sample evaluated in the long sample
  for ( const double x : short_sample )
  { u_short += long_sample.rank ( x ) ; }
  
  const double n1 = static_cast<double> ( e1.size() ) ;
  const double n2 = static_cast<double> ( e2.size() ) ;
  
  // Derive the alternative statistic using the exact algebraic identity: U1 + U2 = n1 * n2
  const double u_long = n1 * n2 - u_short ;
  
  const double u_mean = ( n1 * n2 ) / 2.0;

  return std::abs ( std::min ( u_short , u_long ) - u_mean ) ;
  
}
// ===========================================================================
/* Symmetrized version of Savage two-sample test statistic (Exponential Scores Test).
 *
 *  Calculates the Savage test statistic for two unweighted empirical cumulative
 *  distribution functions (ECDF). The Savage test is a rank-based nonparametric test
 *  utilizing exponential scores (expected values of order statistics from an exponential
 *  distribution), making it particularly powerful for detecting differences in the
 *  right tails of distributions.
 *
 *  @param e1 First ECDF sample
 *  @param e2 Second ECDF sample
 *  @return Savage test statistic value for the first sample
 *
 *  @note Savage, I.R. (1956). "Contributions to the theory of rank order statistics — the two-sample case".
 *        The Annals of Mathematical Statistics, 27(3): 590–615.
 *  @note Hájek, J., Šidák, Z. (1967). "Theory of Rank Tests". Academic Press.
 */
// ===========================================================================
double Ostap::Math::savage
( const Ostap::Math::ECDF& e1 ,
  const Ostap::Math::ECDF& e2 ) 
{
  if ( !e1.ok() || !e2.ok() ) { return 0.0 ; }
  
  // Merge both samples into a pooled ordered ECDF to determine global ranks
  const Ostap::Math::ECDF e_pool = e1.merge ( e2 ) ;
  const std::size_t       N      = e_pool.size () ;
  if ( 0 == N ) { return 0.0 ; }
  
  // Precompute Savage (exponential) scores for all ranks from 1 to N:
  // E_r = \sum_{j=r}^{N} (1 / j)
  std::vector<long double> savage_scores ( N + 1 , 0.0L ) ;
  long double running_sum = 0.0 ;
  for ( std::size_t j = N ; 1 <= j ; --j )
  {
    running_sum += 1.0L / static_cast<long double>( j ) ;
    savage_scores [ j ] = running_sum ;
  }
  
  // Accumulate Savage scores for elements belonging to the first sample (e1)
  long double s1 = 0.0 ;
  for ( const double x : e1 )
  {
    const std::size_t r = e_pool.rank ( x ) ;
    if ( 1 <= r && r <= N ) { s1 += savage_scores [ r ] ; }
  }

  // Accumulate Savage scores for elements belonging to the second sample (e2)
  long double s2 = 0.0 ;
  for ( const double x : e2 )
  {
    const std::size_t r = e_pool.rank ( x ) ;
    if ( 1 <= r && r <= N ) { s2 += savage_scores [ r ] ; }
  }

  const long double n1 = static_cast<long double>( e1.size() ) ;
  const long double n2 = static_cast<long double>( e2.size() ) ;
  
  const long double nn = n1 + n2 ;
  const long double f1 = n1 / nn ;
  const long double f2 = n2 / nn ; 
  
  // symmetrized Savage score 
  return std::abs ( f2 * s1 - f1 * s2 ) ;
}
// ===========================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
