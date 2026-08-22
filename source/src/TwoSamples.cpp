// ============================================================================
// Include files
// ============================================================================
// STD& STL 
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
 *  Implement Twi-Samples (weighted) Test
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
    double cdf1 { 0.0 } ;
    double cdf2 { 0.0 } ;
    double H    { 0.0 } ;
    double dH   { 0.0 } ;
  } ;
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
      pt.cdf1 = cnt1 / N1 ;
      pt.cdf2 = cnt2 / N2 ;
      pt.dH   = static_cast<double>( c1 + c2 ) / W_tot ;
      pt.H    = static_cast<double>( cnt1 + cnt2 ) / W_tot ;
      
      result.push_back ( pt ) ;
    }
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
      pt.cdf1 = acc_w1 / sum_w1 ;
      pt.cdf2 = acc_w2 / sum_w2 ;
      pt.dH   = ( step_w1 + step_w2 ) / W_tot ;
      pt.H    = ( acc_w1 + acc_w2 ) / W_tot ;

      result.push_back ( pt ) ;
    }
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
  //
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

#undef IMPLEMENT_TWOSAMPLE_NO_EPS
#undef IMPLEMENT_TWOSAMPLE_EPS

    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
