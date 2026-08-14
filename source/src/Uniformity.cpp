// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD& STL 
// ============================================================================
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Uniformity.h"
// ============================================================================
// Local
// ============================================================================
#include "local_math.h"
// ============================================================================
/** @file
 *  Implementation file
 *  @see Analysis::Ustat
 *  @date 2021-09-27 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 */
// ============================================================================
/* Calculates the standardized t-value (z-score) using Greenwood's Spacing Statistic 
 * to test if the u-values are Uniformly distributed on [0, 1].
 *
 * Ideal for Distance-to-Nearest-Neighbor (DNN) Goodness-of-Fit tests, 
 * as it has maximum power against local clustering and density flaws.
 *
 * @param uvalues Vector of u-values in range [0, 1].
 * @return Standardized t-value (~ N(0,1) under H0). Values near 0 indicate good fit,
 *         large positive values (|t| >> 2) indicate rejection of uniformity.
 */
// ============================================================================
double Ostap::Math::greenwood_t_value
( std::vector<double> uvalues )
{
  
  const size_t N = uvalues.size();
  if ( N < 2 ) return 0.0;
  
  // 1. Sort u-values in ascending order: u_(1) <= u_(2) <= ... <= u_(N)
  std::stable_sort ( uvalues.begin() , uvalues.end() ) ;
  
  // 2. Compute the sum of squared spacings: S = sum(D_i^2)
  // Spacings are intervals between consecutive points: D_i = u_(i) - u_(i-1)
  double sum_sq_spacings = 0.0;
  double prev            = 0.0; // Left boundary (u_0 = 0.0)
  
  for ( const double u : uvalues )
  {
    const double d = u - prev;
    sum_sq_spacings += d * d;
    prev = u;
  }
  
  // Account for the final spacing to the right boundary (u_{N+1} = 1.0)
  const double last_d   = 1.0 - prev;
  sum_sq_spacings      += last_d * last_d;
  
  // 3. Theoretical Mean and Variance of Greenwood's S under H0 ~ U(0,1)
  const double n_dbl  = static_cast<double> ( N );
  const double mean_S = 2.0 / ( n_dbl + 2.0 );
  const double var_S  = ( 4.0 * n_dbl ) / ( std::pow ( n_dbl + 2.0, 2 ) *  ( n_dbl + 3.0 ) ) ;
  
  // 4. Return standardized t-value (z-score)
  return ( sum_sq_spacings - mean_S ) / std::sqrt ( var_S ) ;
}
// ============================================================================
/* Calculates standardized t-value using the Cramér–von Mises statistic (W^2).
 * Excellent all-round Goodness-of-Fit test for general shape deviations.
 *
 * @param uvalues Vector of u-values in range [0, 1].
 * @return Standardized t-value (~ N(0,1) under H0).
 */
// ============================================================================
double Ostap::Math::cramer_von_mises_t_value
( std::vector<double> uvalues )
{
  const size_t N = uvalues.size();
  if (N < 2) return 0.0;
  
  // 1. Sort u-values in ascending order
  std::stable_sort(uvalues.begin(), uvalues.end());
  
  // 2. Compute Cramér–von Mises W^2 statistic
  const double n_dbl = static_cast<double>(N);
  double sum_sq_diff = 0.0;
  
  for ( size_t i = 0 ; i < N; ++i)
  {
    // Step midpoint: e_i = (2*i + 1) / (2*N) for 0-indexed loop
    const double e_i  = ( 2.0 * static_cast<double> ( i ) + 1.0 ) / ( 2.0 * n_dbl );
    const double diff = uvalues [ i ] - e_i;
    sum_sq_diff += diff * diff;
  }
  
  const double W2 = sum_sq_diff + ( 1.0 / ( 12.0 * n_dbl ) );
  
  // 3. Theoretical Mean and Variance under H0 ~ U(0,1)
  const double mean_W2 = 1.0 / 6.0;
  const double var_W2 = (4.0 * n_dbl - 3.0) / (180.0 * n_dbl * n_dbl);
  
  // 4. Return standardized t-value
  return (W2 - mean_W2) / std::sqrt(var_W2);
}
// ============================================================================
/* Calculate standardized t-value using Logarithmic Tail Mean (-E[ln(u) + ln(1-u)]).
 * Extremely sensitive to tail behavior (u -> 0 or u -> 1).
 * Fast O(N) performance — does NOT require sorting.
 *
 * @param uvalues Vector of u-values in range [0, 1].
 * @param delta Safety margin to avoid log(0).
 * @return Standardized t-value (~ N(0,1) under H0).
 */
// ============================================================================
double Ostap::Math::logarithmic_tail_t_value
( const std::vector<double>& uvalues ,
  const double               delta   )
{
  
  const size_t N = uvalues.size();
  if ( N < 2 ) { return 0.0; } 
  
  const double n_dbl = static_cast<double>(N);
  double sum_log = 0.0;
  
  // 1. Accumulate -ln(u) - ln(1 - u)
  for ( const double u_raw : uvalues )
  {
    // Clamp value to safe interval [delta, 1 - delta]
    const double u = std::clamp ( u_raw , delta , 1.0 - delta ) ;
    sum_log += std::log(u) + std::log(1.0 - u);
  }
  
  const double mean_log_stat = -sum_log / n_dbl;
  
  // 2. Theoretical Mean = 2.0, Variance = (2*pi^2/3 - 4) / N approx 2.5797 / N
  const double mean_H0 =  2.0 ;
  const double var_H0 = ( 2.0 * s_pi * s_pi / 3.0 - 4.0) / n_dbl;
  
  // 3. Return standardized t-value
  return ( mean_log_stat - mean_H0 ) / std::sqrt ( var_H0 ) ;
}
// ============================================================================
/* Calculates standardized t-value (Stephens-adjusted) for Kolmogorov-Smirnov test.
 * Measures maximum absolute deviation between empirical and theoretical CDFs.
 *
 * @param uvalues Vector of u-values in range [0, 1].
 * @return Standardized KS statistic (~ 0 for uniform, > 1.358 for p < 0.05).
 */
// ============================================================================
double Ostap::Math::kolmogorov_smirnov_t_value
( std::vector<double> uvalues )
{
  const size_t N = uvalues.size();
  if ( N < 2 ) { return 0.0 ; } 
  
  // 1. Sort u-values
  std::stable_sort ( uvalues.begin() , uvalues.end() );
  
  const double n_dbl = static_cast<double>(N);
  double max_d = 0.0;
  
  // 2. Compute max positive and negative deviations: D+ and D-
  for ( size_t i = 0; i < N; ++i )
  {
    //
    const double i_curr  = static_cast<double> ( i + 1 ) ;
    const double d_plus  = ( i_curr / n_dbl ) - uvalues [ i ] ;
    const double d_minus = uvalues[i] - ( ( i_curr - 1.0 ) / n_dbl ) ;
    //
    if ( d_plus  > max_d ) { max_d = d_plus  ; }
    if ( d_minus > max_d ) { max_d = d_minus ; }
  }
  //
  // 3. Apply Stephens (1970) modification factor for standardized T-value
  const double sqrt_n  = std::sqrt ( n_dbl ) ;
  const double t_value = max_d * ( sqrt_n + 0.12 + 0.11 / sqrt_n ) ;
  //
  return t_value;
}

// ============================================================================
//                                                                     The END 
// ============================================================================
