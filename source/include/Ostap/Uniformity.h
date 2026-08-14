// $Id$
// ============================================================================
#ifndef OSTAP_UNIFORMITY_H 
#define OSTAP_UNIFORMITY_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <vector>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** Calculate the standardized t-value (z-score) using Greenwood's Spacing Statistic 
     * to test if the u-values are Uniformly distributed on [0, 1].
     *
     * Ideal for Distance-to-Nearest-Neighbor (DNN) Goodness-of-Fit tests, 
     * as it has maximum power against local clustering and density flaws.
     *
     * @see Greenwood, P. E. (1946).
     *      "The Statistical Inspection of Series of Consecutive Random Occurrence"
     *      Annals of Mathematical Statistics, 17(2), 257–259.
     * 
     * @param uvalues Vector of u-values in range [0, 1].
     * @return Standardized t-value (~ N(0,1) under H0). Values near 0 indicate good fit,
     *         large positive values (|t| >> 2) indicate rejection of uniformity.
     */
    double greenwood_t_value
    ( std::vector<double> uvalues ) ;
    // ========================================================================
    /** Calculate standardized t-value using the Cramér–von Mises statistic (W^2).
     * Excellent all-round Goodness-of-Fit test for general shape deviations.
     *
     * @param uvalues Vector of u-values in range [0, 1].
     * @return Standardized t-value (~ N(0,1) under H0).
     */
    // ============================================================================
    double cramer_von_mises_t_value
    ( std::vector<double> uvalues ) ;
    // ============================================================================
    /** Calculate standardized t-value using Logarithmic Tail Mean (-E[ln(u) + ln(1-u)]).
     * Extremely sensitive to tail behavior (u -> 0 or u -> 1).
     * Fast O(N) performance — does NOT require sorting.
     *
     * @param uvalues Vector of u-values in range [0, 1].
     * @param delta Safety margin to avoid log(0).
     * @return Standardized t-value (~ N(0,1) under H0).
     */
    double logarithmic_tail_t_value
      ( const std::vector<double>& uvalues       ,
        const double               delta = 1.e-8 ) ; 
    // =======================================================================
    /** Calculate standardized t-value (Stephens-adjusted) for Kolmogorov-Smirnov test.
     *  Measures maximum absolute deviation between empirical and theoretical CDFs.
     *  @param uvalues Vector of u-values in range [0, 1].
     *  @return Standardized KS statistic (~ 0 for uniform, > 1.358 for p < 0.05).
     */
    double kolmogorov_smirnov_t_value
    ( std::vector<double> uvalues ) ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_UNIFORMITY_H
// ============================================================================
