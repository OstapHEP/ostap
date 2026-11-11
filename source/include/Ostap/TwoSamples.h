// ============================================================================
#ifndef OSTAP_TWOSAMPLES_H
#define OSTAP_TWOSAMPLES_H 1
// ============================================================================
#include "Ostap/ECDF.h"
// ============================================================================
/** @file Ostap/TwoSamples.h
 *  @brief Collection of non-parametric two-sample goodness-of-fit test statistics.
 *
 *  This header provides non-parametric statistics comparing two empirical 
 *  cumulative distribution functions (unweighted ECDF or weighted WECDF).
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2026-02-22
 */
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    // Forward declarations 
    // ========================================================================
    class  ECDF ; 
    class WECDF ;    
    // ========================================================================
    /** @brief Kolmogorov-Smirnov two-sample statistic (\f$D\f$).
     *
     *  Calculates the supremum difference between two empirical CDFs:
     *  \f[ D = \sup_x |F_1(x) - F_2(x)| \f]
     *
     *  @see https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
     *  @note Kolmogorov, A.N. (1933). "Sulla determinazione empirica di una legge di distribuzione". 
     *        Giornale dell'Istituto Italiano degli Attuari, 4: 83–91.
     *  @note Smirnov, N.V. (1939). "On the estimation of the discrepancy between empirical curves of distribution". 
     *        Rec. Math. (Matematicheskii Sbornik) N.S., 6(48): 3–26.
     */
    double kolmogorov_smirnov ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double kolmogorov_smirnov ( const WECDF& e1 , const WECDF& e2 ) ;
    double kolmogorov_smirnov ( const ECDF&  e1 , const WECDF& e2 ) ;
    double kolmogorov_smirnov ( const WECDF& e1 , const ECDF&  e2 ) ;
    // ========================================================================
    /** @brief Kuiper's two-sample test statistic (\f$V\f$).
     *
     *  A variant of the Kolmogorov-Smirnov test invariant under cyclic transformations:
     *  \f[ V = D^+ + D^- = \sup_x (F_1(x) - F_2(x)) + \sup_x (F_2(x) - F_1(x)) \f]
     *
     *  @see https://en.wikipedia.org/wiki/Kuiper%27s_test
     *  @note Kuiper, N.H. (1960). "Tests concerning random points on a circle". 
     *        Proceedings of the Koninklijke Nederlandse Akademie van Wetenschappen, Series A, 63: 38–47.
     */
    double kuiper             ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double kuiper             ( const WECDF& e1 , const WECDF& e2 ) ;
    double kuiper             ( const ECDF&  e1 , const WECDF& e2 ) ;
    double kuiper             ( const WECDF& e1 , const ECDF&  e2 ) ;
    // ========================================================================
    /** @brief Cramér-von Mises two-sample statistic (\f$T\f$ or \f$W^2\f$).
     *
     *  Integrated squared difference between two empirical CDFs with respect to pooled CDF \f$H(x)\f$:
     *  \f[ T = \int_{-\infty}^{\infty} (F_1(x) - F_2(x))^2 dH(x) \f]
     *
     *  @see https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93von_Mises_criterion
     *  @note Cramér, H. (1928). "On the composition of elementary errors". Skandinavisk Aktuarietidskrift, 11: 13–74.
     *  @note von Mises, R.E. (1931). "Wahrscheinlichkeitsrechnung und ihre Anwendung in der Statistik und theoretischen Physik".
     */
    double cramer_von_mises   ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double cramer_von_mises   ( const WECDF& e1 , const WECDF& e2 ) ;
    double cramer_von_mises   ( const ECDF&  e1 , const WECDF& e2 ) ;
    double cramer_von_mises   ( const WECDF& e1 , const ECDF&  e2 ) ;
    // ========================================================================
    /** @brief Anderson-Darling two-sample statistic (\f$A^2\f$).
     *
     *  Weighted integrated squared difference giving higher weight to distribution tails:
     *  \f[ A^2 = \int_{-\infty}^{\infty} \frac{(F_1(x) - F_2(x))^2}{H(x)(1 - H(x))} dH(x) \f]
     *
     *  @param eps Threshold to clip \f$H(x)\f$ into \f$[\epsilon, 1-\epsilon]\f$ avoiding division by zero.
     *
     *  @see https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test
     *  @note Anderson, T.W., Darling, D.A. (1952). "Asymptotic theory of certain 'goodness of fit' criteria based on stochastic processes". 
     *        The Annals of Mathematical Statistics, 23(2): 193–212.
     */
    double anderson_darling   ( const ECDF&  e1 , const ECDF&  e2 , const double eps = 1e-10 ) ;
    double anderson_darling   ( const WECDF& e1 , const WECDF& e2 , const double eps = 1e-10 ) ;
    double anderson_darling   ( const ECDF&  e1 , const WECDF& e2 , const double eps = 1e-10 ) ;
    double anderson_darling   ( const WECDF& e1 , const ECDF&  e2 , const double eps = 1e-10 ) ;
    // ========================================================================
    /** @brief Berk-Jones two-sample statistic (\f$R_{BJ}\f$).
     *
     *  Likelihood ratio-based statistic defined as maximum pointwise Kullback-Leibler divergence:
     *  \f[ R_{BJ} = \sup_x \max \left( K(F_1(x), H(x)), K(F_2(x), H(x)) \right) \f]
     *
     *  @param eps Threshold to clip CDF values into \f$[\epsilon, 1-\epsilon]\f$ avoiding \f$\ln(0)\f$.
     *
     *  @see https://en.wikipedia.org/wiki/Berk%E2%80%93Jones_statistic
     *  @note Berk, R.H., Jones, L.L. (1979). "Goodness-of-fit statistics based on maximum likelihood ratios". 
     *        Z. Wahrscheinlichkeitstheorie verw. Gebiete, 47(1): 47–59.
     */
    double berk_jones         ( const ECDF&  e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    double berk_jones         ( const WECDF& e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double berk_jones         ( const ECDF&  e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double berk_jones         ( const WECDF& e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    // ========================================================================
    /** @brief Zhang's \f$Z_A\f$ two-sample test statistic.
     *
     *  Information-theoretic integral analogue of Anderson-Darling test based on likelihood ratios:
     *  \f[ Z_A = \int_{-\infty}^{\infty} \left[ K(F_1(x), H(x)) + K(F_2(x), H(x)) \right] dH(x) \f]
     *
     *  @param eps Threshold to clip CDF values into \f$[\epsilon, 1-\epsilon]\f$.
     *
     *  @note Zhang, J. (2002). "Powerful Goodness-of-Fit Tests Based on Likelihood Ratio". 
     *        Journal of the Royal Statistical Society, Series B, 64(3): 689–694.
     */
    double ZA                 ( const ECDF&  e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    double ZA                 ( const WECDF& e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZA                 ( const ECDF&  e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZA                 ( const WECDF& e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    // ========================================================================
    /** @brief Zhang's \f$Z_C\f$ two-sample test statistic.
     *
     *  Information-theoretic analogue of Cramér-von Mises test. Mathematically equivalent 
     *  to the two-sample Anderson-Darling statistic.
     *
     *  @param eps Threshold to clip CDF values into \f$[\epsilon, 1-\epsilon]\f$.
     *
     *  @note Zhang, J. (2006). "Powerful two-sample tests based on the likelihood ratio". 
     *        Computational Statistics & Data Analysis, 50(1): 187–199.
     */
    double ZC                 ( const ECDF&  e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    double ZC                 ( const WECDF& e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZC                 ( const ECDF&  e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZC                 ( const WECDF& e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    // ========================================================================
    /** @brief Zhang's \f$Z_K\f$ two-sample test statistic.
     *
     *  Information-theoretic supremum statistic based on maximum Kullback-Leibler divergence:
     *  \f[ Z_K = \sup_x \left[ K(F_1(x), H(x)) + K(F_2(x), H(x)) \right] \f]
     *
     *  @param eps Threshold to clip CDF values into \f$[\epsilon, 1-\epsilon]\f$.
     *
     *  @note Zhang, J. (2002). "Powerful Goodness-of-Fit Tests Based on Likelihood Ratio". 
     *        Journal of the Royal Statistical Society, Series B, 64(3): 689–694.
     */
    double ZK                 ( const ECDF&  e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    double ZK                 ( const WECDF& e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZK                 ( const ECDF&  e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZK                 ( const WECDF& e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    // ========================================================================
    /** @brief First Wasserstein distance (\f$W_1\f$ / Earth Mover's Distance).
     *
     *  Calculates the 1st Wasserstein distance between two empirical distributions:
     *  \f[ W_1(F_1, F_2) = \int_{-\infty}^{\infty} |F_1(x) - F_2(x)| \, dx \f]
     *
     *  @see https://en.wikipedia.org/wiki/Wasserstein_metric
     */
    double wasserstein1       ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double wasserstein1       ( const WECDF& e1 , const WECDF& e2 ) ;
    double wasserstein1       ( const ECDF&  e1 , const WECDF& e2 ) ;
    double wasserstein1       ( const WECDF& e1 , const ECDF&  e2 ) ;
    // ========================================================================
    /** @brief Second Wasserstein distance (\f$W_2\f$).
     *
     *  Calculates the square root of the integrated squared differences between empirical CDFs:
     *  \f[ W_2(F_1, F_2) = \left( \int_{-\infty}^{\infty} (F_1(x) - F_2(x))^2 \, dx \right)^{1/2} \f]
     *
     *  @see https://en.wikipedia.org/wiki/Wasserstein_metric
     */
    double wasserstein2       ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double wasserstein2       ( const WECDF& e1 , const WECDF& e2 ) ;
    double wasserstein2       ( const ECDF&  e1 , const WECDF& e2 ) ;
    double wasserstein2       ( const WECDF& e1 , const ECDF&  e2 ) ;
    // ========================================================================
    /** @brief Rosenblatt-type two-sample test statistic (\f$R\f$).
     *
     *  Calculates cumulative squared discrepancies between two empirical CDFs 
     *  evaluated on the pooled sample points:
     *  \f[ R = \sum_i (F_1(x_i) - F_2(x_i))^2 \f]
     *
     *  @see Rosenblatt, M. (1952)
     */
    double rosenblatt         ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double rosenblatt         ( const WECDF& e1 , const WECDF& e2 ) ;
    double rosenblatt         ( const ECDF&  e1 , const WECDF& e2 ) ;
    double rosenblatt         ( const WECDF& e1 , const ECDF&  e2 ) ;
    // ========================================================================
    /** @brief Watson's two-sample test statistic (\f$U^2\f$).
     *
     *  A Cramér-von Mises type variant invariant under cyclic transformations:
     *  \f[ U^2 = \int_{-\infty}^{\infty} \left( (F_1(x) - F_2(x)) - \bar{D} \right)^2 dH(x) \f]
     *
     *  @see https://en.wikipedia.org/wiki/Watson_test
     *  @note Watson, G.S. (1961, 1962).
     */
    double watson_u2          ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double watson_u2          ( const WECDF& e1 , const WECDF& e2 ) ;
    double watson_u2          ( const ECDF&  e1 , const WECDF& e2 ) ;
    double watson_u2          ( const WECDF& e1 , const ECDF&  e2 ) ;    
    // ========================================================================

    // ========================================================================
    // Rank-based estimators (for unweighted data only)
    // ========================================================================
    
    // ========================================================================
    /**  Mann-Whitney U two-sample test statistic.
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
    double mann_whitney ( const ECDF& e1 , const ECDF& e2 ) ;
    
    // ========================================================================
    /** Symmetrized version of Savage two-sample test statistic (Exponential Scores Test).
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
    // ========================================================================
    double savage ( const ECDF& e1 , const ECDF& e2 ) ; 

    // ========================================================================
  } //                                        The end  of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_TWOSAMPLES_H
// ============================================================================
//                                                                      The END 
// ============================================================================
