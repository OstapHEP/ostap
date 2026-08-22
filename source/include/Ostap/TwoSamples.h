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
     *  @note Pettitt, A.N. (1976). "A two-sample Anderson-Darling rank statistic". Biometrika, 63(1): 161–168.
     */
    double kolmogorov_smirnov ( const ECDF&  e1 , const ECDF&  e2 ) ;
    double kolmogorov_smirnov ( const WECDF& e1 , const WECDF& e2 ) ;
    double kolmogorov_smirnov ( const ECDF&  e1 , const WECDF& e2 ) ;
    double kolmogorov_smirnov ( const WECDF& e1 , const ECDF&  e2 ) ;
    // ========================================================================
    /** @brief Kuiper's two-sample test statistic (\f$V\f$).
     *
     *  A variant of the Kolmogorov-Smirnov test invariant under cyclic transformations of the independent variable:
     *  \f[ V = D^+ + D^- = \sup_x (F_1(x) - F_2(x)) + \sup_x (F_2(x) - F_1(x)) \f]
     *
     *  @see https://en.wikipedia.org/wiki/Kuiper%27s_test
     *  @note Kuiper, N.H. (1960). "Tests concerning random points on a circle". 
     *        Proceedings of the Koninklijke Nederlandse Akademie van Wetenschappen, Series A, 63: 38–47.
     *        DOI: 10.1016/S1385-7258(60)50006-0
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
     *  @note Anderson, T.W. (1962). "On the distribution of the two-sample Cramer-von Mises criterion". 
     *        The Annals of Mathematical Statistics, 33(3): 1148–1159. DOI: 10.1214/aoms/1177704477
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
     *        The Annals of Mathematical Statistics, 23(2): 193–212. DOI: 10.1214/aoms/1177729437
     *  @note Pettitt, A.N. (1976). "A two-sample Anderson-Darling rank statistic". 
     *        Biometrika, 63(1): 161–168. DOI: 10.1093/biomet/63.1.161
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
     *        Z. Wahrscheinlichkeitstheorie verw. Gebiete, 47(1): 47–59. DOI: 10.1007/BF00533250
     *  @note Jager, L., Wellner, J.A. (2007). "Goodness-of-fit tests via phi-divergences". 
     *        The Annals of Statistics, 35(5): 2018–2053. DOI: 10.1214/009053607000000244
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
     *        Journal of the Royal Statistical Society, Series B, 64(3): 689–694. DOI: 10.1111/1467-9868.00353
     *  @note Zhang, J. (2006). "Powerful two-sample tests based on the likelihood ratio". 
     *        Computational Statistics & Data Analysis, 50(1): 187–199. DOI: 10.1016/j.csda.2004.07.023
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
     *        Computational Statistics & Data Analysis, 50(1): 187–199. DOI: 10.1016/j.csda.2004.07.023
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
     *        Journal of the Royal Statistical Society, Series B, 64(3): 689–694. DOI: 10.1111/1467-9868.00353
     *  @note Zhang, J. (2006). "Powerful two-sample tests based on the likelihood ratio". 
     *        Computational Statistics & Data Analysis, 50(1): 187–199. DOI: 10.1016/j.csda.2004.07.023
     */
    double ZK                 ( const ECDF&  e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    double ZK                 ( const WECDF& e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZK                 ( const ECDF&  e1 , const WECDF& e2 , const double eps = 1e-12 ) ;
    double ZK                 ( const WECDF& e1 , const ECDF&  e2 , const double eps = 1e-12 ) ;
    // ========================================================================
  } //                                        The end  of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_TWOSAMPLES_H
// ============================================================================
//                                                                      The END 
// ============================================================================
