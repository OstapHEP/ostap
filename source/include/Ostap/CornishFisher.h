// ===========================================================================
#ifndef OSTAP_CORNISH_FISHER_H
#define OSTAP_CORNISH_FISHER_H 1
// ===========================================================================
namespace Ostap
{
    // =======================================================================
    namespace Math
    {
        // ===================================================================
        /** An asymptotic Cornish-Fisher expansion  
         * of the quantile function in terms of cumulants  
         * 
         * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
         * 
         * @param p  probablity, \f$ 0 < p < 1 \f$
         * @param mu mean value (1st cumulant)
         * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
         * @return approximate quantile for probabilty p   
         */
        double cornish_fisher 
        ( const double p        , 
          const double mu       , 
          const double sigma    ) ;
        // ===================================================================
        /** An asymptotic Cornish-Fisher expansion  
         * of the quantile function in terms of cumulants  
         * 
         * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
         * 
         * @param p  probablity, \f$ 0 < p < 1 \f$
         * @param mu mean value (1st cumulant)
         * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
         * @param skewness  skewness  , 3rd cumulant 
         * @return approximate quantile for probabilty p   
         */
        double cornish_fisher 
        ( const double p        , 
          const double mu       , 
          const double sigma    , 
          const double skewness ) ;
        // ===================================================================
        /** An asymptotic Cornish-Fisher expansion  
         * of the quantile function in terms of cumulants  
         * 
         * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
         * 
         * @param p  probablity, \f$ 0 < p < 1 \f$
         * @param mu mean value (1st cumulant)
         * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
         * @param skewness  skewness  , 3rd cumulant 
         * @param kurtosis   excess kurtosis, 4thcumulant  
         * @return approximate quantile for probabilty p   
         */
        double cornish_fisher 
        ( const double p        , 
          const double mu       , 
          const double sigma    , 
          const double skewness , 
          const double kurtosis ) ;
        // ===================================================================
        /** An asymptotic Cornish-Fisher expansion  
         * of the quantile function in terms of cumulants  
         * 
         * @see https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion
         * 
         * @param p  probablity, \f$ 0 < p < 1 \f$
         * @param mu mean value (1st cumulant)
         * @param sigma positive square root fromm variance (square root fromm 2nd cumulant)
         * @param skewness  skewness  , 3rd cumulant 
         * @param kurtosis  excess kurtosis, 4th cumulant
         * @param kappa5    5th cumulant   
         * @return approximate quantile for probabilty p   
         */
        double cornish_fisher 
        ( const double p        , 
          const double mu       , 
          const double sigma    , 
          const double skewness , 
          const double kurtosis , 
          const double kappa5   ) ;
        // ===================================================================
    } //                                      The end of namespace Ostap::Math
    // =======================================================================
} //                                                The end of namespace Ostap
//============================================================================ 
#endif // OSTAP_CORNISH_FISHER_H
// ===========================================================================
//                                                                     The END  
// ===========================================================================