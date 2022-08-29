// ============================================================================
#ifndef OSTAP_QMATH_H 
#define OSTAP_QMATH_H 1
// ============================================================================
// Include files
// ============================================================================
/** @file Ostap/qMath.h
 *  Collection of functions related qto Tsallis statistics 
 *  @see https://en.wikipedia.org/wiki/Tsallis_statistics
 *  @see Umarov, Sabir; Tsallis, Constantino; Steinberg, Stanly (2008). 
 *      "On a q-Central Limit Theorem Consistent with Nonextensive 
 *      Statistical Mechanics" Milan J. Math. Birkhauser Verlag. 76: 307–328. 
 *  @see doi:10.1007/s00032-008-0087-y. S2CID 55967725..
 *  @see https://doi.org/10.1007%2Fs00032-008-0087-y
 *  @date 2022-08-28 
 *  @author anya Belyaev Ivan/Belyaev@itep.ru
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ======================================================================== 
    /** q-sum of two variables in Tsallis statistics 
     *  \f$ x + \oplus_q y = x + y + (1-q)xy\f$ 
     */
    double tsallis_qsum
    ( const double x     ,
      const double y     , 
      const double q = 1 ) ;
    // ========================================================================
    /** q-subtraction of two variables in Tsallis statistics 
     *  \f$ x + \ominus_q y = \frac{x-y}{1+(1-q)y}\f$ 
     */
    double tsallis_qsubtraction  
    ( const double x     ,
      const double y     , 
      const double q = 1 ) ;
    // =======================================================================
    /** q-product of two variables in Tsallis statistics 
     *  \f$ x + \otimes_q y = \left{ x^{1-1} + y^{1-q} -1 \right]_+^{\frac{1}{1-q}}\f$ 
     */
    double tsallis_qproduct 
    ( const double x     ,
      const double y     , 
      const double q = 1 ) ;
    // ======================================================================== 
    /** q-division of two variables in Tsallis statistics 
     *  \f$ x + \oslash_q y = \left{ x^{1-1} - y^{1-q} +1 \right]_+^{\frac{1}{1-q}}\f$ 
     */
    double tsallis_qdivision     
    ( const double x     ,
      const double y     , 
      const double q = 1 ) ;
    // ======================================================================== 
    /** q-exponent in Tsallis statistics 
     *  \f$ e_q(x) = \left[1+(1-q)x\right]_+^{\frac{1}{1-q}\f$ 
     */
    double tsallis_qexp         
    ( const double x     ,
      const double q = 1 ) ;
    // ======================================================================== 
    /** q-logarithm in Tsallis statistics 
     *  \f$ \log_q(x) = \frac{x^{1-q}-1}{1-q}\f$ 
     */
    double tsallis_qlog  
    ( const double x     ,
      const double q = 1 ) ;
    // ======================================================================== 
    /** unnormalized q-gaussian in Tsallis statistics
     *  \f$ G_q(x, \beta , q) = 
     *   e_q( - \left| \beta\right| x^2 ) \f$ 
     *  - for \f$ q=1\f$, it correspoind to Gaussian 
     *  - for \f$ q<1\f$, it it a function with the finite support 
     *        \f$ \left[ -\frac{1}{\sqrt{\beta(1-q)}}, \frac{1}{\sqrt{\beta(1-q)}}\right] \f$  
     *  - for \f$ 1<q\f$, it it a variant of generalized Student't distribution
     *  - for \f$ q=2\f$, it it a Cauchy distribution 
     *  - for \f$ 3\le q\f$, it cannot be normalzed 
     *  @see Ostap::Math::tsallis_qexp 
     */
    double tsallis_qgaussian_u
    ( const double x        ,
      const double beta     , 
      const double q    = 1 ) ;
    // ========================================================================
    /** normalized q-gaussian in Tsallis statistics for \f$ q < 3 \f$  ) 
     *  \f$ G_q(x, \beta,q) = 
     *   \frac{\sqrt{ \left| beta\right|}}{C_q} e_q( - \left| \beta\right| x^2 ) \f$ 
     *  - for \f$ q=1\f$, it correspoind to Gaussian 
     *  - for \f$ q<1\f$, it it a function with the finite support  
     *        \f$ \left[ -\frac{1}{\sqrt{\beta(1-q)}}, \frac{1}{\sqrt{\beta(1-q)}}\right] \f$  
     *  - for \f$ 1<q<3\f$, it it a variant of generalized Student't distribution
     *  - for \f$ q=2\f$, it it a Cauchy distribution 
     *  @see Ostap::Math::tsallis_qexp 
     */
    double tsallis_qgaussian
    ( const double x     ,
      const double beta  , 
      const double q     ) ;
    // ========================================================================
    /** normalized q-gaussian in Tsallis statistics for \f$ q < 3 \f$ 
     *  \f$ G_q(x,\mu. \sigma,q) = \frac{1}{2} G_q ( x-\mu , \frac{1}{2}, q \f$ 
     *  - for \f$ q<1\f$, it it a function with a finite support  
     *        \f$ \left[ \mu--\sqrt{\frac{2}{(1-q)}}, \mu+\sigma\sqrt{\frac{2}{(1-q)}}\right] \f$  
     *  - for \f$ q=1\f$, it is a Gaussian
     *  - for \f$ q=2\f$, it it a Cauchy distribution 
     *  - for \f$ 1<q<3\f$, it it a variant of generalized Student't distribution
     */
    double tsallis_qgaussian
    ( const double x     ,
      const double mu    , 
      const double sigma , 
      const double q     ) ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_QMATH_H
// ============================================================================
