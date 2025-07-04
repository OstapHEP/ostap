// ============================================================================
#ifndef OSTAP_QUANTILE_TYPES_H 
#define OSTAP_QUANTILE_TYPES_H 1
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace QuantileTypes
  {
    // ========================================================================
    /** @enum HydnmanFanType
     * Hyndman-Fan Taxomomy of quantile estimators
     * @see https://en.wikipedia.org/wiki/Quantile
     * @see https://doi.org/10.2307%2F2684934
     */
    enum HyndmanFanType
      {
        One = 1 ,
        Two     ,
        Three   ,
        Four    ,
        Five    ,
        Six     ,
        Seven   ,
        Eight   ,
        Nine    , 
      } ;
    // =========================================================================
    /** @class ABQuantileType
     *  Helper class for implementation of ABQuantile 
     *  @see Oastap::Math::ABQuantile 
     *  @see scipy.stats.mstats
     *
     *  Typical values for alphap, betap are:
     * 
     * - (0,1) : p(k) = k/n : linear interpolation of cdf (R type 4)
     * - (.5,.5) : p(k) = (k - 1/2.)/n : piecewise linear function (R type 5)
     * - (0,0) : p(k) = k/(n+1) : (R type 6)
     * - (1,1) : p(k) = (k-1)/(n-1): p(k) = mode[F(x[k])]. (R type 7, R default)
     * - (1/3,1/3): p(k) = (k-1/3)/(n+1/3): Then p(k) ~ median[F(x[k])].
     *   The resulting quantile estimates are approximately median-unbiased 
     *   regardless of the distribution of x. (R type 8)
     * - (3/8,3/8): p(k) = (k-3/8)/(n+1/4): Blom.
     *   The resulting quantile estimates are approximately 
     *   unbiased if x is normally distributed (R type 9)
     * - 0(.4,.4) : approximately quantile unbiased (Cunnane)
     * - (.35,.35): APL, used with PWM     
     */
    class ABQuantileType 
    {
    public :
      // ========================================================================
      /// constructor 
      ABQuantileType
        ( const double alpha = 0.4 ,
          const double beta  = 0.4 ) ;
      // =======================================================================
      /// get alpha 
      inline double alpha () const { return m_alpha ; }
      /// get beta 
      inline double beta  () const { return m_beta  ; }
      // =======================================================================
    public :      
      // =======================================================================
      /// get m(p)
      inline double m  ( const double p ) const
      { return m_alpha + p *  ( 1 - m_alpha - m_beta ) ; }
      // =======================================================================
    private:
      // =======================================================================
      /// alpha-parameter 
      double m_alpha  { 0.4 } ; // alpha-parameter 
      /// beta-parameter 
      double m_beta   { 0.4 } ; // beta-parameter       
      // =======================================================================
    } ;    
    // =========================================================================
    /** @class  HarrellDavisType 
     *  Lightweight marker to use HerrellDavis quantile estimatiro 
     *  @see Ostap::Math::HerrelDavis 
     */
    class HarrellDavisType
    {
    public :
      // ======================================================================
      HarrellDavisType () ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  P2QuantileType 
     *  Lightweight marker to use (approximate) P2-quanitle estimator 
     *  @see Ostap::Math::Quantile 
     */
    class P2QuantileType
    {
    public :
      // ======================================================================
      P2QuantileType () ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                The end of namespace Ostap::QuantileTypes
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
#endif // OSTAP_QUANTILE_TYPES_H
// ============================================================================
//                                                                      The END 
// ============================================================================
