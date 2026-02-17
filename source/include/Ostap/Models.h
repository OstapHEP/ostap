// ============================================================================
#ifndef OSTAP_MODELS_H
#define OSTAP_MODELS_H 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/NSphere.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein1D.h"
#include "Ostap/Workspace.h"
#include "Ostap/PhaseSpace.h"
#include "Ostap/BSpline.h"
#include "Ostap/MoreMath.h"
// ============================================================================
/** @file Ostap/Models.h
 *  Set of useful continiod PDFs/models
 *
 *  Models defined for \f$\left[ +\infty, -\infty \right]\f$: 
 *  @see Ostap::Math::Gumbel
 * 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Gumbel
     *  Gumbel distribution
     *  @see https://en.wikipedia.org/wiki/Gumbel_distribution
     *  \f$   G(x;\mu,\beta) = \frac{1}{\left|\beta\right|} e^{ -(z+e^{-z})}\f$,
     *  where \f$z= \frac{x-\mu}{\beta}\f$
     *  
     *  Very useful important case:
     *  If  x is distributed according  to \f$ f(x) \propto e^{-\tau x} \f$, 
     *  then z, \f$ z  =   log(x) \f$, is distributed according to 
     *  \f$ F(z) = G(x, -\log(\tau), 1 ) \f$ 
     * 
     *  As a   result, if    x is a sum of exponential components 
     *  with different slopes, the transforomation \f$ z=log(x) \f$ will convert each 
     *  exponential components into bump-like structure
     */  
    class Gumbel
    {
    public:
      // ======================================================================
      /** constructor  from all parameters 
       *  @param mu location, bias parameter 
       *  @param beta scale parameter 
       */
      Gumbel 
      ( const double mu   = 0 , 
        const double beta = 1 );
      // ======================================================================
    public: // primary getters 
      // ======================================================================
      inline double mu       () const { return m_mu      ; }
      inline double peak     () const { return   mu   () ; }
      inline double location () const { return   mu   () ; }
      inline double beta     () const { return m_beta    ; }
      inline double scale    () const { return   beta () ; }
      // ======================================================================
    public: // settetrs 
      // ======================================================================
      bool setMu       ( const double value ) ;
      bool setBeta     ( const double value ) ;
      // ======================================================================
      bool setPeak     ( const double value ) { return setMu   ( value ) ; }
      bool setLocation ( const double value ) { return setMu   ( value ) ; }
      bool setScale    ( const double value ) { return setBeta ( value ) ; }
      // ======================================================================
    public: // derived quantities 
      // ======================================================================
      double mean        () const ;
      double median      () const ;
      double variance    () const ;
      double sigma       () const ;
      double skewness    () const ;
      // ======================================================================
      inline double mode        () const { return mu       () ; }
      inline double rms         () const { return sigma    () ; }
      inline double dispersion  () const { return variance () ; }
      inline double sigma2      () const { return variance () ; }
      inline double kurtosis    () const { return 12.0 / 5    ; }
      // ======================================================================
    public: // the main block 
      // ======================================================================
      /// get a value for the function 
      inline double operator() ( const double x )  const { return pdf ( x ) ; }
      inline double evaluate   ( const double x )  const { return pdf ( x ) ; }
      /// get a value for the function      
      double pdf  ( const double x ) const ;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      double cdf      ( const double x ) const ;
      /// get the integral
      double integral () const { return 1 ; }
      /// get the integral between low and high limits
      double integral
      ( const double low  ,
	const double high ) const ;
      // ======================================================================
      /** quantile function
       *  @parameter p  probability \f$ 0 < p < 1 \f$
       */
      double quantile ( const double p ) const ;	
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag ()   const ;
      // ======================================================================
    private:      
      // ======================================================================
      /// mode, location, bias parameter 
      double m_mu   ; // mode, location, bias parameter 
      /// scale parameter 
      double m_beta ; // scale parameter
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GammaDist
     *  Gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  - we have added also shift-parameter 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  GammaDist
    {
    public:
      // ======================================================================
      /** constructor from scale & shape parameters
       *  param k       \f$k\f$      parameter (shape)
       *  param theta   \f$\theta\f$ parameter (scale)
       *  param shift   shift parameter
       */
      GammaDist 
      ( const double k     = 2 ,   // shape parameter
        const double theta = 1 ,   // scale parameter
	const double shift = 0 ) ; // shift parametet 
      // ======================================================================
    public:
      // ======================================================================
      double pdf        ( const double x ) const ;
      /// calculate gamma distribution shape
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      inline double evaluate   ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // variables
      // ======================================================================
      inline double k          () const  { return m_k       ; }
      inline double theta      () const  { return m_theta   ; }
      inline double shift      () const  { return m_shift   ; }
      // ======================================================================
      inline double alpha      () const  { return m_k       ; }
      inline double beta       () const  { return 1/m_theta ; }
      // ======================================================================
    public:
      // ======================================================================
      /// mean value 
      inline double mean       () const  { return m_shift + m_k * m_theta  ; }
      /// variance 
      inline double variance   () const  { return dispersion ()            ; }
      /// variance 
      inline double dispersion () const  { return m_k * m_theta * m_theta  ; }
      /// kurtosis 
      inline double kurtosis   () const  { return 6 / alpha  ()            ; }
      // 
      double        rms        () const  ;
      double        skewness   () const  ;
      double        mode       () const  ; 
      // ======================================================================
      inline double xmin       () const { return m_shift  ; }
      // ======================================================================      
    public:
      // ======================================================================
      /** effective \f$ \chi^2 \f$-parameters
       *  If   \f$ Q  \sim \chi^2(\nu)\f$  and c is a positive constant,
       *  than \f$ cQ \sim \Gamma (k = \nu/2, \theta = 2c) \f$
       *  @attention for shift!=0  NaN is returned
       */
      double nu () const ;
      /** effective \f$ \chi^2 \f$-parameters
       *  If   \f$ Q  \sim \chi^2(\nu)\f$  and c is a positive constant,
       *  than \f$ cQ \sim \Gamma (k = \nu/2, \theta = 2c) \f$
       *  @attention for shift!=0  NaN is returned
       */
      double c  () const ; 
      // ======================================================================
    public:
      // ======================================================================
      bool setK     ( const double value  ) ;
      bool setTheta ( const double value  ) ;
      bool setShift ( const double value  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral 
      ( const double low  ,
        const double high ) const ;
      /// get CDF 
      double cdf 
      ( const double x    ) const ;
      // ======================================================================
    public: // quantiles
      // ======================================================================
      /// calculate the quantile   (0<p<1)
      double quantile ( const double p ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// shape
      double m_k      { 2 } ; // shape
      /// scale
      double m_theta  { 1 } ; // scale
      /// shift
      double m_shift  { 0 } ; // shift 
      // ======================================================================
    private:
      // ======================================================================
      /// auxillary intermediate parameter  -log Gamma ( k ) 
      mutable double m_lgk { 0 } ; 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GramCharlierA4
     *  Gram-Charlier type A approximation
     *  http://en.wikipedia.org/wiki/Edgeworth_series
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-06-13
     */
    class  GramCharlierA
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param mean   the mean value for distribution
       *  @param sigma  the sigma
       *  @param kappa3 the standartized 3rd cumulant
       *  @param kappa4 the standartized 4th cumulant
       */
      GramCharlierA  
      ( const double mean   = 0 ,
        const double sigma  = 1 ,
        const double kappa3 = 1 ,
        const double kappa4 = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Gram-Charlier type A approximation
      double pdf         ( const double x ) const ;
      /// evaluate Gram-Charlier type A approximation
      inline double operator () ( const double x ) const { return pdf ( x ) ; }
      inline double evaluate    ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double   mean   () const { return m_mean          ; }
      inline double   m0     () const { return   mean   ()     ; }
      inline double   peak   () const { return   mean   ()     ; }
      inline double   sigma  () const { return m_sigma         ; }
      inline double   kappa3 () const { return m_kappa3        ; }
      inline double   kappa4 () const { return m_kappa4        ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool setM0      ( const double value ) ;
      bool setMean    ( const double value ) { return setM0   ( value ) ; }
      bool setPeak    ( const double value ) { return setM0   ( value ) ; }
      bool setMass    ( const double value ) { return setPeak ( value ) ; }
      //
      bool setSigma   ( const double value ) ;
      bool setKappa3  ( const double value ) ;
      bool setKappa4  ( const double value ) ;
      // ======================================================================
    public: //
      // ======================================================================
      /// get (possibly truncated) integral
      double integral () const ;
      /// get integral between low and high
      double integral
      ( const double low ,
	const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// mean value
      double              m_mean   ;           //         mean
      /// rms
      double              m_sigma  ;           //          rms
      /// the standartized 3rd cumulant
      double              m_kappa3 ;           // the standartized 3rd cumulant
      /// the standartized 4th cumulant
      double              m_kappa4 ;           // the standartized 4th cumulant
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LogGammaDist
     *  Distribution for log(x) where x has gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  LogGammaDist
    {
    public:
      // ======================================================================
      /** constructor from scale & shape parameters
       *  param k      \f$k\f$ parameter (shape)
       *  param theta  \f$\theta\f$ parameter (scale)
       */
      LogGammaDist
      ( const double k     = 2 ,   // shape parameter
        const double theta = 1 ,   // scale parameter
	const double shift = 0 ) ; // shift parameter
      /// destructor
      virtual ~LogGammaDist() ;  // desctructor
      // ======================================================================
    public:
      // ======================================================================
      /// calculate log-gamma distribution shape
      virtual double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // variables
      // ======================================================================
      inline double k          () const  { return m_gamma.k     () ; }
      inline double theta      () const  { return m_gamma.theta () ; }
      inline double shift      () const  { return m_gamma.shift () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline double gamma_mean       () const { return m_gamma.mean       () ; }
      inline double gamma_dispersion () const { return m_gamma.dispersion () ; }
      inline double gamma_rms        () const { return m_gamma.rms        () ; }
      inline double gamma_skewness   () const { return m_gamma.skewness   () ; }
      inline double gamma_kurtosis   () const { return m_gamma.kurtosis   () ; }
      // ======================================================================
    public:
      // ======================================================================
      /** effective \f$\chi^2\f$-parameters
       *  If   \f$ Q  \sim \chi^2(\nu)\f$  and c is a positive constant,
       *  than \f$ cQ \sim \Gamma (k = \nu/2, \theta = 2c) \f$
       */
      inline double nu () const { return m_gamma.nu () ; }
      inline double c  () const { return m_gamma.c  () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying gamma distribution
      inline const GammaDist& gamma() const { return m_gamma ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setK     ( const double value  ) { return m_gamma.setK     ( value ) ; }
      bool setTheta ( const double value  ) { return m_gamma.setTheta ( value ) ; }
      bool setShift ( const double value  ) { return m_gamma.setShift ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      virtual double integral
      ( const double low  ,
	const double high ) const ;
      // ======================================================================
    public: // quantiles
      // ======================================================================
      /// calculate the quantile   (0<p<1)
      double quantile ( const double p ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    protected:
      // ======================================================================
      /// helper gamma distribution
      GammaDist m_gamma ;  // helper gamma distribution
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Log10GammaDist
     *  Distribution for log10(x) where x has gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  Log10GammaDist : public LogGammaDist
    {
    public:
      // ======================================================================
      /** constructor from scale & shape parameters
       *  param k      \f$k\f$ parameter (shape)
       *  param theta  \f$\theta\f$ parameter (scale)
       */
      Log10GammaDist 
      ( const double k     = 2 ,   // shape parameter
        const double theta = 1 ,   // scale parameter
        const double shift = 0 ) ; // shift parameter
      /// destructor
      virtual ~Log10GammaDist() ;  // destructor
      // ======================================================================
    public:
      // ======================================================================
      /// calculate log-gamma distribution shape
      double operator() ( const double x    ) const  override;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral
      ( const double low  ,
	const double high ) const override;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    public: // quantiles
      // ======================================================================
      /// calculate the quantile   (0<p<1)
      double quantile ( const double p ) const ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenGammaDist
     *  Generalized Gamma-distribution with additional shift parameter
     *  http://en.wikipedia.org/wiki/Generalized_gamma_distribution
     *  special cases :
     *   - p == 1      : Gamma  distribution
     *   - p == k      : Weibull distribution
     *   - p == k == 1 : Exponential distribution
     *   - p == k == 2 : Rayleigh    distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  GenGammaDist
    {
    public:
      // ======================================================================
      /** constructor
       *  param k     \f$k\f$ parameter      (shape)
       *  param theta \f$\theta\f$ parameter (scale)
       *  param p     \f$p\f$ parameter
       *  param low   bias
       */
      GenGammaDist
      ( const double k     = 2 ,
        const double theta = 1 ,
        const double p     = 1 , // 1 corresponds to gamma distribution
        const double low   = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate gamma distribution shape
      double pdf        ( const double x ) const ;
      /// calculate gamma distribution shape
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      inline double evaluate   ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // getters
      // ======================================================================
      inline double k          () const  { return m_k                          ; }
      inline double theta      () const  { return m_theta                      ; }
      inline double p          () const  { return m_p                          ; }
      inline double low        () const  { return m_low                        ; }
      // ======================================================================
    public:
      // ======================================================================
      /// Wikipedia notations
      inline double a          () const { return theta () ; }
      inline double d          () const { return k     () ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      inline double mean       () const  { return m_k * m_theta +   low ()     ; }
      inline double dispersion () const  { return m_k * m_theta * m_theta      ; }
      inline double variance   () const  { return dispersion ()                ; }
      inline double sigma      () const  { return std::sqrt ( dispersion ()  ) ; }
      inline double rms        () const  { return std::sqrt ( dispersion ()  ) ; }
      inline double skewness   () const  { return 2.0 / std::sqrt ( m_k )      ; }
      // ======================================================================
      // get the mode 
      double        mode   () const ; 
      // =====================================================================
      inline double xmin () const { return low () ; }
      // ======================================================================
    public:  // setters
      // ======================================================================
      bool setK     ( const double value  ) ;
      bool setTheta ( const double value  ) ;
      bool setP     ( const double value  ) ;
      bool setLow   ( const double value  ) ;
      // ======================================================================
    public:
      // ======================================================================
      double cdf      ( const double x ) const ;
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// shape
      double m_k      ; // shape
      /// scale
      double m_theta  ; // scale
      /// parameter
      double m_p      ; // parameter
      /// shift
      double m_low    ; // shift
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Amoroso
     *  Another view on generalized gamma distribtion
     *  http://arxiv.org/pdf/1005.3274
     *  @see Ostap::Math::GenGammaDist
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  Amoroso
    {
    public:
      // ======================================================================
      /** constructor
       *  param theta \f$\theta\f$-parameter
       *  param alpha \f$\alpha\f$-parameter (>0)
       *  param beta  \f$\beta\f$-parameter
       *  param a     a-parameter
       *  Note that   \f$\alpha\beta\f$ is equal to k-parameter
       */
      Amoroso
      ( const double theta = 1 ,
        const double alpha = 1 ,
        const double beta  = 1 ,
        const double a     = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Amoroso distribtion
      double pdf               ( const double x ) const ;
      /// evaluate Amoroso distribtion
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:  // direct getters
      // ======================================================================
      double a     () const { return m_a     ; }
      double theta () const { return m_theta ; }
      double alpha () const { return m_alpha ; }
      double beta  () const { return m_beta  ; }
      // ======================================================================
    public:  // derived getters
      // ======================================================================
      double d      () const { return alpha () * beta () ; }
      double k      () const { return alpha () * beta () ; }
      double p      () const { return            beta () ; }
      // ======================================================================
    public:  // helper getters
      // ======================================================================
      double theta2 () const { return m_theta * m_theta  ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool setA      ( const double value ) ;
      bool setTheta  ( const double value ) ;
      bool setAlpha  ( const double value ) ;
      bool setBeta   ( const double value ) ;
      bool setP      ( const double value ) { return setBeta ( value ) ; }
      // ======================================================================
    public: // general properties
      // ======================================================================
      double mode       () const ;
      double mean       () const ;
      double variance   () const ;
      double dispersion () const { return variance () ; }
      double sigma2     () const { return variance () ; }
      double sigma      () const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral () const ;
      double cdf      ( const double x    ) const ;
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_a     ;
      double m_theta ;
      double m_alpha ;
      double m_beta  ;
      // ======================================================================
    };
    // ========================================================================
    /** @class LogGamma
     *  - http://arxiv.org/pdf/1005.3274
     *  - Prentice, R. L. (1974). A log gamma model and its maximum likelikhood
     *                            estimation. Biometrika 61, 539
     *  - Johnson, N. L., Kotz, S., and Balakrishnan, N. (1995). Continuous
     *            univariate distributions, 2nd ed. Vol. 2. Wiley, New York.
     *  - Bartlett, M. S. and G., K. M. (1946). The statistical analysis of
     *                  variance-heterogeneity and the logarithmic transformation.
     *                 J. Roy. Statist. Soc. Suppl. 8, 1, 128.
     *
     *  @attention Do not confuse with Ostap::Math::LogGammaDist
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  LogGamma
    {
    public:
      // ======================================================================
      /** constructor from scale & shape parameters
       *  param nu      \f$\nu\f$ parameter      (location)
       *  param lambda  \f$\lambda\f$ parameter
       *  param alpha   \f$\alpha\f$ parameter    (>0)
       */
      LogGamma
      ( const double nu     = 0 ,   // shape parameter
        const double lambda = 1 ,   // scale parameter
        const double alpha  = 1 ) ; // scale parameter
      // ======================================================================
    public:
      // ======================================================================
      /// calculate log-gamma shape
      double pdf        ( const double x ) const ;
      /// calculate log-gamma shape
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getter s
      // ======================================================================
      inline double nu         () const  { return m_nu     ; }
      inline double lambda     () const  { return m_lambda ; }
      inline double lambd      () const  { return m_lambda ; }
      inline double alpha      () const  { return m_alpha  ; }
      // ======================================================================
    public: // general properties
      // ======================================================================
      double mean       () const ;
      double mode       () const ;
      double variance   () const ;
      double dispersion () const { return variance () ; }
      double sigma2     () const { return variance () ; }
      double sigma      () const ;
      double skewness   () const ;
      double kurtosis   () const ;
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setNu     ( const double value  ) ;
      bool setLambda ( const double value  ) ;
      bool setAlpha  ( const double value  ) ;
      // ======================================================================
    public:
      // ======================================================================
      double cdf ( const double x ) const ;
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double  m_nu      ;
      double  m_lambda  ;
      double  m_alpha   ;
      // ======================================================================
    } ;    
    // ========================================================================
    /** @class Beta
     *  http://en.wikipedia.org/wiki/Beta_prime_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2025-08-10
     */
    class  Beta
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param alpha \f$\alpha\f$-parameter
       *  @param beta  \f$\beta\f$-parameter
       *  @param scale  scale-parameter
       *  @param shift  shift-parameter
       */
      Beta
      ( const double alpha = 1 ,
        const double beta  = 1 ,
        const double scale = 1 ,
        const double shift = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate beta-distributions
      double        evaluate  ( const double x ) const ;
      /// evaluate beta'-distributions
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate beta'-distributions
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double alpha () const { return m_alpha ; }
      double beta  () const { return m_beta  ; }
      double scale () const { return m_scale ; }
      double shift () const { return m_shift ; }
      // ======================================================================
    public: // general properties
      // ======================================================================
      /// mean 
      double        mean       () const ;
      ///  mode 
      double        mode       () const ;
      ///  median 
      double        median     () const ;
      /// Variance 
      double        variance   () const ;
      /// RMS == sqrt (variance)
      double        rms        () const ;
      /// skewness 
      double        skewness   () const ;
      /// (excess kurtosis)
      double        kurtosis   () const ;
      /// dispersion == variance
      inline double dispersion () const { return variance () ; }
      /// sigma == RMS 
      inline double sigma      () const { return rms      () ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setAlpha ( const double value ) ;
      bool   setBeta  ( const double value ) ;
      bool   setScale ( const double value ) ;
      bool   setShift ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ()                    const ;
      double cdf      ( const double x    ) const ;
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public: // internal <--> external tsranformation 
      // ======================================================================
      /// internal t to external x 
      inline double x ( const double t ) const { return m_shift + t * m_scale     ; }
      /// externa x to internal t 
      inline double t ( const double x ) const { return ( x - m_shift ) / m_scale ; }
      // ======================================================================
    public: // support 
      // ======================================================================
      /// minimal x 
      inline double xmin () const { return m_shift           ; } 
      /// maximal x 
      inline double xmax () const { return m_shift + m_scale ; } 
      // ======================================================================
    public: // quantile 
      // ======================================================================
      /// quantile \f$ 0 \le p \le1 \f$
      double quantile ( const double p  ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_alpha { 1 } ;
      double m_beta  { 5 } ;
      double m_scale { 1 } ;
      double m_shift { 0 } ;
      // ======================================================================
    } ;    
    // ========================================================================
    /** @class BetaPrime
     *  http://en.wikipedia.org/wiki/Beta_prime_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  BetaPrime
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param alpha \f$\alpha\f$-parameter
       *  @param beta  \f$\beta\f$-parameter
       *  @param scale  scale-parameter
       *  @param shift  shift-parameter
       */
      BetaPrime
      ( const double alpha = 1 ,
        const double beta  = 5 ,
        const double scale = 1 ,
        const double shift = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate beta'-distributions
      double pdf        ( const double x ) const ;
      /// evaluate beta'-distributions
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double alpha () const { return m_alpha ; }
      double beta  () const { return m_beta  ; }
      double scale () const { return m_scale ; }
      double shift () const { return m_shift ; }
      // ======================================================================
    public: // general properties
      // ======================================================================
      double mean       () const ;
      double mode       () const ;
      double variance   () const ;
      double dispersion () const { return variance () ; }
      double sigma2     () const { return variance () ; }
      double sigma      () const ;
      double skewness   () const ;
      double kurtosis   () const ;
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setAlpha ( const double value ) ;
      bool   setBeta  ( const double value ) ;
      bool   setScale ( const double value ) ;
      bool   setShift ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ()                    const ;
      double cdf      ( const double x    ) const ;
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_alpha { 1 } ;
      double m_beta  { 5 } ;
      double m_scale { 1 } ;
      double m_shift { 0 } ;
      // ======================================================================
    private:
      // ======================================================================
      /// auxillary intermediate parameter
      double m_aux { 1 } ;                  // auxillary intermediate parameter
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FDistribution
     *  Speficic re-parameterization of Beta-prime distribution
     *  @see https://en.wikipedia.org/wiki/F-distribution
     *  http://en.wikipedia.org/wiki/Beta_prime_distribution
     */
    class FDistribution
    {
      // ======================================================================
    public :
      // ======================================================================
      /** constructor from all parameters 
       *  @param d1  d1-parameter d1>0
       *  @param d2  d2-parameter d2>0
       *  @param scale  (additional) scale-parameter
       *  @param shift  (additional) shift-parameter
       */
      FDistribution
      ( const double d1    = 1 ,
	const double d2    = 1 ,
        const double scale = 1 ,
        const double shift = 0 ) ;
      // =====================================================================
    public:
      // =====================================================================
      /// evaluate F-distribution
      double pdf        ( const double x ) const ;
      /// evaluate F-distribution
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get parameter d1 
      inline double d1    () const { return 2 * m_betap.alpha () ; } 
      /// get parameter d2 
      inline double d2    () const { return 2 * m_betap.beta  () ; }
      /// scale 
      inline double scale () const { return m_scale ; }
      /// shift 
      inline double shift () const { return m_shift ; }
      // ======================================================================
    public: // general properties
      // ======================================================================
      double mean       () const ;
      double mode       () const ;
      double variance   () const ;
      double dispersion () const { return variance () ; }
      double sigma2     () const { return variance () ; }
      double sigma      () const ;
      // ======================================================================
      inline double skewness   () const { return m_betap.skewness () ; } 
      inline double kurtosis   () const { return m_betap.kurtosis () ; } 
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set parameter d1 
      inline bool setD1    ( const double value ) { return m_betap.setAlpha ( 0.5 * value ) ; } 
      /// set parameter d1 
      inline bool setD2    ( const double value ) { return m_betap.setBeta  ( 0.5 * value ) ; }
      /// global scale 
      bool        setScale ( const double value ) ;
      /// global shift 
      bool        setShift ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ()                    const ;
      double cdf      ( const double x    ) const ;
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// beta' distribution
      Ostap::Math::BetaPrime m_betap { 1 , 1 } ;
      /// global scale 
      double                 m_scale { 1 } ;
      /// global shift 
      double                 m_shift { 0 } ;
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class GenBetaPrime
     *  Generalized BetaPrime distribution 
     *  http://en.wikipedia.org/wiki/Beta_prime_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class GenBetaPrime
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param alpha \f$\alpha\f$-parameter
       *  @param beta  \f$\beta\f$-parameter
       *  @param scale  scale-parameter
       *  @param shift  shift-parameter
       */
      GenBetaPrime
      ( const double alpha = 1 ,
        const double beta  = 5 ,
	const double p     = 1 ,
	const double q     = 1 ,	
        const double scale = 1 ,
        const double shift = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate beta'-distributions
      double pdf        ( const double x ) const ;
      /// evaluate beta'-distributions
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double alpha      () const { return m_alpha ; }
      double beta       () const { return m_beta  ; }
      double p          () const { return m_p     ; }
      double q          () const { return m_q     ; }
      double scale      () const { return m_scale ; }
      double shift      () const { return m_shift ; }
      // ======================================================================
    public: // general properties
      // ======================================================================
      double mean       () const ;
      double mode       () const ;
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setAlpha  ( const double value ) ;
      bool   setBeta   ( const double value ) ;
      bool   setP      ( const double value ) ;
      bool   setQ      ( const double value ) ;
      bool   setScale  ( const double value ) ;
      bool   setShift  ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ()                    const ;
      double cdf      ( const double x    ) const ;
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_alpha { 1 } ;
      double m_beta  { 5 } ;
      double m_p     { 1 } ;
      double m_q     { 1 } ;
      double m_scale { 1 } ;
      double m_shift { 0 } ;
      // ======================================================================
    private:
      // ======================================================================
      /// auxillary intermediate parameter
      double                 m_aux { 1 }     ; // auxillary intermediate parameter
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================
    } ;    
    // ========================================================================
    /** @class GenBeta
     *  Generalized Beta distribution (the most general form)
     *  - \f$ c \equiv \sin^2 \frac{\pi\gamma}{2} \f$ to ensure \f$ 0 \le c \le 1 \f$ 
     * @see https://en.wikipedia.org/wiki/Generalized_beta_distribution
     */
    class GenBeta
    {
    public :
      // =======================================================================
      /// the most general form
      GenBeta
      ( const double a     = 1 , // shape 
	const double b     = 1 , // scale 
	const double gamma = 0 , // c = sin^2 gamma
	const double p     = 1 , // shape 
	const double q     = 1 , // shape 
	const double shift = 0 ) ;
      // =======================================================================
    public :
      // =======================================================================
      /// evaluate generalized Beta distribution 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate generalized Beta distribution 
      double        evaluate   ( const double x ) const ;
      // =======================================================================
    public :
      // =======================================================================
      inline double a            () const { return m_a     ; } 
      inline double b            () const { return m_b     ; } 
      inline double c            () const { return m_c     ; } 
      inline double gamma        () const { return m_gamma ; } 
      inline double p            () const { return m_p     ; } 
      inline double q            () const { return m_q     ; } 
      inline double shift        () const { return m_shift ; } 
      inline double scale        () const { return b ()    ; }
      // =======================================================================
    public:
      // =======================================================================
      /// finite range ? (beta-like or beta'-like?)
      inline bool   finite_range () const { return 0 < m_a && !m_c1 ; }
      /// minimal x 
      double        xmin         () const ;
      /** maximal x
       *  @attention it can be infinite!
       */
      double        xmax         () const ;
      // =======================================================================
    public :
      // =======================================================================
      bool setA      ( const double value  ) ;
      bool setB      ( const double value  ) ;
      bool setGamma  ( const double value  ) ;
      bool setP      ( const double value  ) ;
      bool setQ      ( const double value  ) ;
      bool setShift  ( const double value  ) ;
      bool setAB
      ( const double valuep ,
	const double valueq ) ;		       
      bool setPQ
      ( const double valuea ,
	const double valueb ) ;		       
      // =======================================================================
    public :
      // =======================================================================
      /// get the integral 
      double integral () const ;
      /// get the integral 
      double integral
      ( const double low  ,
	const double high ) const ; 
      /// get the CDF 
      double cdf 
      ( const double x    ) const ; 
      // =======================================================================
    public :
      // ======================================================================
      /// unique tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      ///  parameter a 
      double m_a         {  1 } ; // parameter a
      ///  parameter b 
      double m_b         {  1 } ; // parameter b
      ///  parameter c
      double m_c         { -1 } ; // parameter c
      ///  parameter p
      double m_p         {  1 } ; // parameter p
      ///  parameter q
      double m_q         {  1 } ; // parameter q
      ///  shift/location parameter 
      double m_shift     {  0 } ; // shift/location parameter q
      /// parameter gamma
      double m_gamma     {  0 } ; // parameter gamma
      // ======================================================================
    private: 
      // ======================================================================
      ///  \f$ \log \Beta Beta ( p , q ) 
      double m_logBpq    {  0    } ; //  \f$ \log \Beta Beta ( p , q )
      /// \f$ \log b \f$
      double m_logb      {  0    } ;
      /// \f$ c = 1 ?\f$
      bool   m_c1        { false } ;
      /// \f$ \f$
      double m_ac        { -1    } ;  // (1-c)**(1/a) 
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================      
    } ;    
    // ========================================================================
    /** @class Landau
     *  http://en.wikipedia.org/wiki/Landau_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  Landau 
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param scale scale-parameter
       *  @param shift shift-parameter
       */
      Landau
      ( const double scale = 1 ,
        const double shift = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate beta'-distributions
      double        pdf        ( const double x ) const ;
      /// evaluate beta'-distributions
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      inline double scale () const { return m_scale ; }
      inline double shift () const { return m_shift ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setScale ( const double value ) ;
      bool   setShift ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double cdf
      ( const double x    ) const ;
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_scale ;
      double m_shift ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Weibull 
     *  3-parameter  Weibull distribution 
     *  \f$ f(x,\lambda,k,x_0) = \frac{k}{\lambda}  y^{k-1} e^{-y^k}\f$, where 
     *  \f$ y \equiv \frac{x-x_0}{\lambda}\f$
     *  @see https://en.wikipedia.org/wiki/Weibull_distribution
     *  @@date 2018-02-26
     */
    class Weibull
    {
    public:
      // ======================================================================
      /** constructor from all parameters    
       *  @param scale the scale parameter "lambda"  >0 
       *  @param shape the shape parameter "k"       >0
       *  @param shift the shift parameter "x0"
       */
      Weibull
      ( const double scale = 1 ,  
        const double shape = 1 ,  
        const double shift = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Weibull-distributions
      double         pdf        ( const double x ) const ;
      /// evaluate Weibull-distributions
      inline  double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: //  direct & derived getters 
      // ======================================================================
      inline double scale () const { return m_scale ; }
      inline double shape () const { return m_shape ; }
      inline double shift () const { return m_shift ; }
      // ======================================================================
      inline double lambd () const { return scale () ; }
      inline double k     () const { return shape () ; }
      inline double x0    () const { return shift () ; }
      inline double xmin  () const { return x0    () ; }      
      // ======================================================================
    public: //  setters  
      // ======================================================================
      bool setScale  ( const double value ) ;
      bool setShape  ( const double value ) ;
      bool setShift  ( const double value ) ;
      // ======================================================================
      inline bool setLambda ( const double value ) { return setScale ( value ) ; }
      inline bool setK      ( const double value ) { return setShape ( value ) ; }
      inline bool setX0     ( const double value ) { return setShift ( value ) ; }
      // ======================================================================
    public: // various derived  quantities 
      // ======================================================================
      /// the mean value 
      double mean     () const ;
      /// the mode 
      double mode     () const ;
      /// the median
      double median   () const ;
      /// variance
      double variance () const ;
      /// rms 
      double rms      () const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double cdf      ( const double x    ) const ;
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      ///  scale 
      double m_scale  ; //  scale 
      ///    shape  
      double m_shape  ; // shape
      ///  shift 
      double m_shift  ; // shift 
      // ======================================================================
    };    
    // ========================================================================
    /** @class Rice
     *  Rice distribution 
     *  \f$ f(x; \nu , \varsigma) = 
     *  \frac{\delta x}{\varsigma^2} \mathrm{e}^{-\frac{ \delta x^2+\nu^2}{2\varsigma^2} } 
     *   I_0 (\frac{\delta x\nu}{\varsigma^2}) \f$, 
     *   where  \f$ \delta x = x - \x_0\f$ and 
     *    - \f$ x\ge x_0\f$  
     *    - \f$ \nu \ge 0 \f$ 
     *    - \f$ \varsigma \ge 0 \f$ 
     *  @see https://en.wikipedia.org/wiki/Rice_distribution
     */
    class Rice
    {
    public: 
      // ======================================================================
      /// contructor with all parameters
      Rice
      ( const double nu       = 0 , 
        const double varsigma = 1 ,
        const double shift    = 0 ) ;
      // ======================================================================
    public: // 
      // ======================================================================
      /// evaluate the function
      double        evaluate   ( const double x ) const ;
      /// evaluate the function
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate the function
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: // setters & getters 
      // ======================================================================
      /// parameter nu 
      inline double nu       () const { return m_nu       ; }
      /// parameter varsigma 
      inline double varsigma () const { return m_varsigma ; }
      // parameger shift 
      inline double shift    () const { return m_shift    ; }
      /// minimal x 
      inline double xmin     () const { return m_shift    ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setNu       ( const double value ) ;
      bool setVarsigma ( const double value ) ;
      bool setShift    ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral   () const ;
      /// get the integral between low and high
      double integral
      ( const double low  ,
	const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// mean value 
      double mean       () const ;
      /// variance 
      double variance   () const ;
      /// dispersion 
      double dispersion () const { return variance () ; }
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_nu       { 0 } ;
      double m_varsigma { 1 } ;
      double m_shift    { 0 } ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace {} ; // workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenInvGauss
     *  Generalised Inverse Gaussian distribution using
     *  \f$ (\theta,\eta) \f$ parameterisation  
     *  - |f$ \theta = \sqrt{ab}\$ 
     *  - |f$ \eta   = \sqrt{\frac{b}{a}}\$ 
     *  @see https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution
     */
    class GenInvGauss 
    {
      // ======================================================================
    public:
      // ======================================================================
      GenInvGauss 
      ( const double theta = 1 , 
        const double eta   = 1 , 
        const double p     = 1 , 
        const double shift = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      // evaluat the function 
      double evaluate ( const double x ) const ;
      // evaluate the function 
      inline double operator() ( const double x ) const 
      { return evaluate ( x ) ; }
      // evaluate the function 
      inline double pdf        ( const double x ) const 
      { return evaluate ( x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// parameter theta 
      inline double theta () const { return m_theta ; }
      /// parameter eta 
      inline double eta   () const { return m_eta   ; }  
      /// parameter p 
      inline double p     () const { return m_p     ; }  
      /// parameter shif
      inline double shift () const { return m_shift ; }  
      // ======================================================================
    public: // derived parameters
      // ======================================================================      
      // parameter a (derived) 
      double a () const ;
      // parameter b (derived) 
      double b () const ;
      /// minimal x 
      inline double xmin () const { return m_shift ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setTheta  ( const double value ) ;
      bool setEta    ( const double value ) ;
      bool setP      ( const double value ) ;
      bool setShift  ( const double value ) ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean value 
      double mean     () const ;
      /// mode  
      double mode     () const ;
      /// variance  
      double variance () const ;
      /// RMS 
      double rms      () const ;
      /// disparsion 
      inline double dispersion () const { return variance () ; }
      /// sigma 
      inline double sigma      () const { return rms      () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral   () const ;
      /// get the integral between low and high
      double integral  
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag() const ;
      // ======================================================================
    private:      
      // ======================================================================
      /// parameter theta 
      double m_theta      {  1 } ; // parameter theta 
      /// parameter eta 
      double m_eta        {  1 } ; // parameter eta 
      /// p-parameter
      double m_p          {  1 } ; // p-parameter
      /// shift parameter
      double m_shift      {  1 } ; // shift parameter
      /// cache 
      double m_iKp_theta  { -1 } ; // cache 
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace {} ; //
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class IrwinHall
     *  Irwin-Hall distribution 
     *  @see https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution
     */
    class IrwinHall 
    {
    public: 
      // ======================================================================
      /** constructor from n-parameter
       *  @param n  n-parameter (shape) 
       */
      IrwinHall ( const unsigned short n = 2 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double evaluate ( const double x )  const ;
      /// evaluate the function 
      inline double operator () ( const double x )  const 
      { return evaluate ( x ) ;  }
      /// evaluate the function 
      inline double pdf         ( const double x )  const 
      { return evaluate ( x ) ;  }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// n-parameter 
      inline unsigned short n    () const { return m_bspline.degree () + 1 ; } ;
      // ======================================================================
    public: // 
      // ======================================================================
      inline double         xmin () const { return m_bspline.xmin   () ; }
      inline double         xmax () const { return m_bspline.xmax   () ; }      
      // ======================================================================
    public: // some statistical properties 
      // ======================================================================
      /// get mode 
      inline double mode       () const { return 0.5 * n ()  ; }
      /// get mean 
      inline double mean       () const { return 0.5 * n ()  ; }
      /// get median
      inline double median     () const { return 0.5 * n ()  ; }
      /// get variance 
      inline double variance   () const { return n () / 12.0 ; }
      /// get doispersion 
      inline double dispersion () const { return n () / 12.0 ; }
      /// get rms 
      double rms        () const ; 
      /// get skewness 
      double skewness   () const { return 0 ; }
      /// get kurtosis 
      double kurtosis   () const ; 
      // ======================================================================
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual bspliine 
      Ostap::Math::BSpline m_bspline ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Bates 
     *  Bates distribution 
     *  @see https://en.wikipedia.org/wiki/Bates_distribution
     *  Essentially it is a scaled version of Irwin-Hall distribution 
     *  @see https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution
     *  @see Ostap::Math::IrwinHall 
     */
    class Bates
    {
    public:
      // ======================================================================
      /** constructor from n-parameter
       *  @param n  n-parameter (shape) 
       */
      Bates  ( const unsigned short n = 2 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double evaluate ( const double x )  const ;
      /// evaluate the function 
      inline double operator () ( const double x )  const 
      { return evaluate ( x ) ;  }
      /// evaluate the function 
      inline double pdf         ( const double x )  const 
      { return evaluate ( x ) ;  }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// n-parameter 
      inline unsigned short n    () const { return m_ih.n () ; } ;
      // ======================================================================
    public: // 
      // ======================================================================
      inline double         xmin () const { return 0 ; }
      inline double         xmax () const { return 1 ; }      
      // ======================================================================
    public: // some statistical properties 
      // ====================================================================== 
      /// get mode 
      inline double mode       () const { return 0.5 ; }
      /// get mean 
      inline double mean       () const { return 0.5 ; }
      /// get median 
      inline double median     () const { return 0.5 ; }
      /// get variance 
      double        variance   () const ; 
      /// get rms 
      double        rms        () const ;
      /// get kurtosis 
      double        kurtosis   () const ; 
      /// get dispersion 
      inline double dispersion () const { return variance () ; }
      /// get skewness 
      inline double skewness   () const { return 0           ; }
      // ======================================================================
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual Irwin-Hall distribution  
      Ostap::Math::IrwinHall m_ih ;
      // ======================================================================
    };
    // ========================================================================
    /** @class BatesShape 
     *  Modified Bates distribution such that it has mean of \f$ \mu \f$
     *  and rms of \f$ \sigma \f$, \f$ n>0\f$ is just a shape parameters 
     *  @see https://en.wikipedia.org/wiki/Bates_distribution
     *  Essentially it is a scaled version of Irwin-Hall distribution 
     *  @see https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution
     *  @see Ostap::Math::Bates  
     *  @see Ostap::Math::IrwinHall 
     */
    class BatesShape 
    {
    public:
      // ======================================================================
      /** constructor from n-parameter
       *  @param n     n-parameter (shape) 
       *  @param mu    mu-parameter (location)
       *  @param sigma sigma-parameter (scale/width)
       *  - \f$ \mu = \frac{1}{2} \f$ and \f$ \sigma = \frac{1}{\sqrt{12n} \f$ 
       *     give the original Bates shape 
       */
      BatesShape 
      ( const double         mu    = 0 , 
        const double         sigma = 1 , 
        const unsigned short n     = 2 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double evaluate ( const double x )  const ;
      /// evaluate the function 
      inline double operator () ( const double x )  const 
      { return evaluate ( x ) ;  }
      /// evaluate the function 
      inline double pdf         ( const double x )  const 
      { return evaluate ( x ) ;  }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// mu-parameter 
      inline double         mu    () const { return m_mu         ; }
      /// sigma-parameter 
      inline double         sigma () const { return m_sigma      ; }
      /// n-parameter 
      inline unsigned short n     () const { return m_bates.n () ; } 
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set parameter mu 
      bool setMu     ( const double value ) ;
      /// set parameter sigma 
      bool setSigma  ( const double value ) ;
      // ======================================================================
    public: // 
      // ======================================================================
      /// get x for given t 
      inline double x ( const double t ) const 
      { return m_mu + ( t - 0.5 ) * m_sigma * m_sq12n ; }
      /// get t for given x 
      inline double t ( const double x ) const  
      { return ( x - m_mu ) / ( m_sq12n * m_sigma ) + 0.5 ; }
      // ======================================================================
    public: // 
      // ======================================================================
      /// minimal x 
      double         xmin () const ; 
      /// maximal x 
      double         xmax () const ; 
      // ======================================================================
    public: // some statistical properties 
      // ====================================================================== 
      /// get mode 
      inline double mode       () const { return m_mu ; }
      /// get mean 
      inline double mean       () const { return m_mu ; }
      /// get median 
      inline double median     () const { return m_mu ; }
      /// get variance 
      inline double variance   () const { return m_sigma * m_sigma ; }
      /// get dispersion 
      inline double dispersion () const { return variance () ; }
      /// get rms 
      inline double rms        () const { return m_sigma     ; }
      /// get skewness 
      inline double skewness   () const { return 0           ; }
      /// get kurtosis 
      double kurtosis   () const ; 
      // ======================================================================
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// helper constant \f$ \sqrt{12n}\f$ 
      double             m_sq12n  ; // sqrt(12n)
      /// parameter mu 
      double             m_mu     ; // parameter mu 
      /// parameter sigma 
      double             m_sigma  ; // parameter sigma 
      /// the actual Bates distribution  
      Ostap::Math::Bates m_bates  ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenPareto
     *  Generalized Pareto Distribution
     *  @see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
     */
    class GenPareto
    {
      // ======================================================================
    public:
      // ======================================================================
      GenPareto
      ( const double mu    = 0 , 
        const double scale = 1 ,
        const double shape = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      // evaluate function 
      double evaluate   ( const double x ) const ;
      // evaluate function 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      inline double mu    () const { return m_mu    ; }
      inline double scale () const { return m_scale ; }
      inline double shape () const { return m_shape ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool        setMu    ( const double value ) ;
      bool        setScale ( const double value ) ;
      bool        setShape ( const double value ) ;
      inline bool setXi    ( const double value ) { return setShape ( value ) ; }
      // ====================================================================== 
    public:
      // ====================================================================== 
      inline double xi    () const { return m_shape ; }
      /// xmin 
      inline double xmin  () const { return m_mu    ; }
      /// xmax: can be infinite 
      double        xmax  () const ;
      // ======================================================================      
    public: //properties 
      // ====================================================================== 
      /// mean value, defined only for shape<1 
      double mean       () const ;
      /// median value,
      double median     () const ;
      /// mode  value,
      double mode       () const { return m_mu         ; }
      /// variance, defined only for shape<0.5 
      double variance   () const ;
      /// dispersion, defined only for shape<0.5 
      double dispersion () const { return variance ()  ; }
      /// rms,  defined only for shape<0.5 
      double rms        () const ;
      /// skewness  defined only for shape<1/3 
      double skewness   () const ;
      /// (ex)skewness  defined only for shape<1/4 
      double kurtosis   () const ;  
      // ======================================================================      
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      /// get cdf 
      double cdf ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mu    { 0 } ;
      double m_scale { 1 } ;
      double m_shape { 0 } ;
      // ======================================================================      
    } ;  
    // ========================================================================
    /** @class ExGenPareto
     *  reparameterized Expontiated Generalized Pareto Distribution
     *  @see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
     *  \f[ f(x;\mu,\sigma,\x0 = 
     *  \left\{ \begin{array}{ll}
     *   e^z\left( 1+ \xi e^z)^{-frac{1}{\xi}-1} & \text{for}~\xi\ne0 \ \
     *   e^{z-e^z} & \text{for}~\xi=0 \                                 \
     *  \end{array}\right. \f]
     *  - where \f$ z = \frac{x-mu}{\sigma}\f$ 
     */
    class ExGenPareto
    {
      // ======================================================================
    public:
      // ======================================================================
      ExGenPareto
      ( const double mu    = 0 , 
        const double scale = 1 ,
        const double shape = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      // evaluate function 
      double evaluate   ( const double x ) const ;
      // evaluate function 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      inline double mu    () const { return m_mu    ; }
      inline double scale () const { return m_scale ; }
      inline double shape () const { return m_shape ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu    ( const double value ) ;
      bool setScale ( const double value ) ;
      bool setShape ( const double value ) ;
      bool setXi    ( const double value ) { return setShape ( value ) ; }
      // ====================================================================== 
    public:
      // ====================================================================== 
      /// get xi (same as shape) 
      inline double xi  () const { return m_shape ; }
      /// mode of distribution 
      double mode       () const ;
      /// mean value
      double mean       () const ;
      /// variance 
      double variance   () const ;
      /// rms 
      double rms        () const ;
      /// dispersion
      inline double dispersion () const { return variance () ; }
      // ======================================================================      
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mu         { 0 } ;
      double m_scale      { 1 } ;
      double m_shape      { 0 } ;
      // ======================================================================
      /// value of log ( abs ( xi )  ) for xi!=0 
      double m_alog_xi    { 0 } ; // value of log(abs(xi)) for xi!=0
      /// bias for mean 
      double m_bias_mean  { 0 } ; // bias for for mean 
      /// bias for variance 
      double m_bias_var   { 0 } ; // bias for for variance
      // ======================================================================      
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ; // workspace
      // ======================================================================
    } ;
    // ========================================================================
    /* @class Benini
     * A modified version of Benini distribution 
     * @see https://en.wikipedia.org/wiki/Benini_distribution
     * Parameters 
     *  - shift  \f$ \mu  \f$ 
     *  - scale  \f$ s    \f$ 
     *  - shape parameters \f$ p_i \f$
     * 
     * Cumulative distribution function is defined for \f$x>\mu\f$ as 
     * \f[ F(x) = 1 - \mathrm{e}^{ - \sum_i \left| p_i \right| \Delta^i } \f]
     * where \f$ \Delta = \log \frac{x-\mu}{s} \f$ 
     *
     *  For standard Benini one has only linear and quadratic terms 
     *  and shift is zero \f$ \mu=0\f$ 
     */
    class Benini 
    {
      // ======================================================================
    public:
      // ======================================================================
      Benini
      ( const std::vector<double>& pars      , // shape parameters               
        const double               scale = 1 ,               // scale parameter 
        const double               shift = 0 ) ;             // shift parameter 
      // ======================================================================
      Benini
      ( const double               scale = 1 ,              // scale parameter 
        const double               shift = 0 ) ;            // shift parameter 
      // ======================================================================
      Benini
      ( const unsigned short       n         ,              // number of shape parameters 
        const double               scale     ,              // scale parameter 
        const double               shift     ) ;            // shift parameter 
      // ======================================================================
      /// two shape parameters: alpha and beta 
      Benini
      ( const double               alpha     , 
	const double               beta      ,	
	const double               scale     ,              // scale parameter 
        const double               shift     ) ;            // shift parameter
      // ======================================================================
      /// three shape parameters: alpha, beta & gamma 
      Benini
      ( const double               alpha     , 
	const double               beta      , 	
	const double               gamma     , 	
	const double               scale     ,              // scale parameter 
        const double               shift     ) ;            // shift parameter
      // ======================================================================
      /// four shape parameters: alpha, beta, gamma & delta 
      Benini
      ( const double               alpha     , 
	const double               beta      , 	
	const double               gamma     , 	
	const double               delta     , 	
	const double               scale     ,              // scale parameter 
        const double               shift     ) ;            // shift parameter
      // ======================================================================
    public:
      // ======================================================================
      // evaluate function 
      double        evaluate   ( const double x ) const ;
      // evaluate function 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      /// shift parameter 
      inline double shift  () const { return m_shift ; }
      /// scale parameter 
      inline double scale  () const { return m_scale ; }
      /// shape parameter
      inline double par    ( const unsigned short i ) const 
      { return i < m_pars.size() ? m_pars [ i ] : 0.0 ; }
      /// all shape parametgers 
      inline const std::vector<double>& pars ()       const { return m_pars ; }
      /// linear     term 
      inline double alpha  () const { return par ( 0 ) ; }
      /// quadrattic term 
      inline double beta   () const { return par ( 1 ) ; }
      /// cubic      term 
      inline double gamma  () const { return par ( 2 ) ; }
      /// quartic    term 
      inline double delta  () const { return par ( 3 ) ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setScale ( const double value ) ;
      bool setShift ( const double value ) ;
      bool _setPar  ( const unsigned short i , const double value ) ;
      bool setPar   ( const unsigned short i , const double value ) 
      { return i < m_pars.size() ? _setPar ( i , value ) : false ; }
      /// linear term  (if present) 
      bool setAlpha ( const double value ) { return _setPar ( 0 , value ) ; }
      /// quadratic term 
      bool setBeta  ( const double value ) { return _setPar ( 1 , value ) ; }
      /// cubic term    (i fpresent)
      bool setGamma ( const double value ) { return _setPar ( 2 , value ) ; }
      /// quartic term  (if present) 
      bool setDelta ( const double value ) { return _setPar ( 3 , value ) ; }
      /// set paramters
      bool setPars  ( const std::vector<double>& values ) ;
      // ====================================================================== 
    public:
      // ====================================================================== 
      /// minimal x 
      inline double xmin () const { return m_shift + m_scale ; }
      // ====================================================================== 
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      /// get the CDF 
      double cdf ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// vector of parameters 
      std::vector<double>  m_pars  ; // vector of parameters 
      /// vector of parameters 
      std::vector<double>  m_pars2 ; // auxillary vector of parameters 
      /// scale parameter 
      double m_scale { 1 } ; // scale parameter 
      /// shift parameter 
      double m_shift { 0 } ; // shift parameter 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GEV
     *  Generalized extreme value distribution 
     *  https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
     */
    class GEV
    {
      // ======================================================================
    public:
      // ======================================================================
      GEV
      ( const double mu    = 0 , 
        const double scale = 1 ,
        const double shape = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      // evaluate function 
      double evaluate   ( const double x ) const ;
      // evaluate function 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      inline double mu    () const { return m_mu    ; }
      inline double scale () const { return m_scale ; }
      inline double shape () const { return m_shape ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu    ( const double value ) ;
      bool setScale ( const double value ) ;
      bool setShape ( const double value ) ;
      bool setXi    ( const double value ) { return setShape ( value ) ; }
      // ====================================================================== 
    public:
      // ====================================================================== 
      /// get xi (same as shape) 
      inline double xi  () const { return m_shape ; }
      /// mode of distribution 
      double mode       () const ;
      /// median
      double median     () const ;
      // ======================================================================      
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      /// get the CDF 
      double cdf ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mu         { 0 } ;
      double m_scale      { 1 } ;
      double m_shape      { 0 } ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class MPERT
     *  Modified PERT distribution 
     *  @see https://en.wikipedia.org/wiki/PERT_distribution
     *  @see https://www.vosesoftware.com/riskwiki/ModifiedPERTdistribution.php
     */
    class MPERT 
    {
    public: 
      // ======================================================================
      /// constructor from two parameters
      MPERT
      ( const double xmin  = 0   , 
        const double xmax  = 1   ) ;
      /// constructor from four parameters
      MPERT
      ( const double xmin        , 
        const double xmax        , 
        const double xi          , // mode parameter 
        const double gamma = 4   ) ;
      // ======================================================================      
    public:
      // ======================================================================      
      /// get the value of MPERT distribution 
      double evaluate ( const double x ) const ;
      /// get the value of MPERT distribution 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// get the value of MPERT distribution 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================      
    public: // getters 
      // ======================================================================      
      /// xmin 
      inline double xmin  () const { return m_xmin  ; }
      /// xmax 
      inline double xmax  () const { return m_xmax  ; }
      /// xi - mode parameter 
      inline double xi    () const { return m_xi    ; }
      /// gamma/shape 
      inline double gamma () const { return m_gamma ; }
      // ======================================================================      
    public: // shape parameters 
      // ======================================================================      
      /// shape parameter alpha1 
      inline double alpha1 () const { return m_alpha1 ; }
      /// shape parameter alpha2 
      inline double alpha2 () const { return m_alpha2 ; }
      // ======================================================================      
    public: // setters 
      // ======================================================================      
      /// set xi parameter/mode 
      bool setXi    ( const double value ) ;
      /// set gamma 
      bool setGamma ( const double gamma ) ;
      // ======================================================================      
    public: // derived quantities 
      // ======================================================================      
      /// mode 
      inline double mode  () const { return m_xi  ; }
      /// mean value/mu  
      inline double mu    () const 
      { return ( m_xmin + m_gamma * m_xi + m_xmax ) / ( m_gamma + 2 ) ; }  
      /// mean value/mu 
      inline double mean () const { return mu  () ; }
      /// variance 
      double variance () const ;
      /// dispersion 
      inline double dispersion () const { return variance () ; }
      /// rms 
      double rms      () const ;
      /// skewness 
      double skewness () const ;
      /// (excess) kurtosis
      double kurtosis () const ;
      // ======================================================================      
    public: // integrals 
      // ======================================================================
      /// get the intergral between xmin and xmax 
      double integral () const ;
      /// get the integral 
      double integral
      ( const double xlow , 
        const double xhigh ) const ;
      /// get the CDF value 
      double cdf ( const double x ) const ;
      // ======================================================================      
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// low edge 
      double m_xmin   { 0   } ; // low edge 
      /// high edge 
      double m_xmax   { 1   } ; // high edge 
      /// xi mode parameter 
      double m_xi     { 0.5 } ; // xi - mode parameter 
      /// gamma/shape  
      double m_gamma  { 4   } ; // gamma/shape
      // ======================================================================
      /// shape parameter alpha1 
      double m_alpha1 { -1 } ; // shape parameter alpha2 
      /// shape parameter alpha1 
      double m_alpha2 { -1 } ; // shape parameter alpha2 
      /// normalization constant
      double m_N      { -1 } ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FisherZ
     *  Fisher's Z-distribution with additional location-scale parameters 
     *  @see https://en.wikipedia.org/wiki/Fisher%27s_z-distribution
     */
    class FisherZ 
    {
    public: 
      // ======================================================================
      /// constructor from parameters
      FisherZ
      ( const double mu     = 0 , 
	const double d1     = 5 , 
        const double d2     = 5 , 
	const double scale  = 1 ) ;
      // ======================================================================      
    public:
      // ======================================================================      
      /// get the value of Fisher-Z distribution 
      double evaluate ( const double x ) const ;
      /// get the value of Fisher-Z distribution 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// get the value of Fisher-Z distribution 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================      
    public: // getters 
      // ======================================================================      
      /// mu 
      inline double mu    () const { return m_mu    ; }
      /// scale 
      inline double scale () const { return m_scale ; }
      /// d1 - shape parameter 
      inline double d1    () const { return m_d1    ; }
      /// d2 - shape parameter 
      inline double d2    () const { return m_d2    ; }
      // ======================================================================      
    public: // setters 
      // ======================================================================      
      /// set mu     parameter/mode 
      bool setMu    ( const double value ) ;
      /// set scale  parameter 
      bool setScale ( const double value ) ;
      /// set d1    parameter
      bool setD1    ( const double value ) ;
      /// set d2    parameter
      bool setD2    ( const double value ) ;
      // ======================================================================      
    public: // integrals 
      // ======================================================================
      /// get the integral 
      double integral () const ;
      /// get the intergral between xmin and xmax 
      double integral
      ( const double xlow , 
        const double xhigh ) const ;
      // ======================================================================      
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// get normalization constant (withoh scale!)
      double norm () const ; // get normalization constant (without scale!)
      // ======================================================================
    private:
      // ======================================================================
      /// mu-parameter 
      double m_mu    { 0 } ; // mu-parameter
      /// scale-parameter 
      double m_scale { 1 } ; // scale-parameter
      /// d1/shape parameter 
      double m_d1    { 5 } ; // d1-shape parameter
      /// d2/shape parameter 
      double m_d2    { 5 } ; // d2-shape parameter
      // ======================================================================
    private:
      // ======================================================================
      /// normalization constant
      double m_C     { -1 } ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BirnbaumSaunders
     *  Birnbaum-Saunders distribution 
     *  @see https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution
     *  \f[ f(x;\mu, \beta,\gamma) = 
     *   \frac{ z + z^{-1}}{2\gamma(x-\mu)}\phi( \frac{1}{\gamma}(z-z^{-1}) \f]
     *  where
     *   - \f$ z=\sqrt{\frac{x-\mu}{\beta}}\f$
     *   - \f$ \phi\f$ is Gaussian PDF 
     */
    class BirnbaumSaunders
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param mu   location parameter 
       *  @param beta  scale parameter 
       *  @param gamma shape parameter 
       */
      BirnbaumSaunders
      ( const double mu    = 0 ,   // location 
	const double beta  = 1 ,   // scale
	const double gamma = 1 ) ; // shape 
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of B-S distribution 
      double evaluate ( const double x ) const ;
      /// get the value of B-S distribution 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// get the value of B-S distribution 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================      
    public: // getters 
      // ======================================================================
      /// get location parameter 
      inline double mu    () const { return m_mu ; } 
      /// get scale parameter 
      inline double beta  () const { return m_beta     ; } 
      /// get scale parameter 
      inline double scale () const { return m_beta     ; } 
      /// get shape paramteer 
      inline double gamma () const { return m_gamma    ; }
      /// get shape paramteer 
      inline double shape () const { return m_gamma    ; }
      /// get alpha parameter
      inline double alpha () const { return   gamma () ; } 
      /// xmin
      inline double xmin  () const { return   mu ()    ; }
      // ======================================================================      
    public: // getters 
      // ======================================================================
      double mean     () const ;
      double variance () const ;
      double rms      () const ;
      double skewness () const ;
      double kurtosis () const ;
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu    ( const double value ) ;
      bool setBeta  ( const double value ) ;
      bool setGamma ( const double value ) ;
      bool setScale ( const double value ) { return setBeta  ( value ) ; }
      bool setShape ( const double value ) { return setGamma ( value ) ; }
      bool setAlpha ( const double value ) { return setGamma ( value ) ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// get the integral 
      double integral () const ;
      /// get CDF
      double  cdf ( const double x ) const ; 
      /// get the intergral between xmin and xmax 
      double integral
      ( const double xlow , 
        const double xhigh ) const ;
      // ======================================================================      
    public:
      // ======================================================================
      /// get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// location parameter
      double m_mu    { 0 } ; // location parameter
      /// scale parameter
      double m_beta  { 1 } ; // scale parameter
      /// shape parameter
      double m_gamma { 1 } ; // shape  parameter
      // ======================================================================
    };
    // ========================================================================
    /** @class Freshet distribution
     *  @see https://en.wikipedia.org/wiki/Fr%C3%A9chet_distribution
     */
    class Frechet
    {
    public :
      // ======================================================================
      Frechet
      ( const double alpha = 1 ,   // shape
	const double scale = 1 ,   // scale 
	const double shift = 0 ) ; // shift 
      // ======================================================================
    public: 
      // ======================================================================
      /// evaluate Frechet distribution
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Frechet distribution
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Frechet distribution
      double evaluate ( const double x ) const ;
      // ======================================================================
    public :
      // ====================================================================== 
      /// shape parameter 
      inline double alpha () const { return m_alpha ; }
      /// scale parameter 
      inline double scale () const { return m_scale ; }
      /// shift parameter 
      inline double shift () const { return m_shift ; }
      // ======================================================================
      /// shape parameter 
      inline double shape () const { return alpha () ; }
      /// shift/bias 
      inline double xmin  () const { return shift () ; } 
      // ======================================================================
    public : 
      // ======================================================================
      bool         setAlpha ( const double value ) ;
      bool         setScale ( const double value ) ;
      bool         setShift ( const double value ) ;
      inline  bool setShape ( const double value ) { return setAlpha ( value ) ; }
      // ======================================================================
    public :
      // ======================================================================
      /// mean 
      double mean     () const ;
      /// median
      double median   () const ;
      /// mode 
      double mode     () const ;
      /// variance 
      double variance () const ;
      /// RMS 
      double rms      () const ;
      /// skewness 
      double skewness () const ;
      /// excess kurtotis 
      double kurtosis () const ;
      // ======================================================================
    public :
      // ======================================================================
      /** get quantile 
       *  @param p probability \f$ 0 \le p < 1 \f$
       */
      double quantile ( const double p ) const  ;
      /// integral 
      double integral () const ;
      /// integral 
      double integral
      ( const double a ,
	const double b ) const ;
      /// CDF
      double cdf
      ( const double x ) const ;
      // ======================================================================
    public :
      // ======================================================================
      /// uniquie tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// shape parameter alpha 
      double m_alpha { 1 } ; // shape parameter alpha 
      /// scale parameter  
      double m_scale { 1 } ; // scale parameter  
      /// shift parameter  
      double m_shift { 0 } ; // shift parameter  
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Dagum
     *  Dagum distribution (with bias parameter)
     *  @see https://en.wikipedia.org/wiki/Dagum_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     */
    class Dagum
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param p shape parameter \f$ 0 < p \f$
       *  @param a shape parameter \f$ 0 < a \f$
       *  @param b scale parameter \f$ 0 < b \f$
       *  @param shift shift parameter        
       */
      Dagum
      ( const double p     = 1 ,
	const double a     = 1 ,
	const double b     = 1 ,
	const double shift = 0 ) ;       	
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate Dagum distribution 
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Dagum distribution 
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Dagum distribution 
      double        evaluate   ( const double x ) const ;
      // ======================================================================
    public :
      // ======================================================================
      /// shape parameter p 
      inline double  p     () const { return m_p      ; } 
      /// shape parameter a 
      inline double  a     () const { return m_a      ; } 
      /// scale parameter b
      inline double  b     () const { return m_b      ; } 
      /// shift parameter 
      inline double  shift () const { return m_shift  ; } 
      /// shift parameter 
      inline double  xmin  () const { return shift () ; } 
      // ======================================================================
    public:
      // ======================================================================
      /// set shape parameter p 
      bool setP     ( const double value ) ;
      /// set shape parameter a
      bool setA     ( const double value ) ;
      /// set scale parameter b
      bool setB     ( const double value ) ;
      /// set shift parameter 
      bool setShift ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// mean value
      double mean     () const ;
      /// mode value
      double mode     () const ;
      /// median value
      double median   () const ;
      /// variance value
      double variance () const ;
      /// RMS value
      double rms      () const ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral () const ;
      /// get the integral 
      double integral
      ( const double low  ,
	const double high ) const ;
      /// get CDF
      double cdf
      ( const double x    ) const ;
      /// get quantile \f$ 0 < p < 1 \f$
      double quantile
      ( const double u ) const ;
      // ======================================================================
    public :
      // ======================================================================
      /// uniquie tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// shape parameter p
      double m_p     { 1 } ; // shape parameter p 
      /// shape parameter a
      double m_a     { 1 } ; // shape parameter a
      /// scale parameter b
      double m_b     { 1 } ; // scale parameter a
      /// shift parameter 
      double m_shift { 0 } ; // shift parameter a
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class BenktanderI
     *  Variant of Benktander type 1 distribution
     *  @see https://en.wikipedia.org/wiki/Benktander_type_I_distribution
     *  - \f$ z = \frac{x-x_0}{\sigma} + 1 \f$ for \f$ 0\le x \f$
     *  - \f$ b = p \frac{a(a+1)}{2}\f$ , where \f$ 0 < p \le 1 \f$
     *  - \f$ p = 1 / \sqrt{ r^2 + 1 } \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     */
    class BenktanderI
    {      
      // ======================================================================
    public : 
      // ======================================================================
      /// constructor from all parameters
      BenktanderI
      ( const double a     = 1 ,
	const double r     = 0 , 
	const double scale = 1 ,
	const double shift = 0 ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate Benktander Type I distribution
      inline double operator() ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate Benktander Type I distribution
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Benktander Type I distribution
      double        evaluate   ( const double x ) const ;
      // ======================================================================
    public :
      // ======================================================================
      inline double a     () const { return m_a      ; } 
      inline double r     () const { return m_r      ; } 
      inline double scale () const { return m_scale  ; } 
      inline double shift () const { return m_shift  ; }
      // =====================================================================
      /// helper parameter "p"
      inline double p     () const { return m_p      ; }
      /// original parameter "b" 
      inline double b     () const { return m_p * m_a * ( m_a + 1 ) / 2 ; }
      /// minimal x 
      inline double xmin  () const { return shift () ; }
      // =====================================================================
    public :
      // =====================================================================
      bool setA     ( const double value ) ;
      bool setR     ( const double value ) ;
      bool setScale ( const double value ) ;
      bool setShift ( const double value ) ;
      // =====================================================================
    public: 
      // =====================================================================
      /// mean value 
      double mean     () const ;
      /// variance value 
      double variance () const ;
      /// RMS value 
      double rms      () const ;
      // =====================================================================
    public : 
      // =====================================================================
      /// integral 
      double integral () const ;
      /// integral 
      double integral
      ( const double low  ,
	const double high ) const ;
      /// CDF  
      double cdf 
      ( const double x    ) const ;
      // =====================================================================
    public :
      // ======================================================================
      /// unique tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// set the p from delta 
      void setP ( const double delta ) ;
      // ======================================================================
    private:
      // ======================================================================
      ///  parameter a 
      double m_a         {  1 } ; // parameter a
      ///  parameter r 
      double m_r         {  0 } ; // parameter r 
      /// scale parameter 
      double m_scale     {  1 }  ; // scale parameter 
      /// shift parameter 
      double m_shift     {  0 }  ; // shift parameter
      /// parameter p 
      double m_p         { -1 } ; // parameter p 
      // ======================================================================
    } ; 
    // ========================================================================
    /** @class BenktanderII
     *  Variant of Benktander type II distribution
     *  @see https://en.wikipedia.org/wiki/Benktander_type_II_distribution
     *  - \f$ z = \frac{x-x_0}{\sigma} + 1 \f$ for \f$ 0\le x \f$
     *  - \f$ b = 1 / \sqrt{ r^1 + 1 }  \f$ 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     */
    class BenktanderII
    {      
      // ======================================================================
    public : 
      // ======================================================================
      /// constructor from all parameters
      BenktanderII
      ( const double a     = 1 ,
	const double r     = 0 , 
	const double scale = 1 ,
	const double shift = 0 ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate Benktander Type II distribution
      inline double operator() ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate Benktander Type II distribution
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Benktander Type II distribution
      double        evaluate   ( const double x ) const ;
      // ======================================================================
    public :
      // ======================================================================
      inline double a     () const { return m_a     ; } 
      inline double r     () const { return m_r     ; } 
      inline double scale () const { return m_scale ; } 
      inline double shift () const { return m_shift ; }
      /// original parameter "b" 
      inline double b     () const { return m_b      ; }
      /// minimal x 
      inline double xmin  () const { return shift () ; }
      // =====================================================================
    public :
      // =====================================================================
      bool setA     ( const double value ) ;
      bool setR     ( const double value ) ;
      bool setScale ( const double value ) ;
      bool setShift ( const double value ) ;
      // =====================================================================
    public: 
      // =====================================================================
      /// mode 
      double mode     () const ;
      // =====================================================================
    public : 
      // =====================================================================
      /// integral 
      double integral () const ;
      /// integral 
      double integral
      ( const double low  ,
	const double high ) const ;
      /// CDF  
      double cdf 
      ( const double x    ) const ;
      // =====================================================================
    public :
      // ======================================================================
      /// unique tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      ///  parameter a 
      double m_a         {  1 } ; // parameter a
      ///  parameter r             
      double m_r         {  0 } ; // parameter r
      /// scale parameter 
      double m_scale     {  1 }  ; // scale parameter 
      /// shift parameter 
      double m_shift     {  0 }  ; // shift parameter
      ///  parameter b       
      double m_b         { -1 }  ; // parameter b
      // ======================================================================
    } ; 
    // ========================================================================
    /** @class LogNormal
     *  Log-normal distribution
     *  @see https://en.wikipedia.org/wiki/Log-normal_distribution
     *  - We add here an shift parameter
     *  - and use "mu = log(scale)"
     *  - and we use "sigma" as shape parameter
     */
    class LogNormal
    {
    public:
      // ======================================================================
      /// constructor 
      LogNormal
      ( const double shape = 1 ,
	const double scale = 1 ,
	const double shift = 0 ) ; 
      // ======================================================================
    public : 
      // ======================================================================
      /// evaluate log-normal function
      double        evaluate    ( const double x ) const ;
      /// evaluate log-normal function
      inline double operator () ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate log-normal function
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; } 
      // ======================================================================
    public : 
      // ======================================================================
      inline double shape () const { return m_shape ; }
      inline double scale () const { return m_scale ; }
      inline double shift () const { return m_shift ; }
      // ======================================================================
      inline double xmin  () const { return m_shift ; }
      // ======================================================================
    public :  // "canonical parameters"
      // ======================================================================
      /// canonical mu 
      double        canonical_mu    () const ;
      /// canonical sigma 
      inline double canonical_sigma () const { return m_shape ; } 
      // ======================================================================
    public:
      // =====================================================================
      bool setShape ( const double value ) ;
      bool setScale ( const double value ) ;
      bool setShift ( const double value ) ;
      // =====================================================================
    public : 
      // =====================================================================
      /// integral 
      double integral () const ;
      /// integral 
      double integral
      ( const double low  ,
	const double high ) const ;
      /// CDF  
      double cdf 
      ( const double x    ) const ;
      /// quantile  function \f$ 0 < p < 1 \f$ 
      double quantile 
      ( const double p    ) const ;
      // =====================================================================
    public :
      // ======================================================================
      /// unique tag
      std::size_t tag () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// mean value 
      double mean     () const ;
      /// the mode  
      double mode     () const ;
      /// median value 
      double median   () const ;
      /// variance 
      double variance () const ;
      /// rms  
      double rms      () const ;
      /// skewness       
      double skewness () const ; 
      /// (excess) kurtosis
      double kurtosis () const ;      
      // ======================================================================
    private :
      // ======================================================================
      /// shape parameter 
      double m_shape { 1 } ; // shape parameter 
      /// scale parameter 
      double m_scale { 1 } ; // scale parameter 
      // shift parameter
      double m_shift { 0 } ; // shift parameter      
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ExpoLog
     *  Exponential-logarithmic distribution
     *  @see https://en.wikipedia.org/wiki/Exponential-logarithmic_distribution
     *  - We have added a shift parameter
     *  - to ensure \f$ 0 < p < 1 \f$  we use
     *    \f$ p = \frac{1}{2}\left[ 1 + \tanh \psi \right] \f$ 
     */
    class ExpoLog
    {
      // ======================================================================
    public :
      // ======================================================================
      ExpoLog
      ( const double beta  = 1 ,   // scale 
	const double psi   = 0 ,   // related to p 
	const double shift = 0 ) ; // shift
      // ======================================================================
    public : 
      // ======================================================================
      /// evaluate expo-log function
      double        evaluate    ( const double x ) const ;
      /// evaluate expo-log function
      inline double operator () ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate expo-log function
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; } 
      // ======================================================================
    public : 
      // ======================================================================
      inline double beta  () const { return m_beta  ; }
      inline double psi   () const { return m_psi   ; }
      inline double shift () const { return m_shift ; }
      // ======================================================================
      /// original p-parameter 0 < p < 1 
      inline double p     () const { return m_p     ; } 
      inline double xmin  () const { return m_shift ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setBeta  ( const double value ) ;
      bool setPsi   ( const double value ) ;
      bool setShift ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// intergal 
      double integral ()                    const ;
      /// intergal 
      double integral
      ( const double low  ,
        const double high ) const ;
      /// get CDF 
      double cdf      ( const double x    ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// exponential parameter beta 
      double m_beta  { 1   } ; // exponential parameter beta 
      /// parameter psi 
      double m_psi   { 0   } ; // parameter psi 
      /// shift parameter
      double m_shift { 0   } ; // shift parameter
      /// parameter p
      double m_p     { 0.5 } ; // parameter p 
      /// parameter p
      double m_logp  { 100 } ; // log ( p )
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Davis
     * @see https://en.wikipedia.org/wiki/Davis_distribution
     */
    class Davis
    {
    public:
      // ======================================================================
      /// conastructor 
      Davis
      ( const double b  = 1 ,
	const double n  = 5 ,
	const double mu = 0 ) ;	
      // ======================================================================
    public : 
      // ======================================================================
      /// evaluate davis function
      double        evaluate    ( const double x ) const ;
      /// evaluate davis function
      inline double operator () ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate davis  function
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; } 
      // ======================================================================
    public : 
      // ======================================================================
      inline double b    () const { return m_b  ; }
      inline double n    () const { return m_n  ; }
      inline double mu   () const { return m_mu ; }
      inline double xmin () const { return m_mu ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setB  ( const double value ) ;
      bool setN  ( const double value ) ;
      bool setMu ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// mean-value 
      double mean     () const ;
      /// variance 
      double variance () const ;
      /// rms 
      double rms      () const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// intergal 
      double integral ()    const ;
      /// intergal 
      double integral
      ( const double low  ,
        const double high ) const ;
      /// CDF 
      double cdf 
      ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// parameter b 
      double m_b  { 1 } ; // parameter b
      /// parameter n 
      double m_n  { 5 } ; // parameter n
      /// parameter mu 
      double m_mu { 0 } ; // parameter mu
      // ======================================================================
    private:
      // ======================================================================
      /// zeta ( n     ) 
      double m_z0  { 0  } ; // zeta ( n     ) 
      /// zeta ( n - 1 ) 
      double m_z1  { 0  } ; // zeta ( n - 1 ) 
      /// zeta ( n - 2  ) 
      double m_z2  { 0  } ; // zeta ( n - 2 ) 
      /// normailzation
      double m_C   { -1 } ; // normalzation 
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;      
    // ========================================================================
    /** @class Kumaraswami
     *  Kumaraswami distribution with scale and shift
     *  @see https://en.wikipedia.org/wiki/Kumaraswamy_distribution
     */
    class Kumaraswami
    {
      // ======================================================================
    public :
      // ======================================================================
      Kumaraswami
      ( const double a     = 1 , // 0<a 
	const double b     = 1 , // 0<b 
	const double scale = 1 , // 0<scale
	const double shift = 0 ) ;      	
      // ======================================================================
    public : 
      // ======================================================================
      /// evaluate Kumaraswami function
      double        evaluate    ( const double x ) const ;
      /// evaluate Kumaraswami function
      inline double operator () ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate Kumaraswami function
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; } 
      // ======================================================================
    public : 
      // ======================================================================
      inline double a      () const { return m_a     ; }
      inline double b      () const { return m_b     ; }
      inline double scale  () const { return m_scale ; }
      inline double shift  () const { return m_shift ; }
      // ======================================================================
    public :
      // ======================================================================
      inline double xmin  () const { return m_shift           ; }
      inline double xmax  () const { return m_shift + m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setA     ( const double value ) ;
      bool setB     ( const double value ) ;
      bool setScale ( const double value ) ;
      bool setShift ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// mean-value 
      double mean     () const ;
      /// median-value 
      double median   () const ;
      /// mode-value 
      double mode     () const ;
      /// variance 
      double variance () const ;
      /// rms 
      double rms      () const ;
      /// skewness
      double skewness () const ;
      /// (excess) kurtosis
      double kurtosis () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// raw moment for standartized Kumaraswami
      double moment ( const unsigned int n ) const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// intergal 
      double integral ()    const ;
      /// intergal 
      double integral
      ( const double low  ,
        const double high ) const ;
      /// CDF 
      double cdf 
      ( const double x ) const ;
      /// quantile  function \f$ 0 < p < 1 \f$ 
      double quantile 
      ( const double p    ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// parameter a 
      double m_a      { 1 } ; // parameter a 
      /// parameter b
      double m_b      { 1 } ; // parameter b 
      /// scale parameter
      double m_scale  { 1 } ; // scale parameter 
      /// shift parameter
      double m_shift  { 0 } ; // shift  parameter 
      // ======================================================================
    };
    // ========================================================================
    /** @class InverseGamma
     *  Inverse Gamma distribution (with shift)
     *  @see https://en.wikipedia.org/wiki/Inverse-gamma_distribution
     */
    class InverseGamma
    {
      // ======================================================================
    public:
      // ======================================================================
      /// full constructot 
      InverseGamma
      ( const double alpha = 8 ,
	const double beta  = 1 ,
	const double shift = 0 ) ;	      
      // ======================================================================
    public : 
      // ======================================================================
      /// evaluate Inverse-Gamma function
      double        evaluate    ( const double x ) const ;
      /// evaluate Inverse-Gamma function
      inline double operator () ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate Inverse-Gamma function
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; } 
      // ======================================================================
    public : 
      // ======================================================================
      inline double alpha  () const { return m_alpha ; }
      inline double beta   () const { return m_beta  ; }
      inline double shift  () const { return m_shift ; }
      // ======================================================================
    public :
      // ======================================================================
      inline double xmin  () const { return m_shift ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setAlpha ( const double value ) ;
      bool setBeta  ( const double value ) ;
      bool setShift ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// mean
      double mean     () const ;
      /// variance 
      double variance () const ;
      /// RMS 
      double rms      () const ;
      /// skewness 
      double skewness () const ;
      /// (excess) kurtosis
      double kurtosis () const ;      
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// intergal 
      double integral ()    const ;
      /// intergal 
      double integral
      ( const double low  ,
        const double high ) const ;
      /// CDF 
      double cdf 
      ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// parameter alpha (shape) 
      double m_alpha { 8 } ; // parameter alpha
      /// parameter beta  (scale) 
      double m_beta  { 1 } ; // parameter b 
      /// shift parameter
      double m_shift { 0 } ; // shift  parameter 
      // ======================================================================
    };
    // ========================================================================    
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_MODELS_H
// ============================================================================
