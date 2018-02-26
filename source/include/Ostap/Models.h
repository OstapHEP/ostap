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
// OStap
// ============================================================================
#include "Ostap/NSphere.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein1D.h"
#include "Ostap/Workspace.h"
#include "Ostap/PhaseSpace.h"
// ============================================================================
/** @file Ostap/Models.h
 *
 *  set of useful models
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
     *  If  x is distributed accrofing  to \f$ f(x) \propto e^{-\tau x} \f$, 
     *  then z, \f$ z  =   log(x) \f$, is distributed accoring to 
     *  \f$ F(z) = G(x, -\log(\tau), 1 ) \f$ 
     * 
     *  As a   result, if    x is distributes as sum of exponentilcomponents 
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
      Gumbel ( const double mu   = 0 , 
               const double beta = 1 );
      // ======================================================================
    public: // primary getters 
      // ======================================================================
      double mu       () const { return m_mu      ; }
      double peak     () const { return   mu   () ; }
      double location () const { return   mu   () ; }
      double beta     () const { return m_beta    ; }
      double scale    () const { return   beta () ; }
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
      double mode        () const { return mu() ; }
      double variance    () const ;
      double dispersion  () const { return variance () ; }
      double sigma2      () const { return variance () ; }
      double sigma       () const ;
      double skewness    () const ;
      double kurtosis    () const { return  12.0 / 5 ; }
      // ======================================================================
    public: // the main block 
      // ======================================================================
      /// get a value for the function 
      double operator() ( const double x )  const { return pdf ( x ) ; }
      /// get a value for the function      
      double pdf  ( const double x ) const ;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      double cdf      ( const double x ) const ;
      /// get the integral
      double integral () const { return 1 ; }
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
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
      GramCharlierA  ( const double mean   = 0 ,
                       const double sigma  = 1 ,
                       const double kappa3 = 1 ,
                       const double kappa4 = 1 ) ;
      /// destructor
      ~GramCharlierA () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Gram-Charlier type A approximation
      double pdf         ( const double x ) const ;
      /// evaluate Gram-Charlier type A approximation
      double operator () ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      double   mean   () const { return m_mean          ; }
      double   m0     () const { return   mean   ()     ; }
      double   peak   () const { return   mean   ()     ; }
      double   sigma  () const { return m_sigma         ; }
      double   kappa3 () const { return m_kappa3        ; }
      double   kappa4 () const { return m_kappa4        ; }
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
      double integral ( const double low ,
                        const double high ) const ;
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
    /** @class PhaseSpacePol
     *  simple function to represent the product of N-body phase space
     *  and positive polynomial
     *  @see Ostap::Math::PhaseSpaceNL
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpacePol final
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from thresholds and number of particles
       *  @param threshold_L the low-mass  threshold
       *  @param threshold_H the high-mass threshold
       *  @param l           how many particles we consider
       *  @param n           total number of particles ( n>l!)
       *  @param N           degree of polynomial
       */
      PhaseSpacePol ( const double         threshold_L =  0 ,
                      const double         threshold_H = 10 ,
                      const unsigned short l           =  2 ,
                      const unsigned short n           =  3 ,
                      const unsigned short N           =  1 ) ; // degree of polynomial
      // =====================================================================
      /** constructor from phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       */
      PhaseSpacePol ( const PhaseSpaceNL&  ps      ,
                      const unsigned short N  =  1 ) ; // degree of polynomial
      // ======================================================================
      /** constructor from phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       */
      PhaseSpacePol ( const PhaseSpaceNL&  ps      ,
                      const unsigned short N       ,
                      const double         xlow    ,
                      const double         xhigh   ) ;
      // ======================================================================
      /// destructor 
      ~PhaseSpacePol() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N/L-body modulated phase space
      double evaluate    ( const double x ) const ;
      /// evaluate N/L-body modulated phase space
      double operator () ( const double x ) const  
      { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL& phasespace () const { return m_phasespace ; }
      const Ostap::Math::Positive&     polynom    () const { return m_positive   ; }
      const Ostap::Math::Positive&     positive   () const { return m_positive   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_positive.npars () ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value )
      { return m_positive.setPar ( k , value ) ; }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned short k ) const
      { return m_positive.par ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the phase space
      Ostap::Math::PhaseSpaceNL   m_phasespace ; // the phase space
      Ostap::Math::Positive       m_positive   ; // the positive polynom
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace      m_workspace  ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================

   // ========================================================================
    /** @class GammaDist
     *  Gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  GammaDist
    {
    public:
      // ======================================================================
      /** constructor form scale & shape parameters
       *  param k      \f$k\f$ parameter (shape)
       *  param theta  \f$\theta\f$ parameter (scale)
       */
      GammaDist ( const double k     = 2 ,   // shape parameter
                  const double theta = 1 ) ; // scale parameter
      /// desctructor
      ~GammaDist() ;  // destructor
      // ======================================================================
    public:
      // ======================================================================
      double pdf        ( const double x ) const ;
      /// calculate gamma distribution shape
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // variables
      // ======================================================================
      double k          () const  { return m_k                          ; }
      double theta      () const  { return m_theta                      ; }
      // ======================================================================
    public:
      // ======================================================================
      double mean       () const  { return m_k * m_theta                ; }
      double dispersion () const  { return m_k * m_theta * m_theta      ; }
      double variance   () const  { return dispersion ()                ; }
      double sigma      () const  ;
      double skewness   () const  ;
      // ======================================================================
    public:
      // ======================================================================
      /** effective $\chi^2\f$-parameters
       *  If   \f$ Q  \sim \chi^2(\nu)\f$  and c is a positive constant,
       *  than \f$ cQ \sim \Gamma (k = \nu/2, \theta = 2c) \f$
       */
      double nu () const { return 2   * k     () ; }
      double c  () const { return 0.5 * theta () ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setK     ( const double value  ) ;
      bool setTheta ( const double value  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public: // quantiles
      // ======================================================================
      /// calculate the quantile   (0<p<1)
      double quantile ( const double p ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// shape
      double m_k      ; // shape
      /// scale
      double m_theta  ; // scale
      // ======================================================================
    private:
      // ======================================================================
      /// auxillary intermediate parameter
      mutable double m_aux ; // auxillary intermediate parameter
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
      LogGammaDist ( const double k     = 2 ,   // shape parameter
                     const double theta = 1 ) ; // scale parameter
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
      double k          () const  { return m_gamma.k     () ; }
      double theta      () const  { return m_gamma.theta () ; }
      // ======================================================================
    public:
      // ======================================================================
      double mean       () const  { return m_gamma.mean       () ; }
      double dispersion () const  { return m_gamma.dispersion () ; }
      double sigma      () const  { return m_gamma.sigma      () ; }
      double skewness   () const  { return m_gamma.skewness   () ; }
      // ======================================================================
    public:
      // ======================================================================
      /** effective $\chi^2\f$-parameters
       *  If   \f$ Q  \sim \chi^2(\nu)\f$  and c is a positive constant,
       *  than \f$ cQ \sim \Gamma (k = \nu/2, \theta = 2c) \f$
       */
      double nu () const { return m_gamma.nu () ; }
      double c  () const { return m_gamma.c  () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying gamma distribution
      const GammaDist& gamma() const { return m_gamma ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setK     ( const double value  ) { return m_gamma.setK     ( value ) ; }
      bool setTheta ( const double value  ) { return m_gamma.setTheta ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      virtual double integral ( const double low  ,
                                const double high ) const ;
      // ======================================================================
    public: // quantiles
      // ======================================================================
      /// calculate the quantile   (0<p<1)
      double quantile ( const double p ) const ;
      // ======================================================================
    private:
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
      /** constructor form scale & shape parameters
       *  param k      \f$k\f$ parameter (shape)
       *  param theta  \f$\theta\f$ parameter (scale)
       */
      Log10GammaDist ( const double k     = 2 ,   // shape parameter
                       const double theta = 1 ) ; // scale parameter
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
      double integral   ( const double low  ,
                          const double high ) const override;
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
      GenGammaDist ( const double k     = 2 ,
                     const double theta = 1 ,
                     const double p     = 1 , // 1 corresponds to gamma distribution
                     const double low   = 0 ) ;
      /// desctructor
      ~GenGammaDist() ;  // desctructor
      // ======================================================================
    public:
      // ======================================================================
      /// calculate gamma distribution shape
      double pdf        ( const double x ) const ;
      /// calculate gamma distribution shape
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // getters
      // ======================================================================
      double k          () const  { return m_k                          ; }
      double theta      () const  { return m_theta                      ; }
      double p          () const  { return m_p                          ; }
      double low        () const  { return m_low                        ; }
      // ======================================================================
    public:
      // ======================================================================
      /// Wikipedia notations
      double a          () const { return theta () ; }
      double d          () const { return k     () ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      double mean       () const  { return m_k * m_theta +   low ()     ; }
      double dispersion () const  { return m_k * m_theta * m_theta      ; }
      double variance   () const  { return dispersion ()                ; }
      double sigma      () const  { return std::sqrt ( dispersion ()  ) ; }
      double skewness   () const  { return 2.0 / std::sqrt ( m_k )      ; }
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
      double integral ( const double low  ,
                        const double high ) const ;
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
      Amoroso ( const double theta = 1 ,
                const double alpha = 1 ,
                const double beta  = 1 ,
                const double a     = 0 ) ;
      /// destructor
      ~Amoroso () ;  // desctructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Amoroso distribtion
      double pdf        ( const double x ) const ;
      /// evaluate Amoroso distribtion
      double operator() ( const double x ) const { return pdf ( x ) ; }
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
      double integral ( const double low  ,
                        const double high ) const ;
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
     *  dot not mix with Ostap::Math::LogGammaDist
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
      LogGamma ( const double nu     = 0 ,   // shape parameter
                 const double lambda = 1 ,   // scale parameter
                 const double alpha  = 1 ) ; // scale parameter
      /// destructor
      ~LogGamma () ;  // desctructor
      // ======================================================================
    public:
      // ======================================================================
      /// calculate log-gamma shape
      double pdf        ( const double x ) const ;
      /// calculate log-gamma shape
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getter s
      // ======================================================================
      double nu         () const  { return m_nu     ; }
      double lambda     () const  { return m_lambda ; }
      double alpha      () const  { return m_alpha  ; }
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
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double  m_nu      ;
      double  m_lambda  ;
      double  m_alpha   ;
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
      BetaPrime ( const double alpha = 3 ,
                  const double beta  = 3 ,
                  const double scale = 1 ,
                  const double shift = 0 ) ;
      /// destructor
      ~BetaPrime () ;
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
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_alpha ;
      double m_beta  ;
      double m_scale ;
      double m_shift ;
      // ======================================================================
    private:
      // ======================================================================
      /// auxillary intermediate parameter
      double m_aux ; // auxillary intermediate parameter
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
      Landau ( const double scale = 1 ,
               const double shift = 0 ) ;
      /// destructor
      ~Landau () ;
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
      double scale () const { return m_scale ; }
      double shift () const { return m_shift ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setScale ( const double value ) ;
      bool   setShift ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double cdf      ( const double x    ) const ;
      double integral ( const double low  ,
                        const double high ) const ;
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
    class   Weibull
    {
    public:
      // ======================================================================
      /** constructor from all parameters    
       *  @param scale the scale parameter "lambda"  >0 
       *  @param shape the shape parameter "k"       >0
       *  @param shift the shift parameter "x0"
       */
      Weibull (  const double scale = 1 ,  
                 const double shape = 1 ,  
                 const double shift = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Weibull-distributions
      double pdf        ( const double x ) const ;
      /// evaluate Weibull-distributions
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: //  direct & derived getters 
      // ======================================================================
      double scale () const { return m_scale ; }
      double shape () const { return m_shape ; }
      double shift () const { return m_shift ; }
      // ======================================================================
      double lambd () const { return scale () ; }
      double k     () const { return shape () ; }
      double x0    () const { return shift () ; }      
      // ======================================================================
    public: //  setters  
      // ======================================================================
      bool setScale  ( const double value ) ;
      bool setShape  ( const double value ) ;
      bool setShift  ( const double value ) ;
      // ======================================================================
      bool setLambda ( const double value ) { return setScale ( value ) ; }
      bool setK      ( const double value ) { return setShape ( value ) ; }
      bool setX0     ( const double value ) { return setShift ( value ) ; }
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
      double integral ( const double low  ,
                        const double high ) const ;
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
    /** @class Argus
     *  http://en.wikipedia.org/wiki/ARGUS_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  Argus
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param shape shape parameter 
       *  @param high  high parameter 
       *  @param low   low parameter 
       */
      Argus  ( const double shape  = 1   ,
               const double high   = 1   ,
               const double low    = 0   ) ;
      /// destructor
      ~Argus () ;
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
      double shape  () const { return m_shape  ; }
      double low    () const { return m_low    ; }
      double high   () const { return m_high   ; }
      // ======================================================================
    protected:
      // ======================================================================
      double  y_ ( const double x ) const
      { return ( x - m_low  ) / ( m_high - m_low ) ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setHigh  ( const double value ) ;
      bool   setLow   ( const double value ) ;
      bool   setShape ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double cdf      ( const double x    ) const ;
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_shape ;
      double m_high  ;
      double m_low   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ExpoPositive
     *  useful function for parameterizing smooth background:
     *  product of the exponential and positive polinonmial
     *  @see Ostap::Math::Positive
     */
    class  ExpoPositive 
    {
    public:
      // ======================================================================
      /// constructor from the order
      ExpoPositive ( const unsigned short       N     =  0 ,
                     const double               tau   =  0 , // exponent
                     const double               xmin  =  0 ,
                     const double               xmax  =  1 ) ;
      // ======================================================================
      /// constructor from N phases
      ExpoPositive ( const std::vector<double>& pars       ,
                     const double               tau   =  0 , // exponent
                     const double               xmin  =  0 ,
                     const double               xmax  =  1 ) ;

      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get exponential
      double tau    () const { return m_tau ;}
      /// set new value for the exponent
      bool   setTau ( const  double value ) ;
      /// get number of polinomial parameters
      std::size_t npars () const { return 1 + m_positive.npars() ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value )
      { return
          m_positive.npars() == k ?
          setTau            (     value ) :
          m_positive.setPar ( k , value ) ; }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      {
        return
          m_positive.npars() == k ?
          tau            (   )    :
          m_positive.par ( k )    ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the parameter value
      double  parameter ( const unsigned short k ) const { return par ( k ) ; }
      /// get lower edge
      double xmin () const { return m_positive.xmin () ; }
      /// get upper edge
      double xmax () const { return m_positive.xmax () ; }
      /// transform variables
      double x ( const double t ) const { return m_positive. x ( t )  ; }
      double t ( const double x ) const { return m_positive. t ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying positive function
      const Ostap::Math::Positive&  positive  () const { return m_positive  ; }
      // ======================================================================
    public:
      // ======================================================================
      double integral ( const double low , const double high ) const ;
      double integral () const { return integral ( xmin() , xmax() ) ; }
      // ======================================================================
    private:
      // ======================================================================
      Ostap::Math::Positive m_positive  ;
      double                m_tau       ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Sigmoid
     *  Sigmoid function, modulated by the positive polynomial
     *  \f$ f(x) = ( 1 + tanh( \alpha ( x  - x_0) ) \times P_{pos} (x) \f$
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  Sigmoid 
    {
    public:
      // ============================================================
      /// constructor from polynom and parameters "alpha" and "x0"
      Sigmoid
        ( const Ostap::Math::Positive& poly      ,
          const double                 alpha = 0 ,
          const double                 x0    = 0 ) ;
      /// constructor from polynom and parameter "alpha"
      Sigmoid
        ( const unsigned short             N = 0 ,
          const double                 xmin  = 0 ,
          const double                 xmax  = 1 ,
          const double                 alpha = 0 ,
          const double                 x0    = 0 ) ;
      /// constructor from polynom and parameter "alpha"
      Sigmoid
        ( const std::vector<double>&   pars       ,
          const double                 xmin  =  0 ,
          const double                 xmax  =  1 ,
          const double                 alpha =  0 ,
          const double                 x0    =  0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double alpha      () const { return m_alpha ; }
      // set new value of alpha
      bool setAlpha     ( const double value ) ;
      // ======================================================================
      double x0         () const { return m_x0    ; }
      // set new value of alpha
      bool setX0        ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const { return m_positive.xmin () ; }
      double xmax () const { return m_positive.xmax () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return 2 + m_positive.npars () ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value )
      { return
          m_positive.npars()     == k ? setAlpha ( value ) :
          m_positive.npars() + 1 == k ? setX0    ( value ) :
          m_positive.setPar (  k  , value ) ; }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { return
          m_positive.npars()     == k ? alpha () :
          m_positive.npars() + 1 == k ? x0    () : m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned short k ) const { return par ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Positive& positive() const { return m_positive ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const ;
      /// get the integral between low and high
      double integral   ( const double low  ,
                          const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      Ostap::Math::Positive  m_positive ;
      double                 m_alpha    ;
      double                 m_x0       ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace for integration
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    };
    // ========================================================================
    /** @class TwoExpos
     *  simple difference of two exponents
     *  \f$ f \propto
     *        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} =
     *        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  TwoExpos 
    {
    public:
      // ======================================================================
      TwoExpos ( const double alpha = 1 ,
                 const double delta = 1 ,
                 const double x0    = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get alpha
      double alpha () const { return m_alpha ; }
      /// get delta
      double delta () const { return m_delta ; }
      /// get x0
      double x0    () const { return m_x0    ; }
      // ======================================================================
    public:
      // ======================================================================
      /// slope for the first  exponent
      double a1         () const { return m_alpha           ; }
      /// slope for the second exponent
      double a2         () const { return m_alpha + m_delta ; }
      /// mean-value (for -inf,+inf) interval
      double mean       () const ;
      /// mode
      double mode       () const ;
      /// variance
      double variance   () const ;
      /// dispersion
      double dispersion () const { return variance () ; }
      /// sigma
      double sigma      () const ;
      // get normalization constant
      double norm       () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// slope for the first  exponent
      double tau1       () const { return -a1() ; }
      /// slope for the second exponent
      double tau2       () const { return -a2() ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setAlpha ( const double value ) ;
      bool setDelta ( const double value ) ;
      bool setX0    ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between -inf and +inf
      double integral    () const ;
      /// get the integral between low and high
      double integral    ( const double low  ,
                           const double high ) const ;
      /// get the derivative at given value
      double derivative  ( const double x    ) const ;
      /// get the second at given value
      double derivative2 ( const double x    ) const ;
      /// get the Nth derivative at given value
      double derivative  ( const double   x  ,
                           const unsigned N  ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_alpha ;
      double m_delta ;
      double m_x0    ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class TwoExpoPositive
     *  simple difference of two exponents modulated with positive polynomials
     *  @see TwoExpos
     *  @see Positive
     *  @see ExpoPositive
     *  \f$ f(x) = e_2(x) * p_n(x) \f$, where
     *  \f$ e_2(x) \propto
     *        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} =
     *        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
     *  and $p_2(s)$ is positive polynomial function
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-03-28
     */
    class  TwoExpoPositive 
    {
    public:
      // ======================================================================
      TwoExpoPositive
        ( const unsigned short N = 1 ,
          const double alpha     = 1 ,
          const double delta     = 1 ,
          const double x0        = 0 ,
          const double xmin      = 0 ,
          const double xmax      = 1 ) ;
      // ======================================================================
      TwoExpoPositive
        ( const std::vector<double>& pars ,
          const double alpha  = 1 ,
          const double delta  = 1 ,
          const double x0     = 0 ,
          const double xmin   = 0 ,
          const double xmax   = 1 ) ;
      // ======================================================================
      TwoExpoPositive
        ( const Positive& poly    ,
          const double alpha  = 1 ,
          const double delta  = 1 ,
          const double x0     = 0 ) ;
      // ======================================================================
      TwoExpoPositive
        ( const Positive& poly   ,
          const TwoExpos& expos  ) ;
      // ======================================================================
      TwoExpoPositive
        ( const TwoExpos& expos  ,
          const Positive& poly   ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of polinomial parameters
      std::size_t npars () const { return 3 + m_positive.npars() ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value )
      { return
          m_positive.npars()     == k ? setAlpha ( value ) :
          m_positive.npars() + 1 == k ? setDelta ( value ) :
          m_positive.npars() + 2 == k ? setX0    ( value ) :
          m_positive.setPar ( k , value ) ; }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { return
          m_positive.npars()     == k ? alpha () :
          m_positive.npars() + 1 == k ? delta () :
          m_positive.npars() + 2 == k ? x0    () : m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned short k ) const { return par ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin () const { return m_positive.xmin () ; }
      /// get upper edge
      double xmax () const { return m_positive.xmax () ; }
      /// transform variables
      double x ( const double t ) const { return m_positive. x ( t )  ; }
      double t ( const double x ) const { return m_positive. t ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get alpha
      double alpha () const { return m_2exp.alpha () ; }
      /// get delta
      double delta () const { return m_2exp.delta () ; }
      /// get x0
      double x0    () const { return m_2exp.x0    () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// slope for the first  exponent
      double a1         () const { return m_2exp.a1   () ; }
      /// slope for the second exponent
      double a2         () const { return m_2exp.a2   () ; }
      /// slope for the first  exponent
      double tau1       () const { return m_2exp.tau1 () ; }
      /// slope for the second exponent
      double tau2       () const { return m_2exp.tau2 () ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setAlpha ( const double value ) { return m_2exp.setAlpha ( value ) ; }
      bool setDelta ( const double value ) { return m_2exp.setDelta ( value ) ; }
      bool setX0    ( const double value ) { return m_2exp.setX0    ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral    () const ;
      /// get the integral between low and high
      double integral    ( const double low  ,
                           const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying positive function
      const Ostap::Math::Positive&  positive  () const { return m_positive  ; }
      /// get the underlying exponents
      const Ostap::Math::TwoExpos&  twoexpos  () const { return m_2exp      ; }
      // ======================================================================
    private:
      // ======================================================================
      Ostap::Math::Positive m_positive ;
      Ostap::Math::TwoExpos m_2exp     ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Tsallis
     *
     *  Useful function to describe pT-spectra of particles
     *
     *  - C. Tsallis,
     *  "Possible generalization of Boltzmann-Gibbs statistics,
     *  J. Statist. Phys. 52 (1988) 479.
     *  - C. Tsallis,
     *  Nonextensive statistics: theoretical, experimental and computational
     *  evidences and connections, Braz. J. Phys. 29 (1999) 1.
     *
     *  \f[ \frac{d\sigma}{dp_T} \propto
     *    p_T\times \left( 1 + \frac{E_{kin}}{Tn}\right)^{-n}\f],
     *  where \f$E_{kin} = \sqrt{p_T^2-M^2}-M\f$
     *  is transverse kinetic energy
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-07-11
     */
    class  Tsallis 
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param mass particle mass  (M>0)
       *  @param n    n-parameter    (N>1)
       *  @param T    T-parameter    (T>0)
       */
      Tsallis
        ( const double mass         = 0   ,
          const double n            = 10  ,
          const double T            = 1.1 ) ;
      /// destructor
      ~Tsallis() ;
      // ======================================================================
    public:
      // ======================================================================
      /// get Tsallis PDF
      double pdf ( const double x ) const ;
      /// get Tsallis PDF
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get mass-parameter
      double mass () const { return m_mass  ; } // get mass-parameter
      /// get n-parameter
      double n    () const { return m_n     ; } // get n-parameter
      /// get T-parameter
      double T    () const { return m_T     ; } // get T-parameter
      // ======================================================================
      // aliases
      // ======================================================================
      /// get mass-parameter
      double m    () const { return  mass () ; } // get mass-parameter
      /// get mass-parameter
      double M    () const { return  mass () ; } // get mass-parameter
      /// get n-parameter
      double N    () const { return  n    () ; } // get n-parameter
      // ======================================================================
    public:
      // ======================================================================
      /// update mass-parameter
      bool setMass ( const double value ) ; // update mass-parameter
      bool setM    ( const double value ) { return setMass ( value ) ; }
      /// update n-parameter
      bool setN    ( const double value ) ; // update n-parameter
      /// update T-parameter
      bool setT    ( const double value ) ; // update T-parameter
      // ======================================================================
    public:
      // ======================================================================
      /// get the min-value
      double xmin() const { return 0 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the transverse kinetic energy
      inline double eTkin ( const double x ) const
      { return std::sqrt ( x * x + m_mass * m_mass ) - m_mass ; }
      // ======================================================================
      /// get the transverse mass
      inline double mT   ( const double x ) const
      { return std::sqrt ( x * x + m_mass * m_mass ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral    ( const double low  ,
                           const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mass ;
      double m_n    ;
      double m_T    ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ; // workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class QGSM
     *
     *  Useful function to describe pT-spectra of particles
     *
     * - A. B. Kaidalov and O. I. Piskunova, Z. Phys. C 30 (1986) 145.
     * - O. I. Piskounova, arXiv:1301.6539 [hep-ph];
     * - O. I. Piskounova, arXiv:1405.4398 [hep-ph].
     * - A. A. Bylinkin and O. I. Piskounova,
     *  "Transverse momentum distributions of baryons at LHC energies",
     *  arXiv:1501.07706.
     *
     *  \f[ \frac{d\sigma}{dp_T} \propto
     *  p_T \times \mathrm{e}^{ -b_0 (m_T-m)} \f],
     *  where transverse mass is defined as \f$m_T = \sqrt{p_T^2+m^2}\f$
     *
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-07-11
     */
    class  QGSM
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param mass particle mass        (M>0)
       *  @param b    b-parameter/slope    (b>0)
       */
      QGSM
        ( const double mass         = 0   ,    // partile mass
          const double b            = 1   ) ;  // slope
      /// destructor
      ~QGSM () ;
      // ======================================================================
    public:
      // ======================================================================
      /// define QGSM pdf
      double pdf ( const double x ) const ;
      /// define QGSM pdf
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get mass-parameter
      double mass () const { return m_mass  ; } // get mass-parameter
      /// get n-parameter
      double b    () const { return m_b     ; } // get b-parameter
      // ======================================================================
      // aliases
      // ======================================================================
      /// get mass-parameter
      double m    () const { return  mass () ; } // get mass-parameter
      /// get mass-parameter
      double M    () const { return  mass () ; } // get mass-parameter
      /// get b-parameter
      double B    () const { return  b    () ; } // get n-parameter
      /// get b-parameter
      double B0   () const { return  b    () ; } // get n-parameter
      /// get b-parameter
      double b0   () const { return  b    () ; } // get n-parameter
      // =====================================================================
    public:
      // ======================================================================
      /// update mass-parameter
      bool setMass ( const double value ) ; // update mass-parameter
      bool setM    ( const double value ) { return setMass ( value ) ; }
      /// update b-parameter
      bool setB    ( const double value ) ; // update b-parameter
      // ======================================================================
    public:
      // ======================================================================
      /// get the min-value
      double xmin() const { return 0 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the transverse kinetic energy
      inline double eTkin ( const double x ) const
      { return mT( x ) - m_mass ; }
      // ======================================================================
      /// get the transverse mass
      inline double mT   ( const double x ) const
      { return std::sqrt ( x * x + m_mass * m_mass ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral    ( const double low  ,
                           const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mass ;
      double m_b    ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ; // workspace
      // ======================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_FUNCTIONS_H
// ============================================================================
