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
// ============================================================================
/** @file Ostap/Models.h
 *  set of useful models
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
     *  If  x is distributed accroding  to \f$ f(x) \propto e^{-\tau x} \f$, 
     *  then z, \f$ z  =   log(x) \f$, is distributed accoring to 
     *  \f$ F(z) = G(x, -\log(\tau), 1 ) \f$ 
     * 
     *  As a   result, if    x is distributes as sum of exponential components 
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
    /** @class PhaseSpacePol
     *  Function to represent the product of l/n-body phase space
     *  and positive polynomial
     *  \f[ \Phi_{l,n}^{(N)(x)} \equiv 
     *      \Phi_{l,n}(x;x_{low},x_{high}) P_{N}(x) \f]
     *  where :
     *  -  \f$  \Phi_{l,n}(x;x_{low},x_{high}) \f$  is a phase space of 
     *     l-particles from n-body decay
     *  -  \f$ P_{N}(x) \f$ is a positive polynomial of degree N
     *  
     *  @see Ostap::Math::PhaseSpaceNL
     *  @see Ostap::Math::Positive
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
      PhaseSpacePol
        ( const double         threshold_L =  0 ,
          const double         threshold_H = 10 ,
          const unsigned short l           =  2 ,
          const unsigned short n           =  3 ,
          const unsigned short N           =  1 ) ; // degree of polynomial
      // =====================================================================
      /** constructor from phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       */
      PhaseSpacePol 
        ( const PhaseSpaceNL&  ps      ,
          const unsigned short N  =  1 ) ; // degree of polynomial
      // ======================================================================
      /** constructor from phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       */
      PhaseSpacePol 
        ( const PhaseSpaceNL&  ps      ,
          const unsigned short N       ,
          const double         xlow    ,
          const double         xhigh   ) ;
      // ======================================================================
      /// constructor from phase space and polynomial
      PhaseSpacePol 
        ( const PhaseSpaceNL&          ps  ,
          const Ostap::Math::Positive& pol ) ;
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
      ///  get all parameters 
      // const  std::vector<double>& pars() const { return m_positive.pars() ; }
      // get the order of polynomial 
      unsigned short n () const { return m_positive.degree() ; }
      // ======================================================================
      double xmin () const { return m_positive.xmin() ; }
      double xmax () const { return m_positive.xmax() ; }
      // ======================================================================
    public:
      // ======================================================================
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
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL* operator->() const 
      { return &m_phasespace ; }
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
    /** @class PhaseSpaceLeftExpoPol
     *  Function to represent the product of l-body phase space, 
     *  positive polynomial and the exponential function 
     *  \f[ \Phi_{l}^{(N)(x)} \propto
     *      \Phi_{l}(x;x_{low}) \mathrm{e}^{-\left|\tau\right| x } P_{N}(x) \f]
     *  where :
     *  -  \f$  \Phi_{l}(x;x_{low}) \f$  is a phase space of 
     *     l-particles near the threshold 
     *  -  \f$ P_{N}(x) \f$ is a positive polynomial of degree N
     *  
     *  @see Ostap::Math::PhaseSpaceLeft
     *  @see Ostap::Math::Positive
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2018-10-21
     */
    class PhaseSpaceLeftExpoPol final
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from threshold and number of particles
       *  @param threshold_L the low-mass  threshold
       *  @param l           how many particles we consider
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol 
        ( const double         threshold_L =  0 ,   // low threshold 
          const unsigned short l           =  2 ,   // number of particles 
          const unsigned short N           =  1 ,   // degree of polynomial
          const double         tau         =  0 ,   // the exponent 
          const double         xhigh       =  1 ) ; // high edge 
      // =====================================================================
      /** constructor from threshold and number of particles
       *  @param threshold_L the low-mass  threshold
       *  @param l           how many particles we consider
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xlow        the low  edge 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
        ( const double         threshold_L ,   // low threshold 
          const unsigned short l           ,   // number of particles 
          const unsigned short N           ,   // degree of polynomial
          const double         tau         ,   // the exponent 
          const double         xlow        ,   // low edge 
          const double         xhigh       ) ; // high edge 
      // =====================================================================
      /** constructor from the phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
        ( const PhaseSpaceLeft& ps        ,
          const unsigned short  N     = 1 ,   // degree of polynomial
          const double          tau   = 0 ,   // the exponent 
          const double          xhigh = 1 ) ; // high edge 
      // =========================================================================
      /** constructor from the phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xlow        the low  edge 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
        ( const PhaseSpaceLeft& ps    ,
          const unsigned short  N     ,   // degree of polynomial
          const double          tau   ,   // the exponent 
          const double          xlow  ,   // low edge 
          const double          xhigh ) ; // high edge
      // ======================================================================
      /** constructor from the phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xlow        the low  edge 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
        ( const PhaseSpaceLeft&        ps  ,   // pjase space 
          const Ostap::Math::Positive& pol ,   // polynomial 
          const double                 tau ) ; // the exponent 
      // ======================================================================
      /// destructor 
      ~PhaseSpaceLeftExpoPol() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate modulated phase space
      double evaluate    ( const double x ) const ;
      /// evaluate modulated phase space
      double operator () ( const double x ) const  
      { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceLeft& phasespace () const { return m_phasespace ; }
      const Ostap::Math::Positive&       polynom    () const { return m_positive   ; }
      const Ostap::Math::Positive&       positive   () const { return m_positive   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars   () const { return m_positive.npars () ; }
      /// set k-parameter
      bool setPar         ( const unsigned short k , const double value )
      { return m_positive.setPar ( k , value ) ; }
      /// set k-parameter
      bool setParameter   ( const unsigned short k , const double value )
      { return setPar     ( k , value ) ; }
      /// get the parameter value
      double  par         ( const unsigned short k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter   ( const unsigned short k ) const
      { return m_positive.par ( k ) ; }
      ///  get all parameters 
      // const  std::vector<double>& pars() const { return m_positive.pars() ; }
      // get the order of polynomial 
      unsigned short n   () const { return m_positive.degree() ; }
      /// get the exponent 
      double   tau       () const { return m_tau ; }
      /// get the threshold  
      double   threshold () const { return m_phasespace.threshold() ; }
      // ======================================================================
      double xmin () const { return m_positive.xmin() ; }
      double xmax () const { return m_positive.xmax() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set the new exponent 
      bool setTau   ( const double value ) ;
      /// set the   scale  
      bool setScale ( const double value ) 
      { return m_phasespace.setScale ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
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
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceLeft* operator->() const 
      { return &m_phasespace ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the phase space
      Ostap::Math::PhaseSpaceLeft m_phasespace ; // the phase space
      Ostap::Math::Positive       m_positive   ; // the positive polynom
      double                      m_tau        ; // the exponent
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace      m_workspace {} ; // integration workspace
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
      GammaDist 
      ( const double k     = 2 ,   // shape parameter
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
      double k          () const  { return m_k       ; }
      double theta      () const  { return m_theta   ; }
      // ======================================================================
      double alpha      () const  { return m_k       ; }
      double beta       () const  { return 1/m_theta ; }
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
      /** effective \f$ \chi^2 \f$-parameters
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
      double integral 
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
      LogGammaDist
      ( const double k     = 2 ,   // shape parameter
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
      /** effective \f$\chi^2\f$-parameters
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
      /** constructor form scale & shape parameters
       *  param k      \f$k\f$ parameter (shape)
       *  param theta  \f$\theta\f$ parameter (scale)
       */
      Log10GammaDist 
      ( const double k     = 2 ,   // shape parameter
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
      LogGamma
      ( const double nu     = 0 ,   // shape parameter
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
      double lambd      () const  { return m_lambda ; }
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
      /// destructor
      ~GenBetaPrime () ;
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
      double xmin  () const { return x0    () ; }      
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
      ExpoPositive 
      ( const unsigned short       N     =  0 ,
        const double               tau   =  0 , // exponent
        const double               xmin  =  0 ,
        const double               xmax  =  1 ) ;
      // ======================================================================
      /// constructor from N phases
      ExpoPositive 
      ( const std::vector<double>& pars       ,
        const double               tau   =  0 , // exponent
        const double               xmin  =  0 ,
        const double               xmax  =  1 ) ;
      // ======================================================================
      /// constructor from polynom and exponential 
      ExpoPositive 
      ( const Ostap::Math::Positive& pol        , 
        const double                 tau   =  0 ) ;// exponent
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
      double integral () const { return integral ( xmin() , xmax() ) ; }
      double integral ( const double low , 
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
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
      TwoExpos 
      ( const double alpha = 1 ,
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
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================      
      /// get the derivative at given value
      double derivative  ( const double x    ) const ;
      /// get the second at given value
      double derivative2 ( const double x    ) const ;
      /// get the Nth derivative at given value
      double derivative 
      ( const double   x  ,
        const unsigned N  ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
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
    public:
      // ======================================================================
      // get the tag
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      Ostap::Math::Positive m_positive ;
      Ostap::Math::TwoExpos m_2exp     ;
      // ======================================================================
    } ;
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
      double operator() ( const double x ) const ;
      // ======================================================================
    public: // setters & getters 
      // ======================================================================
      /// parameter nu 
      double nu       () const { return m_nu       ; }
      /// parameter varsigma 
      double varsigma () const { return m_varsigma ; }
      // parameger shift 
      double shift    () const { return m_shift    ; }
      /// minimal x 
      double xmin     () const { return m_shift    ; }
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
      double integral   ( const double low  ,
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
      // evaluat the function 
      inline double operator() ( const double x ) const 
      { return evaluate ( x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// parameter theta 
      double theta () const { return m_theta ; }
      /// parameter eta 
      double eta   () const { return m_eta   ; }  
      /// parameter p 
      double p     () const { return m_p ; }  
      /// parameter shif
      double shift () const { return m_shift ; }  
      // ======================================================================
    public: // derived parameters
      // ======================================================================      
      // parameter a (derived) 
      double a () const ;
      // parameter b (derived) 
      double b () const ;
      /// minimal x 
      double xmin () const { return m_shift ; }
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
      double integral   ( const double low  ,
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
    /** @class Argus 
     *  Slightly modified version of Argus distribution, with 
     *  support in the interval  \f$ \mu - c \le x \le \mu \f$
     *  @see https://en.wikipedia.org/wiki/ARGUS_distribution
     *  @see ARGUS Collaboration, H. Albrecht et al., 
     *      "Measurement of the polarization in the decay B → J/ψK*". 
     *      Physics Letters B. 340 (3): 217–220.
     *  @see doi:10.1016/0370-2693(94)01302-0.
     *  @see https://doi.org/10.1016%2F0370-2693%2894%2901302-0
     */
    class Argus 
    {
      // ======================================================================
    public :
      // ======================================================================
      // constructor for all elements 
      Argus 
      ( const double mu  = 1 , 
        const double  c  = 1 , 
        const double chi = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate function 
      double        evaluate   ( const double x ) const ;
      /// get PDF
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// get PDF
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:  // gettters 
      // ====================================================================== 
      /// parameter mu 
      double mu  () const { return m_mu  ; }
      /// parameter c 
      double c   () const { return m_c   ; }
      /// parameter chi 
      double chi () const { return m_chi ; }        
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set mu parameter
      bool setMu  ( const double value ) ;
      /// set c parameter
      bool setC   ( const double value ) ;
      /// set chi parameter
      bool setChi ( const double value ) ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean of the distribution 
      double mean     () const ;
      /// mode of the distribution 
      double mode     () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// xmin
      double xmin () const { return m_mu - m_c ; }
      /// xmax
      double xmax () const { return m_mu       ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral   () const ;
      /// get the integral between low and high
      double integral  
      ( const double low  ,
        const double high ) const ;
      /// get CDF 
      double cdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private: // helper function
      // ======================================================================
      /** helper function 
       *  \f$ \Psi ( \chi ) = \Phi(\chi )  - \chi \phi  (\chi ) - \frac{1}{2} \f$ 
       */
      double psi ( const double value ) const ;
      // ======================================================================
   private:
      // ======================================================================
      /// parameter mu 
      double   m_mu   {  1 } ; // parameter mu      
      /// parameter c 
      double   m_c    {  1 } ; // parameter c 
      /// parameter chi 
      double   m_chi  {  1 } ; // parameter chi
      /// normalization 
      double   m_norm { -1 } ; // normalization
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenArgus 
     *  Slightly modified version of generalized Argus distribution, with 
     *  support in the interval  \f$ \mu - c \le x \le \mu \f$
     *  @see https://en.wikipedia.org/wiki/ARGUS_distribution
     *  @see ARGUS Collab oration, H. Albrecht et al., 
     *      "Measurement of the polarization in the decay B → J/ψK*". 
     *      Physics Letters B. 340 (3): 217–220.
     *  @see doi:10.1016/0370-2693(94)01302-0.
     *  @see https://doi.org/10.1016%2F0370-2693%2894%2901302-0
     */
    class GenArgus 
    {
      // ======================================================================
    public :
      // ======================================================================
      // constructor for all elements 
      GenArgus 
      ( const double mu  = 1   , 
        const double  c  = 1   ,  
        const double chi = 1   ,
        const double dp  = 1.5 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate function 
      double        evaluate   ( const double x ) const ;
      /// get PDF
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// get PDF
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:  // gettters 
      // ====================================================================== 
      /// parameter mu 
      double mu  () const { return m_mu  ; }
      /// parameter c 
      double c   () const { return m_c   ; }
      /// parameter chi 
      double chi () const { return m_chi ; }        
      /// parameter dp 
      double dp  () const { return m_dp  ; }        
      // ======================================================================
    public:
      // ======================================================================
      /// parameter p 
      double p () const { return m_dp - 1 ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set mu parameter
      bool setMu  ( const double value ) ;
      /// set c parameter
      bool setC   ( const double value ) ;
      /// set chi parameter
      bool setChi ( const double value ) ;
      /// set dp parameter
      bool setDp  ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// xmin
      double xmin () const { return m_mu - m_c ; }
      /// xmax
      double xmax () const { return m_mu       ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral   () const ;
      /// get the integral between low and high
      double integral   ( const double low  ,
                          const double high ) const ;
      /// get CDF 
      double cdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
   private:
      // ======================================================================
      /// parameter mu 
      double   m_mu   {  1 } ; // parameter mu      
      /// parameter c 
      double   m_c    {  1 } ; // parameter c 
      /// parameter chi 
      double   m_chi  {  1 } ; // parameter chi
      /// parameter dp 
      double   m_dp   {  1 } ; // parameter chi
      /// normalization 
      double   m_norm { -1 } ; // normalization
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
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// n-parameter 
      unsigned short n    () const { return m_bspline.degree () + 1 ; } ;
      // ======================================================================
    public: // 
      // ======================================================================
      double         xmin () const { return m_bspline.xmin   () ; }
      double         xmax () const { return m_bspline.xmax   () ; }      
      // ======================================================================
    public: // some statistical properties 
      // ======================================================================
      /// get mode 
      double mode       () const { return 0.5 * n ()  ; }
      /// get mean 
      double mean       () const { return 0.5 * n ()  ; }
      /// get median
      double median     () const { return 0.5 * n ()  ; }
      /// get variance 
      double variance   () const { return n () / 12.0 ; }
      /// get doispersion 
      double dispersion () const { return n () / 12.0 ; }
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
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// n-parameter 
      unsigned short n    () const { return m_ih.n () ; } ;
      // ======================================================================
    public: // 
      // ======================================================================
      double         xmin () const { return 0 ; }
      double         xmax () const { return 1 ; }      
      // ======================================================================
    public: // some statistical properties 
      // ====================================================================== 
      /// get mode 
      double mode       () const { return 0.5 ; }
      /// get mean 
      double mean       () const { return 0.5 ; }
      /// get median 
      double median     () const { return 0.5 ; }
      /// get variance 
      double variance   () const ; 
      /// get dispersion 
      double dispersion () const { return variance () ; }
      /// get rms 
      double rms        () const ;
      /// get skewness 
      double skewness   () const { return 0           ; }
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
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// mu-parameter 
      double         mu    () const { return m_mu         ; }
      /// sigma-parameter 
      double         sigma () const { return m_sigma      ; }
      /// n-parameter 
      unsigned short n     () const { return m_bates.n () ; } 
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
      double mode       () const { return m_mu ; }
      /// get mean 
      double mean       () const { return m_mu ; }
      /// get median 
      double median     () const { return m_mu ; }
      /// get variance 
      double variance   () const { return m_sigma * m_sigma ; }
      /// get dispersion 
      double dispersion () const { return variance () ; }
      /// get rms 
      double rms        () const { return m_sigma     ; }
      /// get skewness 
      double skewness   () const { return 0           ; }
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
     *  where \f$E_{kin} = \sqrt{p_T^2+M^2}-M\f$
     *  is transverse kinetic energy
     *  @see Ostap::Math::Tsallis2 
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
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
      ( const double mass         = 1   ,
        const double n            = 10  ,
        const double T            = 1.1 ) ;
      /// destructor
      ~Tsallis() ;
      // ======================================================================
    public:
      // ======================================================================
      /** evaluate Tsallis functuon
       *  @param pt transverse momentum of the particle 
       */
      double evaluate   ( const double pt ) const ;
      /** evaluate Tsallis functuon
       *  @param pt transverse momentum of the particle 
       */
      double operator() ( const double pt ) const { return evaluate ( pt ) ; }
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
      /// q-parameter for Tsallis enthropy 
      double q    () const { return m_n / ( m_n - 1 ) ; }
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
      inline double eTkin ( const double pt ) const
      { return mT ( pt )  - m_mass ; }
      // ======================================================================
      /// get the transverse mass
      inline double mT   ( const double pt ) const
      { return std::hypot ( pt ,  m_mass ) ; }
      // ======================================================================
    public:
      // ======================================================================
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
      double B_0   () const { return  b    () ; } // get n-parameter
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
      inline double eTkin ( const double pt ) const
      { return mT ( pt ) - m_mass ; }
      // ======================================================================
      /// get the transverse mass
      inline double mT   ( const double pt ) const
      { return std::hypot ( pt ,  m_mass  ) ; }
      // ======================================================================
    public:
      // ======================================================================
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
    /** @class Hagedorn 
     *  simple function to des ribe pT spectra of particles 
     *  @see R.Hagedorn, "Multiplicities, p_T distributions and the 
     *       expected hadron \to Quark - Gluon Phase Transition", 
     *       Riv.Nuovo Cim. 6N10 (1983) 1-50
     *  @see https://doi.org/10.1007/BF02740917 
     *  @see https://inspirehep.net/literature/193590
     *  
     *  \f[ f(p_T; m, T) \propto 
     *   p_T \sqrt{p^2_T + m^2} K_1( \beta \sqrt{ p^2_T+m^2} ) \f] 
     *
     *  where \f$ \beta \f$ is inverse temporature 
     *  \f$ \beta = \frac{1}{T} f$ 
     *
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2022-12-06
     */
    class Hagedorn 
    {
    public: 
      // ======================================================================      
      /** Constrcutor for all parameters
       *  @param mass   mas sof th eparticle 
       *  @param beta   inverse temperature
       */
      Hagedorn 
      ( const double mass = 0 , 
        const double beta = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluiate the function 
      double        evaluate    ( const double x ) const ;
      /// ditto
      inline double operator()  ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: /// getters 
      // ======================================================================
      /// get the particle mass 
      inline double mass () const { return m_mass ; }
      /// get the inverse temperature 
      inline double beta () const { return m_beta ; }      
      // ======================================================================
    public: /// setters 
      // ======================================================================
      /// set new mass 
      bool setMass ( const double value ) ;
      /// set new inverse temperatoire  
      bool setBeta ( const double beta  ) ;
      // ======================================================================
    public:  /// some statistics 
      // ======================================================================
      /// get the mean value 
      double mean () const ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the transverse mass 
      inline double mT ( const double pT ) const 
      { return std::hypot ( pT , m_mass ) ; }  
      // ======================================================================
    public:
      // ======================================================================
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
      /// particle mass 
      double m_mass { 0 } ; // particle mass 
      /// inverse temperature 
      double m_beta { 1 } ; // inverse temperature 
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ; // workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class HORNSdini 
     *  \f[ f(x;a,\delta, \phi) = 
     *  \frac{3}{2\delta}\left( z \right)^2
     *  \left( \cos^2( \phi + \frac{\pi}{4}) ( 1 + z ) +
     *         \sin^2( \phi + \frac{\pi}{4}) ( 1 - z ) \right) \f]
     *  where  \f$ z = \frac{ x - ( a - \delta ) } { \delta } \f$ 
     *  for \f$ a \le x \le a + 2\delta\$ and zero otherwise 
     *  
     * The first factor accound for two-horn parabolic shape, 
     * and the second factor accouns for the linear correction factor 
     * ("efficiency")
     *
     *  - For the actual use it needs to be convoluted with resolution function 
     * @see 
     */
    class HORNSdini 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor fron all parameters
       *  @param a position of th eleft paraboilic horn
       *  @param delta distance fron left to right parabolic horn 
       *  @param phi  correction parameter ("efficiency")
       */
      HORNSdini
      ( const double a     = 0 , 
        const double delta = 1 ,
        const double phi   = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evalaute the function 
      double evaluate ( const double x )  const ;
      /// evalaute the function 
      inline double operator () ( const double x )  const 
      { return evaluate ( x ) ;  }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// left horn  
      double a     () const { return m_a     ; }
      /// right horn
      double b     () const { return m_a + 2 * m_delta ; }
      /// delta 
      double delta () const { return m_delta ; }
      /// phi 
      double phi   () const { return m_phi  ; }
      // ======================================================================
    public :
      // ======================================================================
      double xmin () const { return a () ; }
      double xmax () const { return b () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setA      ( const double value ) ;
      bool setDelta  ( const double value ) ;
      bool setPhi    ( const double value ) ;
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
      double m_a        { 0 } ;
      double m_delta    { 1 } ;
      double m_phi      { 0 } ;
      double m_cos2_phi { 0 } ;
      double m_sin2_phi { 0 } ;
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class HILLdini 
     *  \f[ f(x;a,\delta, \phi) = 
     *  \frac{3}{2\delta}\left( 1 - z \right)^2
     *  \left( \cos^2( \phi + \frac{\pi}{4}) ( 1 + z ) +
     *         \sin^2( \phi + \frac{\pi}{4}) ( 1 - z ) \right) \f]
     *  where  \f$ z = \frac{ x - ( a - \delta ) } { \delta } \f$ 
     *  for \f$ a \le x \le a + 2\delta\$ and zero otherwise 
     *  
     * The first factor accound for a parabolic shape, 
     * and the second factor accouns for the linear correction factor 
     * ("efficiency")
     *
     *  - For the actual use it needs to be convoluted with resolution function 
     * @see 
     */
    class HILLdini 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor fron all parameters
       *  @param a position of th eleft paraboilic horn
       *  @param delta distance fron left to right parabolic horn 
       *  @param phi  correction parameter ("efficiency")
       */
      HILLdini
      ( const double a     = 0 , 
        const double delta = 1 ,
        const double phi   = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evalaute the function 
      double evaluate ( const double x )  const ;
      /// evalaute the function 
      inline double operator () ( const double x )  const 
      { return evaluate ( x ) ;  }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// left horn  
      double a     () const { return m_a     ; }
      /// right horn
      double b     () const { return m_a + 2 * m_delta ; }
      /// delta 
      double delta () const { return m_delta ; }
      /// phi 
      double phi   () const { return m_phi  ; }
      // ======================================================================
    public :
      // ======================================================================
      double xmin () const { return a () ; }
      double xmax () const { return b () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setA      ( const double value ) ;
      bool setDelta  ( const double value ) ;
      bool setPhi    ( const double value ) ;
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
      double m_a        { 0 } ;
      double m_delta    { 1 } ;
      double m_phi      { 0 } ;
      double m_cos2_phi { 0 } ;
      double m_sin2_phi { 0 } ;
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
      double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      double mu    () const { return m_mu    ; }
      double scale () const { return m_scale ; }
      double shape () const { return m_shape ; }
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
      double xmin  () const { return m_mu    ; }
      double xi    () const { return m_shape ; }
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
      double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      double mu    () const { return m_mu    ; }
      double scale () const { return m_scale ; }
      double shape () const { return m_shape ; }
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
      double xi         () const { return m_shape ; }
      /// mode of distribution 
      double mode       () const ;
      /// mean value
      double mean       () const ;
      /// variance 
      double variance   () const ;
      /// dispersion
      double dispersion () const { return variance () ; }
      /// rms 
      double rms        () const ;
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
     * Cumulative distributuon fnuctiin is defined as
     * \f[ F(x) = 1 - \mathem{e}^{ - \sum_i \left| p_i \right| \Delta^i } \f]
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
      /// two shape parameters: alpha and beta 
      Benini
      ( const double               scale = 1 ,              // scale parameter 
        const double               shift = 0 ) ;            // shift parameter 
      // ======================================================================
      Benini
      ( const unsigned short       n         ,              // number of shape parameters 
        const double               scale     ,              // scale parameter 
        const double               shift     ) ;            // shift parameter 
      // ======================================================================
    public:
      // ======================================================================
      // evaluate function 
      double evaluate   ( const double x ) const ;
      // evaluate function 
      double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      /// shift parameter 
      double shift  () const { return m_shift ; }
      /// scale parameter 
      double scale  () const { return m_scale ; }
      /// shape parameter
      double par    ( const unsigned short i ) const 
      { return i < m_pars.size() ? m_pars [ i ] : 0.0 ; }
      /// all shape parametgers 
      const std::vector<double>& pars ()       const { return m_pars ; }
      /// linear     term 
      double alpha  () const { return par ( 0 ) ; }
      /// quadrattic term 
      double beta   () const { return par ( 1 ) ; }
      /// cubic      term 
      double gamma  () const { return par ( 2 ) ; }
      /// quartic    term 
      double delta  () const { return par ( 3 ) ; }
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
      double xmin () const { return m_shift + m_scale ; }
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
      double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // evaluate function 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : // getters 
      // ======================================================================
      double mu    () const { return m_mu    ; }
      double scale () const { return m_scale ; }
      double shape () const { return m_shape ; }
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
      double xi         () const { return m_shape ; }
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
      inline double pdf        ( const double x ) const 
      { return evaluate ( x ) ; }
      /// get the value of MPERT distribution 
      inline double operator() ( const double x ) const 
      { return evaluate ( x ) ; }
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
    /** @class CutOffGauss 
     *  Useful function for smooth Gaussian cut-off:
     *  \f[ f(x;x_0;\sigma) = \left\{ 
     *    \begin{array}{ll}
     *    1  & \mathrm{for~} x \le x_0  \\
     *    \mathrm{e}^{-\frac{1}{2}\left( \frac{ (x-x_0)^2}{\sigma^2} \right)}
     *       & \mathrm{for~} x >   x_0 
     *    \end{array}\right. \f] 
     */
    class CutOffGauss 
    {
    public:
      // ======================================================================
      /** Constructor from all parameters
       *  @param right dump direction
       *  @param x0    threshold value 
       *  @param sigma sigma  
       */
      CutOffGauss 
      ( const bool   right = true , 
        const double x0    = 0    , 
        const double sigma = 1    ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// the main method 
      double operator() ( const double x ) const ;
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// right ?
      bool right   () const { return m_right ; }
      /// x_0 
      double x0    () const { return m_x0    ; }
      /// sigma
      double sigma () const { return m_sigma ; }      
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// update x_0
      bool setX0    ( const double value ) ;
      /// update sigma 
      bool setSigma ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
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
      bool   m_right ;
      double m_x0    ;
      double m_sigma ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CutOffStudent
     *  Useful function for smooth Student's t=-like (power-law) cut-off:
     *  \f[ f(x;x_0;\sigma) = \left\{ 
     *    \begin{array}{ll}
     *    1  & \mathrm{for~} x \le x_0  \\
     *    \left( \frac{1}{\nu} \left( \frac{(x-x_0)}{\sigma^2} \right)^{ - \frac{\nu+1}{2}} \right) 
     *       & \mathrm{for~} x >   x_0 
     *    \end{array}\right. \f] 
     */
    class CutOffStudent
    {
    public:
      // ======================================================================
      /** Constructor from all parameters
       *  @param right dump direction
       *  @param x0    threshold value 
       *  @param nu    parameter nu 
       *  @param sigma parameter sigma  
       */
      CutOffStudent 
      ( const bool   right = true , 
        const double x0    = 0    , 
        const double n     = 1    ,
        const double sigma = 1    ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// the main method 
      double operator() ( const double x ) const ;
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// right ?
      bool right   () const { return m_right ; }
      /// x_0 
      double x0    () const { return m_x0    ; }
      /// n 
      double nu    () const { return m_nu    ; }      
      /// sigma
      double sigma () const { return m_sigma ; }      
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// update x_0
      bool setX0    ( const double value ) ;
      /// update nu 
      bool setNu    ( const double value ) ;
      /// update sigma 
      bool setSigma ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
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
      bool   m_right ;
      double m_x0    ;
      double m_nu    ;
      double m_sigma ;
      double m_C     ;
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
