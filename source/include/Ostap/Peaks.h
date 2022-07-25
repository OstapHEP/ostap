// ============================================================================
#ifndef OSTAP_PEAKS_H 
#define OSTAP_PEAKS_H 1
// ============================================================================
//  Include  files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
// ============================================================================
/** @file Ostap/Peaks.h
 *  Collection of useful peak-like models
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
    /** @class BifurcatedGauss
     *  simple representation of bifurcated gaussian function
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class  BifurcatedGauss
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param peak    the peak posiion
       *  @param sigmaL  (left  sigma)
       *  @param sigmaR  (right sigma)
       */
      BifurcatedGauss
      ( const double peak   = 0 ,
        const double sigmaL = 1 ,
        const double sigmaR = 1 ) ;
      // ======================================================================
      /// destructor
      ~BifurcatedGauss() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Bifurcated Gaussian
      double evaluate   ( const double x ) const ;
      /// evaluate Bifurcated Gaussian
      double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Bifurcated Gaussian
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      double peak    () const { return m_peak    ; }
      /// peak position
      double m0      () const { return peak()    ; }
      /// left sigma
      double sigmaL  () const { return m_sigmaL  ; }
      /// right sigma
      double sigmaR  () const { return m_sigmaR  ; }
      // ======================================================================
      /// ssigma 
      double sigma   () const { return 0.5  * ( m_sigmaL + m_sigmaR )            ; }
      /// sigma-asymmetry 
      double asym    () const { return 0.5  * ( m_sigmaL - m_sigmaR ) / sigma () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position 
      bool setPeak    ( const double value ) ;
      /// peak position 
      bool setM0      ( const double value ) { return setPeak ( value ) ; }
      /// peak position 
      bool setMass    ( const double value ) { return setPeak ( value ) ; }
      /// left sigma 
      bool setSigmaL  ( const double value ) ;
      /// right sigma 
      bool setSigmaR  ( const double value ) ;
      // ======================================================================
    public: // integrals & CDF
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      ///  get CDF 
      double cdf      ( const double x    ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private: // parameters
      // ======================================================================
      /// the peak position
      double m_peak   ;       //                              the peak position
      /// sigma left
      double m_sigmaL ;       // sigma-left
      /// sigma right
      double m_sigmaR ;       // sigma-right
      // ======================================================================
    } ;
    // ========================================================================
    /** @class DoubleGauss
     *  simple representation of double gaussian function
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class  DoubleGauss 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param peak     the peak position
       *  @param sigma    the sigma for first component
       *  @param fraction the fraction of the first component 
       *  @param scale    the ratio of sigmas for second and first components
       */
      DoubleGauss
      ( const double peak     = 0   ,
        const double sigma    = 1   , 
        const double fraction = 0.9 , 
        const double scale    = 1.1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Bifurcated Gaussian
      double pdf        ( const double x ) const ;
      double operator() ( const double x ) const { return pdf ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      double peak     () const { return m_peak    ; }
      /// peak position
      double mean     () const { return peak ()   ; }
      /// peak position
      double m0       () const { return peak ()   ; }
      /// sigma 
      double sigma    () const { return m_sigma   ; }
      /// sigma-1
      double sigma1   () const { return m_sigma   ; }
      /// sigma-2
      double sigma2   () const { return m_sigma * m_scale ; }
      /// scale 
      double scale    () const { return m_scale    ; }
      /// fraction 
      double fraction () const { return m_fraction ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak positon
      bool setPeak     ( const double value ) ;
      /// peak positon
      bool setM0       ( const double value ) { return setPeak ( value ) ; }
      /// peak positon
      bool setMass     ( const double value ) { return setPeak ( value ) ; }
      /// sigma 
      bool setSigma    ( const double value ) ;
      /// scale 
      bool setScale    ( const double value ) ;
      /// fraction 
      bool setFraction ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the cdf 
      double cdf      ( const double x ) const ;
      /// get the integral
      double integral () const { return 1  ; }  ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private: // parameters
      // ======================================================================
      /// the peak position
      double m_peak     ;       //                              the peak position
      /// sigma
      double m_sigma    ;       // sigma
      /// fraction
      double m_fraction ;       // fraction
      /// scale
      double m_scale    ;       // scale
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  Gauss 
     *  trivial gaussian , just for completeness 
     */
    class Gauss
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param peak    the peak posiion
       *  @param sigmaL  (left  sigma)
       *  @param sigmaR  (right sigma)
       */
      Gauss
      ( const double peak  = 0 ,
        const double sigma = 1 ) ;
      // ======================================================================
      /// destructor
      ~Gauss() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Bifurcated Gaussian
      double evaluate   ( const double x ) const ;
      /// evaluate Bifurcated Gaussian
      double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Bifurcated Gaussian
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      double peak    () const { return m_peak    ; }
      /// peak position
      double m0      () const { return peak()    ; }
      /// sigma
      double sigma   () const { return m_sigma   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position 
      bool setPeak    ( const double value ) ;
      /// peak position 
      bool setM0      ( const double value ) { return setPeak ( value ) ; }
      /// peak position 
      bool setMass    ( const double value ) { return setPeak ( value ) ; }
      /// left sigma 
      bool setSigma   ( const double value ) ;
      // ======================================================================
    public: // integrals & CDF
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      ///  get CDF 
      double cdf      ( const double x    ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private: // parameters
      // ======================================================================
      /// the peak position
      double m_peak   ;       //                              the peak position
      /// sigma 
      double m_sigma  ;       // sigma
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenGaussV1
     *  Simple class that implements the generalized normal distribution v1
     *  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-08-25
     *  @see M. T. Subbotin, “On the Law of Frequency of Error”, Mat. Sb., 31:2 (1923), 296–301
     *  @see http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=sm&paperid=6854&option_lang=eng
     *  @see Nadarajah, Saralees (September 2005). "A generalized normal distribution".
     *       Journal of Applied Statistics. 32 (7): 685–694. doi:10.1080/02664760500079464.
     *  @see https://doi.org/10.1080%2F02664760500079464     
     */
    class  GenGaussV1
    {
    public:
      // ======================================================================
      /** constructor from all agruments
       *  @param mu     location/peak posiiton
       *  @param alpha  "scale" parameter
       *  @param beta   "shape" parameter
       */
      GenGaussV1
      ( const double mu    = 0 ,
        const double alpha = 1 ,
        const double beta  = 2 ) ; // beta=2 correponds to gaussian
      /// desctructor
      ~GenGaussV1() ;
      // ======================================================================
    public: // primary getters
      // ======================================================================
      /// peak position 
      double mu          () const { return m_mu       ; }
      /// peak position 
      double peak        () const { return   mu    () ; }
      /// peak position 
      double location    () const { return   mu    () ; }
      ///  alpha 
      double alpha       () const { return m_alpha    ; }
      ///  scale 
      double scale       () const { return   alpha () ; }
      ///   beta 
      double beta        () const { return m_beta     ; }
      ///  shape 
      double shape       () const { return   beta  () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      ///  set position 
      bool  setMu        ( const double value ) ;
      ///  alpha 
      bool  setAlpha     ( const double value ) ;
      ///  bneta 
      bool  setBeta      ( const double value ) ;
      //
      ///  set position 
      bool  setPeak      ( const double value ) { return setMu    ( value ) ; }
      ///  set position 
      bool  setLocation  ( const double value ) { return setMu    ( value ) ; }
      ///  scale 
      bool  setScale     ( const double value ) { return setAlpha ( value ) ; }
      ///  shape 
      bool  setShape     ( const double value ) { return setBeta  ( value ) ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      /// mean
      double mean        () const { return   mu    () ; }
      /// median
      double median      () const { return   mu    () ; }
      /// mode 
      double mode        () const { return   mu    () ; }
      /// variance 
      double variance    () const ;
      /// dispersion 
      double dispersion  () const { return variance () ; }
      /// sigma^2
      double sigma2      () const { return variance () ; }
      /// sigma 
      double sigma       () const ;
      /// skewness 
      double skewness    () const { return 0 ; }
      /// kurtosis 
      double kurtosis    () const ;
      // ======================================================================
    public :
      // ======================================================================
      /// get pdf
      double operator() ( const double x ) const { return pdf ( x ) ; }
      /// get pdf
      double pdf        ( const double x ) const ;
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
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// location 
      double m_mu     ;  // location
      ///  scale 
      double m_alpha  ;  // scale
      /// shape 
      double m_beta   ;  // shape
      ///  aux 
      double m_gbeta1 ;  // helper parameter
      ///  aux 
      double m_gbeta2 ;  // helper parameter
    } ;
    // ========================================================================
    /** @class GenGaussV2
     *  Simple class that implements the generalized normal distribution v2
     *  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_2
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-08-25
     */
    class  GenGaussV2
    {
    public:
      // ======================================================================
      /** constructor from all agruments
       *  @param xi     location/peak posiiton
       *  @param alpha  "scale" parameter
       *  @param kappa  "shape" parameter
       */
      GenGaussV2
      ( const double xi    = 0 ,
        const double alpha = 1 ,
        const double kappa = 0 ) ; // kappa=0 correponds to gaussian
      /// desctructor
      ~GenGaussV2() ;
      // ======================================================================
    public: // primary getters
      // ======================================================================
      double xi          () const { return m_xi       ; }
      double peak        () const { return   xi    () ; }
      double location    () const { return   xi    () ; }
      double alpha       () const { return m_alpha    ; }
      double scale       () const { return   alpha () ; }
      double kappa       () const { return m_kappa    ; }
      double shape       () const { return   kappa () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool  setXi        ( const double value ) ;
      bool  setAlpha     ( const double value ) ;
      bool  setKappa     ( const double value ) ;
      //
      bool  setPeak      ( const double value ) { return setXi    ( value ) ; }
      bool  setLocation  ( const double value ) { return setXi    ( value ) ; }
      bool  setScale     ( const double value ) { return setAlpha ( value ) ; }
      bool  setShape     ( const double value ) { return setKappa ( value ) ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      double mean        () const ;
      double median      () const { return   xi    () ; }
      //
      double variance    () const ;
      double dispersion  () const { return variance () ; }
      double sigma2      () const { return variance () ; }
      double sigma       () const ;
      //
      double skewness    () const ;
      double kurtosis    () const ;
      // ======================================================================
    public :
      // ======================================================================
      /// get pdf
      double operator() ( const double x ) const { return pdf ( x ) ; }
      double pdf        ( const double x ) const ;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      double cdf        ( const double x ) const ;
      /// get the integral
      double integral   () const { return 1 ; }
      /// get the integral between low and high limits
      double integral   ( const double low  ,
                          const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double  y ( const double x ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_xi      ;  // location
      double m_alpha   ;  // scale
      double m_kappa   ;  // shape
    } ;
    // ========================================================================
    /** @class SkewGauss
     *  Simple class that implements the skew normal distribution
     *  @see http://en.wikipedia.org/wiki/Skew_normal_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-08-25
     */
    class  SkewGauss
    {
    public:
      // ======================================================================
      /** constructor from all agruments
       *  @param xi     location/peak posiiton
       *  @param omega  "scale" parameter
       *  @param alpha  "shape" parameter
       */
      SkewGauss
      ( const double xi    = 0 ,
        const double omega = 1 ,
        const double alpha = 0 ) ; // alpha=0 correponds to gaussian
      /// desctructor
      ~SkewGauss () ;
      // ======================================================================
    public: // primary getters
      // ======================================================================
      double xi          () const { return m_xi       ; }
      double peak        () const { return   xi    () ; }
      double location    () const { return   xi    () ; }
      double omega       () const { return m_omega    ; }
      double scale       () const { return   omega () ; }
      double alpha       () const { return m_alpha    ; }
      double shape       () const { return   alpha () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool  setXi        ( const double value ) ;
      bool  setOmega     ( const double value ) ;
      bool  setAlpha     ( const double value ) ;
      //
      bool  setPeak      ( const double value ) { return setXi    ( value ) ; }
      bool  setLocation  ( const double value ) { return setXi    ( value ) ; }
      bool  setScale     ( const double value ) { return setOmega ( value ) ; }
      bool  setShape     ( const double value ) { return setAlpha ( value ) ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      double mean        () const ;
      double variance    () const ;
      double dispersion  () const { return variance () ; }
      double sigma2      () const { return variance () ; }
      double sigma       () const ;
      double skewness    () const ;
      double kurtosis    () const ;
      // ======================================================================
    public :
      // ======================================================================
      /// get pdf
      double operator() ( const double x ) const { return pdf ( x ) ; }
      double pdf        ( const double x ) const ;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      double cdf        ( const double x ) const ;
      /// get the integral
      double integral   () const { return 1 ; }
      /// get the integral between low and high limits
      double integral   ( const double low  ,
                          const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_xi     ;  // location
      double m_omega  ;  // scale
      double m_alpha  ;  // shape
      // =======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    /** @class ExGauss 
     *  Exponentially modified Gaussian function, EMG
     *  @see https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
     *
     *  It is a distibutiin for the varibale that is a 
     *  sum (or difference for negative \f$ k\f$) 
     *  of a Gaussian and exponential variables: \f$ X \sim Y + sign(k) Z \f$,  
     *  where 
     *  - \f$ Y \sim N(\mu,\sigma) \f$
     *  - \f$ Z \sim  \frac{1}{k\sigma}\mathrm{e}^{-\frac{x}{k\sigma}} \f$ 
     *  
     *  For \f$ k=0\f$ one gets a Gaussian distribution
     *  - \f$ k>0\f$ corresponds to the rigth tail  
     *  - \f$ kM0\f$ corresponds to the left tail  
     *
     *  It can be considered as "single-tail" version of the Normal Laplace distribution:
     *  - \f$ k = 0 \f$ corresponds to Gaussian distribution
     *  - \f$ k > 0 \f$ corresponds to Normal Laplace \f$ NL(\mu,\sigma,0,k)\f$ 
     *  - \f$ k < 0 \f$ corresponds to Normal Laplace \f$ NL(\mu,\sigma,\left|\tau\right|,0)\f$ 
     *
     *  @see Reed, W.J, "The Normal-Laplace Distribution and Its Relatives". 
     *       In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
     *       "Advances in Distribution Theory, Order Statistics, and Inference. 
     *       Statistics for Industry and Technology". Birkhäuser Boston. 
     *  @see https://doi.org/10.1007/0-8176-4487-3_4
     *  @see Ostap::Math::NormalLaplace 
     */
    class ExGauss 
    {
    public:
      // ======================================================================
      /// constructor from all parameters 
      ExGauss
      ( const double mu       = 0 , 
        const double varsigma = 1 , 
        const double k        = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      double evaluate          ( const double x ) const ;
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// parameter mu
      double mu       () const { return m_mu        ; }
      /// parameter varsigma
      double varsigma () const { return m_varsigma  ; }
      /// parameter k 
      double k        () const { return m_k         ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu       ( const double value ) ;
      bool setVarsigma ( const double value ) ;
      bool setK        ( const double value ) ;
      // ======================================================================      
    public:  // integrals
      // ======================================================================
      /// get CDF
      double cdf        ( const double x ) const ;
      /// get the integral
      double integral   () const ;
      /// get the integral between low and high limits
      double integral   ( const double low  ,
                          const double high ) const ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean value 
      double mean        () const ;
      /// variance 
      double variance    () const ;
      /// RMS 
      double rms         () const ;
      /// dispersion 
      double dispersion  () const { return variance () ; }
      /// skewness 
      double skewness    () const ;
      /// kurtosis 
      double kurtosis    () const ;      
      /// get cumulant 
      double cumulant    ( const unsigned short r ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// parameter mu 
      double m_mu       { 0 } ; // parameter mu 
      /// parameter varsigma
      double m_varsigma { 1 } ; // parameter varsigma 
      /// parameter k 
      double m_k        { 0 } ; // parameter k 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class NormalLaplace 
     *  Distribution for a sum of Gaussian and (asymmertric) Laplace variables 
     *  It behaves line core Gaussian with exponential tails 
     *  @see Wiliam J. Reed, "The Normal-Laplace Distribution Relatives", 
     *  October, 2004
     *  @see https://www.math.uvic.ca/faculty/reed/NL.draft.1.pdf
     *  @see Reed, W.J, "The Normal-Laplace Distribution and Its Relatives". 
     *       In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
     *       "Advances in Distribution Theory, Order Statistics, and Inference. 
     *       Statistics for Industry and Technology". Birkhäuser Boston. 
     *  @see https://doi.org/10.1007/0-8176-4487-3_4
     *
     *   \f$ f(x; \mu, \sigma, k_L , k_R ) = 
     *   \frac{1}{\sigma ( k_L + k_R) } 
     *   \phi ( z ) \left( R ( \frac{1}{k_R} - z ) + 
     *                     R ( \frac{1}{k_L} + z ) \right) 
     *   \f$, where
     *   - \f$ k_L,k_R \ge 0 \f$ 
     *   - \f$ z = \frac{x-\mu}{\sigma} \f$ 
     *   - \f$ \phi(z) \f$ is Gaussian PDF  
     *   - \f$  R(x)   \f$ is Mill's ratio 
     */    
    class NormalLaplace 
    {
    public : 
      // ======================================================================
      NormalLaplace 
      ( const double mu       = 0 ,
        const double varsigma = 1 ,
        const double kL       = 0 , 
        const double kR       = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      double evaluate          ( const double x ) const ;
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// parameter mu
      double mu       () const { return m_mu        ; }
      /// parameter varsigma
      double varsigma () const { return m_varsigma  ; }
      /// left  exponential 
      double kL       () const { return m_kL        ; }
      /// right exponential 
      double kR       () const { return m_kR        ; }
      // ======================================================================
    public: // original parameterisation 
      // ======================================================================
      /// parameter alpha 
      double alpha () const ;
      /// parametyer beta  
      double beta  () const  ;      
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu       ( const double value ) ;
      bool setVarsigma ( const double value ) ;
      bool setKL       ( const double value ) ;
      bool setKR       ( const double value ) ;
      // ======================================================================      
    public:  // integrals
      // ======================================================================
      /// get CDF
      double cdf        ( const double x ) const ;
      /// get the integral
      double integral   () const ;
      /// get the integral between low and high limits
      double integral   ( const double low  ,
                          const double high ) const ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean value 
      double mean        () const ;
      /// variance 
      double variance    () const ;
      /// RMS 
      double rms         () const ;
      /// dispersion 
      double dispersion  () const { return variance () ; }
      /// skewness 
      double skewness    () const ;
      /// (excess) kurtosis  
      double kurtosis    () const ;      
      /// get cumulant 
      double cumulant    ( const unsigned short r ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// parameter mu 
      double m_mu       { 0 } ; // parameter mu 
      /// parameter varsigma
      double m_varsigma { 1 } ; // parameter varsigma
      /// left exponential
      double m_kL       { 0 } ; // left exponential
      /// right exponential
      double m_kR       { 0 } ; // right exponential
      // ======================================================================
    };
    // ========================================================================
    /** @class Bukin
     *  ``Bukin-function'', aka "Modified Novosibirsk function"
     *  for description of asymmetric peaks with the exponential tails
     *
     *  @see http://arxiv.org/abs/1107.5751
     *  @see https://doi.org/10.1007/JHEP06(2012)141
     *  @date 2011-04-19
     */
    class  Bukin 
    {
    public :
      // ======================================================================
      /** constructor from all parameters
       *  @param peak  the peak posiion
       *  @param sigma the effective sigma, defined as FWHM/2.35
       *  @param xi    the asymmetry parameter
       *  @param rhoL  the left  tail parameter
       *  @param rhoR  the right tail parameter
       */
      Bukin
      ( const double peak   = 0 ,
        const double sigma  = 1 ,
        const double xi     = 0 ,
        const double rhoL   = 0 ,
        const double rhoR   = 0 ) ;
      // ======================================================================
      /// destructor
      ~Bukin () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Bukin's function
      double pdf        ( const double x ) const ;
      /// evaluate Bukin's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double peak  () const { return m_peak    ; }
      double m0    () const { return   peak () ; }
      double sigma () const { return m_sigma   ; }
      double xi    () const { return m_xi      ; }
      double rho_L () const { return m_rho_L   ; }
      double rho_R () const { return m_rho_R   ; }
      // ======================================================================
      double x1    () const { return m_x1      ; }
      double x2    () const { return m_x2      ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setPeak  ( const double value ) ;
      bool setM0    ( const double value ) { return setPeak ( value ) ; }
      bool setMass  ( const double value ) { return setPeak ( value ) ; }
      bool setSigma ( const double value ) ;
      bool setXi    ( const double value ) ;
      bool setRhoL  ( const double value ) ;
      bool setRhoR  ( const double value ) ;
      bool setRho_L ( const double value ) { return setRhoL ( value ) ; }
      bool setRho_R ( const double value ) { return setRhoR ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private: // parameters
      // ======================================================================
      /// the peak position
      double m_peak    ;      //                              the peak position
      /// the effective resolution, defined as FWHM/2.35
      double m_sigma   ;      // the effective resolution, defined as FWHM/2.35
      /// the asymmetry parameter
      double m_xi      ;      //                        the asymmetry parameter
      /// the left  tail parameter
      double m_rho_L   ;      //                        the left tail parameter
      /// the right tail parameter
      double m_rho_R   ;      //                        the right tail parameter
      // ======================================================================
    private: // internals
      // ======================================================================
      /// A/2 -region : left edge
      double m_x1        ;   // A/2 -region : left edge
      /// A/2 -region : right  edge
      double m_x2        ;   // A/2 -region : right  edge
      //
      /// the first magic constant for the central region
      double m_A         ;   // the first  magic constant for the central region
      /// the second magic constant for the central region
      double m_B2        ;   // the second magic constant for the central region
      //
      /// tails parameters (times  Bukin's  constants)
      double m_L         ;   // left  tail
      double m_R         ;   // right tail
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Novosibirsk
     *  ``Novosibirsk-function'' for description of gaussian with tails
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class  Novosibirsk
    {
    public :
      // ======================================================================
      /** constructor from all parameters
       *  @param m0    the peak posiion
       *  @param sigma the effective sigma
       *  @param tau   the tail paramter
       */
      Novosibirsk
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double tau   = 0 ) ;
      /// destructor
      ~Novosibirsk () ;                                           // destructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Novosibirsk's function
      double pdf        ( const double x ) const ;
      /// evaluate Novosibirsk's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double m0    () const { return m_m0       ; }
      double peak  () const { return   m0    () ; }
      double mass  () const { return   m0    () ; }
      double sigma () const { return m_sigma    ; }
      double tau   () const { return m_tau      ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0    ( const double value ) ;
      bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      bool setMass  ( const double value ) { return setM0 ( value ) ; }
      bool setSigma ( const double value ) ;
      bool setTau   ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private: // recalculate constants
      // ======================================================================
      /// recalculate integral
      void integrate () ; // recalculate integral
      /// get parameter lambda
      void getLambda () ; // get parameter lambda
      // ======================================================================
    private: // parameters
      // ======================================================================
      /// the peak position
      double m_m0      ;      //                              the peak position
      /// the effective resolution
      double m_sigma   ;      //                       the effective resolution
      /// the tail parameter
      double m_tau     ;      //                             the tail parameter
      // ======================================================================
    private: // internals
      // ======================================================================
      /// lambda value
      double m_lambda    ;   // lambda value
      // integration
      double m_integral  ;
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================
    // Crystal Ball & Co
    // ========================================================================
    /** @class CrystalBall
     *  ``Crystal Ball-function'' for description of gaussian with the tail
     *  @see http://en.wikipedia.org/wiki/Crystal_Ball_function
     *
     *  for \f$\alpha>0\f$
     *
     *  \f[ f(x;\alpha,n,x_0,\sigma) = \frac{1}{ \sqrt{2\pi\sigma^2} } \left\{
     *  \begin{array}{ll}
     *  \mathrm{e}^{-\frac{1}{2}\left(\frac{x-x_0}{\sigma}\right)^2}
     *  & \text{for}~\frac{x-x_0}\ge-\alpha\sigma \\
     *  \mathrm{- \frac{\alpha^2}{2}} \times
     *  \left(  \frac{n+1}{ n+1 - \alpha^2 - \left|\alpha\right|\frac{x-x_0}{\sigma}}\right)^{n+1}
     *  & \text{for}~\frac{x-x_0}\le-\alpha\sigma
     *  \end{array}
     *  \right.\f]
     *
     * where
     *
     * \f[ C = \frac{n+1}{\left|\alpha\right|\times \frac{1}{n} \times \mathrm{e}^{-\frac{\alpha^2}{2}}}  \f]
     * \f[ B = \sqrt{\frac{\pi}{2}}\left(1+\mathrm{erf}\left(-\frac{\alpha}{\sqrt{2}}\right)\right) \f]
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  CrystalBall
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0       parameter
       *  @param sigma  sigma    parameter
       *  @param alpha  alpha    parameter
       *  @param n      n        parameter (equal for N-1 for "standard" definition)
       */
      CrystalBall
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double alpha = 2 ,
        const double n     = 1 ) ;
      /// destructor
      ~CrystalBall () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      double m0    () const { return m_m0    ; }
      double peak  () const { return   m0 () ; }
      double sigma () const { return m_sigma ; }
      double alpha () const { return m_alpha ; }
      double n     () const { return m_n     ; }
      // ======================================================================
      double aa    () const { return std::abs ( m_alpha ) ; }
      double np1   () const { return n()  + 1 ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool setM0    ( const double value ) ;
      bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      bool setMass  ( const double value ) { return setPeak ( value ) ; }
      bool setSigma ( const double value ) ;
      bool setAlpha ( const double value ) ;
      bool setN     ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get (possibly truncated, if n==0 or alpha=0) integral
      double integral () const ;
      /// get the integral between low and high
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
      /// the peak position
      double m_m0       ;  // the peak position
      /// the peak resolution
      double m_sigma    ;  // the peak resolution
      /// parameter alpha
      double m_alpha    ;  // parameter alpha
      /// parameter n
      double m_n        ;  // parameter n
      // ======================================================================
    private :
      // ======================================================================
      /// helper constants
      double m_A        ;  // exp(-0.5*alpha^2)
      double m_B        ;  // integral over the gaussian part
      double m_C        ;  // integral over the power-law tail
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Needham
     *  The special parametrization by Matthew NEEDHAM of
     *  ``Crystal Ball-function'' suitable for \f$J/\psi/\Upsilon\f$-peaks     
     *  - thanks to Matthew Needham
     *
     *  Recommended constants for \f$J/psi\f$-peak:
     *    -  \f$a_0 =  1.975   \f$
     *    -  \f$a_1 =  0.0011  \f$
     *    -  \f$a_2 = -0.00018 \f$
     *
     *  Recommended constants for \f$\Upsilon\f$-peaks:
     *    -  \f$a_0 =  1.91    \f$
     *    -  \f$a_1 =  0.0017  \f$
     *    -  \f$a_2 = -5.22\times10^{-6} \f$
     *
     *  @see Ostap::Math::CrystalBall
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-05-13
     */
    class  Needham 
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0       parameter
       *  @param sigma  sigma    parameter
       *  @param a0     a0       parameter
       *  @param a1     a1       parameter
       *  @param a2     a2       parameter
       */
      Needham
      ( const double m0    = 3096.0     ,  // for J/psi
        const double sigma =   13.5     ,
        const double a0    =    1.975   ,
        const double a1    =    0.0011  ,
        const double a2    =   -0.00018 ) ;
      /// destructor
      ~Needham() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Needham's function
      double pdf        ( const double x ) const ;
      /// evaluate Needham's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      double m0    () const { return m_cb.m0    () ; }
      double peak  () const { return      m0    () ; }
      double sigma () const { return m_cb.sigma () ; }
      double a0    () const { return m_a0          ; }
      double a1    () const { return m_a1          ; }
      double a2    () const { return m_a2          ; }
      double alpha () const { return a0 () + sigma() * ( a1() + sigma() * a2 () ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool setM0    ( const double value ) { return m_cb.setM0    ( value ) ; }
      bool setPeak  ( const double value ) { return      setM0    ( value ) ; }
      bool setMass  ( const double value ) { return      setPeak  ( value ) ; }
      bool setSigma ( const double value ) { return m_cb.setSigma ( value ) ; }
      bool setA0    ( const double value ) ;
      bool setA1    ( const double value ) ;
      bool setA2    ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get (possibly truncated) integral
      double integral () const { return m_cb.integral() ; }
      /// get integral between low and high
      double integral ( const double low ,
                        const double high ) const
      { return m_cb.integral ( low , high ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the function itself
      Ostap::Math::CrystalBall m_cb ; // the function itself
      /// a0-parameter
      double m_a0       ;  // a0_parameter
      /// a0-parameter
      double m_a1       ;  // a1_parameter
      /// a0-parameter
      double m_a2       ;  // a2_parameter
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallRightSide
     *  ritgh-sided Crystal Ball function
     *  @see CrystalBall
     *  @date 2011-05-25
     */
    class  CrystalBallRightSide 
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0       parameter
       *  @param sigma  sigma    parameter
       *  @param alpha  alpha    parameter
       *  @param n      n        parameter (equal for N-1 for "standard" definition)
       */
      CrystalBallRightSide
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double alpha = 2 ,
        const double n     = 1 ) ;
      /// destructor
      ~CrystalBallRightSide () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      double m0    () const { return m_cb.m0    () ; }
      double peak  () const { return      m0    () ; }
      double sigma () const { return m_cb.sigma () ; }
      double alpha () const { return m_cb.alpha () ; }
      double n     () const { return m_cb.n     () ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool setM0    ( const double value ) { return m_cb.setM0    ( value ) ; }
      bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      bool setMass  ( const double value ) { return setPeak ( value ) ; }
      bool setSigma ( const double value ) { return m_cb.setSigma ( value ) ; }
      bool setAlpha ( const double value ) { return m_cb.setAlpha ( value ) ; }
      bool setN     ( const double value ) { return m_cb.setN     ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get (possibly truncated, if n==0 or alpha=0) integral
      double integral () const ;
      /// get the integral between low and high
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
      /// the actual CB-function:
      Ostap::Math::CrystalBall m_cb ;                 // the actual CB-function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallDoubleSided
     *  ``Crystal Ball-function'' for description of gaussian with the tail
     *  @see CrystalBall
     *  @see CrystalBallRightSide
     *  @date 2011-05-25
     */
    class  CrystalBallDoubleSided
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0      m0          parameter
       *  @param sigma   sigma       parameter
       *  @param alpha_L alpha_L     parameter
       *  @param n_L     n_L         parameter  (N-1 for "standard" definition)
       *  @param alpha_R alpha_R parameter
       *  @param n_R     n_R         parameter  (N-1 for "standard" definition)
       */
      CrystalBallDoubleSided
      ( const double m0      = 1 ,
        const double sigma   = 1 ,
        const double alpha_L = 2 ,
        const double n_L     = 1 ,
        const double alpha_R = 2 ,
        const double n_R     = 1 ) ;
      /// destructor
      ~CrystalBallDoubleSided() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      double m0      () const { return m_m0      ; }
      double peak    () const { return   m0 ()   ; }
      double sigma   () const { return m_sigma   ; }
      double alpha_L () const { return m_alpha_L ; }
      double n_L     () const { return m_n_L     ; }
      double alpha_R () const { return m_alpha_R ; }
      double n_R     () const { return m_n_R     ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool setM0      ( const double value ) ;
      bool setPeak    ( const double value ) { return setM0 ( value ) ; }
      bool setMass    ( const double value ) { return setPeak ( value ) ; }
      bool setSigma   ( const double value ) ;
      bool setAlpha_L ( const double value ) ;
      bool setN_L     ( const double value ) ;
      bool setAlpha_R ( const double value ) ;
      bool setN_R     ( const double value ) ;
      // ======================================================================
    public: //
      // ======================================================================
      /// get (possibly truncated) integral
      double integral () const ;
      /// get integral between low and high
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the peak position
      double m_m0       ;  // the peak position
      /// the peak resolution
      double m_sigma    ;  // the peak resolution
      /// parameter alpha
      double m_alpha_L  ;  // parameter alpha
      /// parameter N
      double m_n_L      ;  // parameter N
      /// parameter alpha_R
      double m_alpha_R  ;  // parameter alpha
      /// parameter N_R
      double m_n_R      ;  // parameter N
      // ======================================================================
    private:
      // ======================================================================
      /// helper constants
      double m_AL       ;  // exp(-0.5*alpha_L^2)
      double m_AR       ;  // exp(-0.5*alpha_R^2)
      double m_B        ;  // integral over the gaussian part
      double m_TL       ;  // integral over the left  power-law tail
      double m_TR       ;  // integral over the right power-law tail
      // ======================================================================
    } ;

    // ========================================================================
    /** @class Apollonios
     *  A modified gaussian with power-law tail on right side
     *  and an exponential tail on low-side
     *
     *  The function is proposed by Diego Martinez Santos
     *  @see https://indico.cern.ch/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=262633
     *  @see http://arxiv.org/abs/1312.5000
     *  Here a bit modified version is used with redefined parameter <code>n</code>
     *  to be coherent with local definitions of Crystal Ball
     *
     *  \f[ f(x;\alpha,n,x_0,\sigma) = \left\{
     *  \begin{array}{ll}
     *  \mathrm{e}^{-\left|b\right|\sqrt{1+(\delta x)^2}} & \text{for}~~\delta x >-a \\
     *  A \times \left( \frac{\left|n\right|+1}{ \left|n\right|+1 - \frac{(a+\delta x)\left|ab\right|}
     *  {\sqrt{1+a^2}} } \right)^{ \left|n\right|+1} & \text{otherwise}
     *  \end{array}
     *  \right. \f]
     *
     * where
     *
     * \f[ \delta x  = \frac{ x - x_0}{\left|\sigma\right|} \f]
     * \f[ A = \mathrm{e}^{-\left|b\right|\sqrt{1+a^2}}     \f]
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date  2013-12-01
     */
    class  Apollonios
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0       parameter
       *  @param sigma  sigma    parameter
       *  @param alpha  alpha    parameter
       *  @param n      n        parameter (equal for N-1 for "standard" definition)
       *  @param b      b        parameter
       */
      Apollonios
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double alpha = 2 ,
        const double n     = 1 ,
        const double b     = 1 ) ;
      /// destructor
      ~Apollonios () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Apollonios's function
      double pdf        ( const double x ) const ;
      /// evaluate Apollonios's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      double m0    () const { return m_m0     ; }
      double peak  () const { return   m0 ()  ; }
      double sigma () const { return m_sigma  ; }
      double alpha () const { return m_alpha  ; }
      double n     () const { return m_n      ; }
      double b     () const { return m_b      ; }
      // ======================================================================
      double a1    () const { return std::sqrt ( 1 + alpha() * alpha() ) ; }
      double aa    () const { return std::abs ( alpha() * b() ) / a1 ()  ; }
      double np1   () const { return n()  + 1 ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool setM0    ( const double value ) ;
      bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      bool setMass  ( const double value ) { return setPeak ( value ) ; }
      bool setSigma ( const double value ) ;
      bool setAlpha ( const double value ) ;
      bool setN     ( const double value ) ;
      bool setB     ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
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
      /// the peak position
      double m_m0       ;  // the peak position
      /// the peak resolution
      double m_sigma    ;  // the peak resolution
      /// parameter alpha
      double m_alpha    ;  // parameter alpha
      /// parameter n
      double m_n        ;  // parameter
      /// parameter b
      double m_b        ;  // parameter n
      /// helper constants
      double m_A        ;  // exp(-0.5*alpha^2)
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Apollonios2
     *  "Bifurcated Apollonios"
     *  A modified gaussian with asymmetric exponential tails on both sides
     *
     *  A convinient reparameterization is applied to keep reduce
     *  the correlations between "sigma"s and "beta"
     *
     *  \f[ f(x;\mu,\sigma_l,\sigma_r,\beta) \propto
     *  \mathrm{e}^{\left|\beta\right|( \left|\beta\right| - \sqrt{ \beta^2+\left(\delta x\right)^2}}
     *  \f]
     *
     * where
     *
     * \f[ \delta x  = \left\{ \begin{array}{ccc}
     *     \frac{x-\mu}{\sigma_l} &  for  & x \le \mu \\
     *     \frac{x-\mu}{\sigma_r} &  for  & x \ge \mu \\
     *     \end{array}
     *     \right.\f]
     *
     *  Large betas corresponds to gaussian
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date  2013-12-01
     */
    class  Apollonios2 
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0      m0        parameter
       *  @param sigmaL  sigmaL    parameter
       *  @param alphaR  alphaR    parameter
       *  @param beta    beta      parameter
       */
      Apollonios2
        ( const double m0      = 0   ,
          const double sigmaL  = 1   ,
          const double alphaR  = 1   ,
          const double beta    = 100 ) ;  // large beta correponds to gaussian
      /// destructor
      ~Apollonios2 () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Apollonios2's function
      double pdf        ( const double x ) const ;
      /// evaluate Apollonios2's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      double m0     () const { return m_m0     ; }
      double peak   () const { return   m0 ()  ; }
      double sigmaL () const { return m_sigmaL ; }
      double sigmaR () const { return m_sigmaR ; }
      double beta   () const { return m_beta   ; }
      // ======================================================================
      double sigma  () const { return 0.5 * ( m_sigmaL + m_sigmaR )           ; }
      double asym   () const { return 0.5 * ( m_sigmaL - m_sigmaR ) / sigma() ; }
      double b2     () const { return m_beta * m_beta ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool setM0     ( const double value ) ;
      bool setPeak   ( const double value ) { return setM0 ( value ) ; }
      bool setMass   ( const double value ) { return setPeak ( value ) ; }
      bool setSigmaL ( const double value ) ;
      bool setSigmaR ( const double value ) ;
      bool setBeta   ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
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
      /// the peak position
      double m_m0       ;  // the peak position
      /// the peak resolution
      double m_sigmaL  ;  // the peak resolution
      /// the peak resolution
      double m_sigmaR  ;  // the peak resolution
      /// parameter beta
      double m_beta    ;  // parameter beta
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    /** @class StudentT
     *  simple function to parameterize the symmetric peak using
     *  Student's ditribution
     *
     *  \f[  f(y) = \frac{1}{\sqrt{\pi n}} \frac { \Gamma( \frac{n+1}{2}) } { \Gamma( \frac{n}{2}  ) }
     *  \left( 1 + \frac{y^2}{n} \right)^{ -\frac{n+1}{2}} \f],
     *  where \f$ y = \frac{x - \mu}{\sigma} \f$
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-01-05
     */
    class  StudentT
    {
    public:
      // ======================================================================
      /** constructor from mass, resolution and "n"-parameter
       *  @param mass  mass
       *  @param sigma width parameter
       *  @param n     n-parameter  ( actually  n=1+|N| )
       */
      StudentT 
      ( const double mass  = 0 ,
        const double sigma = 1 ,
        const double n     = 2 ) ;
      /// destructor
      ~StudentT() ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate StudentT's shape
      double operator() ( const double x ) const{ return pdf ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      // variables
      // ======================================================================
      double M      () const  { return m_M      ; }
      double m0     () const  { return   M   () ; }
      double mass   () const  { return   M   () ; }
      double peak   () const  { return   M   () ; }
      // ======================================================================
      double sigma  () const  { return m_s      ; }
      double s      () const  { return sigma () ; }
      double gamma  () const  { return sigma () ; }
      double width  () const  { return sigma () ; }
      // ======================================================================
      double nu     () const  { return m_n      ; }
      double n      () const  { return m_n      ; }
      // ======================================================================
      bool setM     ( const double value  ) ;
      bool setM0    ( const double value  ) { return setM  ( value ) ; }
      bool setMass  ( const double value  ) { return setM  ( value ) ; }
      bool setPeak  ( const double value  ) { return setM  ( value ) ; }
      // ======================================================================
      bool setSigma ( const double value  ) ;
      bool setS     ( const double value  ) { return setSigma ( value ) ; }
      bool setGamma ( const double value  ) { return setSigma ( value ) ; }
      bool setWidth ( const double value  ) { return setSigma ( value ) ; }
      // ======================================================================
      bool setN     ( const double value  ) ;
      // ======================================================================
    public:
      // ======================================================================
      double pdf    ( const double x ) const ;
      double cdf    ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// mass
      double m_M  ; //
      /// width parameter
      double m_s ; // width parameter
      /// n-parameter
      double m_n ; // n-parameter
      // ======================================================================
    private: // normalization
      // ======================================================================
      double m_norm  ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BifurcatedStudentT
     *  simple function to parameterize the asymmetric peak using
     *  Student's ditribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-01-05
     */
    class  BifurcatedStudentT
    {
    public:
      // ======================================================================
      /** constructor from mass, resolution and "n"-parameter
       *  @param mass   mass
       *  @param sigmaL left width parameter
       *  @param sigmaR right width parameter
       *  @param nL     left n-parameter  ( actually  n=1+|N| )
       *  @param nR     right n-parameter  ( actually  n=1+|N| )
       */
      BifurcatedStudentT ( const double mass   = 0 ,
                           const double sigmaL = 1 ,
                           const double sigmaR = 1 ,
                           const double nL     = 2 ,
                           const double nR     = 2 ) ;
      /// destructor
      ~BifurcatedStudentT() ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate bifurcated StudentT's shape
      double operator() ( const double x ) const{ return pdf ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      // variables
      // ======================================================================
      double M       () const  { return m_M      ; }
      double m0      () const  { return   M   () ; }
      double mass    () const  { return   M   () ; }
      double peak    () const  { return   M   () ; }
      // ======================================================================
      double sigmaL  () const  { return m_sL      ; }
      double sL      () const  { return sigmaL () ; }
      double gammaL  () const  { return sigmaL () ; }
      double widthL  () const  { return sigmaL () ; }
      // ======================================================================
      double sigmaR  () const  { return m_sR      ; }
      double sR      () const  { return sigmaR () ; }
      double gammaR  () const  { return sigmaR () ; }
      double widthR  () const  { return sigmaR () ; }
      // ======================================================================
      double nuL     () const  { return m_nL      ; }
      double nL      () const  { return m_nL      ; }
      // =========== ===========================================================
      double nuR     () const  { return m_nR      ; }
      double nR      () const  { return m_nR      ; }
      // ======================================================================
      bool setM      ( const double value  ) ;
      bool setM0     ( const double value  ) { return setM  ( value ) ; }
      bool setMass   ( const double value  ) { return setM  ( value ) ; }
      bool setPeak   ( const double value  ) { return setM  ( value ) ; }
      // ======================================================================
      bool setSigmaL ( const double value  ) ;
      bool setSL     ( const double value  ) { return setSigmaL ( value ) ; }
      bool setGammaL ( const double value  ) { return setSigmaL ( value ) ; }
      bool setWidthL ( const double value  ) { return setSigmaL ( value ) ; }
      // ======================================================================
      bool setSigmaR ( const double value  ) ;
      bool setSR     ( const double value  ) { return setSigmaR ( value ) ; }
      bool setGammaR ( const double value  ) { return setSigmaR ( value ) ; }
      bool setWidthR ( const double value  ) { return setSigmaR ( value ) ; }
      // ======================================================================
      bool setNL     ( const double value  ) ;
      bool setNR     ( const double value  ) ;
      // ======================================================================
    public:
      // ======================================================================
      double pdf    ( const double x ) const ;
      double cdf    ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// mass
      double m_M  ; //
      /// width parameter
      double m_sL ; // width parameter
      /// width parameter
      double m_sR ; // width parameter
      /// nL-parameter
      double m_nL ; // n-parameter
      /// nR-parameter
      double m_nR ; // n-parameter
      // ======================================================================
    private: // normalization
      // ======================================================================
      double m_normL  ;
      double m_normR  ;
      // ======================================================================
    } ;
    // ========================================================================
    
    // ========================================================================
    /**  @class PearsonIV 
     *   Pearson Type IV distribution  
     *   \f$ f(x;\mu, n, \kappa) = 
     *   C \left( 1 + y^{2}\right)^{-(\frac{1}{2}+n)}
     *   \mathrm{e}^{ -\kappa \atan y }}\f$, where 
     *   - \f$  y = \frac{x-\mu}{\sigma}\f$,
     *   - \f$ 0 < n \f$  
     *  @see https://en.wikipedia.org/wiki/Pearson_distribution
     *  For $\kappa=0\f$ one gets Student's t-distribution
     *  @see J. Heinrich, "A guide to the Pearson Type IV distribution", 
     *       CDF/MEMO/STATISTICS/PUBLIC/6820, 2004 
     *  @see http://www-cdf.fnal.gov/physics/statistics/notes/cdf6820_pearson4.pdf
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2033-07-10
     */     
    class PearsonIV 
    {
    public: 
      // ========================================================================
      /** constructor from all parameters 
       *  @param mu    location parameter 
       *  @param sigma width/scale parameter 
       *  @param n     n-parameter 
       *  @param kappa asymmetry parameter 
       */
      PearsonIV
      ( const double mu    = 0 , 
        const double sigma = 1 , 
        const double n     = 2 , 
        const double kappa = 0 );
      // ======================================================================
    public:
      // ======================================================================
      /// get value of the function 
      double evaluate           ( const double x ) const ;
      /// get value of the function 
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      /// get value of the function 
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================      
    public: // getters 
      // ======================================================================
      /// location parameter 
      double mu       () const { return m_mu       ; } ;
      /// width/scale parameter 
      double varsigma () const { return m_varsigma ; } ;
      /// n-parameter 
      double n        () const { return m_n        ; } ;
      /// 
      double kappa    () const { return m_kappa    ; } ;
      // ======================================================================      
    public : // derived parameters 
      // ======================================================================      
      /// parameteter m
      inline double m  () const { return m_n + 0.5       ; }
      /// parameteter nu
      inline double nu () const { return m_kappa         ; }
      /// parameter r 
      inline double r  () const { return 2 * ( m () -1 ) ; }
      /// parameter a  
      inline double a  () const { return  m_varsigma     ; }
      // ======================================================================      
    public: // setters 
      // ======================================================================
      bool setMu       ( const double value ) ;
      bool setVarsigma ( const double value ) ;
      bool setN        ( const double value ) ;
      bool setKappa    ( const double value ) ;
      // ======================================================================      
    public: // properties  
      // ======================================================================
      /// mode 
      double mode            () const ; // mode of the distribution 
      /// mean value of the distribution  (for m>1)
      double mean            () const ; // mean value of the distribution 
      /// variance                        (for m>1.5)
      double variance        () const ; // variance of the distribution 
      /// RMS  
      double rms             () const ; // RMS  
      /// skewness (for m>2) 
      double skewness        () const ; // skewness 
      /// (excess) kurtosis  (for m > 5/2)
      double kurtosis        () const ; // (excess) kurtosis 
      /// (central) moment 
      double moment          ( const unsigned short k ) const ;
      /// beta1 parameter of Pearson family (m>2) 
      double beta1           () const ; // beta1 parameter of Pearson family 
      /// beta2 parameter of Pearson family (m>5/2)
      double beta2           () const ; // beta2 parameter of Pearson family 
      /** distance between two infection points:
       *  distance between two points with \f$ f^{\prime\prime}=0\f$.
       *  the twp points are equidstance fro mthe mode 
       */
      double infection_width () const ;
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
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// location parameter 
      double   m_mu       {  0 } ; // location parameter 
      /// width/scale  parameter 
      double   m_varsigma {  1 } ; // width/scale parameter 
      /// n-parameter 
      double   m_n        {  1 } ; // n-parameter 
      /// asymmetry parameter 
      double   m_kappa    {  0 } ; // asymmetry parameter
      /// normalization 
      double   m_C        { -1 } ; // normalization factor 
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================        
    } ;
    
      



    // ========================================================================
    /** @class SinhAsinh
     *
     *  Jones, M. C.; Pewsey, A. (2009).
     *  "Sinh-arcsinh distributions". Biometrika 96 (4): 761.
     *  doi:10.1093/biomet/asp053
     *  http://oro.open.ac.uk/22510
     *
     *  Location & scale  parameters are the
     *  usual representation of the family of
     *  distributions
     *  - \f$\epsilon\f$ parameter control the skewness
     *  - \f$\delta\f$   parameter control the kurtosis
     *  Normal distribtion reappears as \f$\epsilon=0\f$
     *  and \f$\delta=1\f$
     *  The heavy tails correspond to \f$\delta<1\f$,
     *  light tails correpond to \f$\delta>1\f$
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-08-02
     */
    class  SinhAsinh
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param location \f$\mu\f$-parameter       \f$-\infty<\mu<+\infty\f$
       *  @param scale    \f$\sigma\f$-parameter    \f$0<\sigma\f$
       *  @param epsilon  \f$\epsilon\f$-parameter  \f$-\infty<\epsilon<+\infty\f$
       *  @param delta    \f$\delta\f$-parameter    \f$0<\epsilon<+\infty\f$
       */
      SinhAsinh  ( const double location  = 1   ,
                   const double scale     = 1   ,
                   const double epsilon   = 0   ,
                   const double delta     = 1   ) ;
      /// destructor
      ~SinhAsinh() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate sinhasinh-distributions
      double pdf        ( const double x ) const ;
      /// evaluate sinhasinh-distributions
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double location () const { return mu    () ; }
      double scale    () const { return sigma () ; }
      // ======================================================================
      double mu       () const { return m_mu      ; }
      double sigma    () const { return m_sigma   ; }
      double epsilon  () const { return m_epsilon ; }
      double delta    () const { return m_delta   ; }
      // ======================================================================
      double peak     () const { return m_mu      ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setLocation ( const double value ) { return setMu    ( value ) ; }
      bool setScale    ( const double value ) { return setSigma ( value ) ; }
      bool setMu       ( const double value ) ;
      bool setSigma    ( const double value ) ;
      bool setEpsilon  ( const double value ) ;
      bool setDelta    ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double cdf      ( const double x    ) const ;
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mu      ;
      double m_sigma   ;
      double m_epsilon ;
      double m_delta   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class JohnsonSU
     *
     *  Johnson, N. L. (1949)
     *  "Systems of frequency curves generated by methods of translation"
     *  Biometrika 36: 149–176 JSTOR 2332539
     *  @see https://en.wikipedia.org/wiki/Johnson_SU_distribution
     *
     *  When variable \f$x\f$ follows Johnson-SU distribution,
     *  the variable
     *  \f$ z = \gamma + \delta \sinh^{-1}\frac{ x - \xi}{\lambda} \f$
     *  follows normal distribtion with mean 0 and sigma 1.
     *
     *  Note:
     *  Symmetric case of JonhsonSU distribution is
     *  recovered by \f$\delta\rightarrow0\f$ for
     *  "sinh-asinh" distribution, see
     *  Jones, M. C.; Pewsey, A. (2009).
     *  "Sinh-arcsinh distributions". Biometrika 96 (4): 761.
     *  doi:10.1093/biomet/asp053
     *  http://oro.open.ac.uk/22510
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2015-07-11
     */
    class  JohnsonSU
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param xi     \f$\xi\f$-parameter       \f$-\infty<\xi<+\infty\f$
       *  @param lambda \f$\lambda\f$-parameter   \f$   0<\lambda<+\infty\f$
       *  @param delta  \f$\delta\f$-parameter    \f$   0<\delta<+\infty\f$
       *  @param gamma  \f$\gamma\f$-parameter    \f$-\infty<\epsilon<+\infty\f$
       */
      JohnsonSU  ( const double xi      = 0 ,   // related to location
                   const double lambda  = 1 ,   // related to variance
                   const double delta   = 1 ,   // shape
                   const double gamma   = 0 ) ; // shape
      /// destructor
      ~JohnsonSU () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate JohnsonSU-distributions
      double pdf        ( const double x ) const ;
      /// evaluate JohnsonSU-distributions
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double xi       () const { return m_xi       ; }
      double lam      () const { return m_lambda   ; }
      double lambda   () const { return m_lambda   ; }
      double delta    () const { return m_delta    ; }
      double gamma    () const { return m_gamma    ; }
      // ======================================================================
    public:
      // ======================================================================
      double mean       () const ;
      double variance   () const ;
      double dispersion () const { return variance () ; }
      double sigma      () const { return std::sqrt ( variance() ) ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setXi       ( const double value ) ;
      bool setLambda   ( const double value ) ;
      bool setDelta    ( const double value ) ;
      bool setGamma    ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double cdf      ( const double x    ) const ;
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_xi      ;
      double m_lambda  ;
      double m_delta   ;
      double m_gamma   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Atlas
     *  Modified gaussian function
     *  \f$  f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\delta x/2}}}{2})\f$,
     *  where \f$\delta x = \left| x - \mu \right|/\sigma\f$
     *  Function is taken from http://arxiv.org/abs/arXiv:1507.07099
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-08-21
     */
    class  Atlas 
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param mean  \f$\mu\f$-parameter
       *  @param sigma \f$\sigma\f$-parameter
       */
      Atlas   ( const double mean   = 0  ,
                const double sigma  = 1  ) ;
      /// destructor
      ~Atlas () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate atlas function
      double pdf        ( const double x ) const ;
      /// evaluate atlas function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double mean   () const { return m_mean   ; }
      /// get parameters "sigma"
      double sigma  () const { return m_sigma  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      double peak     () const { return mean () ; }
      /// get mode
      double mode     () const { return mean () ; }
      /// get variance:  good numerical approximation
      double variance () const ;
      /// get rms :  good numerical approximation
      double rms      () const ;
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setMean  ( const double value ) ;
      bool   setSigma ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ( const double low  ,
                        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// parameteter "mu", mean, mode
      double m_mean  ; // parameter mu , mean , mode
      /// parameter   "sigma"
      double m_sigma ; // parameter sigma
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Sech
     *  Hyperbolic secant distribution or "inverse-cosh" distribution
     *
     *  The hyperbolic secant distribution shares many properties with the
     *  standard normal distribution:
     *  - it is symmetric with unit variance and zero mean,
     *    median and mode
     *  -its pdf is proportional to its characteristic function.
     *
     *  However, the hyperbolic secant distribution is leptokurtic;
     *  that is, it has a more acute peak near its mean, and heavier tails,
     *  compared with the standard normal distribution.
     *
     *  \f$ f(x,\mu,\sigma) \propto \frac{1}{2} \mathrm{sech} ( \frac{\pi}{2}\frac{x-\mu}{\sigma} )\f$
     *  @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-04-25
     */
    class  Sech 
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param mean  \f$\mu\f$-parameter
       *  @param sigma \f$\sigma\f$-parameter
       */
      Sech   ( const double mean   = 0  ,
               const double sigma  = 1  ) ;
      /// destructor
      ~Sech () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate sech function
      double pdf        ( const double x ) const ;
      /// evaluate sech function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double mean   () const { return m_mean    ; }
      double peak   () const { return   mean () ; }
      /// get parameters "sigma"
      double sigma  () const { return m_sigma   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get mode
      double mode     () const { return mean()            ; }
      /// get variance
      double variance () const { return m_sigma * m_sigma ; }
      /// get rms
      double rms      () const { return m_sigma           ; }
      /// get skewness
      double skewness () const { return 0 ; }
      /// get kurtosis
      double kurtosis () const { return 2 ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setMean  ( const double value ) ;
      bool   setSigma ( const double value ) ;
      // ======================================================================
    public: // quantile (0<p<1)
      // ======================================================================
      /// get quantile (0<p<1)
      double quantile ( const double p ) const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get integral from low to high
      double integral ( const double low  ,
                        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const ;
      /// evaluate atlas function
      double cdf      ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// parameteter "mu", mean, mode
      double m_mean  ; // parameter mu,mean,mode
      /// parameter   "sigma"
      double m_sigma ; // parameter sigma
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Logistic
     *  aka "Sech-square"
     *  \f$ f(x;\mu;s) = \frac{1}{4s}sech^2\left(\frac{x-\mu}{2s}\right)\f$,
     *  where
     *  \f$  s = \sigma \frac{\sqrt{3}}{\pi}\f$
     *  @see https://en.wikipedia.org/wiki/Logistic_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-14
     */
    class  Logistic
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param mean  \f$\mu  \f$-parameter
       *  @param sigma \f$sigma\f$-parameter
       */
      Logistic  ( const double mean  = 0  ,
                  const double sigma = 1  ) ;
      /// destructor
      ~Logistic () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate sech function
      double pdf        ( const double x ) const ;
      /// evaluate sech function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double mean   () const { return m_mean   ; }
      /// get parameters "sigma"
      double sigma  () const { return m_sigma  ; }
      /// get "s"
      double s      () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get peak
      double peak     () const { return peak () ; }
      /// get mode
      double mode     () const { return mean () ; }
      /// get median
      double median   () const { return mean () ; }
      /// get variance
      double variance () const { return m_sigma * m_sigma ; }
      /// get rms
      double rms      () const { return m_sigma ; }
      /// get skewness
      double skewness () const { return   0 ; }
      /// get kurtosis
      double kurtosis () const { return 1.2 ; }
      // ======================================================================
    public: // quanilies
      // ======================================================================
      /// quantile fnuiction  (0<p<1)
      double quantile ( const double p ) const ;
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setMean  ( const double value ) ;
      bool   setSigma ( const double value ) ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get integral from low to high
      double integral ( const double low  ,
                        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const ;
      /// evaluate Logistc CDF function
      double cdf      ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// parameteter "mu", mean, mode
      double m_mean  ; // parameter mu,mean,mode
      /// parameter   "sigma"
      double m_sigma ; // parameter sigma
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Losev 
     *  ``Losev distribution'' - asymmetric variant of Hyperbolic secant/Sech-function
     *  \f[ f(x;\mu,\alpha,\beta) \equiv 
     *   \frac{A}{\mathrm{e}^{-\left|\alpha\right| (x-\mu)} + 
     *                         \mathrm{e}^{\left|\beta\right|(x-mu)}}, \f]
     *  where \f$ A = \frac{\left|\alpha\right|+\left|\beta\right|}{\pi}.
     *  \sin \frac{\pi\left| \beta\right| }{\left|\alpha\right|+\left|\beta\right|}\f$. 
     *  - Leptokurtic distribution with exponential tails 
     *  @see Losev, A., "A new lineshape for fitting x‐ray photoelectron peaks", 
     *           Surf. Interface Anal., 14: 845-849. doi:10.1002/sia.740141207
     *  @see  https://doi.org/10.1002/sia.740141207
     *  @see  https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
     */
    class Losev 
    {
    public:
      // ======================================================================
      /** constructor from positive parameters alpha and beta 
       *  @param mean  \f$\mu\f$-parameter 
       *  @param alpha \f$\alpha\f$-parameter 
       *  @param beta  \f$\beta\f$-parameter 
       */ 
      Losev ( const double mu    = 0 , 
              const double alpha = 1 , 
              const double beta  = 1 ) ;        
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double operator() ( const double x ) const { return pdf ( x ) ; }
      /// evaluate the function 
      double pdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// parameter mu 
      double mu    () const { return m_mu    ; }
      /// parameter alpha
      double alpha () const { return m_alpha ; }
      /// parameter beta 
      double beta  () const { return m_beta  ; }      
      // ======================================================================
    public:
      // ======================================================================
      bool setMu    ( const double mu ) ;
      bool setAlpha ( const double mu ) ;
      bool setBeta  ( const double mu ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// the mode of the distribution 
      double mode () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================      
    public :
      // ======================================================================
      /// get the integral \f$ \int_{-\infty}^{+\infty} f(x) dx \f$
      double integral () const { return 1 ; }
      /** get the integral between low and high values 
       *  \f$ \int_{low}^{high}f(x) dx\f$
       */
      double integral ( const double low  , 
                        const double high ) const ;
      // ======================================================================      
    private :
      // ======================================================================
      /// parameteter "mu"
      double m_mu     { 0 }  ; // parameter
      /// left exponent 
      double  m_alpha { 1 } ; // left exponent 
      /// right exponent 
      double  m_beta  { 1 } ; // right exponent 
      // ======================================================================
    private:
      // ======================================================================
      /// normalization 
      mutable double m_norm { -1 } ; // normalization 
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class  Slash 
     *  ``Slash''-distribution -  symmetric peak with very heavy tail
     *  @see https://en.wikipedia.org/wiki/Slash_distribution
     *  Tails arew so heavy that moments (e.g. variance) do not exist 
     */
    class Slash 
    {
    public :
      // ======================================================================
      /** Constructor from location and mean 
       *  @param mu location 
       *  @param scale the scale, scale>0
       */
      Slash ( const double mu    = 0 ,   // location 
              const double scale = 1 ) ; // scale ;      
      /// destructor
      ~Slash() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate slash function
      double pdf        ( const double x ) const ;
      /// evaluate slash function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double mu       () const { return m_mu    ; }
      /// get parameters "scale"
      double scale    () const { return m_scale ; }
      // ======================================================================
    public: //   derived getters  
      // ======================================================================
      double mean     () const { return mu ()   ; }
      double peak     () const { return mu ()   ; }
      double location () const { return mu ()   ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool setMu       ( const double value ) ;
      bool setScale    ( const double value ) ;
      // ======================================================================      
      bool setMean     ( const double value ) { return setMu ( value ) ; }
      bool setLocation ( const double value ) { return setMu ( value ) ; }
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get integral from low to high
      double integral ( const double low  ,
                        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const { return 1 ; }
      /// evaluate sslash CDF function
      double cdf      ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      ///  peak location
      double m_mu    ; //  peak location
      /// the scale
      double m_scale ; // the scale 
      // ======================================================================      
    } ;  
    // ========================================================================
    /** @class AsymmetricLaplace
     *  @see https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution
     *  Here we use the ``inversed'' slopes
     *  \f$ f(x) \propto e^{ \pm\frac {x-\mu} { \lambda_{L,R}}} \f$ 
     */ 
    class AsymmetricLaplace
    {
    public:
      // ======================================================================
      /** constructor from all parameters 
       *  @param mu  peak location
       *  @param lambdaL ``left''  exponential slope  (lambdaL>0)
       *  @param lambdaR ``right'' exponential slope  (lambdaR>0)
       */
      AsymmetricLaplace ( const double mu      = 0 ,   // location 
                          const double lambdaL = 1 ,   // left  exponential slope 
                          const double lambdaR = 1 ) ; // right exponential slope 
      ///  destructor 
      ~AsymmetricLaplace() ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate asymmetic laplace function
      double pdf        ( const double x ) const ;
      /// evaluate asymmetric laplace function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      double mu       () const { return m_mu      ; }
      /// get ``left'' lambda parameter 
      double lambdaL  () const { return m_lambdaL ; }
      /// get ``right'' lambda parameter 
      double lambdaR  () const { return m_lambdaR ; }
      // ======================================================================
    public: //   derived getters  
      // ======================================================================
      double mean     () const { return mu ()   ; }
      double peak     () const { return mu ()   ; }
      double location () const { return mu ()   ; }
      // ======================================================================
    public: // the standard parameterization (slopes are inverse)
      // ======================================================================
      double lambda   () const { return 1.0 / std::sqrt ( m_lambdaL * m_lambdaR ) ; }
      /// get the ``asymmetry''   0<k<+inf 
      double k        () const { return       std::sqrt ( m_lambdaR / m_lambdaL ) ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setMu       ( const double value ) ;
      bool   setLambdaL  ( const double value ) ;
      bool   setLambdaR  ( const double value ) ;
      // ======================================================================      
      bool   setMean     ( const double value ) { return setMu ( value ) ; }
      bool   setLocation ( const double value ) { return setMu ( value ) ; }
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get integral from low to high
      double integral ( const double low  ,
                        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const {  return  1 ; }
      /// evaluate CDF function
      double cdf      ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      ///  peak location
      double m_mu      ; //  peak location
      /// ``left''  exponential slope 
      double m_lambdaL ; /// ``left''  exponential slope 
      /// ``right''  exponential slope 
      double m_lambdaR ; /// ``right''  exponential slope      
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  RaisingCosine 
     *  "Raising cosine" distribution
     *  \f$ f(x,\mu,s) = \frac{1}{2s}   \left( 1   +\cos \pi y \right)  \f$, 
     *  where \f$  y  \equiv = \frac{x-\mu}{s}\f$ 
     *  @see https://en.wikipedia.org/wiki/Raised_cosine_distribution
     */
    class RaisingCosine 
    {
    public:
      // ======================================================================
      /** constructor with all arguments 
       *  @param mu  the mean/mode/median of the distribution
       *  @param s   the width-parameteter
       */
      RaisingCosine ( const double mu = 0 , 
                      const double s  = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate raising cosine distribution
      double pdf (  const double x ) const ;
      /// evaluate raising cosine distribution
      double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // main getters 
      // ======================================================================
      double mu    () const { return  m_mu ; }
      double s     () const { return  m_s  ; }
      // ======================================================================
    public: // derived getters 
      // ======================================================================
      double location () const { return mu () ; }
      double peak     () const { return mu () ; }
      double scale    () const { return s  () ; }
      // ======================================================================
    public: //  setters 
      // ======================================================================
      bool setMu       ( const double value ) ;
      bool setS        ( const double value ) ;
      bool setScale    ( const double value ) { return setS  ( value ) ; }
      bool setLocation ( const double value ) { return setMu ( value ) ; }      
      bool setMean     ( const double value ) { return setMu ( value ) ; }      
      // ======================================================================
    public: // derived getters & stats 
      // ======================================================================
      /// mean 
      double mean     () const {  return  mu() ; }
      /// mode 
      double mode     () const {  return  mu() ; }
      ///  median
      double median   () const {  return  mu() ; }
      ///  variance  
      double variance () const ;
      ///  rms   
      double rms      () const ;
      /// skewness
      double skewness () const { return 0 ; }
      /// kurtosis
      double kurtosis () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get CDF 
      double cdf      ( const double x ) const ;
      /// evaluate the integral
      double integral ( const double low  , 
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// mean/mode/median
      double m_mu ; // mean/mode/median
      /// width-parameter   
      double m_s  ; // width-parameter   
      // ======================================================================      
    } ;
    // ========================================================================    
    /** @class QGaussian  
     *  q-Gaussian distribution:
     *  \f$ f(x) = \frac{ \sqrt{\beta}}{C_q} e_q (-\beta (x-\mu)^2)  \f$, 
     *  where  \f$ e_q (x) = \left( 1 + (1-q)x\right)^{\frac{1}{1-q}}\f$ 
     *  @see https://en.wikipedia.org/wiki/Q-Gaussian_distribution
     *  If is equal to 
     *  - scaled version of Student' t-distribution for 1<q<3
     *  - Gaussian distribution for q = 1 
     *  - has finite  support for q<1 
     *  Here we use \f$ \beta = \frac{1}{2\sigma^2}\f$
     */
    class QGaussian  
    {
    public:
      // ======================================================================
      /** constructor from all arguments            
       *  @param mean  the mean/mode/location of the peak 
       *  @param q     q-value   (q<3, for q>3 q = 6-q)
       *  @param scale 
       */
      QGaussian ( const double mean  = 0 ,   // mean/mode/location 
                  const double q     = 1 ,   //  q-parameter 
                  const double scale = 1 ) ; // scale/sigma
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate  pdf  for q-Gaussian distribution
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  for q-Gaussian distribution
      double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public : // primary getters 
      // ======================================================================
      double mean   () const { return  m_mean  ; }
      double q      () const { return  m_q     ; }
      double scale  () const { return  m_scale ; }
      // ======================================================================
    public : // derived getters 
      // ======================================================================
      double peak     () const { return mean  () ; }
      double mu       () const { return mean  () ; }
      double mode     () const { return mean  () ; }
      double median   () const { return mean  () ; }
      double location () const { return mean  () ; }
      double sigma    () const { return scale () ; }
      /// get the original beta 
      double beta     () const { return 0.5 / ( m_scale * m_scale ) ; }
      // ======================================================================
    public : // primay getters 
      // ======================================================================
      // set mean 
      bool setMean  ( const double value ) ;
      // set q :   if q>3, q=6-q
      bool setQ     ( const double value ) ;
      // set scale
      bool setScale ( const double value ) ;
      // ======================================================================
    public : // derived setters 
      // ======================================================================
      bool setMu       ( const double value ) { return setMean  ( value )  ; }
      bool setLocation ( const double value ) { return setMean  ( value )  ; }
      bool setSigma    ( const double value ) { return setScale ( value )  ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral () const { return 1 ; }
      /// get the integral 
      double integral ( const double low  , 
                        const double high ) const ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// mean/mode/location 
      double m_mean  ; // mean/mode/location 
      /// q-value 
      double m_q     ; // q-value 
      /// scale/sigma 
      double m_scale ; // scale/sigma
      // ======================================================================
    private:
      // ======================================================================
      /// get C_q constant 
      double m_cq ; // get C_q constant 
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;     
    // ========================================================================
    /** @class Hyperbolic 
     *  Hyperbolic disribtion
     *  @see  https://en.wikipedia.org/wiki/Hyperbolic_distribution
     *  @see  Barndorff-Nielsen, Ole, 
     *    "Exponentially decreasing distributions for the logarithm of particle size". 
     *    Proceedings of the Royal Society of London. Series A, 
     *    Mathematical and Physical Sciences. 
     *    The Royal Society. 353 (1674): 401–409
     *     doi:10.1098/rspa.1977.0041. JSTOR 79167.
     * 
     *  \f[  f(x;\mu, \beta, \delta, \gamma) = 
     *  \frac{\gamma}{2\alpha\delta K_1(\delta \gamma)}
     *  \mathrm{e}^{ -\sqrt{ \alpha^2\delta^2 + \alpha^2 (x-\mu)^2 } + \beta ( x - \mu)}
     *  \f]
     *  where 
     *  - \f$ \alpha^2 = \beta^2\f + \gamma^2$
     *  - \f$ K_1\f$ is a modified Bessel function of the second kind 
     *  
     * In the code we adopt parameterisation in terms of
     *  - location parameter      \f$\mu\f$
     *  - parameter               \f$\sigma \gt  0 \f$, related to the width;
     *  - dimensionless parameter \f$\kappa\f$,         related to the asymmetry;
     *  - dimensionless parameter \f$\zeta   \ge 0 \f$, related to the shape 
     *
     * The parameters are defined as:
     * \f[\begin{array}{lcl}
     *     \sigma^2 & \equiv & \gamma^{-2} \zeta \frac{K_2(\zeta)}{\zetaK_1(zeta)} \\
     *     \kappa   & \equiv & \frac{\beta}{\sigma} \\
     *     \zeta    & \equiv & \delta \gamma \end{array} \f]
     * - For \f$ \beta=0 (\kappa=0)\f$,  \f$\sigma^2\f$ is a variance of the distribution.
     * - Large values of \f$\zeta\f$ distribtionhas small kurtosis 
     * - For small \f$ \zeta \f$ distribution shows kurtosis of 3 
     *
     * The inverse transformation is:
     * \f[ \begin{array}{lcl}
     *     \beta    & = & \frac{\kappa}{\sigma}            \\
     *     \delta   & = & \frac{\zeta}{\gamma}             \\
     *     \gamma   & = & \frac{\sqrt{A^*(\zeta)}}{\sigma} \\
     *     \alpha   & = & \sqrt { \beta^2 + \gamma^2} \end{array} \f]
     * 
     * where \f$ A^{*}(\zeta) = \frac{\zeta K^*_2(\zeta)}{K^*_1(zeta)} \f$. 
     * It is largely inspired by NIM A764 (2014) 150, arXiv:1312.5000, 
     * but has much better properties when \f$ \zeta \rigtarrow 0 \f$ 
     * @see D. Martinez Santos and F. Dupertuis,
     *         "Mass distributions marginalized over per-event errors",
     *          NIM A764 (2014) 150, arXiv:1312.5000
     *          DOI: 10.1016/j.nima.2014.06.081",
     *
     *  The final form of the distribution is 
     *  \f[  f(x;\mu,\sigma,\zeta,\kappa) = 
     *     \frac{ A^{*2}(\zeta) } { 2 \sigma \sqrt{\kappa^2+A^{*2}(\zeta)} \zeta K^*_1(\zeta) } 
     *     \mathrm{e}^{\zeta - \sqrt{ (\kappa^2+A^{*2}(\zeta)) \left(
     *     \frac{\zeta^2}{A^{*2}(\zeta)} +\left( \frac{x-\mu}{\sigma}\right)^2  \right) } } 
     *  \f]
     *  where \f$ K^*_n(x)\f$ is a scaled modified Bessel functon to the second kind 
     *   \f$ K^*_n(x) = \mathrm{e}^{x}K_1(x) \f$ 
     *
     *  In all expressions \f$ \left| \sigma \right|\f$ and 
     *  \f$ \left| \zeta \right|\f$ are used instead of \f$\sigma\f$ and \f$\zeta\f$ 
     * 
     *  - For the finite values of \f$ \zeta\f$ it has 
     *    Gaussian-like core and exponential tails 
     * 
     *  Useful subclasses: 
     *  - \f$ \zeta\rigtharrow+\infty, \kappa=0\f$    : Gaussian function  
     *  - \f$ \zeta\rigtharrow+\infty, \kappa\neq0\f$ : shifted Gaussian function  
     *  - \f$ \zeta\rigtharrow+0 , \kappa=0\f$        : symmetric Laplace distribution
     *  - \f$ \zeta\rigtharrow+0 , \kappa\neq0\f$     : asymmetric Laplace distribution  
     */
    class Hyperbolic 
    {
    public :
      // =======================================================================
      /** constructor from mu, sigma, zeta and kappa 
       *  @param mu    related to location 
       *  @param beta  related to asymmetry
       *  @param sigma related to width 
       *  @param zeta  related to what ?
       */
      Hyperbolic ( const double mu     = 0 ,   // related to location 
                   const double sigma  = 1 ,   // related to withs  
                   const double zeta   = 1 ,   // related to shape 
                   const double kappa  = 0 ) ; // related to asymmetry 
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate  pdf  for Hyperbolic distribution
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  for Hyperbolic distribution
      double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// location parameter 
      double mu       () const { return m_mu    ; } // location parameter 
      /// sigma-parameter 
      double sigma    () const { return m_sigma ; } // sigma-parameter 
      /// squared sigma parameter 
      double sigma2   () const { return m_sigma * m_sigma  ; } 
      /// zeta-parameter 
      double zeta     () const { return m_zeta  ; } // zeta-parameter 
      /// squared zeta parameters
      double zeta2    () const { return m_zeta  * m_zeta   ; }
      /// asymmetry parameter kappa 
      double kappa    () const { return m_kappa ; } // asymmetry parameter kappa 
      /// squared asymmetry parameter
      double kappa2   () const { return m_kappa * m_kappa ; }
      // ======================================================================
    public : // original parameters 
      // ======================================================================
      /// alpha parameter 
      double alpha    () const { return std::hypot ( beta () , gamma () ) ; }
      /// squared alpha 
      double alpha2   () const { return beta2 () + gamma2 ()       ; }      
      /// beta parameter 
      double beta     () const { return m_kappa / m_sigma ; } // beta-parameter 
      /// squared beta parameter 
      double beta2    () const { return std:: pow ( beta  () , 2 ) ; }
      /// gamma parameter 
      double gamma    () const { return m_AL    / m_sigma          ; }
      /// squared gamma parameter  
      double gamma2   () const { return std::pow  ( gamma () , 2 ) ; }
      /// delta parameter  
      double delta    () const { return m_zeta * m_sigma / m_AL    ; }
      /// squared delta parameter  
      double delta2   () const { return std::pow  ( delta () , 2 ) ; }
      //  =====================================================================
    public : // features 
      // ======================================================================
      double location    () const { return mu ()   ; }
      /// get the actual mode of the distribution
      double mode        () const ; 
      /// get mean value 
      double mean        () const ;
      /// get variance 
      double variance    () const ;
      /// get dispersion 
      double dispersion  () const { return variance () ; }
      /// get RMS 
      double rms         () const { return std::sqrt ( variance () ) ; }
      /// get RMS 
      double RMS         () const { return rms ()      ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu       ( const double value ) ;
      bool setLocation ( const double value ) { return setMu ( value ) ; }
      bool setSigma    ( const double value ) ;
      bool setZeta     ( const double value ) ;
      bool setKappa    ( const double value ) ;
      // ======================================================================
    public: // extended setter 
      // ======================================================================
      /** set "standard" parameters 
       *  @param mu     mu-parameter, location 
       *  @param beta   beta-parameter, asymmetry 
       *  @param gamma  gamma-parameter, shape 
       *  @param delta  delta-parameter, scale 
       *  \f$ \delta \ge 0, left| \beta \right|<  \alpha \f$ 
       */ 
      bool setStandard
      ( const double mu    ,
        const double beta  , 
        const double gamma , 
        const double delta ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral () const { return 1 ; }
      /// get the integral 
      double integral ( const double low  , 
                        const double high ) const ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ; // get the unique tag 
      // ======================================================================
    private :
      // ======================================================================
      /// location parameter 
      double m_mu    { 0 } ;  // location 
      /// scale/width parameter 
      double m_sigma { 1 } ;  // width parameter 
      /// shape parameter
      double m_zeta  { 1 } ;  // shape parameter
      /// asymmetry parameter 
      double m_kappa { 0 } ;  // asymmetry parameter
      // ======================================================================
    private:  // helper parameters 
      // ======================================================================
      /** "constant" that relates sigma and gamma 
       *  \f$ A \equiv \frac{\zeta K^*_2(\zeta)}{K^*_1(\zeta)} \f$, 
       *  where \f$ K^*_n(x) \f$ is scaled modified irregular Bessel function 
       *  \f$ K^*_1(x) = \mathrmf{e}^{x}K_1(x)\f$
       *  - it behaves nicely when \f$\zeta\rightarrow 0 \f$
       *  - it behaves nicely when \f$\zeta\rightarrow +\infty \f$ 
       */
      double m_AL { -1 } ; // constant that relates sigma and gamma for fixed zeta 
      // ======================================================================
      /** helper normalization constant 
       *  \f$ k_1 = \frac{1}{\zeta K^*_1 ( \zeta ) }\f$, 
       *  where \f$ K^*_1(x) \f$ is scaled modified irregular Bessel function 
       *  \f$ K^*_1(x) = \mathrmf{e}^{x}K_1(x)\f$
       *  - it behaves nicely when \f$\zeta\rightarrow 0 \f$
       *  - it behaves nicely when \f$\zeta\rightarrow +\infty \f$ 
       */
      double m_N { -1 } ; // helper normalization constant 
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ; // integration workspace
      // ======================================================================      
    } ;
    // ========================================================================    
    /** @class GenHyperbolic
     *  Generalized Hyperbolic distribution
     *  @see https://en.wikipedia.org/wiki/Generalised_hyperbolic_distribution 
     *  
     *  \f[ f(x;\lambda, \alpha,\beta,\gamma,\delta,\mu)= 
     *  \frac{  (\gamma/\delta)^{\lambda} }{ \sqrt{2\pi} K_{\lambda}(\delta\gamma) }
     *   \mathrm{e}^{\beta (x-\mu)}
     *   \frac{  K_{\lambda -1/2} ( \alpha \sqrt{ \delta^2 + (x-\mu)^2} ) }
     *   { (  \sqrt{ \delta^2 + (x-\mu)^{2} }  /\alpha)^{1/2-\lambda} } \f]
     * where 
     *  - $\alpha=\sqrt{\beta^2+\gamma^2}$
     *
     * In the code we adopt parameterisation in terms of
     *  - location parameter      \f$\mu\f$
     *  - shape parameter         \f$\lambda\f$
     *  - parameter               \f$\sigma \gt  0 \f$, related to the width;
     *  - dimensionless parameter \f$\kappa\f$,         related to the asymmetry;
     *  - dimensionless parameter \f$\zeta   \ge 0 \f$, related to the shape 
     *  
     * The parameters are defined as:
     * \f[\begin{array}{lcl}
     *     \sigma^2 & \equiv & \gamma^{-2} \zeta \frac{K_{\lambda+1}(\zeta)}{\zetaK_{\lambda}(zeta)} \\
     *     \kappa   & \equiv & \frac{\beta}{\sigma} \\
     *     \zeta    & \equiv & \delta \gamma \end{array} \f]
     * - For \f$ \beta=0 (\kappa=0)\f$,  \f$\sigma^2\f$ is a variance of the distribution.
     * - Large values of \f$\zeta\f$ distribtionhas small kurtosis 
     * - For small \f$\zeta\f$ distribution shows kurtosis of 3 
     *
     * The inverse transformation is:
     * \f[ \begin{array}{lcl}
     *     \beta    & = & \frac{\kappa}{\sigma}                      \\
     *     \delta   & = & \frac{\zeta}{\gamma}                       \\
     *     \gamma   & = & \frac{\sqrt{A_{\lambda}^*(\zeta)}}{\sigma} \\
     *     \alpha   & = & \sqrt { \beta^2 + \gamma^2} \end{array} \f]
     *
     * In general it has exponential tails for \f$ \lambda >0 \f$ and Gaussian core.
     * For negative \f$ \lambda \f$ tails are more heavy..
     *
     * Special cases:
     *  - \f$ \lampba = 1 \f$    : hyperbolic  
     *  - \f$ \lampba = -1/2 \f$ : NIG   (normal inverse Gaussian)
     *  - \f$ \lambda = 0  \f$   : hyperbola 
     *  - \f$ \lambda = 1/2  \f$ : hyperboloid 
     * 
     * @see Ostap::Math::Hyperbolic
     * Useful subclasses 
     *  - \f$ \lambda=1\f$ : Hyperbolic distribution  
     *  - \f$ \lambda=-\frac{1}{2}\f$ : Normal Inverse Gaussian distributtion
     *  - \f$ \lambda=-\frac{n}{2}, \zeta\rightarrow+0\f$ : Student's t-distibution 
     *  - \f$ \lambda \rightarrow \pm\infty, \kappa=0\f$ : Gaussian distribution 
     *  - \f$ \zeta \rightarrow +\infty, \kappa=0\f$ : Gaussian distribution 
     *
     */
    class GenHyperbolic
    {
    public:
      // ======================================================================
      /** constructor from mu, sigma, zeta and kappa 
       *  @param mu     related to location 
       *  @param sigma  related to width 
       *  @param kappa  related to asymmetry
       *  @param zeta   related to kurtosis 
       *  @param lambda shape-related    
       */
      GenHyperbolic 
      ( const double mu     = 0 ,  // related to location 
        const double sigma  = 1 ,  // related to width 
        const double zeta   = 1 ,  // related to kurtosis
        const double kappa  = 0 ,  // related to asymmetry 
        const double lambda = 1 ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate  pdf  for Generalised Hyperbolic distribution
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  for Generalised Hyperbolic distribution
      inline double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // accessors
      // ======================================================================
      /// location parameter 
      inline double mu       () const { return m_mu     ; } // location parameter
      /// location parameter 
      inline double location () const { return m_mu     ; } // location parameter 
      /// sigma parameter 
      inline double sigma    () const { return m_sigma  ; } // sigma parameter 
      /// squared sigma parameter 
      inline double sigma2   () const { return m_sigma * m_sigma  ; } 
      /// asymmetry parameter 
      inline double kappa    () const { return m_kappa  ; } // kappa-parameter 
      /// squared asymmetry parameter
      inline double kappa2   () const { return m_kappa * m_kappa ; }
      /// zeta parameters
      inline double zeta     () const { return m_zeta   ; } // zeta-parameter 
      /// squared zeta parameters
      inline double zeta2    () const { return m_zeta  * m_zeta   ; }
      /// shape parameter 
      inline double lambda   () const { return m_lambda ; } // lambda-parameter 
      /// shape parameter 
      inline double lambd    () const { return m_lambda ; } // lambda-parameter  
      // ======================================================================
    public : // original parameters 
      // ======================================================================
      /// alpha parameter 
      inline double alpha    () const { return std::hypot ( beta () , gamma () ) ; }
      /// squared alpha 
      inline  double alpha2   () const { return beta2 () + gamma2 ()       ; }      
      /// beta parameter 
      inline double beta     () const { return m_kappa / m_sigma ; } // beta-parameter 
      /// squared beta parameter 
      inline double beta2    () const { return std:: pow ( beta  () , 2 ) ; }
      /// gamma parameter 
      inline double gamma    () const { return m_AL    / m_sigma          ; }
      /// squared gamma parameter  
      inline double gamma2   () const { return std::pow  ( gamma () , 2 ) ; }
      /// delta parameter  
      inline double delta    () const { return m_zeta * m_sigma / m_AL    ; }
      /// squared delta parameter  
      inline double delta2   () const { return std::pow  ( delta () , 2 ) ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool        setMu       ( const double value ) ;
      bool        setSigma    ( const double value ) ;
      bool        setKappa    ( const double value ) ;
      bool        setZeta     ( const double value ) ;
      bool        setLambda   ( const double value ) ;
      inline bool setLocation ( const double value ) { return setMu     ( value ) ; }
      inline bool setLambd    ( const double value ) { return setLambda ( value ) ; }
      // ======================================================================
    public: // extended setter 
      // ======================================================================
      /** set "standard" parameters 
       *  @param mu     mu-parameter, location 
       *  @param beta   beta-parameter, asymmetry 
       *  @param gamma  gamma-parameter, shape 
       *  @param delta  delta-parameter, scale 
       *  @param lambda lambda-parameter, shape 
       *  
       *  \f$ \alpha = \sqrt{\beta^2 + \gamma^2} \f$ 
       *
       *  \f[ \begin{array}{ll} 
       *    \delta \ge 0, left| \beta \right|<  \alpha &  if~\lambda > 0 \\ 
       *    \delta >   0, left| \beta \right|<  \alpha &  if~\lambda = 0 \\ 
       *    \delta >   0, left| \beta \right|\le\alpha &  if~\lambda < 0 
       *  \end{array}\f] 
       */ 
      bool setStandard
      ( const double mu     , 
        const double beta   , 
        const double gamma  , 
        const double delta  ,
        const double lambda ) ;
      // ======================================================================
    public : // features 
      // ======================================================================
      /// get mean value 
      double mean        () const ;
      /// get variance 
      double variance    () const ;
      /// get dispersion 
      inline double dispersion  () const { return variance () ; }
      /// get RMS 
      inline double rms         () const { return std::sqrt ( variance () ) ; }
      /// get RMS 
      inline double  RMS        () const { return rms ()      ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      inline double integral () const { return 1 ; }
      /// get the integral 
      double        integral
      ( const double low  , 
        const double high ) const ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ; // get the unique tag 
      // ======================================================================
    private: // main parametters 
      // ======================================================================
      /// mu - location parameter 
      double  m_mu    { 0 } ;    // mu : location parameter 
      /// sigma width parameter 
      double  m_sigma { 1 } ;    // sigma - width parameter 
      /// zeta  - related to kurtosis parameter 
      double  m_zeta  { 1 } ;    // zeta  - related to kurtosis parameter 
      /// kappa - asymmetry parameter 
      double  m_kappa { 0 } ;    // kappa - asymmetry parameter 
      /// lambda - shape parameter 
      double m_lambda { 1 } ;    // lambda - shape parameter 
      // ======================================================================
    private : // helper "constants"
      // ======================================================================
      /** "constant" that relates sigma and gamma 
       *  \f$ A \equiv \frac{\zeta K^*_{\lambda+1}(\zeta)}{K^*_{\lambda}(\zeta)} \f$, 
       *  where \f$ K^*_{\nu}(x) \f$ is scaled modified irregular Bessel function 
       *  \f$ K^*_{\nu}(x) = \mathrmf{e}^{x}K_{\nu}(x)\f$
       *  - it behaves nicely when \f$\zeta\rightarrow 0 \f$
       *  - it behaves nicely when \f$\zeta\rightarrow +\infty \f$ 
       */
      double  m_AL { -1 } ; // helper constant
      // = ====================================================================
      /** normalization constant 
       *  \f$ \frac{1}{ \zeta^{\lambda} K_{\lambda}(zeta)}\f$,
       *  where  \f$ K^*_{\nu}(z) \f$ is modified Bessel function 
       */
      double  m_N  { -1 } ; // normalization      
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ; // integration workspace
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Das
     *  Simple gaussian function with exponential tails.
     *  It corresponds to <code>ExpGaussExp</code> function, 
     *  \f[ 
     *    f (x ; \mu, \sigma, k_L, k_R ) = \frac{1}{\sqrt{2\pi}\sigma}
     *   \left\{ \begin{array}[lcl}
     *  \mathrm{e}^{  \frac{k_L^2}{2} + k_L\left(\frac{x-mu}{\sigma}\right) }
     *   & \mathrm{for}  &  \left(\frac{x-\mu}{\sigma}\right) < -k_L \\   
     *  \mathrm{e}^{ \frac{1}{s} \left( \frac{x-\mu}{\sigma}\right)^2}
     *   & \mathrm{for}  &  -k_L < \left(\frac{x-\mu}{\sigma}\right) < k_R \\    
     *  \mathrm{e}^{  \frac{k_R^2}{2} - k_R\left(\frac{x-mu}{\sigma}\right) }
     *   & \mathrm{for}  &  \left(\frac{x-\mu}{\sigma}\right)> k_R   
     *  \end{array} \right. \f]
     *  - \f$ k_L \ge 0\f$
     *  - \f$ k_R \ge 0\f$
     *
     *  @see Souvik Das, "A simple alternative to Crystall Ball fnuction"
     *                   arXiv:1603.08591  [hep-ex]
     *  @see https://arxiv.org/abs/1603.08591
     *  @attention - the function is not normalized! 
     *  Function was used in 
     *  @see CMS collaboration, V.Khachatryan, 
     *       "Search for resonant pair production of Higgs bosons decaying 
     *        to two bottom quark\textendash{}antiquark pairs 
     *        in proton-proton collisions at 8 TeV}",
     *        Phys. Lett. B749 (2015) 560 
     * @see https://arxiv.org/abs/1503.04114 
     * @see https://doi.org/10.1016/j.physletb.2015.08.047 
     * - Gaussian function is restored when \f$k_L,k_R \rigtharrow +\infty\f$ 
     */
    class Das
    {
      // ======================================================================
    public:
      // ====================================================================== 
      /** constructor with full parameters 
       *  @param mu peak location 
       *  @param sigma sigma for Gaussian Core 
       *  @param kL    left tail parameter 
       *  @param kR    right tail parameter 
       */
      Das 
      ( const double mu     = 0 ,    // location parameter 
        const double sigma  = 1 ,    // width parameter 
        const double kL     = 2 ,    // left tails 
        const double kR     = 2 ) ;  // right tail 
      // ======================================================================      
    public :
      // ======================================================================
      /// evaluate  pdf  
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  
      double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get location parameter 
      inline double mu       () const { return m_mu    ; }
      /// get width parameter 
      inline double sigma    () const { return m_sigma ; }
      /// get left tail 
      inline double kL       () const { return m_kL    ; }
      /// get right tail 
      inline double kR       () const { return m_kR    ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set location parameter 
      bool        setMu       ( const double value ) ;
      /// set width parameter 
      bool        setSigma    ( const double value ) ;
      /// set left tail 
      bool        setKL       ( const double value ) ;
      /// set right tail 
      bool        setKR       ( const double value ) ;
      /// set location parameter 
      inline bool setLocation ( const double value ) { return setMu ( value ) ; }
      // ======================================================================
    public: // derived  
      // ======================================================================
      /// get location parameter 
      inline double location () const { return mu () ; }
      /// get mode  
      inline double mode     () const { return mu () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ; // get the unique tag 
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral () const ;
      /// get the integral 
      double integral
      ( const double low  , 
        const double high ) const ;      
      // ======================================================================
    private:
      // ======================================================================
      /// location parameter 
      double m_mu    { 0 } ; // location parameter 
      /// width  parameter 
      double m_sigma { 1 } ; // width parameter 
      /// left tail 
      double m_kL    { 2 } ; // left tail
      /// right tail 
      double m_kR    { 2 } ; // right tail      
      // ======================================================================
    };
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PEAKS_H
// ============================================================================
