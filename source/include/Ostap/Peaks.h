// ============================================================================
#ifndef OSTAP_PEAKS_H 
#define OSTAP_PEAKS_H 1
// ============================================================================
//  Include  files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
#include "Ostap/Tails.h"
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
       *  @param sigma   the peak width 
       */
      Gauss
      ( const double peak  = 0 ,
        const double sigma = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Gaussian
      double        evaluate   ( const double x ) const ;
      /// evaluate Gaussian
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Gaussian
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      double peak    () const { return m_peak    ; }
      double m0      () const { return m_peak    ; }
      double mu      () const { return m_peak    ; }
      double mode    () const { return m_peak    ; }
      double mass    () const { return m_peak    ; }
      /// sigma
      double sigma   () const { return m_sigma   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position 
      bool setPeak    ( const double value ) ;
      /// left sigma 
      bool setSigma   ( const double value ) ;
      // ======================================================================
      /// peak position 
      inline bool setM0   ( const double value ) { return setPeak ( value ) ; }
      inline bool setMu   ( const double value ) { return setPeak ( value ) ; }
      inline bool setMode ( const double value ) { return setPeak ( value ) ; }
      inline bool setMass ( const double value ) { return setPeak ( value ) ; }
      // ======================================================================
    public: // integrals & CDF
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral 
      ( const double low  ,
        const double high ) const ;
      ///  get CDF 
      double cdf  ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the logarithmic derivative
       * \f$ \frac{ f6\prime}{f}  \f$
       */  
      double  dFoF ( const  double x ) const ; 
      // ======================================================================
      // get normalized version of the variable 
      inline double t ( const double x ) const 
      { return ( x - m_peak ) / m_sigma ;}
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
    /** @class BifurcatedGauss
     *  Aka a split normal distribution 
     *  @see https://en.wikipedia.org/wiki/Split_normal_distribution
     *  simple representation of bifurcated gaussian function/
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
      ( const double peak   ,
        const double sigmaL ,
        const double sigmaR ) ;
      // =====================================================================
      /** constructor from all parameters
       *  @param peak    the peak posiion
       *  @param siggma  average sigma 
       */
      BifurcatedGauss
      ( const double peak   =  0 ,
        const double sigma  =  1 ) ;
      // =====================================================================
      /** constructor from all parameters
       *  @param peak    the peak posiion
       *  @param gauss   Gaussian peak 
       */
      BifurcatedGauss
      ( const Gauss& gauss ) ; 
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Bifurcated Gaussian
      double evaluate   ( const double x ) const ;
      /// evaluate Bifurcated Gaussian
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate Bifurcated Gaussian
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      inline double peak      () const { return m_peak   ; }
      inline double m0        () const { return m_peak   ; }
      inline double mu        () const { return m_peak   ; }
      inline double mode      () const { return m_peak   ; }
      inline double mass      () const { return m_peak   ; }
      /// left sigma
      inline double sigmaL    () const { return m_sigmaL ; }
      /// right sigma
      inline double sigmaR    () const { return m_sigmaR ; }
      /// left sigma squared 
      inline double sigmaL2   () const { return m_sigmaL * m_sigmaL ; }
      /// right sigma squared 
      inline double sigmaR2   () const { return m_sigmaR * m_sigmaR ; }
      // ======================================================================
    public:
      // ======================================================================      
      /// average sigma 
      inline double sigma     () const { return 0.5  * ( m_sigmaL + m_sigmaR ) ; }
      /// sigma-asymmetry 
      inline double asymmetry () const { return m_kappa ; }
      /// sigma-asymmetry       
      inline double asym      () const { return m_kappa ; }
      /// sigma-asymmetry 
      inline double kappa     () const { return m_kappa ; }
      /** sigma-asymmetry:
       *  \f$ \kappa  \equiv \tanh \psi \f$ 
       */ 
      inline double psi       () const { return m_psi   ; } 
      // ======================================================================
      /** set asymmetry keeping average sigma untouched
       *  \f$ \left| \kappa \right| < 1 \f$ 
       */
      bool setKappa ( const double value ) ;
      // ======================================================================
      /** set asymmetry keeping average sigma untouched
       */
      bool setPsi   ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// peak position 
      bool setPeak    ( const double value ) ;
      /// set left  sigma (& keep right sigma)
      bool setSigmaL  ( const double value ) ;
      /// set right sigma (&keep left   sigma) 
      bool setSigmaR  ( const double value ) ;
      /// set both sigmas simultaneously (main method)
      bool setSigma   
      ( const double valueL ,
        const double valueR ) ;
      /// set both sigmas simultaneously
      inline bool setSigma   
      ( const double value ) { return setSigma ( value , value ) ; }
      // ======================================================================
      inline bool setM0   ( const double value ) { return setPeak ( value ) ; } 
      inline bool setMu   ( const double value ) { return setPeak ( value ) ; } 
      inline bool setMode ( const double value ) { return setPeak ( value ) ; } 
      inline bool setMass ( const double value ) { return setPeak ( value ) ; } 
      // ======================================================================
    public: // integrals & CDF
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral
      ( const double low  ,
        const double high ) const ;
      ///  get CDF 
      double cdf      ( const double x    ) const ;
      // ======================================================================
    public: // logarithmic derivative 
      // ======================================================================
      /** log-derivative \f$ \frac{ f^\prime}{f}  \f$
       *  - Useful to attach the tail to ensure 
       *  the continuity of the function and the 1st derivatibve 
       */
      double dFoF ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ;
      // ======================================================================
    private: // parameters
      // ======================================================================
      /// the peak position
      double m_peak   { 0 } ;       //                              the peak position
      /// sigma left
      double m_sigmaL { 1 } ;       // sigma-left
      /// sigma right
      double m_sigmaR { 1 } ;       // sigma-right
      /// asymmetry 
      double m_kappa  { 0 } ; 
      /// psi 
      double m_psi    { 0 } ; 
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
      double        pdf        ( const double x ) const ;
      inline double operator() ( const double x ) const { return pdf ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      inline double peak     () const { return m_peak    ; }
      inline double mean     () const { return m_peak    ; }
      inline double m0       () const { return m_peak    ; }
      inline double mu       () const { return m_peak    ; }
      inline double mode     () const { return m_peak    ; }      
      inline double mass     () const { return m_peak    ; }      
      /// sigma 
      inline double sigma    () const { return m_sigma   ; }
      /// sigma-1
      inline double sigma1   () const { return m_sigma   ; }
      /// sigma-2
      inline double sigma2   () const { return m_sigma * m_scale ; }
      /// scale 
      inline double scale    () const { return m_scale    ; }
      /// fraction 
      inline double fraction () const { return m_fraction ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak positon
      bool setPeak     ( const double value ) ;
      /// sigma 
      bool setSigma    ( const double value ) ;
      /// scale 
      bool setScale    ( const double value ) ;
      /// fraction 
      bool setFraction ( const double value ) ;
      // ======================================================================
      inline bool setM0   ( const double value ) { return setPeak ( value ) ; }
      inline bool setMu   ( const double value ) { return setPeak ( value ) ; }
      inline bool setMode ( const double value ) { return setPeak ( value ) ; }
      inline bool setMass ( const double value ) { return setPeak ( value ) ; }
      // ======================================================================      
    public:
      // ======================================================================
      /// get the cdf 
      double cdf      ( const double x ) const ;
      /// get the integral
      double integral () const { return 1  ; }  ;
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
      // ======================================================================
    public :
      // ======================================================================
      /// get pdf
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      /// get pdf
      double        pdf        ( const double x ) const ;
      // ======================================================================
    public: // primary getters
      // ======================================================================
      /// peak position 
      inline double mu          () const { return m_mu       ; }
      /// peak position 
      inline double peak        () const { return m_mu       ; }
      /// peak position 
      inline double location    () const { return m_mu       ; }
      ///  alpha 
      inline double alpha       () const { return m_alpha    ; }
      ///  scale 
      inline double scale       () const { return   alpha () ; }
      ///   beta 
      inline double beta        () const { return m_beta     ; }
      ///  shape 
      inline double shape       () const { return   beta  () ; }
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
      inline bool  setPeak      ( const double value ) { return setMu    ( value ) ; }
      ///  set position 
      inline bool  setLocation  ( const double value ) { return setMu    ( value ) ; }
      ///  scale 
      inline bool  setScale     ( const double value ) { return setAlpha ( value ) ; }
      ///  shape 
      inline bool  setShape     ( const double value ) { return setBeta  ( value ) ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      /// mean
      inline double mean   () const { return m_mu ; }
      /// median
      inline double median () const { return m_mu ; }
      /// mode 
      inline double mode   () const { return m_mu ; }
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
      // ======================================================================
    public :
      // ======================================================================
      /// get PDF
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      /// get PDF
      double        pdf        ( const double x ) const ;
      // ======================================================================
    public: // primary getters
      // ======================================================================
      inline double xi          () const { return m_xi       ; }
      inline double peak        () const { return   xi    () ; }
      inline double location    () const { return   xi    () ; }
      inline double alpha       () const { return m_alpha    ; }
      inline double scale       () const { return   alpha () ; }
      inline double kappa       () const { return m_kappa    ; }
      inline double shape       () const { return   kappa () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool  setXi        ( const double value ) ;
      bool  setAlpha     ( const double value ) ;
      bool  setKappa     ( const double value ) ;
      //
      inline bool  setPeak      ( const double value ) { return setXi    ( value ) ; }
      inline bool  setLocation  ( const double value ) { return setXi    ( value ) ; }
      inline bool  setScale     ( const double value ) { return setAlpha ( value ) ; }
      inline bool  setShape     ( const double value ) { return setKappa ( value ) ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      /// get the mean   of distribution s
      double mean        () const ;
      /// get the median of distribution s
      double median      () const { return   xi    () ; }
      /// get the mode   of distribution s
      double mode        () const ;
      //
      double variance    () const ;
      double dispersion  () const { return variance () ; }
      double sigma2      () const { return variance () ; }
      double sigma       () const ;
      //
      double skewness    () const ;
      double kurtosis    () const ;
      // ======================================================================
      // ======================================================================
    public:  // integrals
      // ======================================================================
      /// get the CDF 
      double cdf        ( const double x ) const ;
      /// get the integral
      double integral   () const { return 1 ; }
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
      double  y ( const double x ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_xi      ;  // location
      double m_alpha   ;  // scale
      double m_kappa   ;  // shape
      // ======================================================================
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
      // ======================================================================
    public :
      // ======================================================================
      /// get pdf
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      inline double evaluate   ( const double x ) const { return pdf ( x ) ; }
      double        pdf        ( const double x ) const ;
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
      double mode        () const ;
      double variance    () const ;
      double dispersion  () const { return variance () ; }
      double sigma2      () const { return variance () ; }
      double sigma       () const ;
      double skewness    () const ;
      double kurtosis    () const ;
      // ======================================================================
      /// approximate mode 
      double approximate_mode () const ;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      double cdf        ( const double x ) const ;
      /// get the integral
      double integral   () const { return 1 ; }
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
     *  It is a distibution for the variable that is a 
     *  sum (or difference for negative \f$ k\f$) 
     *  of a Gaussian and exponential variables: \f$ X \sim Y + sign(k) Z \f$,  
     *  where 
     *  - \f$ Y \sim N(\mu,\sigma) \f$
     *  - \f$ Z \sim  \frac{1}{k\sigma}\mathrm{e}^{-\frac{x}{k\sigma}} \f$ 
     *  
     *  For \f$ k=0\f$ one gets a Gaussian distribution
     *  - \f$ k>0\f$ corresponds to the rigth tail  
     *  - \f$ k<0\f$ corresponds to the left tail  
     *
     *  It can be considered as "single-tail" version of the Normal Laplace distribution:
     *  - \f$ k = 0 \f$ corresponds to Gaussian distribution
     *  - \f$ k > 0 \f$ corresponds to Normal Laplace \f$ NL(\mu,\sigma,0,k)\f$ 
     *  - \f$ k < 0 \f$ corresponds to Normal Laplace \f$ NL(\mu,\sigma,\left|k\right|,0)\f$ 
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
      inline  double mu       () const { return m_mu        ; }
      /// parameter varsigma
      inline double varsigma () const { return m_varsigma  ; }
      /// parameter k 
      inline double k        () const { return m_k         ; }
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
      double cdf 
      ( const double x ) const ;
      /// get the integral
      double integral   () const ;
      /// get the integral between low and high limits
      double integral  
      ( const double low  ,
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
      /// mode 
      double mode        () const ;
      // ======================================================================
    public:
      // ======================================================================
      // difference between mode and mu 
      double delta       () const { return m_varsigma * m_mk ; }
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
      /// mode-related (cache) parameter 
      double m_mk       { 0 } ;
      // ======================================================================
    } ;
    // ========================================================================
    /* @class ExGauss2 
     * Reparameterization of Exponentiallymodified Gaussian distribution
     * using the mode as parameters
     * @see Ostap::Math::ExGauss 
     */
    class ExGauss2 
    {
    public:
      // ======================================================================
      /// constructor from all parameters 
      ExGauss2
      ( const double mode     = 0 , // the mode 
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
      double mu       () const { return m_emg.mode     () ; }
      /// parameter varsigma
      double varsigma () const { return m_emg.varsigma () ; }
      /// parameter k 
      double k        () const { return m_emg.k        () ; }
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
      double cdf ( const double x ) const ; 
      /// get the integral
      double integral   ()          const ; 
      /// get the integral between low and high limits
      double integral  
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean value 
      double mean        () const { return m_emg.mean     () ; }
      /// variance 
      double variance    () const { return m_emg.variance () ; }
      /// RMS 
      double rms         () const { return m_emg.rms      () ; }
      /// dispersion 
      double dispersion  () const { return m_emg.variance () ; }
      /// skewness 
      double skewness    () const { return m_emg.skewness () ; }
      /// kurtosis 
      double kurtosis    () const { return m_emg.kurtosis () ; }      
      /// get cumulant 
      double cumulant    ( const unsigned short r ) const 
      { return m_emg.cumulant ( r ) ; }
      /// mode 
      double mode        () const { return m_emg.mode     () ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the function 
      Ostap::Math::ExGauss  m_emg {} ; // the function      
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bukin2
     *  Varinant of \f$ f_{2B1}\f$ function
     *  Essentially it is a sum of two ExGauss2 functionn with the same mode.
     *  It is more flexible that NormalLaplace (2 more parameters) 
     *  \f[ f(x; \mu , \sigma_A, \sigma_B , k_A , k_B , \phi ) = 
     *   \sin^2 ( \phi + \frac{\pi}{4} ) \times E ( \mu , \sigma_A , k_A  ) + 
     *   \cos^2 ( ]phi + \frac{\pi}{4} ) \times E ( \mu , \sigma_B , k_B ) \f]
     *  = where \f$ E \f$  stands for exponential modified Gaussian functions
     *  parameterized with the mode parameter..
     *  @see  A.Bukin, "Fitting function for asymmetric peaks}",
     *                 https://arxiv.org/abs/0711.4449
     *  @see https://arxiv.org/abs/0711.4449
     *  @see Ostap::Math::ExGauss2
     *  @see Ostap::Math::ExGauss
     *  @see Ostap::Math::NormalLaplace
     */
    class Bukin2
    {
    public:
      // ======================================================================
      /// constructor with all parameters
      Bukin2
      ( const double mu        =  0.0 ,   // the mode 
        const double varsigmaA =  1.0 ,
        const double varsigmaB =  1.0 ,
        const double kA        = -1.0 ,   // left taiil 
        const double kB        = +1.0 ,   // right tail
        const double phi       =  0.0 ) ; // equal fractionsn 
      // ======================================================================
    public:
      // ======================================================================
      double        evaluate   ( const double x  ) const ;
      inline double operator() ( const double x  ) const { return evaluate ( x ) ; } 
      inline double pdf        ( const double x  ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: // getters  
      // ======================================================================
      inline double mu        () const { return m_A.mu       () ; }
      inline double varsigmaA () const { return m_A.varsigma () ; }
      inline double varsigmaB () const { return m_B.varsigma () ; }
      inline double kA        () const { return m_A.k        () ; }
      inline double kB        () const { return m_B.k        () ; }
      inline double phi       () const { return m_phi           ; }
      // ======================================================================
    public: // derived parameters 
      // ======================================================================
      /// fraction of "A" componet 
      inline double fA        ()  const { return m_fA ; }
      /// fraction of "B" componet 
      inline double fB        ()  const { return m_fB ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool        setMu         ( const double value  ) ;
      bool        setPhi        ( const double value  ) ;
      inline bool setVarsigmaA  ( const double value  ) { return m_A.setVarsigma ( value ) ; }
      inline bool setVarsigmaB  ( const double value  ) { return m_B.setVarsigma ( value ) ; }
      inline bool setKA         ( const double value  ) { return m_A.setK        ( value ) ; }
      inline bool setKB         ( const double value  ) { return m_B.setK        ( value ) ; }
      // ======================================================================
    public:  // integrals
      // ======================================================================
      /// get CDF
      double cdf ( const double x ) const ; 
      /// get the integral
      double integral   ()          const ; 
      /// get the integral between low and high limits
      double integral  
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean value 
      double mean        () const ;
      /// mode 
      double mode        () const { return mu () ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// "A" component
      Ostap::Math::ExGauss2 m_A   {}      ; // "A" component      
      /// "B" component
      Ostap::Math::ExGauss2 m_B   {}      ; // "B" component
      /// angle phi
      double                m_phi { 0 }   ; // angle phi
      /// fraction of "A" component
      double                m_fA  { 0.5 } ; // fraction of "A" component
      /// fraction of "B" component
      double                m_fB  { 0.5 } ; // fraction of "B" component      
      // ======================================================================
    };
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
      inline double mu       () const { return m_mu        ; }
      /// parameter varsigma
      inline double varsigma () const { return m_varsigma  ; }
      /// left  exponential 
      inline double kL       () const { return m_kL        ; }
      /// right exponential 
      inline double kR       () const { return m_kR        ; }
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
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean value 
      double mean        () const ;
      /// mode value 
      double mode        () const ;
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
    public:
      // ======================================================================
      /// evaluate Bukin's function
      double pdf        ( const double x ) const ;
      /// evaluate Bukin's function
      double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      inline double peak  () const { return m_peak    ; }
      inline double m0    () const { return m_peak    ; }
      inline double mu    () const { return m_peak    ; }
      inline double mass  () const { return m_peak    ; }
      inline double mode  () const { return m_peak    ; }
      inline double sigma () const { return m_sigma   ; }
      inline double xi    () const { return m_xi      ; }
      inline double rho_L () const { return m_rho_L   ; }
      inline double rho_R () const { return m_rho_R   ; }
      // ======================================================================
      inline double x1    () const { return m_x1      ; }
      inline double x2    () const { return m_x2      ; }
      // ======================================================================
    public:
      // ======================================================================
      bool        setPeak  ( const double value ) ;
      bool        setSigma ( const double value ) ;
      bool        setXi    ( const double value ) ;
      bool        setRhoL  ( const double value ) ;
      bool        setRhoR  ( const double value ) ;
      inline bool setM0    ( const double value ) { return setPeak ( value ) ; }
      inline bool setMu    ( const double value ) { return setPeak ( value ) ; }
      inline bool setMode  ( const double value ) { return setPeak ( value ) ; }
      inline bool setMass  ( const double value ) { return setPeak ( value ) ; }
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
     *  Novosibirsk-function for description of gaussian with tails
     *  @see H.Ikeda et al., 'A detailed test of the CsI(Tl) calorimeter 
     *      for BELLE with photon beams of energy between 20MeV and 5.4 GeV',
     *       Nucl. Instrum. Meth. A441, (2000) 401.
     *  @see DOI: 10.1016/S0168-9002(99)00992-4
     *  @see https://inspirehep.net/literature/508223 
     *  @see https://doi.org/10.1016/S0168-9002(99)00992-4
     *
     *  \f$ f(x;\mu,\sigma,\tau) = \frac{1}{\sqrt{2\pi}\sigma}
     *  \mathrm{e}^{  -\frac{1}{2} \frac { \log^2 \left( 1 + \Lambda \tau \delta \right) }{\tau^2} 
     *                -\frac{\tau^2}{2} } \f$
     *  where 
     *  - \f$ \delta  = \frac{ x - \mu}{\sigma}\f$ 
     *  - \f$ \Lambda = \frac{  \sinh{ \tau \sqrt{\log 4}} }{\tau\sqrt{\log 4 }}\f$ 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class  Novosibirsk
    {
    public :
      // ======================================================================
      /** constructor from all parameters
       *  @param m0    the peak position
       *  @param sigma the effective sigma
       *  @param tau   the tail paramter
       */
      Novosibirsk
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double tau   = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Novosibirsk's function
      double pdf        ( const double x ) const ;
      /// evaluate Novosibirsk's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      inline double evaluate   ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public:
      // ====================================================================== 
      inline double m0    () const { return m_m0       ; } 
      inline double mu    () const { return m_m0       ; }
      inline double peak  () const { return m_m0       ; }
      inline double mass  () const { return m_m0       ; }
      inline double sigma () const { return m_sigma    ; }
      inline double tau   () const { return m_tau      ; }
      // ======================================================================
      double        mode  () const ;
      // ======================================================================
    public:
      // ======================================================================
      bool        setM0    ( const double value ) ;
      bool        setSigma ( const double value ) ;
      bool        setTau   ( const double value ) ;
      // ======================================================================
      inline bool setMu    ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass  ( const double value ) { return setM0 ( value ) ; }
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
    private: // parameters
      // ======================================================================
      /// the peak position
      double m_m0      { 0 } ; //                            the peak position
      /// the effective resolution
      double m_sigma   { 1 } ; //                     the effective resolution
      /// the tail parameter
      double m_tau     { 0 } ; //                           the tail parameter
      // ======================================================================
    private: // internals
      // ======================================================================
      /// lambda value
      double m_lambda  { 1 } ; // lambda value
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    // Crystal Ball & friends 
    // ========================================================================

    // ========================================================================
    /** @class CrystalBall
     *  ``Crystal Ball-function'' for description of gaussian with the tail
     *  @see http://en.wikipedia.org/wiki/Crystal_Ball_function
     *
     *  @see http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-r-255.pdf
     *  @see http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-r-236.pdf
     *  @see http://inspirehep.net/record/230779/files/f31-86-02.pdf
     *
     *  - J. E. Gaiser, Appendix-F Charmonium Spectroscopy from Radiative Decays of the J/Psi and Psi-Prime, 
     *    Ph.D. Thesis, SLAC-R-255 (1982)
     *  - M. J. Oreglia, A Study of the Reactions psi prime --> gamma gamma psi, 
     *    Ph.D. Thesis, SLAC-R-236 (1980), Appendix D.
     *  - T. Skwarnicki, A study of the radiative CASCADE transitions between the Upsilon-Prime and Upsilon resonances, 
     *    Ph.D Thesis, DESY F31-86-02(1986)
     *  
     *  However, here we adopt a bit different normalization:
     *
     *  \f[ f(x;\alpha,n,x_0,\sigma) \propto \frac{1}{ \sqrt{2\pi\sigma^2} } \left\{
     *  \begin{array}{ll}
     *  \mathrm{e}^{ -\frac{1}{2}\delta^2} 
     *  & \text{for}~\frac{x-x_0}{\sigma}\ge-\left|\alpha\right| \  \
     *  A \times \left( B - \delta )^{-N}}
     *  & \text{for}~\frac{x-x_0}{\sigma}\le-\left|\alpha\right| 
     *  \end{array}
     *  \right.\f]
     *  where :
     *  - \f$ \delta \equiv \frac{x-m_0}{\sigma}\f$ 
     *  - \f$ A \equiv \left( \frac{N}{\left| \alpha\right|}\right)^N \mathrm{e}^{-\frac{1}{2}\alpha^2}} \f$ 
     *  - \f$ B \equiv \frac{N}{\left|\alpha|} - \left|\alpha\right| \f$ 
     *  - \f$ N \equiv \sqrt { 1 + n^2} \f$ 
     *
     *  The function "almost normalized". It allowed to avoid the pathological 
     *  situations with \f$ \alpha \rigtharrow 0\f$ and \f$ M \le 1 f\$
     *
     *  @see https://en.wikipedia.org/wiki/Crystal_Ball_function
     *
     *  @attention: unlike the ff in Tomas'z theses, this function is NOT
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  CrystalBall
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0       parameter
       *  @param sigma  sigma    parameter
       *  @param alpha  alpha    parameter
       *  @param n      n        (external) parameter (not the same as *internal* N)
       */
      CrystalBall
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double alpha = 2 ,
        const double n     = 1 ) ;
      // ======================================================================
      /** constructor from gaussian and tail parameeter 
       */
      CrystalBall 
      ( const Ostap::Math::Gauss& core      , 
	const double              alpha = 2 , 
	const double              n     = 1 ) ;
      // ======================================================================
      /** constructor from gaussian and tail 
       *  @parameter core Gaussian function 
       *  @parameter tail tail     function 
       */
      CrystalBall 
      ( const Ostap::Math::Gauss& core ,
	const Ostap::Math::Tail&  tail ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double        pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0     () const { return m_core.m0     () ; }
      inline double mu     () const { return m_core.m0     () ; }
      inline double peak   () const { return m_core.m0     () ; }
      inline double sigma  () const { return m_core.sigma  () ; }
      inline double alpha  () const { return m_tail.alpha  () ; }
      inline double n      () const { return m_tail.n      () ; }      
      /// Internal N-parameter
      inline double N      () const { return m_tail.N      () ; } // internal N parameter
      /// squared alpha 
      inline double alpha2 () const { return m_tail.alpha2 () ; }
      // ======================================================================
    public:            
      // ======================================================================
      /// mode of the distribution 
      inline double mode () const { return m_core.mode () ; }
      /// the point where Gaussian meets power-law
      inline double xL   () const { return m0 () - alpha () * sigma()  ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline bool setM0    ( const double value ) { return m_core.setM0    ( value ) ; }
      inline bool setSigma ( const double value ) { return m_core.setSigma ( value ) ; }
      inline bool setN     ( const double value ) { return m_tail.setN     ( value ) ; } 
      inline bool setAlpha ( const double value ) { return m_tail.setAlpha ( value ) ; }
      // ======================================================================
      inline bool setMu    ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode  ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass  ( const double value ) { return setM0 ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral
      ( const double low ,
        const double high ) const ;
      // ======================================================================
      /** get the integral from the negative to positive infinity 
       *  @attention +infinity is returned for <code>n=0(N=1)</code>
       */
      double integral () const ;
      // ======================================================================
    public: // components 
      // ======================================================================
      /// get the Gaussian core 
      const Ostap::Math::Gauss&    core       () const { return m_core ; }
      /// get the Gaussian core 
      const Ostap::Math::Gauss&    gauss      () const { return m_core ; }
      /// get left tail
      const Ostap::Math::LeftTail& tail       () const { return m_tail ; }
      /// get left tail
      const Ostap::Math::LeftTail& tail_left  () const { return m_tail ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tail, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// core Gaussian 
      Ostap::Math::Gauss    m_core  {} ; // core gaussian
      /// (left) Tail
      Ostap::Math::LeftTail m_tail  {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Needham
     *  The special parametrization by Matthew NEEDHAM of
     *  `Crystal Ball-function' suitable for \f$J/\psi/\Upsilon\f$-peaks     
     *  - thanks to Matthew Needham
     *
     *  - alpha is parameterized as function of sigma 
     *  \f$ \alpha(\sigma) = \sqrt { \alpha_{min}^2 + ( c_0\frac{ (\sigma/c_1)^{c_2}}{ 1 + (\sigma/c_1)^{c_2} } )^2\f$ 
     *
     *  @attention For majority of physics cases <code>n</code> 
     *             can be fixed <code>n=0</code> (corresponds to <code>N=1</code>
     *
     *  @attention parameter \f$ c_1 \f$ is inverse with respect to the original 
     *             Matt's code
     *
     *  Reasonable values:
     *  - for \f$ c_0 \f$ :  \f$ 2.0 \le c_0 \le 3.0 \f$ 
     *  - for \f$ c_1 \f$ :  \f$ c_2 \approx O(\sigma) \f$ 
     *  - for \f$ c_2 \f$ :  \f$ c_2 \approx O(10) \f$   
     *  - for \f$ \alpha_{min} \f$ : \f$ \alpha_{min} \approx (0.01) \ll 1 \f$,
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
       *  @param c0     c0       parameter (should be between 1 and 4, a bit less than 3 ) 
       *  @param c1     c1       parameter (should be close to `sigma`)
       *  @param c2     c2       parameter (should be in excess of 2) 
       *  @param n      n        parameter
       *  @param amin   amin     parameter
       */
      Needham
      ( const double m0    = 3096.0  ,   // for J/psi
        const double sigma =   13.5  ,
        const double c0    =    2.5  ,   // (should be between 1 and 4, a bit less than 3 ) 
        const double c1    =   13.5  ,   // (should be sclose to `sigma`)
        const double c2    =   10    ,   // (should be in excess of 2) 
        const double n     =    0    ,   // note that it is different from internal N!
	const double amin  =    0.01 ) ; // cut-off parameter for alpha      
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Needham's function
      double        pdf        ( const double x ) const ;
      /// evaluate Needham's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0    () const { return m_cb.m0    () ; }
      inline double mu    () const { return m_cb.mu    () ; }
      inline double peak  () const { return m_cb.peak  () ; }
      inline double mode  () const { return m_cb.mode  () ; }
      inline double sigma () const { return m_cb.sigma () ; }
      inline double c0    () const { return m_c0          ; }
      inline double c1    () const { return m_c1          ; }
      inline double c2    () const { return m_c2          ; }
      inline double alpha () const { return m_cb.alpha () ; }
      inline double n     () const { return m_cb.n     () ; }
      inline double N     () const { return m_cb.N     () ; }      
      /// the point where Gaussian meets power-law
      inline double xL    () const { return m_cb.xL    () ; }
      // ======================================================================
    public: // show alpha as function of sigma 
      // ======================================================================
      /// alpha as function of sigma 
      double alpha ( const double sigma ) const ;
      /// minimal/cut-off value of alpha
      inline double amin      () const { return m_amin ; }
      /// minimal/cut-off value of alpha
      inline double alpha_min () const { return m_amin ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline bool setM0    ( const double value ) { return m_cb.setM0    ( value ) ; }
      inline bool setMu    ( const double value ) { return m_cb.setMu    ( value ) ; }
      inline bool setPeak  ( const double value ) { return m_cb.setPeak  ( value ) ; }
      inline bool setMode  ( const double value ) { return m_cb.setMode  ( value ) ; }
      inline bool setMass  ( const double value ) { return m_cb.setMass  ( value ) ; }
      inline bool setN     ( const double value ) { return m_cb.setN     ( value ) ; }
      // ====================================================================
    public:
      // ====================================================================
      // setting sigma causses some change in alpha ...
      bool setSigma ( const double value ) ;
      // =====================================================================
      /// set all three c-values together
      bool setC
      ( const double c0 ,
        const double c1 ,
        const double c2 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get integral between low and high
      inline double integral
      ( const double low ,
        const double high ) const
      { return m_cb.integral ( low , high ) ; }
      // ======================================================================
      /** get the integral from the the negative to positive infinity 
       *  @attention +infinity is returned for <code>n=0(N=1)</code>
       */      
      inline double integral () const { return m_cb.integral () ; }
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tail, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       */
      inline double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const
      { return m_cb.non_gaussian ( xlow , xhigh ) ; } 
      // ======================================================================
    public: // components 
      // ======================================================================
      /// get the Gaussian core 
      const Ostap::Math::Gauss&    core      () const { return m_cb.core () ; }
      /// get the Gaussian core 
      const Ostap::Math::Gauss&    gauss     () const { return m_cb.core () ; }
      /// get left tail
      const Ostap::Math::LeftTail& tail      () const { return m_cb.tail () ; }
      /// get left tail
      const Ostap::Math::LeftTail& tail_left () const { return m_cb.tail () ; }
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
      /// c0-parameter
      double m_c0                   ;  // c0_parameter
      /// c0-parameter
      double m_c1                   ;  // c1_parameter
      /// c0-parameter
      double m_c2                   ;  // c2_parameter
      /// alpha min parameter
      double m_amin          { 0.01 } ; // amin-parameter
      // ======================================================================
    } ; //                                The end of class Ostap::Math::Needham 
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
       *  @param n      n        (external) parameter (not the same as *internal* N)
       */
      CrystalBallRightSide 
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double alpha = 2 ,
        const double n     = 1 ) ;
      // ======================================================================
      /** constructor from gaussian and tail parameeter 
       */
      CrystalBallRightSide
      ( const Ostap::Math::Gauss& core      , 
	const double              alpha = 2 , 
	const double              n     = 1 ) ;
      // ======================================================================
      /** constructor from gaussian and tail 
       *  @parameter core Gaussian function 
       *  @parameter tail tail     function 
       */
      CrystalBallRightSide  
      ( const Ostap::Math::Gauss& core ,
	const Ostap::Math::Tail&  tail ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double        pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0     () const { return m_core.m0     () ; }
      inline double mu     () const { return m_core.m0     () ; }
      inline double peak   () const { return m_core.m0     () ; }
      inline double sigma  () const { return m_core.sigma  () ; }
      inline double alpha  () const { return m_tail.alpha  () ; }
      inline double n      () const { return m_tail.n      () ; }      
      /// Internal N-parameter
      inline double N      () const { return m_tail.N      () ; } // internal N parameter
      /// squared alpha 
      inline double alpha2 () const { return m_tail.alpha2 () ; }
      // ======================================================================
    public:            
      // ======================================================================
      /// mode of the distribution 
      inline double mode () const { return m_core.mode () ; }
      /// the point where Gaussian meets power-law
      inline double xR   () const { return m_core.m0 () + alpha () * sigma()  ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline bool setM0     ( const double value ) { return m_core.setM0    ( value ) ; }
      inline bool setSigma  ( const double value ) { return m_core.setSigma ( value ) ; }
      inline bool setN      ( const double value ) { return m_tail.setN     ( value ) ; } 
      inline bool setAlpha  ( const double value ) { return m_tail.setAlpha ( value ) ; } 
      // ======================================================================
      inline bool setMu    ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode  ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass  ( const double value ) { return setM0 ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral
      ( const double low ,
        const double high ) const ;
      // ======================================================================
      /** get the integral from the the negative to positive infinity 
       *  @attention +infinity is returned for <code>n=0(N=1)</code>
       */      
      double integral () const ;
      // ======================================================================
    public: // components 
      // ======================================================================
      /// get the Gaussian core 
      const Ostap::Math::Gauss&     core       () const { return m_core ; }
      /// get the Gaussian core 
      const Ostap::Math::Gauss&     gauss      () const { return m_core ; }
      /// get left tail
      const Ostap::Math::RightTail& tail       () const { return m_tail ; }
      /// get left tail
      const Ostap::Math::RightTail& tail_right () const { return m_tail ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tail, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// core Gaussian 
      Ostap::Math::Gauss     m_core  {} ; // core gaussian
      /// (right) Tail
      Ostap::Math::RightTail m_tail  {} ;
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
       *  @param m0     m0          parameter
       *  @param sigma  sigma       parameter
       *  @param alphaL alpha_L     parameter
       *  @param nL     n_L         parameter  (N-1 for "standard" definition)
       *  @param alphaR alpha_R parameter
       *  @param nR     n_R         parameter  (N-1 for "standard" definition)
       */
      CrystalBallDoubleSided
      ( const double m0      = 1 ,
        const double sigma   = 1 ,
        const double alphaL  = 2 ,
        const double nL      = 1 ,
        const double alphaR  = 2 ,
        const double nR      = 1 ) ;
      // ======================================================================
      /** constructor from all parameters
       *  @param core core Gaussian 
       *  @param alphaL  alpha_L     parameter
       *  @param nL      n_L         parameter  (N-1 for "standard" definition)
       *  @param alphaR  alpha_R parameter
       *  @param nR      n_R         parameter  (N-1 for "standard" definition)
       */
      CrystalBallDoubleSided
      ( const Ostap::Math::Gauss& core , 
        const double              alpha_L = 2 ,
        const double              n_L     = 1 ,
        const double              alpha_R = 2 ,
        const double              n_R     = 1 ) ;
      // ========================================================================
      /** constructor from all components  
       *  @param core core Gaussian 
       *  @param left  left  tail
       *  @param right right tail
       */
      CrystalBallDoubleSided
      ( const Ostap::Math::Gauss&     core  ,
        const Ostap::Math::LeftTail&  left  ,
	      const Ostap::Math::RightTail& right ) ;
      // ========================================================================
      /** constructor from all components  
       *  @param cb CrystalBall function (left tail) 
       *  @param right right tail
       */
      CrystalBallDoubleSided
      ( const Ostap::Math::CrystalBall& cb    , 
	      const Ostap::Math::RightTail  & right ) ;
      // =======================================================================
      /** constructor from all components  
       *  @param cb CrystalBallRightSide function (right tail) 
       *  @param left  left  tail
       */
      CrystalBallDoubleSided
      ( const Ostap::Math::CrystalBallRightSide& cb   , 
	      const Ostap::Math::LeftTail            & left ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double        pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0      () const { return m_core.m0      () ; }
      inline double mu      () const { return m_core.mu      () ; }
      inline double peak    () const { return m_core.peak    () ; }
      inline double mode    () const { return m_core.mode    () ; }
      inline double mass    () const { return m_core.mass    () ; }
      //
      inline double sigma   () const { return m_core .sigma  () ; }
      //
      inline double alphaL  () const { return m_left .alpha  () ; }
      inline double alphaR  () const { return m_right.alpha  () ; }      
      inline double nL      () const { return m_left .n      () ; }
      inline double nR      () const { return m_right.n      () ; }
      inline double NL      () const { return m_left .N      () ; }
      inline double NR      () const { return m_right.N      () ; }
      /// squared alphaL
      inline double alphaL2 () const { return m_left .alpha2 () ; }
      /// squared alphaR
      inline double alphaR2 () const { return m_right.alpha2 () ; }      
      // =====================================================================
    public:
      // =====================================================================
      /// the point where Gaussian meets Power-Law 
      inline double xL      () const { return m0 () - alphaL() * sigma () ; }
      /// the point where Gaussian meets Power-Law 
      inline double xR      () const { return m0 () + alphaR() * sigma () ; }
      // =====================================================================
    public: // trivial accessors
      // ======================================================================
      inline bool setM0     ( const double value ) { return m_core  .setM0    ( value ) ; }       
      inline bool setSigma  ( const double value ) { return m_core  .setSigma ( value ) ; } 
      inline bool setNL     ( const double value ) { return m_left  .setN     ( value ) ; }
      inline bool setNR     ( const double value ) { return m_right .setN     ( value ) ; }
      // ======================================================================
      inline bool setAlphaL ( const double value ) { return m_left  .setAlpha ( value ) ; }
      inline bool setAlphaR ( const double value ) { return m_right .setAlpha ( value ) ; }
      // ======================================================================
      inline bool setMu     ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass   ( const double value ) { return setM0 ( value ) ; }
      // ======================================================================
      inline bool setN
      ( const double valueL , 
	      const double valueR ) 
      {
	      const bool updatedL = m_left .setN ( valueL ) ;
	      const bool updatedR = m_right.setN ( valueR ) ;
	      return updatedL || updatedR ;
      }
      // ======================================================================      
      inline bool setAlpha
      ( const double valueL , 
	      const double valueR ) 
      {
	      const bool updatedL = m_left .setAlpha ( valueL ) ;
	      const bool updatedR = m_right.setAlpha ( valueR ) ;
	      return updatedL || updatedR ;
      }      
      // ======================================================================
      /// set both alpha-parameters 
      inline bool setAlpha ( const double value ) { return setAlpha ( value , value ) ; } 
      /// set both n-parameters 
      inline bool setN     ( const double value ) { return setN     ( value , value ) ; } 
      // ======================================================================      
    public: //
      // ======================================================================
      /// get integral between low and high
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
      /** get the integral from the the negative to positive infinity 
       *  @attention +infinity is returned for <code>n=0(N=1)</code>
       */      
      double integral () const ;
      // ======================================================================
    public: // components 
      // ======================================================================
      /// get the Gaussian core 
      const Ostap::Math::Gauss&     core       () const { return m_core  ; }
      /// get the Gaussian core 
      const Ostap::Math::Gauss&     gauss      () const { return m_core  ; }
      /// get left  tail
      const Ostap::Math::LeftTail&  tail_left  () const { return m_left  ; }
      /// get right tail
      const Ostap::Math::RightTail& tail_right () const { return m_right ; }
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       */
      double non_gaussian 
      ( const double xlow  ,
	      const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// core Gaussian 
      Ostap::Math::Gauss     m_core  {} ; // core gaussian
      /// (left) Tail
      Ostap::Math::LeftTail  m_left  {} ;
      /// (right) Tail
      Ostap::Math::RightTail m_right {} ;
      // =======================================================================
    } ;
    // ========================================================================
    /** @Class CrystalBallA 
     *  Variant of Crystal Ball function with asymmetric/bifurcated  core 
     *  @see Ostap::Math::CrystalBall 
     *  @see Ostap::Math::BifurcatedGauss 
     *  @see Ostap::Math::LeftTail 
     */
    class  CrystalBallA 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0       parameter
       *  @param sigmaL left sigma   parameter
       *  @param sigmaR right sigma  parameter
       *  @param alpha  alpha    parameter
       *  @param n      n        (external) parameter (not the same as *internal* N)
       */
      CrystalBallA 
      ( const double m0     = 0 ,
        const double sigmaL = 1 ,
        const double sigmaR = 1 ,
        const double alpha  = 2 ,
        const double n      = 1 ) ;
      // ======================================================================
      /** constructor from gaussian and tail parameeter 
       */
      CrystalBallA 
      ( const Ostap::Math::BifurcatedGauss& core      , 
	const double                        alpha = 2 , 
	const double                        n     = 1 ) ;
      // ======================================================================
      /** constructor from gaussian and tail 
       *  @parameter core bifurcated Gaussian function 
       *  @parameter tail tail     function 
       */
      CrystalBallA 
      ( const Ostap::Math::BifurcatedGauss& core ,
	const Ostap::Math::Tail&            tail ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double        pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0     () const { return m_core.m0     () ; }
      inline double mu     () const { return m_core.m0     () ; }
      inline double peak   () const { return m_core.m0     () ; }
      inline double sigmaL () const { return m_core.sigmaL () ; }
      inline double sigmaR () const { return m_core.sigmaR () ; }
      inline double alpha  () const { return m_tail.alpha  () ; }
      inline double n      () const { return m_tail.n      () ; }      
      /// Internal N-parameter
      inline double N      () const { return m_tail.N      () ; } // internal N parameter
      /// squared alpha 
      inline double alpha2 () const { return m_tail.alpha2 () ; }
      // ======================================================================
    public:            
      // ======================================================================
      /// mode of the distribution 
      inline double mode () const { return m_core.mode () ; }
      /// the point where Gaussian meets power-law
      inline double xL   () const { return m_core.m0 () - m_tail.alpha () * m_core.sigmaL ()  ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================== ==============================
      inline bool setM0     ( const double value  ) { return m_core.setM0     ( value ) ; }
      inline bool setSigmaL ( const double value  ) { return m_core.setSigmaL ( value ) ; }
      inline bool setSigmaR ( const double value  ) { return m_core.setSigmaR ( value ) ; }
      inline bool setSigma  ( const double valueL ,
			                        const double valueR ) { return m_core.setSigma  ( valueL , valueR  ) ; }
      inline bool setSigma  ( const double value  ) { return m_core.setSigma  ( value ) ; }
      inline bool setN      ( const double value  ) { return m_tail.setN      ( value ) ; } 
      inline bool setAlpha  ( const double value  ) { return m_tail.setAlpha  ( value ) ; } 
      // ======================================================================
      inline bool setMu    ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak  ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode  ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass  ( const double value ) { return setM0 ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral
      ( const double low ,
        const double high ) const ;
      // ======================================================================
      /** get the integral from the negative to positive infinity 
       *  @attention +infinity is returned for <code>n=0(N=1)</code>
       */
      double integral () const ;
      // ======================================================================
    public: // components 
      // ======================================================================
      /// get the bifurcated Gaussian core 
      const Ostap::Math::BifurcatedGauss& core      () const { return m_core ; }
      /// get left tail
      const Ostap::Math::LeftTail&        tail      () const { return m_tail ; }
      /// get left tail
      const Ostap::Math::LeftTail&        tail_left () const { return m_tail ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tail, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       */
      double non_gaussian 
      ( const double xlow  ,
	      const double xhigh ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// core Gaussian 
      Ostap::Math::BifurcatedGauss m_core  {} ; // core gaussian
      /// (left) Tail
      Ostap::Math::LeftTail        m_tail  {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallDoubleSidedA
     *  ``Crystal Ball-function'' for description of asymmetrui gauussia core with tails 
     *  @see CrystalBall
     *  @see CrystalBallRightSide
     *  @date 2011-05-25
     */
    class  CrystalBallDoubleSidedA
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0          parameter
       *  @param sigmaL left sigma  parameter
       *  @param sigmaR right sigma parameter
       *  @param alphaL alpha_L     parameter
       *  @param nL     n_L         parameter  (N-1 for "standard" definition)
       *  @param alphaR alpha_R parameter
       *  @param nR     n_R         parameter  (N-1 for "standard" definition)
       */
      CrystalBallDoubleSidedA
      ( const double m0      = 1 ,
        const double sigmaL  = 1 ,
        const double sigmaR  = 1 ,
        const double alphaL  = 2 ,
        const double nL      = 1 ,
        const double alphaR  = 2 ,
        const double nR      = 1 ) ;
      // ======================================================================
      /** constructor from all parameters
       *  @param core core Gaussian 
       *  @param alphaL  alpha_L     parameter
       *  @param nL      n_L         parameter  (N-1 for "standard" definition)
       *  @param alphaR  alpha_R parameter
       *  @param nR      n_R         parameter  (N-1 for "standard" definition)
       */
      CrystalBallDoubleSidedA
      ( const Ostap::Math::BifurcatedGauss& core , 
        const double                        alpha_L = 2 ,
        const double                        n_L     = 1 ,
        const double                        alpha_R = 2 ,
        const double                        n_R     = 1 ) ;
      // ========================================================================
      /** constructor from all components  
       *  @param core core Gaussian 
       *  @param left  left  tail
       *  @param right right tail
       */
      CrystalBallDoubleSidedA
      ( const Ostap::Math::BifurcatedGauss& core  ,
        const Ostap::Math::LeftTail&        left  ,
	      const Ostap::Math::RightTail&       right ) ;
      // ========================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double        pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0      () const { return m_core.m0      () ; }
      inline double mu      () const { return m_core.mu      () ; }
      inline double peak    () const { return m_core.peak    () ; }
      inline double mode    () const { return m_core.mode    () ; }
      inline double mass    () const { return m_core.mass    () ; }
      //
      inline double sigmaL  () const { return m_core .sigmaL () ; }
      inline double sigmaR  () const { return m_core .sigmaR () ; }
      //
      inline double alphaL  () const { return m_left .alpha  () ; }
      inline double alphaR  () const { return m_right.alpha  () ; }      
      inline double nL      () const { return m_left .n      () ; }
      inline double nR      () const { return m_right.n      () ; }
      inline double NL      () const { return m_left .N      () ; }
      inline double NR      () const { return m_right.N      () ; }
      /// squared alphaL
      inline double alphaL2 () const { return m_left .alpha2 () ; }
      /// squared alphaR
      inline double alphaR2 () const { return m_right.alpha2 () ; }      
      // =====================================================================
    public:
      // =====================================================================
      /// the point where Gaussian meets Power-Law 
      inline double xL      () const { return m_core.m0 () - m_left .alpha() * m_core.sigmaL () ; }
      /// the point where Gaussian meets Power-Law 
      inline double xR      () const { return m_core.m0 () + m_right.alpha() * m_core.sigmaR () ; }
      // =====================================================================
    public: // trivial accessors
      // ======================================================================
      inline bool setM0     ( const double value ) { return m_core  .setM0     ( value ) ; }       
      inline bool setNL     ( const double value ) { return m_left  .setN      ( value ) ; }
      inline bool setNR     ( const double value ) { return m_right .setN      ( value ) ; }
      // ======================================================================
      inline bool setSigmaL ( const double value ) { return m_core  .setSigmaL ( value ) ; }
      inline bool setSigmaR ( const double value ) { return m_core  .setSigmaR ( value ) ; }
      inline bool setAlphaL ( const double value ) { return m_left  .setAlpha  ( value ) ; }
      inline bool setAlphaR ( const double value ) { return m_right .setAlpha  ( value ) ; }
      // ======================================================================
      inline bool setMu     ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass   ( const double value ) { return setM0 ( value ) ; }
      // ======================================================================      
      inline bool setN
      ( const double valueL , 
	const double valueR ) 
      {
	const bool updatedL = m_left .setN ( valueL ) ;
	const bool updatedR = m_right.setN ( valueR ) ;
	return updatedL || updatedR ;
      }
      // ======================================================================      
      inline bool setAlpha
      ( const double valueL , 
	const double valueR ) 
      {
	const bool updatedL = m_left .setAlpha ( valueL ) ;
	const bool updatedR = m_right.setAlpha ( valueR ) ;
	return updatedL || updatedR ;
      }      
      /// set both alpha-parameters 
      inline bool setAlpha ( const double value )
      { return setAlpha ( value , value ) ; } 
      /// set both n-parameters 
      inline bool setN     ( const double value )
      { return setN     ( value , value ) ; }
      /// set both sigmas 
      inline bool setSigma
      ( const double valueL ,
	const double valueR )
      { return m_core  .setSigma ( valueL , valueR ) ; }
      /// set both sigmas
      inline bool setSigma  ( const double value  )
      { return m_core  .setSigma ( value ) ; }
      // ======================================================================
    public: //
      // ======================================================================
      /// get integral between low and high
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
      /** get the integral from the the negative to positive infinity 
       *  @attention +infinity is returned for <code>n=0(N=1)</code>
       */      
      double integral () const ;
      // ======================================================================
    public: // components 
      // ======================================================================
      /// get the Gaussian core 
      const Ostap::Math::BifurcatedGauss& core       () const { return m_core  ; }
      /// get left  tail
      const Ostap::Math::LeftTail&        tail_left  () const { return m_left  ; }
      /// get right tail
      const Ostap::Math::RightTail&       tail_right () const { return m_right ; }
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// core Bifurcatd Gaussian 
      Ostap::Math::BifurcatedGauss m_core  {} ; // core gaussian
      /// (left) Tail
      Ostap::Math::LeftTail        m_left  {} ;
      /// (right) Tail
      Ostap::Math::RightTail       m_right {} ;
      // =======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallDoubleSidedE
     *  `Crystal Ball-like function' : 
     *  - asymmetric core 
     *  - left powerlaw tail
     *  - right exponential tail
     *  @see CrystalBall
     *  @see CrystalBallRightSide
     *  @date 2011-05-25
     */
    class  CrystalBallDoubleSidedE
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0          parameter
       *  @param sigmaL left sigma  parameter
       *  @param sigmaR right sigma parameter
       *  @param alphaL alpha_L     parameter
       *  @param nL     n_L         parameter
       *  @param alphaR alpha_R     parameter
       */
      CrystalBallDoubleSidedE
      ( const double m0      = 1 ,
        const double sigmaL  = 1 ,
        const double sigmaR  = 1 ,
        const double alphaL  = 2 ,
        const double nL      = 1 ,
        const double alphaR  = 2 ) ; 
      // ======================================================================
      /** constructor from all parameters
       *  @param core core Gaussian 
       *  @param alphaL  alpha_L     parameter
       *  @param nL      n_L         parameter 
       *  @param alphaR  alpha_R parameter
       */
      CrystalBallDoubleSidedE
      ( const Ostap::Math::BifurcatedGauss& core        , 
        const double                        alpha_L = 2 ,
        const double                        n_L     = 1 ,
        const double                        alpha_R = 2 ) ;
      // ========================================================================
      /** constructor from all components  
       *  @param core core Gaussian 
       *  @param left  left  tail
       *  @param right right tail
       */
      CrystalBallDoubleSidedE
      ( const Ostap::Math::BifurcatedGauss& core  ,
        const Ostap::Math::LeftTail&        left  ,
	const Ostap::Math::RightExpTail&    right ) ;
      // ========================================================================
    public:
      // ======================================================================
      /// evaluate CrystalBall's function
      double        pdf        ( const double x ) const ;
      /// evaluate CrystalBall's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0      () const { return m_core.m0      () ; }
      inline double mu      () const { return m_core.mu      () ; }
      inline double peak    () const { return m_core.peak    () ; }
      inline double mode    () const { return m_core.mode    () ; }
      inline double mass    () const { return m_core.mass    () ; }
      //
      inline double sigmaL  () const { return m_core .sigmaL () ; }
      inline double sigmaR  () const { return m_core .sigmaR () ; }
      //
      inline double alphaL  () const { return m_left .alpha  () ; }
      inline double alphaR  () const { return m_right.alpha  () ; }      
      inline double nL      () const { return m_left .n      () ; }
      inline double NL      () const { return m_left .N      () ; }
      /// squared alphaL
      inline double alphaL2 () const { return m_left .alpha2 () ; }
      /// squared alphaR
      inline double alphaR2 () const { return m_right.alpha2 () ; }      
      // =====================================================================
    public:
      // =====================================================================
      /// the point where Gaussian meets Power-Law 
      inline double xL      () const { return m_core.m0 () - m_left .alpha() * m_core.sigmaL () ; }
      /// the point where Gaussian meets Power-Law 
      inline double xR      () const { return m_core.m0 () + m_right.alpha() * m_core.sigmaR () ; }
      // =====================================================================
    public: // trivial accessors
      // ======================================================================
      inline bool setM0     ( const double value ) { return m_core  .setM0     ( value ) ; }       
      inline bool setNL     ( const double value ) { return m_left  .setN      ( value ) ; }
      // ======================================================================
      inline bool setSigmaL ( const double value ) { return m_core  .setSigmaL ( value ) ; }
      inline bool setSigmaR ( const double value ) { return m_core  .setSigmaR ( value ) ; }
      inline bool setAlphaL ( const double value ) { return m_left  .setAlpha  ( value ) ; }
      inline bool setAlphaR ( const double value ) { return m_right .setAlpha  ( value ) ; }
      // ======================================================================
      inline bool setMu     ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass   ( const double value ) { return setM0 ( value ) ; }
      /// set n
      inline bool setN      ( const double value ) { return m_left.setN ( value ) ; } 
      /// set both alphas 
      inline bool setAlpha
      ( const double valueL , 
	const double valueR ) 
      {
	const bool updatedL = m_left .setAlpha ( valueL ) ;
	const bool updatedR = m_right.setAlpha ( valueR ) ;
	return updatedL || updatedR ;
      }      
      /// set both alpha-parameters 
      inline bool setAlpha ( const double value )
      { return setAlpha ( value , value ) ; } 
      /// set both sigmas 
      inline bool setSigma
      ( const double valueL ,
	const double valueR )
      { return m_core  .setSigma ( valueL , valueR ) ; }
      /// set both sigmas 
      inline bool setSigma  ( const double value  )
      { return m_core  .setSigma ( value ) ; }
      // ======================================================================
    public: //
      // ======================================================================
      /// get integral between low and high
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
      /** get the integral from the the negative to positive infinity 
       *  @attention +infinity is returned for <code>n=0(N=1)</code>
       */      
      double integral () const ;
      // ======================================================================
    public: // components 
      // ======================================================================
      /// get the Gaussian core 
      const Ostap::Math::BifurcatedGauss& core       () const { return m_core  ; }
      /// get left  tail
      const Ostap::Math::LeftTail&        tail_left  () const { return m_left  ; }
      /// get right tail
      const Ostap::Math::RightExpTail&    tail_right () const { return m_right ; }
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// core Bifurcatd Gaussian 
      Ostap::Math::BifurcatedGauss m_core  {} ; // core gaussian
      /// (left) Tail
      Ostap::Math::LeftTail        m_left  {} ;
      /// (right) Tail
      Ostap::Math::RightExpTail    m_right {} ;
      // =======================================================================
    } ;
    // ========================================================================
    /** @class Apollonios
     *  "Bifurcated Apollonios"
     *  A modified gaussian with asymmetric exponential tails on both sides
     *
     *  \f[  f(x; \mu, \sigma_L, \sigma_R, \beta) \propto \mathrm{e}^{\beta^\prime\ledt(\beta - \sqrt{\beta^2+\delta^2 }\right) } \f]
     *  where
     *  - \f$ \beta_2 = \sqrt{ 2 + \beta^2}\f$
     *  - \f$ \delta  = \frac{x-\mu}{\sigma_L}\f$ for \f$ x<\mu \f$, otherwise \f$ \delta = \frac{x-mu}{\sigma_R}\f$
     *
     *  Well defined limits are: 
     *  - \f$ \beta \rightharrow +\infty \f$ : Bifucated Gaussian with \f$\sigma_{L,R}\f$ 
     *  - \f$ \beta \rightarrow 0        \f$ : Asymmetric Laplace with slopes of \f$\sqrt{2}\sigma_{L,R}\f$ 
     *
     *  The function is largely inspired by:
     *  @see https://indico.cern.ch/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=262633
     *  @see http://arxiv.org/abs/1312.5000
     *
     *  The adopted parameterization
     *  - more generic, allowing intrinsic asymmetries  
     *  - have well-defined limits  
     *  - RMS only sightly depends on \f$ \beta \f$ in region \f$ \beta \sim O(1) \f$ 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date  2013-12-01
     */
    class  Apollonios
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0      m0        parameter
       *  @param sigmaL  sigmaL    parameter
       *  @param alphaR  alphaR    parameter
       *  @param beta    beta      parameter
       */
      Apollonios
      ( const double m0      = 0 ,
	const double sigmaL  = 1 ,
	const double alphaR  = 1 ,
	const double beta    = 1 ) ;  // large beta correponds to gaussian
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Apollonios's function
      double        pdf        ( const double x ) const ;
      /// evaluate Apollonios's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0     () const { return m_m0     ; }
      inline double mu     () const { return m_m0     ; }
      inline double peak   () const { return m_m0     ; }
      inline double mode   () const { return m_m0     ; }
      inline double sigmaL () const { return m_sigmaL ; }
      inline double sigmaR () const { return m_sigmaR ; }
      inline double beta   () const { return m_beta   ; }
      // ======================================================================
      inline double sigma  () const { return 0.5 * ( m_sigmaL + m_sigmaR )           ; }
      inline double asym   () const { return 0.5 * ( m_sigmaL - m_sigmaR ) / sigma() ; }
      inline double beta2  () const { return m_beta * m_beta ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      bool        setM0     ( const double value ) ;
      bool        setSigmaL ( const double value ) ;
      bool        setSigmaR ( const double value ) ;
      bool        setBeta   ( const double value ) ;
      // ======================================================================
      inline bool setMu     ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode   ( const double value ) { return setM0 ( value ) ; }
      // ======================================================================
      /// set both sigmas simultaneously 
      bool        setSigma
      ( const double valueL ,
	const double valueR ) ;
      /// set both sigmas simultaneously 
      inline bool setSigma ( const double value )
      { return setSigma ( value , value ) ; }
      // ======================================================================      
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public: // logarithmic derivative 
      // ======================================================================
      /** log-derivative \f$ \frac{ f^\prime}{f}  \f$
       *  Useful to attach the radiative tail to ensure 
       *  the continuity of the function and the 1st derivatibve 
       */
      double dFoF ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the peak position
      double m_m0      { 0 } ; // the peak position
      /// the peak resolution
      double m_sigmaL  { 1 } ; // the peak resolution
      /// the peak resolution
      double m_sigmaR  { 1 } ; // the peak resolution
      /// parameter beta
      double m_beta    { 1 } ; // parameter beta
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ApolloniosL
     *  An Apollonious core with power-law tail on left side
     *  @see  Ostap::Math::Apollonios
     *  @see  Ostap::Math::LeftTail 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date  2013-12-01
     */
    class  ApolloniosL
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0     m0       parameter
       *  @param sigmaL sigmaL   left-sigma parameter
       *  @param sigmaR sigmaR   right-sigma parameter
       *  @param beta   beta     parameter
       *  @param alpha  alpha    tails alpha-parameter
       *  @param n      n        tail  n-parameter
       */
      ApolloniosL
      ( const double m0     = 0 ,
        const double sigmaL = 1 ,
        const double sigmaR = 1 ,
        const double beta   = 1 , 
	const double alpha  = 2 , 
        const double n      = 1 ) ;
      // =====================================================================
      /** constructor from two parameters
       *  @param core core component 
       *  @param tail left-tail component 
       */
      ApolloniosL
      ( const Ostap::Math::Apollonios& core ,
	const Ostap::Math::Tail&       tail ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Apollonios's function
      double        pdf        ( const double x ) const ;
      /// evaluate Apollonios's function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline double m0     () const { return m_core.m0     () ; }
      inline double mu     () const { return m_core.mu     () ; }
      inline double peak   () const { return m_core.peak   () ; }
      inline double mode   () const { return m_core.mode   () ; }
      inline double sigmaL () const { return m_core.sigmaL () ; }
      inline double sigmaR () const { return m_core.sigmaR () ; }
      inline double alpha  () const { return m_tail.alpha  () ; }
      inline double n      () const { return m_tail.n      () ; }
      inline double N      () const { return m_tail.n      () ; }
      // ======================================================================
      /// the point where core meats power-law tail 
      inline double xL     () const { return m0 ()  - alpha() * sigmaL() ; } 
      // ======================================================================
    public: // trivial accessors
      // ======================================================================
      inline bool setM0     ( const double value ) { return m_core.setM0     ( value ) ; }
      inline bool setSigmaL ( const double value ) { return m_core.setSigmaL ( value ) ; }
      inline bool setSigmaR ( const double value ) { return m_core.setSigmaR ( value ) ; }
      inline bool setBeta   ( const double value ) { return m_core.setBeta   ( value ) ; }
      inline bool setAlpha  ( const double value ) { return m_tail.setAlpha  ( value ) ; }
      inline bool setN      ( const double value ) { return m_tail.setN      ( value ) ; }
      // ======================================================================
      inline bool setMu     ( const double value ) { return setM0 ( value ) ; }
      inline bool setPeak   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMass   ( const double value ) { return setM0 ( value ) ; }
      inline bool setMode   ( const double value ) { return setM0 ( value ) ; }
      // ======================================================================
      /// set both sigmas simultaneously 
      inline bool setSigma
      ( const double valueL ,
	const double valueR )
      { return m_core.setSigma ( valueL , valueR ) ; } 
      /// set both sigmas simultaneously 
      inline bool setSigma ( const double value )
      { return m_core.setSigma ( value ) ; } 
      // ======================================================================      
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral
      ( const double low ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the core 
      const Ostap::Math::Apollonios& core      () const { return m_core ; }
      /// get the (left) tail 
      const Ostap::Math::LeftTail&   tail      () const { return m_tail ; } 
      /// get the  left  tail 
      const Ostap::Math::LeftTail&   tail_left () const { return m_tail ; } 
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// core 
      Ostap::Math::Apollonios m_core {} ;
      /// the tail
      Ostap::Math::LeftTail   m_tail {} ; 
      // =======================================================================
    } ;
    // ========================================================================
    /** @class StudentT
     *  Simple function to parameterize the symmetric peak using
     *  scale-location version of Student' t-distribution
     *
     *  \f[  f(y|\nu,\mu,\sigma) = 
     *   \frac{1}{\sqrt{\pi \nu}} \frac { \Gamma( \frac{\nu+1}{2}) } { \Gamma( \frac{\nu}{2}  ) }
     *  \left( 1 + \frac{y^2}{\nu} \right)^{ -\frac{\nu+1}{2}} \f],
     *
     *  where \f$ y = \frac{x - \mu}{\sigma} \f$
     *
     *  - since we want to have the finite integral and finite variance,
     *    we use n-parameter, such that \f$ \nu = \nu(n) \ge 2 \f$
     *
     *  @see https://en.wikipedia.org/wiki/Student%27s_t-distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-01-05
     */
    class  StudentT
    {
    public:
      // ======================================================================
      /** constructor from mass, scale and "n"-parameter
       *  @param mass  mass/location 
       *  @param scale scale parameter
       *  @param n     n-parameter  ( actually nu=nu(n) )
       */
      StudentT 
      ( const double mass  = 0 ,
        const double scale = 1 ,
        const double n     = 2 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate StudentT's shape
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      /// calculate StudentT's shape
      inline double evaluate   ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // variables
      // ======================================================================
      inline double M        () const  { return m_M      ; }
      inline double m0       () const  { return m_M      ; }
      inline double mu       () const  { return m_M      ; }
      inline double mass     () const  { return m_M      ; }
      inline double peak     () const  { return m_M      ; }
      inline double mode     () const  { return m_M      ; }
      inline double location () const  { return m_M      ; }
      // ======================================================================
      inline double scale    () const  { return m_scale  ; }
      inline double sigma    () const  { return m_scale  ; }
      inline double tau      () const  { return m_scale  ; }
      inline double gamma    () const  { return m_scale  ; }
      inline double width    () const  { return m_scale  ; }
      // ======================================================================
      inline double n        () const  { return m_n      ; }
      // ======================================================================
      inline double nu       () const  { return m_nu     ; }
      // ======================================================================
    public:
      // ======================================================================
      bool        setM        ( const double value  ) ; // setter for location 
      bool        setScale    ( const double value  ) ; // setter for scale/width 
      bool        setN        ( const double value  ) ; // setter for n
      // ======================================================================
      inline bool setM0       ( const double value  ) { return setM     ( value ) ; }
      inline bool setMu       ( const double value  ) { return setM     ( value ) ; }
      inline bool setMass     ( const double value  ) { return setM     ( value ) ; }
      inline bool setPeak     ( const double value  ) { return setM     ( value ) ; }
      inline bool setMode     ( const double value  ) { return setM     ( value ) ; }
      inline bool setLocation ( const double value  ) { return setM     ( value ) ; }
      // ======================================================================
      inline bool setSigma    ( const double value  ) { return setScale ( value ) ; }
      inline bool setTau      ( const double value  ) { return setScale ( value ) ; }      
      inline bool setGamma    ( const double value  ) { return setGamma ( value ) ; }
      inline bool setWidth    ( const double value  ) { return setWidth ( value ) ; }
      // ======================================================================
    public: // some statistics 
      // ======================================================================
      /// the mean value 
      inline double mean       () const { return m_M ; }
      /// skewness 
      inline double skewness   () const { return 0   ; } 
      /// variance
      double        variance   () const ;
      /// (excess) kurtosis 
      double        kurtosis   () const ;
      /// RMS
      double        rms        () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the PDF 
      double pdf    ( const double x ) const ;
      /// get the PDF 
      double cdf    ( const double x ) const ;
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
    public:
      // ======================================================================
      /// get the expression nu=nu(n) 
      static double nu ( const double n ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// mass
      double m_M     { 0  } ; //
      /// width parameter
      double m_scale { 1  } ; // width parameter
      /// n-parameter
      double m_n     { 0  } ; // n-parameter
      /// nu-parameter
      double m_nu    { 2  } ; // nu-parameter
      // ======================================================================
    private: // normalization
      // ======================================================================
      double m_norm  { -1 } ;
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
       *  @param sigmaL left scale/width parameter
       *  @param sigmaR right scale/width parameter
       *  @param nL     left n-parameter  ( actually  nuL = nu(nL) )
       *  @param nR     right n-parameter  ( actually nuR = nu(nR) )
       */
      BifurcatedStudentT
      ( const double mass   = 0 ,
        const double sigmaL = 1 ,
        const double sigmaR = 1 ,
        const double nL     = 2 ,
        const double nR     = 2 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate bifurcated StudentT's shape
      inline double operator() ( const double x ) const { return pdf ( x )  ; }
      /// calculate bifurcated StudentT's shape
      inline double evaluate   ( const double x ) const { return pdf ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      // variables
      // ======================================================================
      inline double M       () const  { return m_M ; }
      inline double m0      () const  { return m_M ; }
      inline double mu      () const  { return m_M ; }
      inline double mass    () const  { return m_M ; }
      inline double peak    () const  { return m_M ; }
      inline double mode    () const  { return m_M ; }
      // ======================================================================
      inline double sigmaL  () const  { return m_sL      ; }
      inline double sL      () const  { return sigmaL () ; }
      inline double gammaL  () const  { return sigmaL () ; }
      inline double widthL  () const  { return sigmaL () ; }
      // ======================================================================
      inline double sigmaR  () const  { return m_sR      ; }
      inline double sR      () const  { return sigmaR () ; }
      inline double gammaR  () const  { return sigmaR () ; }
      inline double widthR  () const  { return sigmaR () ; }
      // ======================================================================
      inline double nL      () const  { return m_nL      ; }
      inline double nR      () const  { return m_nR      ; }
      // ======================================================================
      inline double nuL     () const  { return m_nuL     ; }
      inline double nuR     () const  { return m_nuR     ; }
      // ======================================================================
    public:
      // ======================================================================
      bool        setM      ( const double value  ) ;
      bool        setSigmaL ( const double value  ) ;
      bool        setSigmaR ( const double value  ) ;
      bool        setNL     ( const double value  ) ; // setter for nL 
      bool        setNR     ( const double value  ) ; // setter for nR 
      // ======================================================================      
      inline bool setM0     ( const double value  ) { return setM  ( value ) ; }
      inline bool setMu     ( const double value  ) { return setM  ( value ) ; }
      inline bool setMass   ( const double value  ) { return setM  ( value ) ; }
      inline bool setPeak   ( const double value  ) { return setM  ( value ) ; }
      inline bool setMode   ( const double value  ) { return setM  ( value ) ; }
      // ======================================================================
      inline bool setSL     ( const double value  ) { return setSigmaL ( value ) ; }
      inline bool setGammaL ( const double value  ) { return setSigmaL ( value ) ; }
      inline bool setWidthL ( const double value  ) { return setSigmaL ( value ) ; }
      // ======================================================================
      inline bool setSR     ( const double value  ) { return setSigmaR ( value ) ; }
      inline bool setGammaR ( const double value  ) { return setSigmaR ( value ) ; }
      inline bool setWidthR ( const double value  ) { return setSigmaR ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// PDF 
      double pdf    ( const double x ) const ;
      /// CDF 
      double cdf    ( const double x ) const ;
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
      /// mass
      double m_M     { 0 } ; //
      /// width parameter
      double m_sL    { 1 } ; // width parameter
      /// width parameter
      double m_sR    { 1 } ; // width parameter
      /// nL-parameter
      double m_nL    { 0 } ; // nL-parameter
      /// nR-parameter
      double m_nR    { 0 } ; // nR-parameter
      // ======================================================================
      /// nuL parameter 
      double m_nuL   { 2 } ;  // nuL parameter 
      /// nuR parameter 
      double m_nuR   { 2 } ; // nuR parameter  
      // ======================================================================
    private: // normalization
      // ======================================================================
      double m_normL { -1 } ;
      double m_normR { -1 } ;
      // ======================================================================
    } ;
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
       *  @param sigma width/scale paramete   ( == a       in Heinrich note ) 
       *  @param n     n-parameter            ( == m - 1/2 in Henrich' note )
       *  @param kappa asymmetry parameter    ( == nu      in Heinrich note )
       */
      PearsonIV
      ( const double mu    = 0 ,  // == mu 
        const double sigma = 1 ,  // == a 
        const double n     = 2 ,  // == m - 1/2 
        const double kappa = 0 ); // == nu 
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
      inline double mu       () const { return m_mu       ; } ;
      /// width/scale parameter 
      inline double varsigma () const { return m_varsigma ; } ;
      /// n-parameter 
      inline double n        () const { return m_n        ; } ;
      /// 
      inline double kappa    () const { return m_kappa    ; } ;
      // ======================================================================      
    public : // derived parameters 
      // ======================================================================      
      /// parameteter m
      inline double m  () const { return m_n + 0.5        ; }
      /// parameteter nu 
      inline double nu () const { return m_kappa          ; }
      /// parameter r 
      inline double r  () const { return 2 * ( m () - 1 ) ; }
      /// parameter a  
      inline double a  () const { return  m_varsigma      ; }
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
       *  the twp points are equidstance from the mode 
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
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-valuw with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
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
      SinhAsinh  
      ( const double location  = 0   ,
        const double scale     = 1   ,
        const double epsilon   = 0   ,
        const double delta     = 1   ) ;
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
      inline double location () const { return mu    () ; }
      inline double scale    () const { return sigma () ; }
      // ======================================================================
      inline double mu       () const { return m_mu      ; }
      inline double sigma    () const { return m_sigma   ; }
      inline double epsilon  () const { return m_epsilon ; }
      inline double delta    () const { return m_delta   ; }
      inline double peak     () const { return m_mu      ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool        setMu       ( const double value ) ;
      bool        setSigma    ( const double value ) ;
      bool        setEpsilon  ( const double value ) ;
      bool        setDelta    ( const double value ) ;
      // ======================================================================
      inline bool setLocation ( const double value ) { return setMu    ( value ) ; }
      inline bool setScale    ( const double value ) { return setSigma ( value ) ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// get median
      double median   () const ;
      /// get the mean for the distgribution
      double mean     () const ;
      /// get the variance for the distribution
      double variance () const ;
      /// get the RMS for the distribution
      double rms      () const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get the CDF 
      double cdf      ( const double x    ) const ;
      /// get the integral 
      double integral ()                    const ; 
      /// get the integral      
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	      const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mu      { 0 } ;
      double m_sigma   { 1 } ;
      double m_epsilon { 0 } ;
      double m_delta   { 1 } ;
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
      JohnsonSU  
      ( const double xi      = 0 ,   // related to location
        const double lambda  = 1 ,   // related to variance
        const double delta   = 1 ,   // shape
        const double gamma   = 0 ) ; // shape
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate JohnsonSU-distributions
      double pdf        ( const double x ) const ;
      /// evaluate JohnsonSU-distributions
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      inline double xi       () const { return m_xi       ; }
      inline double lam      () const { return m_lambda   ; }
      inline double lambda   () const { return m_lambda   ; }
      inline double lambd    () const { return m_lambda   ; }
      inline double delta    () const { return m_delta    ; }
      inline double gamma    () const { return m_gamma    ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setXi       ( const double value ) ;
      bool setLambda   ( const double value ) ;
      bool setDelta    ( const double value ) ;
      bool setGamma    ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      double        mean       () const ;
      double        variance   () const ;
      // ======================================================================      
      inline double dispersion () const { return variance () ; }
      inline double sigma      () const { return std::sqrt ( variance() ) ; }
      inline double rms        () const { return std::sqrt ( variance() ) ; }
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get CDF 
      double cdf      ( const double x    ) const ;
      /// get the integral 
      double integral ()                    const ;
      /// get the integral 
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	      const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_xi      { 0 } ;
      double m_lambda  { 1 } ;
      double m_delta   { 1 } ;
      double m_gamma   { 0 } ;
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
      Atlas   
      ( const double mean   = 0  ,
        const double sigma  = 1  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate atlas function
      double pdf        ( const double x ) const ;
      /// evaluate atlas function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      inline double mean   () const { return m_mean   ; }
      /// get parameters "sigma"
      inline double sigma  () const { return m_sigma  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// peak position
      inline double peak     () const { return mean () ; }
      /// get mode
      inline double mode     () const { return mean () ; }
      // ======================================================================
    public:
      // ======================================================================      
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
      double integral
      ( const double low  ,
        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
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
      Sech 
      ( const double mean   = 0  ,
        const double sigma  = 1  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate sech function
      double pdf        ( const double x ) const ;
      /// evaluate sech function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      inline double mean   () const { return m_mean ; }
      inline double peak   () const { return m_mean ; }
      inline double mu     () const { return m_mean ; }
      inline double m0     () const { return m_mean ; }
      /// get parameters "sigma"
      double sigma  () const { return m_sigma   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get mode
      inline double mode     () const { return m_mean ; }
      /// get variance
      inline double variance () const { return m_sigma * m_sigma ; }
      /// get rms
      inline double rms      () const { return m_sigma           ; }
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
      double integral
      ( const double low  ,
        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const ;
      /// evaluate atlas function
      double cdf      ( const double x ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
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
      Logistic  
      ( const double mean  = 0  ,
        const double sigma = 1  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate sech function
      double pdf        ( const double x ) const ;
      /// evaluate sech function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      inline double mean   () const { return m_mean   ; }
      /// get parameters "sigma"
      inline double sigma  () const { return m_sigma  ; }
      /// get "s"
      double s      () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get peak
      inline double peak     () const { return peak () ; }
      /// get mode
      inline double mode     () const { return mean () ; }
      /// get median
      inline double median   () const { return mean () ; }
      /// get variance
      inline double variance () const { return m_sigma * m_sigma ; }
      /// get rms
      inline double rms      () const { return m_sigma ; }
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
      double integral
      ( const double low  ,
        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const ;
      /// evaluate Logistc CDF function
      double cdf      ( const double x ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
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
    /** @class GenLogisticIV
     *  GeneraliZed Logistic Type IV distribution with location/scale 
     *  https://en.wikipedia.org/wiki/Generalized_logistic_distribution
     *  - Type I   : beta  = 1 
     *  - Type II  : alpha = 1 
     *  - Type III : alpha = beta 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2023-05-09
     */
    class  GenLogisticIV
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param mu    \f$\mu  \f$-parameter
       *  @param sigma \f$sigma\f$-parameter
       *  @param alpha \f$alpha\f$-parameter
       *  @param beta  \f$beta\f$-parameter
       */
      GenLogisticIV
      ( const double mu    = 0 ,   // location 
        const double sigma = 1 ,   // scale 
        const double alpha = 1 ,   // expo tail 
        const double beta  = 1 ) ; // expo tail 
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function
      double evaluate          ( const double x ) const ;
      /// evaluate the function
      inline double pdf        ( const double x ) const { return evaluate ( x ) ;  }
      /// evaluate the function
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      /// get parameter mu 
      inline double mu     () const { return m_mu     ; }
      /// get parameter "sigma"
      inline double sigma  () const { return m_sigma  ; }
      /// get parameter "alpha"
      inline double alpha  () const { return m_alpha  ; }
      /// get parameter "beta"
      inline double beta   () const { return m_beta   ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool   setMu    ( const double value ) ;
      bool   setSigma ( const double value ) ;
      bool   setAlpha ( const double value ) ;
      bool   setBeta  ( const double value ) ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// get the mean value 
      inline double mean       () const { return m_mu ; }
      /// get the variance 
      inline double variance   () const { return m_sigma * m_sigma  ; }
      /// get the dispersion 
      inline double dispersion () const { return variance ()   ; }
      /// get the RMS
      inline double rms        () const { return m_sigma ; }
      /// get the skewness 
      double skewness          () const ;
      /// get the (Excess) kurtosis 
      double kurtosis          () const ;
      /// get the mode
      double mode              () const ;
      /// get cumulant 
      double cumulant          ( const unsigned short k ) const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get integral from low to high
      double integral 
      ( const double low  ,
        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the helper variable y 
      double y ( const double z ) const ;
      /// get z from the y 
      double z ( const double y ) const ;
      /// get the "standard" generalized Type IV Logistic distribution 
      double std_type4 ( const double t ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// parameteter "mu"
      double m_mu          ; // parameter mu,mean,mode
      /// parameter   "sigma"
      double m_sigma       ; // parameter sigma
      /// parameter   "alpha"
      double m_alpha       ; // parameter alpha
      /// parameter   "beta"
      double m_beta        ; // parameter beta
      /// tilda-mu 
      double m_tilda_mu    ;
      /// tilda-s2  
      double m_tilda_s     ;
      ///  normalization : 1 / B(alpha,beta)
      double m_norm { -1 } ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Losev 
     *  `Losev distribution' - asymmetric variant of Hyperbolic secant/Sech-function
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
      Losev
      ( const double mu    = 0 , 
        const double alpha = 1 , 
        const double beta  = 1 ) ;        
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      /// evaluate the function 
      double pdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// parameter mu 
      inline double mu    () const { return m_mu    ; }
      /// parameter alpha
      inline double alpha () const { return m_alpha ; }
      /// parameter beta 
      inline double beta  () const { return m_beta  ; }      
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
      double integral
      ( const double low  , 
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
      Slash
      ( const double mu    = 0 ,   // location 
        const double scale = 1 ) ; // scale ;      
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate slash function
      double pdf        ( const double x ) const ;
      /// evaluate slash function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      inline double mu       () const { return m_mu    ; }
      /// get parameters "scale"
      inline double scale    () const { return m_scale ; }
      // ======================================================================
    public: //   derived getters  
      // ======================================================================
      inline double m0       () const { return m_mu    ; }
      inline double mode     () const { return m_mu    ; }
      inline double mean     () const { return m_mu    ; }
      inline double peak     () const { return m_mu    ; }
      inline double location () const { return m_mu    ; }
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool        setMu       ( const double value ) ;
      bool        setScale    ( const double value ) ;
      // ======================================================================      
      inline bool setM0       ( const double value ) { return setMu ( value ) ; }
      inline bool setMean     ( const double value ) { return setMu ( value ) ; }
      inline bool setMode     ( const double value ) { return setMu ( value ) ; }
      inline bool setPeak     ( const double value ) { return setMu ( value ) ; }
      inline bool setLocation ( const double value ) { return setMu ( value ) ; }
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get integral from low to high
      double integral
      ( const double low  ,
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
     *  @see https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution 
     */ 
    class AsymmetricLaplace
    {
    public:
      // ======================================================================
      /** constructor from all parameters 
       *  @param mu  peak location
       *  @param lambdaL `left'  exponential slope  (lambdaL>0)
       *  @param lambdaR `right' exponential slope  (lambdaR>0)
       */
      AsymmetricLaplace 
      ( const double mu      = 0 ,   // location 
        const double lambdaL = 1 ,   // left  exponential slope 
        const double lambdaR = 1 ) ; // right exponential slope 
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate asymmetic laplace function
      double pdf        ( const double x ) const ;
      /// evaluate asymmetric laplace function
      inline double operator() ( const double x ) const { return pdf ( x ) ; }
      // ======================================================================
    public: // direct getters
      // ======================================================================
      inline double mu       () const { return m_mu      ; }
      /// get `left' lambda parameter 
      inline double lambdaL  () const { return m_lambdaL ; }
      /// get `right' lambda parameter 
      inline double lambdaR  () const { return m_lambdaR ; }
      ///  squared lambda_L
      inline double lambdaL2 () const { return m_lambdaL * m_lambdaL ; }
      ///  squared lambda_R
      inline double lambdaR2 () const { return m_lambdaR * m_lambdaR ; }
      // ======================================================================
    public: //   derived getters  
      // ======================================================================
      inline double mode     () const { return m_mu ; }
      inline double peak     () const { return m_mu ; }
      inline double location () const { return m_mu ; }      
      inline double m0       () const { return m_mu ; }      
      // ======================================================================
    public: // the standard parameterization (slopes are inverse)
      // ======================================================================
      ///  \f$ \lambda ^2 \f$ 
      inline double lambda2  () const { return 1.0 / ( m_lambdaL * m_lambdaR ) ; } 
      inline double Lambda2  () const { return lambda2 () ; } 
      inline double Lambda   () const { return std::sqrt ( lambda2 () ) ; }
      inline double lambda   () const { return Lambda  () ; }
      /// get the `asymmetry' \f$ k^2 \f$ 
      inline double k2       () const { return m_lambdaL / m_lambdaR ; }
      /// get the `asymmetry' \f$ k  \f$ 
      inline double k        () const { return std::sqrt ( k2 () )   ; }      
      // ======================================================================
    public: // direct setters
      // ======================================================================
      bool        setMu       ( const double value ) ;
      bool        setLambdaL  ( const double value ) ;
      bool        setLambdaR  ( const double value ) ;
      // ======================================================================      
      inline bool setM0       ( const double value ) { return setMu ( value ) ; }
      inline bool setMode     ( const double value ) { return setMu ( value ) ; }
      inline bool setPeak     ( const double value ) { return setMu ( value ) ; }
      inline bool setLocation ( const double value ) { return setMu ( value ) ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// get mean 
      double mean     () const ;
      /// get median 
      double median   () const ; 
      /// get variance
      double variance () const ;
      /// get skewness 
      double skewness () const ;
      /// get kurtosis 
      double kurtosis () const ;
      /// get RMS
      inline double rms        () const { return std::sqrt ( variance () ) ; };
      /// get sigma
      inline double sigma      () const { return rms      () ; }
      // get dispersion
      inline double dispersion () const { return variance () ; }
      // ======================================================================
    public: // integrals
      // ======================================================================
      /// get integral from low to high
      double integral
      ( const double low  ,
        const double high ) const ;
      /// integral from -infinity to +infinity
      double integral () const {  return  1 ; }
      /// evaluate CDF function
      double cdf      ( const double x ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      // /** quantify the effect of the tails, difference from Gaussian
      //  *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
      //  * where 
      //  * - \f$ I_{CB} \f$ is integral over Gaussian function 
      //  * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
      //  * - Gaussian is centered at mean-value with sigma = RMS 
      //  */
      // double non_gaussian 
      // ( const double xlow  ,
      // 	const double xhigh ) const ;
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
      RaisingCosine
      ( const double mu = 0 , 
        const double s  = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate raising cosine distribution
      double pdf (  const double x ) const ;
      /// evaluate raising cosine distribution
      inline double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // main getters 
      // ======================================================================
      inline double mu    () const { return  m_mu ; }
      inline double s     () const { return  m_s  ; }
      // ======================================================================
    public: // derived getters 
      // ======================================================================
      inline double location () const { return m_mu ; }
      inline double peak     () const { return m_mu ; }
      inline double m0       () const { return m_mu ; }
      inline double scale    () const { return m_s  ; }
      // ======================================================================
    public: //  setters 
      // ======================================================================
      bool        setMu       ( const double value ) ;
      bool        setS        ( const double value ) ;
      inline bool setScale    ( const double value ) { return setS  ( value ) ; }
      inline bool setLocation ( const double value ) { return setMu ( value ) ; }      
      inline bool setMean     ( const double value ) { return setMu ( value ) ; }      
      // ======================================================================
    public: // derived getters & stats 
      // ======================================================================
      /// mean 
      inline double mean     () const {  return m_mu ; }
      /// mode 
      inline double mode     () const {  return m_mu ; }
      ///  median
      inline double median   () const {  return m_mu ; }
      ///  variance  
      double variance () const ;
      ///  rms   
      double rms      () const ;
      /// kurtosis
      double kurtosis () const ;
      /// skewness
      inline double skewness () const { return 0 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get CDF 
      double cdf      ( const double x ) const ;
      /// evaluate the integral
      double integral ()                 const ;
      /// evaluate the integral
      double integral
      ( const double low  , 
        const double high ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
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
     *  q-Gaussian (Tsallis) distribution:
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
       *  @param scale scale parmaeter/sigma 
       *  @param q     q-value   (q<3, for q>3 q is set to be  6-q)
       */
      QGaussian 
      ( const double mean  = 0 ,   // mean/mode/location 
        const double scale = 1 ,   // scale/sigma
        const double q     = 1 );  // q-parameter/shape 
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate  pdf  for q-Gaussian distribution
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  for q-Gaussian distribution
      inline double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public : // primary getters 
      // ======================================================================
      inline double mean   () const { return  m_mean  ; }
      inline double scale  () const { return  m_scale ; }
      inline double q      () const { return  m_q     ; }
      // ======================================================================
    public : // derived getters 
      // ======================================================================
      inline double peak     () const { return mean  () ; }
      inline double mu       () const { return mean  () ; }
      inline double mode     () const { return mean  () ; }
      inline double median   () const { return mean  () ; }
      inline double location () const { return mean  () ; }
      inline double sigma    () const { return scale () ; }
      /// get the original beta 
      inline double beta     () const { return 0.5 / ( m_scale * m_scale ) ; }
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
      inline bool setMu       ( const double value ) { return setMean  ( value )  ; }
      inline bool setM0       ( const double value ) { return setMean  ( value )  ; }
      inline bool setMode     ( const double value ) { return setMean  ( value )  ; }
      inline bool setPeak     ( const double value ) { return setMean  ( value )  ; }
      inline bool setLocation ( const double value ) { return setMean  ( value )  ; }
      inline bool setSigma    ( const double value ) { return setScale ( value )  ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral () const { return 1 ; }
      /// get the integral 
      double integral 
      ( const double low  , 
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
      /// scale/sigma 
      double m_scale ; // scale/sigma
      /// q-value 
      double m_q     ; // q-value 
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
    /** @class KGaussian  
     *  k-Gaussiand (Kaniadakis) distribution:
     *  @see https://en.wikipedia.org/wiki/Kaniadakis_Gaussian_distribution
     *  Here we use \f$ k = \tanh { \kappa } \f$
     */
    class KGaussian  
    {
    public:
      // ======================================================================
      /** constructor from all arguments            
       *  @param mean  the mean/mode/location of the peak 
       *  @param scale scale parmaeter/sigma 
       *  @param kappa kappa-value
       */
      KGaussian 
      ( const double mean  = 0 ,   // mean/mode/location 
        const double scale = 1 ,   // scale/sigma
        const double kappa = 0 );  // kappa-parameter/shape 
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate  pdf  for k-Gaussian distribution
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  for k-Gaussian distribution
      inline double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public : // primary getters 
      // ======================================================================
      inline double mean   () const { return  m_mean  ; }
      inline double scale  () const { return  m_scale ; }
      inline double k      () const { return  m_k     ; }
      inline double kappa  () const { return  m_kappa ; }
      // ======================================================================
    public : // derived getters 
      // ======================================================================
      inline double peak     () const { return mean  () ; }
      inline double mu       () const { return mean  () ; }
      inline double m0       () const { return mean  () ; }
      inline double mode     () const { return mean  () ; }
      inline double median   () const { return mean  () ; }
      inline double location () const { return mean  () ; }
      inline double sigma    () const { return scale () ; }
      /// get the original beta 
      inline double beta     () const { return 0.5 / ( m_scale * m_scale ) ; }
      // ======================================================================
    public : // primay getters 
      // ======================================================================
      // set mean 
      bool setMean  ( const double value ) ;
      // set kappa 
      bool setKappa ( const double value ) ;
      // set scale
      bool setScale ( const double value ) ;
      // ======================================================================
    public : // derived setters 
      // ======================================================================
      inline bool setMu       ( const double value ) { return setMean  ( value )  ; }
      inline bool setM0       ( const double value ) { return setMean  ( value )  ; }
      inline bool setMode     ( const double value ) { return setMean  ( value )  ; }
      inline bool setPeak     ( const double value ) { return setMean  ( value )  ; }
      inline bool setLocation ( const double value ) { return setMean  ( value )  ; }
      inline bool setSigma    ( const double value ) { return setScale ( value )  ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get variance 
      double variance   () const ;
      /// get dispersion 
      double dispersion () const { return variance   () ; }
      /// get RMS 
      double  rms       () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral () const { return 1 ; }
      /// get the integral 
      double integral 
      ( const double low  , 
        const double high ) const ;      
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	      const double xhigh ) const ;
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
      /// scale/sigma 
      double m_scale ; // scale/sigma
      /// k-value 
      double m_k     ; // k-value 
      /// kappa-value 
      double m_kappa ; // kappa-value 
      // ======================================================================
    private:
      // ======================================================================
      /// get Z_k constant 
      double m_Zk ; // get Z_k constant 
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
      Hyperbolic 
      ( const double mu     = 0 ,   // related to location 
        const double sigma  = 1 ,   // related to withs  
        const double zeta   = 1 ,   // related to shape 
        const double kappa  = 0 ) ; // related to asymmetry 
      // ======================================================================
    public :
      // ======================================================================
      /// evaluate  pdf  for Hyperbolic distribution
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  for Hyperbolic distribution
      inline double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// location parameter 
      inline double mu       () const { return m_mu    ; } // location parameter 
      /// sigma-parameter 
      inline double sigma    () const { return m_sigma ; } // sigma-parameter 
      /// squared sigma parameter 
      inline double sigma2   () const { return m_sigma * m_sigma  ; } 
      /// zeta-parameter 
      inline double zeta     () const { return m_zeta  ; } // zeta-parameter 
      /// squared zeta parameters
      inline double zeta2    () const { return m_zeta  * m_zeta   ; }
      /// asymmetry parameter kappa 
      inline double kappa    () const { return m_kappa ; } // asymmetry parameter kappa 
      /// squared asymmetry parameter
      inline double kappa2   () const { return m_kappa * m_kappa ; }
      // ======================================================================
    public : // original parameters 
      // ======================================================================
      /// alpha parameter 
      inline double alpha    () const { return std::hypot ( beta () , gamma () ) ; }
      /// squared alpha 
      inline double alpha2   () const { return beta2 () + gamma2 ()       ; }      
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
      //  =====================================================================
    public : // features 
      // ======================================================================
      inline double location    () const { return mu ()   ; }
      /// get the actual mode of the distribution
      double mode        () const ; 
      /// get mean value 
      double mean        () const ;
      /// get variance 
      double variance    () const ;
      /// get dispersion 
      inline double dispersion  () const { return variance () ; }
      /// get RMS 
      inline double rms         () const { return std::sqrt ( variance () ) ; }
      /// get RMS 
      inline double RMS         () const { return rms ()      ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool        setMu       ( const double value ) ;
      bool        setSigma    ( const double value ) ;
      bool        setZeta     ( const double value ) ;
      bool        setKappa    ( const double value ) ;
      inline bool setLocation ( const double value ) { return setMu ( value ) ; }
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
      double integral
      ( const double low  , 
        const double high ) const ;      
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
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
      double        mean        () const ;
      /// get variance 
      double        variance    () const ;
      /// get dispersion 
      inline double dispersion  () const { return variance () ; }
      /// get RMS 
      inline double rms         () const { return std::sqrt ( variance () ) ; }
      /// get RMS 
      inline double RMS         () const { return rms ()      ; }
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
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over the function 
       * - Gaussian is centered at mean-value with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	      const double xhigh ) const ;
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
        const double alphaL = 2 ,    // left tail parameter 
        const double alphaR = 2 ) ;  // core asymmetry
      // ======================================================================      
    public :
      // ======================================================================
      /// evaluate  pdf  
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  
      inline double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get location parameter 
      inline double mu       () const { return m_core.mu     () ; }
      /// get width parameter 
      inline double sigma    () const { return m_core.sigma  () ; }
      /// get left tail parameters
      inline double alphaL   () const { return m_left .alpha () ; }
      /// get right tail parameter
      inline double alphaR   () const { return m_right.alpha () ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set location parameter 
      bool        setMu       ( const double value ) { return m_core.setMu    ( value ) ; }
      /// set width parameter 
      bool        setSigma    ( const double value ) { return m_core.setSigma ( value ) ; }
      /// set left tail 
      bool        setAlphaL   ( const double value ) { return m_left.setAlpha  ( value ) ; }
      /// set right tail 
      bool        setAlphaR   ( const double value ) { return m_right.setAlpha ( value ) ; }
      /// set location parameter 
      inline bool setLocation ( const double value ) { return setMu ( value ) ; }
      /// set both K simutaneously :
      bool setAlpha 
      ( const double aL , 
        const double aR )
      {
        const bool updated1 = setAlphaL ( aL ) ;
        const bool updated2 = setAlphaR ( aR ) ; 
        return updated1 && updated2 ;
      }
      // ======================================================================
    public: // derived  
      // ======================================================================
      /// get location parameter 
      inline double location () const { return mu () ; }
      /// get mode  
      inline double mode     () const { return m_core.mode () ; }
      /// left transition point 
      inline double xL () const  { return m_core.mu () - m_core.sigma () * m_left .alpha () ; }
      /// right transition point 
      inline double xR () const  { return m_core.mu () + m_core.sigma () * m_right.alpha () ; }
      // ======================================================================
    public :
      // ======================================================================
      const Ostap::Math::Gauss&        gauss      () const { return m_core  ; }
      const Ostap::Math::Gauss&        core       () const { return m_core  ; }
      const Ostap::Math::LeftExpTail&  tail_left  () const { return m_left  ; }
      const Ostap::Math::RightExpTail& tail_right () const { return m_right ; }
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
      /// core function
      Ostap::Math::Gauss        m_core  {} ; 
      /// left tail 
      Ostap::Math::LeftExpTail  m_left  {} ;
      /// Ritgh tail 
      Ostap::Math::RightExpTail m_right {} ;   
      // ======================================================================
    };
    // ========================================================================
  /** @class ADas
   * Asymmetric version of Das function 
   * - bifurcated Gaussian as core 
   * - left exponential tail 
   * - right exponential tail
   * @see Ostap::Math::BifurcatedGauss
   * @see Ostap::Math::Gauss 
   * @see Ostap::Math::LeftExpTail
   * @see Ostap::Math::RightExpTail 
   */
   class ADas
    {
      // ======================================================================
    public:
      // ====================================================================== 
      /** constructor with full parameters 
       *  @param mu peak location 
       *  @param sigmaL left sigma for Bifurcated Gaussian Core 
       *  @param sigmaR right sigma for bigurcatd Gaussian core
       *  @param alphaL left tail parameter 
       *  @param alphaR right tail parameter 
       */
      ADas 
      ( const double mu     = 0 ,    // location parameter 
        const double sigmaL = 1 ,    // left sigma
        const double sigmaR = 1 ,    // right sigma  
        const double alphaL = 2 ,    // left tail parameter 
        const double alphaR = 2 ) ;  // core asymmetry
      // ======================================================================      
    public :
      // ======================================================================
      /// evaluate  pdf  
      double pdf ( const  double x ) const ;
      /// evaluate  pdf  
      inline double operator() ( const double x ) const { return pdf ( x ) ; }      
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get location parameter 
      inline double mu       () const { return m_core.mu     () ; }
      /// get width parameter 
      inline double sigmaL   () const { return m_core.sigmaL () ; }
      /// get width parameter
      inline double sigmaR   () const { return m_core.sigmaR () ; }
      /// get left tail 
      inline double alphaL   () const { return m_left .alpha () ; }
      /// get right tail 
      inline double alphaR   () const { return m_right.alpha () ; }
      /// sigma-asymmetry
      inline double kappa    () const { return m_core.kappa  () ; }
      /// sigma asymmetry kappa = tanh ( psi )
      inline double psi      () const { return m_core.psi    () ; } 
      /// averagae sigma 
      inline double sigma    () const { return m_core.sigma  () ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set location parameter 
      inline bool  setMu       ( const double value ) { return m_core.setMu      ( value ) ; }
      /// set width parameter 
      inline bool  setSigmaL   ( const double value ) { return m_core.setSigmaL  ( value ) ; }
       /// set width parameter 
      inline  bool setSigmaR   ( const double value ) { return m_core.setSigmaR  ( value ) ; }
      /// set left tail 
      inline bool setAlphaL    ( const double value ) { return m_left .setAlpha  ( value ) ; }
      /// set right tail 
      inline bool setAlphaR    ( const double value ) { return m_right.setAlpha  ( value ) ; }
      /// set location parameter 
      inline bool setLocation  ( const double value ) { return setMu ( value ) ; }
      /// set both K simutaneously :
      bool setAlpha 
      ( const double aL , 
        const double aR )
      {
        const bool updated1 = setAlphaL ( aL ) ;
        const bool updated2 = setAlphaR ( aR ) ; 
        return updated1 && updated2 ;
      }
      // set both sigmas simultaneously
      bool setSigma 
      ( const double valueL , 
        const double valueR ) 
      { return m_core.setSigma ( valueL , valueR ); }
      // ======================================================================
    public: // derived  
      // ======================================================================
      /// get location parameter 
      inline double location () const { return mu () ; }
      /// get mode  
      inline double mode     () const { return m_core.mode () ; }
      /// left transition point 
      inline double xL () const  { return m_core.mu () - m_core.sigmaL () * m_left .alpha () ; }
      /// right transition point 
      inline double xR () const  { return m_core.mu () + m_core.sigmaR () * m_right.alpha () ; }
      // ======================================================================
    public :
      // ======================================================================
      const Ostap::Math::BifurcatedGauss& core       () const { return m_core  ; }
      const Ostap::Math::LeftExpTail&     tail_left  () const { return m_left  ; }
      const Ostap::Math::RightExpTail&    tail_right () const { return m_right ; }
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
      /// core function
      Ostap::Math::BifurcatedGauss m_core  {} ; 
      /// left tail 
      Ostap::Math::LeftExpTail     m_left  {} ;
      /// Ritgh tail 
      Ostap::Math::RightExpTail    m_right {} ;   
      // ======================================================================
    };
    // ========================================================================
    /** @class SkewGenT
     *  Skewed Generalised t-distribution
     *  @see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution
     *  Original function is parameterised in terms of parameters 
     *  - \f$ \mu \$ related to locartion 
     *  - \f$ \sigma \$ related to width/scale 
     *  - \f$ -1 < \lambda < 1 \f$ related to asymmetry/skewness  
     *  - \f$ 0<p, 0<q \f$ related to kutsosis
     *
     *  Mean value is defined if \f$ 1 < pq \f$ 
     *  RMS is defined for \f$ 2 < pq \f$
     * 
     *  In this view here we adopt minor reparameterisation in terms of 
     *  - \f$ 0 < r \f$, such as  \f$  r = \frac{1}{p} 
     *  - \f$ 0< \zeta \f$, such as \f$ pq = \zeta + 4 \f$
     *  - \f$ -\infty < \psi < +\infty \f$, such as \f$ \lambda  = \tanh \psi \f$   
     *
     *  Usage of \f$ \zeta\f$ ensures the existance of the  mean, RMS, skewness & kurtosis
     * 
     *  Special limiting cases:
     *  - \f$ q\rigtharrow +\infty (\zeta \rightarrow +\infty) \f$ 
     *     Generalized Error Distribution 
     *  - \f$ \lambda=0 (\psi = 0)  \f$ Generalized t-distribution 
     *  - \f$ p=2(r=\frac{1}{2}) \f$  Skewed t-distribution 
     *  - \f$ p=1(r=1), q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
     *     Skewed Laplace distribution 
     *  - \f$ \lambda=0, q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
     *     Generalized Error Distribution 
     *  - \f$ p=2(r=\frac{1}{2}), q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
     *     Skewed Normal distribution 
     *  - \f$ \sigma=1, \lambda=0,p=2(r=\frac{1}{2},  q=\frac{n+2}{2} (\alpha=n) \f$
     *     Student's t-distribution 
     *  - \f$ \lambda=0, p=1(r=1), q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
     *     Laplace distribution 
     *  - \f$ \lambda=0, p=2(r=\frac{1}{2}, q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
     *     Skewed Normal distribution 
     *
     *  @see Ostap::Math::SkewGenError 
     *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
     *  @date 2022-01-19
     */
    class SkewGenT
    {
    public:
      // ======================================================================    
      /** constructor with full parameters 
       *  @param mu    related to location 
       *  @param sigma related to RSM/scale/width 
       *  @param psi     related to asymmetry/skewness
       *  @param r      shape parameter 
       *  @param alpha  shape parameter    
       */
      SkewGenT  
      ( const double mu     = 0   ,    // location parameter 
        const double sigma  = 1   ,    // width parameter 
        const double psi    = 0   ,    // related to asymmetry/skewness 
        const double r      = 0.5 ,    // shape parameter 
        const double zeta   = 1   ) ;  // shape parameter 
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the pdf 
      double  pdf ( const double x ) const ;
      /// evaluate the pdf 
      inline double evaluate    ( const double x ) const { return pdf (  x ) ; }
      /// evaluate the pdf 
      inline double operator () ( const double x ) const { return pdf (  x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// location parameter 
      inline double mu    () const { return m_mu    ; }
      /// width/scale parameter 
      inline double sigma () const { return m_sigma ; }
      /// asymmetry/skewness  parameter      
      inline double psi   () const { return m_psi   ; }
      /// shape parameter 
      inline double r     () const { return m_r     ; }
      /// shape parameter 
      inline double zeta  () const { return m_zeta  ; }
      // ======================================================================
      // other parameters 
      // ======================================================================
      /// original lambda parameters 
      inline double lambda  () const { return m_lambda  ; }
      /// original lambda parameters 
      inline double lambda_ () const { return m_lambda  ; }
      /// original lambda parametter 
      inline double lambd   () const { return m_lambda  ; }
      /// original p-parameter 
      inline double p       () const { return 1.0 / m_r ; }
      /// original q-parameter 
      inline double q       () const { return ( m_zeta + 4 ) * m_r ; }
      // ======================================================================
    public : // helper parameters 
      // ======================================================================
      /** helper scale parameter 
       *  \f$ v^{\prime} = \frac{1}{ \sqrt{  (2\lambda^2+1) b_3 - 4 \lambda^2 b_2^2 } } \f$
       */
      double        v_scale () const ;
      /** helper bias parameter 
       *  \f$ m^{\prime} = 2 \sigma \lambda b_2 \f$
       */
      double        m_bias  () const ;
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu    ( const double value ) ;
      bool setSigma ( const double value ) ;
      bool setPsi   ( const double value ) ;
      bool setR     ( const double value ) ;
      bool setZeta  ( const double value ) ;
      // ======================================================================
    public: // some stat quantities 
      // ======================================================================
      /// mean 
      inline double mean       () const { return m_mu    ; }
      /// RMS 
      inline double rms        () const { return m_sigma ; }
      /// variance 
      inline double variance   () const { return m_sigma * m_sigma ; }
      /// dispersion 
      inline double dispersion () const { return variance () ; }
      /// skewness 
      double        skewness   () const ;      
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// integral 
      double integral () const ;
      /// integral from low to high 
      double integral
      ( const double low  , 
        const double high ) const ;
      // // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value  with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ; // get the unique tag 
      // ======================================================================
    private: // calculate helper math constants 
      /// =====================================================================
      /** calculate helper math constants
       *  \f[ \left( \begin{array}{l} 
       *    b_1 \\ b_2 \\ b3
       *   \end{array}\right) = 
       *   \left( \begin{array}{l}
       *    \frac{1}                                    { \Beta( \frac{1}{p} , q )} \\ 
       *   \frac{ \Beta( \frac{2}{p} , q -\frac{1}{p})} { \Beta( \frac{1}{p} , q )} \\
       *   \frac{ \Beta( \frac{3}{p} , q -\frac{3}{p})} { \Beta( \frac{1}{p} , q )} 
       *   \end{array}\right) \f]
       */
      void calc_b 
      ( double& b1 ,    // 1            / B ( 1/p, q ) 
        double& b2 ,    // B(2/p,q-1/p) / B ( 1/p, q) 
        double& b3 ) ;  // B(3/p,q-2/p) / B ( 1/p, q) 
      // ==============================
    private: // true parameters         
      // ======================================================================
      /// location 
      double m_mu     { 0   } ; // location parameter
      /// width/scale 
      double m_sigma  { 1   } ; // width/scale parameter
      /// asymmetry/skewness parameter 
      double m_psi    { 0   } ;
      /// shape parametyer 
      double m_r      { 0.5 } ;
      /// shape parametyer 
      double m_zeta   { 4   } ;
      // ======================================================================
    private: // helper parameters 
      // ======================================================================
      /// original lambda parameter 
      double m_lambda { -100 } ;
      /// =====================================================================
    private: // helper math constants 
      /// =====================================================================
      /// 1 / B ( 1/p, q)  
      double   m_b1   { -100 } ;
      /// B (2/p,q-1/p) / B (1/p,q)
      double   m_b2   { -100 } ;
      /// B (3/p,q-2/p) B (1/p,q)
      double   m_b3   { -100 } ;
      /// =====================================================================      
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================      
    } ;  
    // ========================================================================
    /** @class SkewGenError 
     *  Skewed gheneralised error distribution 
     *  @see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution#Skewed_generalized_error_distribution
     *
     *  The Special  case of Skewed Generaliaed T-distribution 
     *  @see Ostap::Math::SkewGenT 
     * 
     *  Original function is parameterised in terms of parameters 
     *  - \f$ \mu \$ related to location  
     *  - \f$ \sigma \$ related to width/scale 
     *  - \f$ -1 < \lambda < 1 \f$ related to asymmetry/skewness  
     *  - \f$ 0<p \f$ shape parameters 
     *
     *  \f[ f(x;\mu,\sigma,\lambda,p) = 
     *    \frac{p}{2v\sigma\Gamma(1/p)} \mathrm{e}^{ - \Delta^{p}},  
     *   \f]
     *  where 
     *   - \f$ v = \sqrt{ \frac{ \pi \Gamma(1/p)}{  \pi(1+3\lambda^2)\Gamma(3/p) 
     *            -16^{1/p} \lambda^2 \Gamma(1/2+1/p)^2\Gamma(1/p) }  }\f$,
     *   - \f$ \Delta = \frac{\left| \delta x \right|}{v\sigma ( 1+ \lambda \sign \delta x )} \f$
     *   - \f$ \delta x = x - \mu + m \f$
     *   - \f$ m =  2^{2/p} v \sigma \Gamma( 1/2+ 1/p)/\sqrt{\pi}\f$ 
     *
     *  Here we adopt sligth reparameterisation in terms of 
     *  - \f$ -\infty < \psi < +\infty \f$, such as \f$ \lambda  = \tanh \psi \f$   
     *  - r = 1/p 
     * 
     *  special cases: 
     *  - \f$ \psi=0 (\lambda=0), r=1/2)$ corresponds to Gaussian function 
     *  - \f$ \psi=0 (\lambda=0), e=1)$ corresponds to Laplace case 
     *
     *  @see Ostap::Math::SkewGenT 
     *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
     *  @date 2022-01-19
     */
    class SkewGenError
    {
    public:
      // ======================================================================    
      /** constructor with full parameters 
       *  @param mu    related to location 
       *  @param sigma related to RSM/scale/width 
       *  @param psi     related to asymmetry/skewness : lambda = tanh(psi)
       *  @param r      shape parameter 
       *  @param alpha  shape parameter    
       */
      SkewGenError  
      ( const double mu     = 0 ,   // location parameter 
        const double sigma  = 1 ,   // width parameter 
        const double psi    = 0 ,   // asymmetry/skewness parameter 
        const double r      = 2 ) ; // shape parameter 
      // ======================================================================
      public:
      // ======================================================================
      /// evaluate the pdf 
      double  pdf ( const double x ) const ;
      /// evaluate the pdf 
      inline double evaluate    ( const double x ) const { return pdf (  x ) ; }
      /// evaluate the pdf 
      inline double operator () ( const double x ) const { return pdf (  x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// locaton parameter 
      inline double mu    () const { return m_mu    ; }
      /// width/scale parameter 
      inline double sigma () const { return m_sigma ; }
      /// asymmetry/skewness  parameter      
      inline double psi   () const { return m_psi   ; }
      /// shape parameter 
      inline double r     () const { return m_r     ; }
      // ======================================================================
      // other/originam  parameters      
      // ======================================================================
      /// original lambda parametter 
      inline double lambda  () const { return m_lambda ; }
      /// original lambda parametter 
      inline double Lambda  () const { return m_lambda ; }
      /// original lambda parametter 
      inline double lambda_ () const { return m_lambda ; }
      /// original lambda parametter 
      inline double lambd   () const { return m_lambda ; }
      /// original p-parameter 
      inline double p       () const { return m_p      ; }
      // ======================================================================
    public:
      // ======================================================================
      /** helper scale parameter 
       *  \f$ v^{\prime} = \sqrt{ \frac{pi} } { \pi ( 1+3\lambda^2) b_1 - \lambda^2 b_2^2 }  \f$
       */
      double        v_scale () const ;
      // ======================================================================
      /** helper bias parameter 
       *  \f$ m^{\prime} = \lambda \sigma b_2 \f$
       */
      double        m_bias  () const ;
      // ======================================================================
    public: // setters 
      // ======================================================================
      bool setMu    ( const double value ) ;
      bool setSigma ( const double value ) ;
      bool setPsi   ( const double value ) ;
      bool setR     ( const double value ) ;
      // ======================================================================
    public: // some stat quantities 
      // ======================================================================
      /// mean 
      inline double mean       () const { return m_mu    ; }
      /// RMS 
      inline double rms        () const { return m_sigma ; }
      /// variance 
      inline double variance   () const { return m_sigma  * m_sigma ; }
      /// dispersion 
      inline double dispersion () const { return variance () ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// integral 
      double integral () const ;
      /// integral from low to high 
      double integral
      ( const double low  , 
        const double high ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value  with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================      
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ; // get the unique tag 
      // ======================================================================
    private: // calculate helper math constants 
      /// =====================================================================
      /** calculate helper math constants
       *  \f[ \left( \begin{array}{l} 
       *    b_0 \\ b_1 \\ b_2 
       *   \end{array}\right) = 
       *   \left( \begin{array}{l}
       *   \frac{1}{\Gamma(1/p) } \\ 
       *   \frac{\Gamma(3/p)    }{\Gamma^3(1/p)} \\ 
       *   2^{2/p}\frac{\Gamma(1/2+1/p)}{\Gamma(1/p)  } 
       *   \end{array}\right) \f]
       */
      void calc_b 
      ( double& b0 ,    // 1/Gamma(1/p) 
        double& b1 ,    // Gamma(3/p)/Gamma^3(1/p) 
        double& b2 ) ;  // 2^{2/p} Gamma(1/2+1/p)/Gamma(1/p)
      // ==============================
    private: // true parameters         
      // ======================================================================
      /// location 
      double m_mu     { 0   } ; // location parameter
      /// width/scale 
      double m_sigma  { 1   } ; // width/scale parameter
      /// asymmetry/skewness parameter 
      double m_psi    { 0   } ;
      /// shape parametyer 
      double m_r      { 0.5 } ;
      // ======================================================================
    private: // helper parameters 
      // ======================================================================
      /// original lambda parameter 
      double m_lambda { -100 } ;
      /// origfinal p-parameter
      double m_p      { 2    } ;   
      /// =====================================================================
    private: // helper math constants 
      /// =====================================================================
      ///  1/ Gamma(1/p) 
      double   m_b0   { -100 } ; // 1/Gamma(1/p)
      /// Gamma(3/p)/Gamma^3(1/p)
      double   m_b1   { -100 } ; // Gamma(3/p)/Gamma^3(1/p)
      /// 2^{2/p}Gamma(1/2+1/p)/Gamma(1/p)
      double   m_b2   { -100 } ; // 2^{2/p} Gamma(1/2+1/p)/Gamma(1/p)
      /// =====================================================================      
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Meixner
     *  Meixner distribution
     *  @see Grigoletto, M., & Provasi, C. (2008). 
     *       "Simulation and Estimation of the Meixner Distribution". 
     *       Communications in Statistics - Simulation and Computation, 38(1), 58–77. 
     *  @see https://doi.org/10.1080/03610910802395679
     *  @see https://reference.wolfram.com/language/ref/MeixnerDistribution.html
     *
     *  Original distribtion is parameterise with 
     *   - location parameter m 
     *   - scale parameter    a 
     *   - shape parameter b : \f$ - \pi < b < +\pi \f$
     *   - shape parameter f : \f$ 0 < d \f$
     *
     *  Here we use a slight reparameterisation:
     *   - \f$ b   = 2 \atan \psi \f$ 
     *   - \f$ a^2 = \sigma^2 \frac{ \cos b + 1} {d}  \f$ 
     * 
     * Asymptotic:
     *  - \f$  x \rightarrow +\infty\f$ : \f$ f  \sim \left| x \right|^\rho \mathrm{e}^{\sigma_-x}\f$
     *  - \f$  x \rightarrow -\infty\f$ : \f$ f  \sim \left| x \right|^\rho \mathrm{e}^{\sigma_+x}\f$
     *  where 
     *  - \f$  \sigma_+ = \frac{\pi + b }{a} \f$    
     *  - \f$  \sigma_+ = \frac{\pi + b }{a} \f$    
     */ 
    class Meixner
    {
    public :
      // =======================================================================
      Meixner
      ( const double mu    = 0 ,   // location
        const double sigma = 1 ,   // scale
        const double psi   = 0 ,   // b = 2 * atan ( psi )
        const double shape = 1 ) ; // shape 
      // =======================================================================
    public:      
      // =======================================================================
      /// evaluate Meixner function
      inline double operator() ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate Meixner function
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; } 
      /// evaluate Meixner function
      double        evaluate   ( const double x ) const ;
      // =======================================================================
    public: // primary getters 
      // =======================================================================
      inline double mu       () const { return m_mu    ; }
      inline double sigma    () const { return m_sigma ; }
      inline double psi      () const { return m_psi   ; }
      inline double shape    () const { return m_shape ; }
      // =======================================================================
    public: // derived getters 
      // =======================================================================
      inline double a        () const { return m_a      ; }      
      inline double b        () const { return m_b      ; }      
      inline double d        () const { return shape () ; }      
      inline double location () const { return mu    () ; } 
      /// kappa = b/pi:  \f$ -1 < \kapppa < 1 \f$   
      double kappa           () const ;
      // =======================================================================      
    public: // setters 
      // =======================================================================
      /// set mu 
      bool setMu    ( const double value ) ;
      /// set sigma
      bool setSigma ( const double value ) ;
      /// set D 
      bool setShape ( const double value ) ;
      /// set psi 
      bool setPsi   ( const double value ) ;
      // =======================================================================
      /// set location 
      inline bool setLocation ( const double value ) { return setMu    ( value ) ; }  
      /// set D 
      inline bool setD        ( const double value ) { return setShape ( value ) ; }  
      // =======================================================================
    public: 
      // =======================================================================
      /// mean value
      double mean     () const ;
      /// mode value
      double mode     () const ;
      /// variance/dispersion 
      inline double variance () const { return m_sigma * m_sigma ; }
      /// RMS 
      inline double rms      () const { return           m_sigma ; }
      /// get skewness 
      double skewness () const ;
      /// get (excess) kurtosis 
      double kurtosis () const ;      
      // =======================================================================
    public: // integrals 
      // ======================================================================
      /// integral 
      double integral () const ;
      /// integral from low to high 
      double integral
      ( const double low  , 
        const double high ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value  with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================      
    public: // asymptotoc 
      // ======================================================================
      /** Asymptotic :
       *  - \f$ x \rigtharrow +\infty\f$
       *   \f$ f \sim \left| x\right|^{\rho} \mathrm{e}^{-\sigma_+\left|x\right|}\f$
       *  - \f$ x \rigtharrow -\infty\f$
       *   \f$ f \sim \left| x\right| ^{\rho} \mathrm{e}^{-\sigma_-\left|x\right|}\f$
       *  - \f$ \rho     = 2d-1 \f$ 
       *  - \f$ \sigma_+ = \frac{\pi + b}{s}  \f$ 
       *  - \f$ \sigma_- = \frac{\pi - b}{s}  \f$ 
       */
      double rho         () const ;
      double sigma_plus  () const ; 
      double sigma_minus () const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ; // get the unique tag 
      // ======================================================================
    private:
      // ======================================================================
      /// location parameter
      double m_mu    {  0 } ; // location parameter 
      /// sigma parameter
      double m_sigma {  1 } ; // scale parameter
      /// asymmetry/skew parameter
      double m_psi   {  0 } ; // skew/asymmetry parameter
      /// shape parameter
      double m_shape {  1 } ; // shape parameter 
      // =======================================================================
    private :
      // =======================================================================
      /// keep a
      double m_a     {  1 } ; // a-value
      /// keep b
      double m_b     {  0 } ; // b-value
      /// normalization 
      double m_C     { -1 } ; // normalization 
      // =======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================            
    };    
    // ========================================================================
    /// some finite functions 
    // ========================================================================
    /** @class Hat 
     *  Finite smooth functon
     *  \f$ f(x;\mu\sigma) = \frac{C}{\sigma} 
     *   \mathrm{e}^{  - frac{1}{1-y^2} }\f$, where 
     *  \F$ y = \frac{m-\mu}{\sigma}\f$ 
     *  @see Ostap::Math::hat 
     */
    class Hat 
    {
    public:
      // ======================================================================
      /// constructor with location and scale parmaeters 
      Hat
      ( const double mu       = 0 , 
        const double varsigma = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double evaluate ( const double x ) const ;
      /// evaluate the function 
      inline double pdf         ( const double x ) const 
      { return evaluate ( x )  ; }
      /// evaluate the function 
      inline double operator()  ( const double x ) const 
      { return evaluate ( x )  ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get mu 
      inline double mu       () const { return m_mu       ; }
      /// get varsigma 
      inline double varsigma () const { return m_varsigma ; } 
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set mu
      bool setMu        ( const double value ) ;
      /// set varsigma 
      bool setVarsigma  ( const double value ) ;
      // ======================================================================
    public: // some properties 
      // ======================================================================
      /// get the mean     of the distribution 
      inline double mean     () const { return m_mu; }      
      /// get the mode     of the distribution 
      inline double mode     () const { return m_mu; }      
      /// get the median   of the distribution 
      inline double median   () const { return m_mu; }      
      /// get the variance of the distribution 
      double variance () const ;
      /// get the RMS  of the distribution 
      double rms      () const ;
      /// get the skewness 
      double skewness () const { return 0 ; }
      /// get the (excess) kurtosis 
      double kurtosis () const ;
      // ======================================================================
    public: // support 
      // ======================================================================
      /// x-min
      double xmin () const { return m_mu - m_varsigma ; }
      /// x-max
      double xmax () const { return m_mu + m_varsigma ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// integral 
      double integral () const ;
      /// integral from low to high 
      double integral
      ( const double low  , 
        const double high ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value  with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================      
    public: // derivative 
      // ======================================================================
      /// get the value of the derivative 
      double derivative ( const double x ) const ;
      // ======================================================================
    public: // unique tag 
      // ======================================================================
      /// unique tag 
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      /// evaluate the "standard" up function 
      double eval ( const double z )  const ;
      // ======================================================================
    private:
      // ======================================================================
      /// location parametyer 
      double m_mu        { 0 } ;
      /// scale paramewter 
      double m_varsigma  { 0 } ;
      // ====================================================================== 
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Up 
     *  Finite atomic function <code>up</code>, a finite solution  
     *  of the equation 
     *  \f[ f^{\prime(x) = 2 \left( f( 2x+1) - f(2x-1)\right) }\f] with
     *  \f$ f(0) = 1 \f$ 
     *  @see Ostap::Math::up_F 
     */
    class Up 
    {
    public:
      // ======================================================================
      /// constructor with location and scale parmaeters 
      Up 
      ( const double mu       = 0 , 
        const double varsigma = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double evaluate ( const double x ) const ;
      /// evaluate the function 
      inline double pdf         ( const double x ) const 
      { return evaluate ( x )  ; }
      /// evaluate the function 
      inline double operator()  ( const double x ) const 
      { return evaluate ( x )  ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get mu 
      double mu       () const { return m_mu       ; }
      /// get varsigma 
      double varsigma () const { return m_varsigma ; } 
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set mu
      bool setMu        ( const double value ) ;
      /// set varsigma 
      bool setVarsigma  ( const double value ) ;
      // ======================================================================
    public: // some properteis 
      // ======================================================================
      /// get the mean     of the distribution 
      double mean     () const { return m_mu; }      
      /// get the mode     of the distribution 
      double mode     () const { return m_mu; }      
      /// get the median   of the distribution 
      double median   () const { return m_mu; }      
      /// get the variance of the distribution 
      double variance () const ;
      /// get the RMS  of the distribution 
      double rms      () const ;
      /// get the skewness 
      double skewness () const { return 0 ; }
      /// get the (excess) kurtosis 
      double kurtosis () const ;
      // ======================================================================
    public: // support 
      // ======================================================================
      /// x-min
      double xmin () const { return m_mu - m_varsigma ; }
      /// x-max
      double xmax () const { return m_mu + m_varsigma ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// integral 
      double integral () const ;
      /// integral from 
      double integral
      ( const double low  , 
        const double high ) const ;
      // ======================================================================
    public: //
      // ======================================================================
      /** quantify the effect of the tails, difference from Gaussian
       *  \f[ Q = 1 = frac{I_{CB} - I_G}{I_{CB}} \f]
       * where 
       * - \f$ I_{CB} \f$ is integral over Gaussian function 
       * - \f$ I_{G}  \f$ is integral over Crystal Ball function 
       * - Gaussian is centered at mean-value  with sigma = RMS 
       */
      double non_gaussian 
      ( const double xlow  ,
	const double xhigh ) const ;
      // ======================================================================      
    public: // derivative 
      // ======================================================================
      /// get the value of the derivative 
      double derivative ( const double x ) const ;
      // ======================================================================
    public: // unique tag 
      // ======================================================================
      /// unique tag 
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      /// evaluate the "standard" up function 
      double eval ( const double z )  const ;
      // ======================================================================
    private:
      // ======================================================================
      /// location parametyer 
      double m_mu        { 0 } ;
      /// scale paramewter 
      double m_varsigma  { 0 } ;
      // ====================================================================== 
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class FupM 
     *  Finite stomic functon <code>fup_N</code>,  a finite solution
     *  of the equation 
     *  \f[ f^{\prime(x) = 2 \left( f( 2x+1) - f(2x-1)\right) }\f] with
     *  \f$ f(0) = 1 \f$ 
     *  @see Ostap::Math::fupN_F 
     */
    class FupN
    {
    public:
      // ======================================================================
      /// constructor with location and scale parmaeters 
      FupN 
      ( const unsigned short N  = 1 , 
        const double   mu       = 0 , 
        const double   varsigma = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double evaluate ( const double x ) const ;
      /// evaluate the function 
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate the function 
      inline double operator()  ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get N 
      inline unsigned short N        () const { return m_N        ; }
      /// get mu 
      inline double         mu       () const { return m_mu       ; }
      /// get varsigma 
      inline double         varsigma () const { return m_varsigma ; } 
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set mu
      bool setMu        ( const double value ) ;
      /// set varsigma 
      bool setVarsigma  ( const double value ) ;
      // ======================================================================
    public: // some properteis 
      // ======================================================================
      /// get the mean     of the distribution 
      inline double mean     () const { return m_mu ; }      
      /// get the mode     of the distribution 
      inline double mode     () const { return m_mu ; }      
      /// get the median   of the distribution 
      inline double median   () const { return m_mu ; }      
      /// get the skewness 
      inline double skewness () const { return 0    ; }
      // ======================================================================
    public: // support 
      // ======================================================================
      /// x-min
      double xmin () const { return m_mu - 0.5 * ( m_N + 2 ) * m_varsigma ; }
      /// x-max
      double xmax () const { return m_mu + 0.5 * ( m_N + 2 ) * m_varsigma ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// integral 
      double integral () const ;
      /// integral from low to high  
      double integral
      ( const double low  , 
        const double high ) const ;
      // ======================================================================
    public: // unique tag 
      // ======================================================================
      /// unique tag 
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      /// evaluate the "standard" fupN function 
      double eval ( const double z )  const ;
      // ======================================================================
    private:
      // ======================================================================
      /// parameter N 
      unsigned short m_N         { 1 } ;
      /// location parametyer 
      double         m_mu        { 0 } ;
      /// scale paramewter 
      double         m_varsigma  { 1 } ;
      // ====================================================================== 
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ; // integration workspace
      // ======================================================================      
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PEAKS_H
// ============================================================================
