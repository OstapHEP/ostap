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
 *
 *  set of useful peak-like models
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
      double pdf        ( const double x ) const { return evaluate ( x ) ; }
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double peak    () const { return m_peak    ; }
      double m0      () const { return peak()    ; }
      double sigmaL  () const { return m_sigmaL  ; }
      double sigmaR  () const { return m_sigmaR  ; }
      // ======================================================================
      double sigma   () const { return 0.5  * ( m_sigmaL + m_sigmaR )            ; }
      double asym    () const { return 0.5  * ( m_sigmaL - m_sigmaR ) / sigma () ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setPeak    ( const double value ) ;
      bool setM0      ( const double value ) { return setPeak ( value ) ; }
      bool setMass    ( const double value ) { return setPeak ( value ) ; }
      bool setSigmaL  ( const double value ) ;
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
       *  @param peak     the peak posiion
       *  @param sigmaL   the sigma for first component
       *  @param fraction the fraction of first component 
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
      double peak    () const { return m_peak    ; }
      double m0      () const { return peak()    ; }
      double sigma   () const { return m_sigma   ; }
      double sigma1  () const { return m_sigma   ; }
      double sigma2  () const { return m_sigma * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setPeak     ( const double value ) ;
      bool setM0       ( const double value ) { return setPeak ( value ) ; }
      bool setMass     ( const double value ) { return setPeak ( value ) ; }
      bool setSigma    ( const double value ) ;
      bool setScale    ( const double value ) ;
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
      double mu          () const { return m_mu       ; }
      double peak        () const { return   mu    () ; }
      double location    () const { return   mu    () ; }
      double alpha       () const { return m_alpha    ; }
      double scale       () const { return   alpha () ; }
      double beta        () const { return m_beta     ; }
      double shape       () const { return   beta  () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool  setMu        ( const double value ) ;
      bool  setAlpha     ( const double value ) ;
      bool  setBeta      ( const double value ) ;
      //
      bool  setPeak      ( const double value ) { return setMu    ( value ) ; }
      bool  setLocation  ( const double value ) { return setMu    ( value ) ; }
      bool  setScale     ( const double value ) { return setAlpha ( value ) ; }
      bool  setShape     ( const double value ) { return setBeta  ( value ) ; }
      // ======================================================================
    public: // derived getters
      // ======================================================================
      double mean        () const { return   mu    () ; }
      double median      () const { return   mu    () ; }
      double mode        () const { return   mu    () ; }
      //
      double variance    () const ;
      double dispersion  () const { return variance () ; }
      double sigma2      () const { return variance () ; }
      double sigma       () const ;
      //
      double skewness    () const { return 0 ; }
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
      double cdf      ( const double x ) const ;
      /// get the integral
      double integral () const { return 1 ; }
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_mu     ;  // location
      double m_alpha  ;  // scale
      double m_beta   ;  // shape
      double m_gbeta1 ;  // helper parameter
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
    private:
      // ======================================================================
      double m_xi     ;  // location
      double m_omega  ;  // scale
      double m_alpha  ;  // shape
      // =======================================================================
    } ;
    // ========================================================================
    /** @class Bukin
     *  ``Bukin-function'', aka "Modified Novosibirsk function"
     *  for description of asymmetric peaks with the exponential tails
     *
     *  @see http://arxiv.org/abs/1107.5751
     *  @see http://dx.doi.org/10.1007/JHEP06(2012)141
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
      /// helper constants
      double m_AL       ;  // exp(-0.5*alpha_L^2)
      double m_AR       ;  // exp(-0.5*alpha_R^2)
      double m_B        ;  // integral over the gaussian part
      double m_TL       ;  // integral over the left  power-law tail
      double m_TR       ;  // integral over the right power-law tail
      // ======================================================================
    } ;

    // ========================================================================
    /** @class Apolonios
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
    class  Apolonios
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
      Apolonios
      ( const double m0    = 0 ,
        const double sigma = 1 ,
        const double alpha = 2 ,
        const double n     = 1 ,
        const double b     = 1 ) ;
      /// destructor
      ~Apolonios () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Apolonios's function
      double pdf        ( const double x ) const ;
      /// evaluate Apolonios's function
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
    /** @class Apolonios2
     *  "Bifurcated Apolonios"
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
    class  Apolonios2 
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0      m0        parameter
       *  @param sigmaL  sigmaL    parameter
       *  @param alphaR  alphaR    parameter
       *  @param beta    beta      parameter
       */
      Apolonios2
        ( const double m0      = 0   ,
          const double sigmaL  = 1   ,
          const double alphaR  = 1   ,
          const double beta    = 100 ) ;  // large beta correponds to gaussian
      /// destructor
      ~Apolonios2 () ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate Apolonios2's function
      double pdf        ( const double x ) const ;
      /// evaluate Apolonios2's function
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
      StudentT ( const double mass  = 0 ,
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
     *  recovere by \f$\delta\rightarrow0\f$ for
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
      /// get mode
      double mode     () const { return mean() ; }
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
    private:
      // ======================================================================
      /// parameteter "mu", mean, mode
      double m_mean  ; // parameter mu,mean,mode
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
      double mean   () const { return m_mean   ; }
      /// get parameters "sigma"
      double sigma  () const { return m_sigma  ; }
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
    private:
      // ======================================================================
      /// parameteter "mu", mean, mode
      double m_mean  ; // parameter mu,mean,mode
      /// parameter   "sigma"
      double m_sigma ; // parameter sigma
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  Slash 
     *  ``Slash''-distribution -  symmetric peak with veyr heavy tail
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
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PEAKS_H
// ============================================================================
