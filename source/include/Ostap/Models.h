// $Id$
// ============================================================================
#ifndef OSTAP_MODELS_H
#define OSTAP_MODELS_H 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <functional>
#include <vector>
#include <complex>
// ============================================================================
// OStap
// ============================================================================
#include "Ostap/NSphere.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Polynomials.h"
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
    /** @class BifurcatedGauss
     *  simple representation of bifurcated gaussian function
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2011-04-19
     */
    class  BifurcatedGauss : public std::unary_function<double,double>
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
      double pdf        ( const double x ) const ;
      double operator() ( const double x ) const { return pdf ( x )  ; }
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
      double m_peak   ;       //                              the peak position
      /// sigma left
      double m_sigmaL ;       // sigma-left
      /// sigma right
      double m_sigmaR ;       // sigma-right
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenGaussV1
     *  Simple class that implements the generalized normal distribution v1
     *  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-08-25
     */
    class  GenGaussV1: public std::unary_function<double,double>
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
    class  GenGaussV2: public std::unary_function<double,double>
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
    // // ========================================================================
    // /** @class SkewGauss
    //  *  Simple class that implements the skew normal distribution
    //  *  @see http://en.wikipedia.org/wiki/Skew_normal_distribution
    //  *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    //  *  @date 2013-08-25
    //  */
    // class  SkewGauss: public std::unary_function<double,double>
    // {
    // public:
    //   // ======================================================================
    //   /** constructor from all agruments
    //    *  @param xi     location/peak posiiton
    //    *  @param omega  "scale" parameter
    //    *  @param alpha  "shape" parameter
    //    */
    //   SkewGauss
    //   ( const double xi    = 0 ,
    //     const double omega = 1 ,
    //     const double alpha = 0 ) ; // alpha=0 correponds to gaussian
    //   /// desctructor
    //   ~SkewGauss () ;
    //   // ======================================================================
    // public: // primary getters
    //   // ======================================================================
    //   double xi          () const { return m_xi       ; }
    //   double peak        () const { return   xi    () ; }
    //   double location    () const { return   xi    () ; }
    //   double omega       () const { return m_omega    ; }
    //   double scale       () const { return   omega () ; }
    //   double alpha       () const { return m_alpha    ; }
    //   double shape       () const { return   alpha () ; }
    //   // ======================================================================
    // public: // setters
    //   // ======================================================================
    //   bool  setXi        ( const double value ) ;
    //   bool  setOmega     ( const double value ) ;
    //   bool  setAlpha     ( const double value ) ;
    //   //
    //   bool  setPeak      ( const double value ) { return setXi    ( value ) ; }
    //   bool  setLocation  ( const double value ) { return setXi    ( value ) ; }
    //   bool  setScale     ( const double value ) { return setOmega ( value ) ; }
    //   bool  setShape     ( const double value ) { return setAlpha ( value ) ; }
    //   // ======================================================================
    // public: // derived getters
    //   // ======================================================================
    //   double mean        () const ;
    //   double mode        () const ;
    //   double variance    () const ;
    //   double dispersion  () const { return variance () ; }
    //   double sigma2      () const { return variance () ; }
    //   double sigma       () const ;
    //   double skewness    () const ;
    //   double kurtosis    () const ;
    //   // ======================================================================
    // public :
    //   // ======================================================================
    //   /// get pdf
    //   double operator() ( const double x ) const { return pdf ( x ) ; }
    //   double pdf        ( const double x ) const ;
    //   // ======================================================================
    // public:  // integrals
    //   // ======================================================================
    //   double cdf        ( const double x ) const ;
    //   /// get the integral
    //   double integral   () const { return 1 ; }
    //   /// get the integral between low and high limits
    //   double integral   ( const double low  ,
    //                       const double high ) const ;
    //   // ======================================================================
    // private:
    //   // ======================================================================
    //   double m_xi     ;  // location
    //   double m_omega  ;  // scale
    //   double m_alpha  ;  // shape
    //   // =======================================================================
    // } ;
    // ========================================================================
    /** @class WorkSpace
     *  helper utility to keep the integration workspace fro GSL integration
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2011-12-03
     */
    class WorkSpace
    {
    public:
      // ======================================================================
      /// constructor
      WorkSpace () ;
      /// (fictive) copy constructor
      WorkSpace ( const WorkSpace& right );
      /// destructor
      ~WorkSpace () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integration workspace
      void* workspace () const ;               // get the integrtaion workspace
      // ======================================================================
    public:
      // ======================================================================
      /// (fictive) assignement operator
      WorkSpace& operator= ( const WorkSpace& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual GSL-workspace
      mutable void*  m_workspace ;  /// the actual GSL-workspace
      // ======================================================================
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
    class  Bukin : public std::unary_function<double,double>
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
    class  Novosibirsk : public std::unary_function<double,double>
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
     * \f[ C = \frac{n+1}{\left|\alpha\right|\times \frac{1}{n} \times \mathrm{e}^{-\frac{\alpha^2}{2}}  \f]
     * \f[ B = \sqrt{\frac{\pi}{2}}\left(1+\mathrm{erf}\left(-\frac{\alpha}{\sqrt{2}}\right)\right) \f]
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  CrystalBall : public std::unary_function<double,double>
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
     *  @thank Matthew Needham
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
    class  Needham : public std::unary_function<double,double>
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
    class  CrystalBallRightSide : public std::unary_function<double,double>
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
      : public std::unary_function<double,double>
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
    class  Apolonios : public std::unary_function<double,double>
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
     *     \frac{x-\mu}{\sigma_l} & \text{for} & x \le \mu \\
     *     \frac{x-\mu}{\sigma_r} & \text{for} & x \gt \mu \\
     *     \end{array}
     *     \right.\f]
     *
     *  Large betas corresponds to gaussian
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date  2013-12-01
     */
    class  Apolonios2 : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0      m0        parameter
       *  @param sigmaL  sigmaL    parameter
       *  @param sigmaR  sigmaR    parameter
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
    /** @class GramCharlierA4
     *  Gram-Charlier type A approximation
     *  http://en.wikipedia.org/wiki/Edgeworth_series
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-06-13
     */
    class  GramCharlierA
      : public std::unary_function<double,double>
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
    /** @class PhaseSpace2
     *  simple function to represent two-body phase space
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpace2
      : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from two masses
      PhaseSpace2 ( const double m1 = 0 ,
                    const double m2 = 1 ) ;
      /// deststructor
      ~PhaseSpace2 () ;                                         // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate 2-body phase space
      double operator () ( const double x ) const ;
      /// integral
      double integral    ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the momentum at center of mass
      double                q_  ( const double x ) const ;
      /// ditto but as complex
      std::complex<double>  q1_ ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m1      () const { return m_m1 ; }
      double m2      () const { return m_m2 ; }
      double lowEdge () const { return m1() + m2() ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the first mass
      double m_m1 ; // the first  mass
      /// the second mass
      double m_m2 ; // the second mass
      // ======================================================================
    public :
      // ======================================================================
      /** calculate the triangle function
       *  \f[ \lambda ( a , b, c ) = a^2 + b^2 + c^2 - 2ab - 2bc - 2ca \f]
       *  @param a parameter a
       *  @param b parameter b
       *  @param c parameter b
       *  @return the value of triangle function
       */
      static double triangle
      ( const double a ,
        const double b ,
        const double c ) ;
      // ======================================================================
      /** calculate the particle momentum in rest frame
       *  \f[ q = \frac{1}{2}\frac{ \lambda^{\frac{1}{2}}
       *        \left( m^2 , m_1^2, m_2^2 \right) }{ m }\f],
       *  @param m the mass
       *  @param m1 the mass of the first particle
       *  @param m2 the mass of the second particle
       *  @return the momentum in rest frame (physical values only)
       */
      static double  q
      ( const double m  ,
        const double m1 ,
        const double m2 ) ;
      // ======================================================================
      /** calculate the particle momentum in rest frame
       *  @param m the mass
       *  @param m1 the mass of the first particle
       *  @param m2 the mass of the second particle
       *  @return the momentum in rest frame  (imaginary for non-physical branch)
       */
      static std::complex<double> q1
      ( const double m  ,
        const double m1 ,
        const double m2 ) ;
      // ======================================================================
      /** calculate the phase space for   m -> m1 + m2
       *  \f[ \Phi = \frac{1}{8\pi} \left( \frac{ \lambda^{\frac{1}{2}}
       *       \left( m^2 , m_1^2, m_2^2 \right) }{ m^2 }\right)^{2L+1}\f],
       *  where \f$\lambda\f$ is a triangle function
       *  @param m the mass
       *  @param m1 the mass of the first particle
       *  @param m2 the mass of the second particle
       *  @return two-body phase space
       */
      static double phasespace
      ( const double         m      ,
        const double         m1     ,
        const double         m2     ,
        const unsigned short L  = 0 ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpace3
     *  simple function to represent three-body phase space
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpace3
      : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from three masses
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       *  @param l1 the angular momentum between 1st and 2nd particle
       *  @param l2 the angular momentum between the pair and 3rd particle
       */
      PhaseSpace3 ( const double         m1 = 0 ,
                    const double         m2 = 1 ,
                    const double         m3 = 2 ,
                    const unsigned short l1 = 0 ,
                    const unsigned short l2 = 0 ) ;
      /// deststructor
      ~PhaseSpace3 () ;                                         // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate 3-body phase space
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double lowEdge () const { return m_m1 + m_m2 + m_m3 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// helper phase space ("23L")
      double ps2_aux ( const double m12 ) const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the mass of the first particle
      double         m_m1 ; // the mass of the first particle
      /// the mass of the second particle
      double         m_m2 ; // the mass of the second particle
      /// the mass of the third particle
      double         m_m3 ; // the mass of the third  particle
      /// the orbital momentum of the first pair
      unsigned short m_l1 ; // the orbital momentum of the first pair
      /// the orbital momentum between the pair and the third particle
      unsigned short m_l2 ; // the orbital momentum between the pair and the third particle
      // ======================================================================
    private:
      // ======================================================================
      /// the temporary mass
      mutable double m_tmp ; /// the temporary mass
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace  ;    // integration workspace
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace2 ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceLeft
     *  simple function to represent N-body phase space near left-threshold
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpaceLeft
      : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from threshold and number of particles
      PhaseSpaceLeft ( const double         threshold = 0 ,
                       const unsigned short num       = 2 ) ;
      /// constructor from list of masses
      PhaseSpaceLeft ( const std::vector<double>& masses ) ;
      /// deststructor
      ~PhaseSpaceLeft() ;                                       // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N-body phase space near left threhsold
      double operator () ( const double x    ) const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      bool setThreshold ( const double x ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// the threshold
      double         m_threshold ; // the threshold
      /// number of particles
      unsigned short m_num       ; // number of particles
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceRight
     *  simple function to represent N/L-body phase space near right-threshold
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpaceRight
      : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from threshold and number of particles
      PhaseSpaceRight ( const double         threshold = 10 ,
                        const unsigned short l         = 2  ,
                        const unsigned short n         = 3  ) ;
      /// deststructor
      ~PhaseSpaceRight () ;                                     // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N/L-body phase space near right  threhsold
      double operator () ( const double x ) const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      bool setThreshold ( const double x ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// the threshold
      double         m_threshold ; // the threshold
      /// number of particles
      unsigned short m_N         ; // number of particles
      /// number of particles
      unsigned short m_L         ; // number of particles
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceNL
     *  simple function to represent the approximation for
     *  the mass distribution of L-particles from N-body
     *  phase space decay
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpaceNL
      : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from thresholds and number of particles
       *  @param threshold_L the low-mass  threshold
       *  @param threshold_H the high-mass threshold
       *  @param L           how many particles we consider
       *  @param N           total number of particles ( N>L!)
       */
      PhaseSpaceNL ( const double         threshold_L =  0 ,
                     const double         threshold_H = 10 ,
                     const unsigned short l           =  2 ,
                     const unsigned short n           =  3 ) ;
      /// destructor
      ~PhaseSpaceNL () ;                                     // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N/L-body phase space
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double          lowEdge () const { return m_threshold1 ; }
      double         highEdge () const { return m_threshold2 ; }
      unsigned short       L  () const { return m_L ; }
      unsigned short       N  () const { return m_N ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set the thresholds
      bool setThresholds
      ( const double mn ,
        const double mx ) ;
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
      /// the threshold
      double         m_threshold1 ; // the threshold
      double         m_threshold2 ; // the threshold
      /// number of particles
      unsigned short m_N          ; // number of particles
      /// number of particles
      unsigned short m_L          ; // number of particles
      /// normalization
      double m_norm               ; // normalization
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
    class  PhaseSpacePol
      : public std::unary_function<double,double>
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
      /// destructor
      ~PhaseSpacePol () ;                                     // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N/L-body modulated phase space
      double operator () ( const double x ) const ;
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
    /** @class PhaseSpace23L
     *  simple function to represent the phase
     *   space of 2 particles from 3-body decays:
     *   \f$ f \propto q^{2\ell+1}p^{2L+1}\f$, where
     *     \f$\ell\f$ is the orbital momentum of the pair of particles,
     *    and \f$L\f$ is the orbital momentum between the pair and
     *    the third particle.
     *   E.g. taking \f$\ell=0, L=1\f$, one can get the S-wave contribution for
     *   \f$\pi^+\pi^-\f$-mass from \f$B^0\rightarrowJ/\psi\pi^+\pi^-\f$ decay.
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2012-04-01
     */
    class  PhaseSpace23L
      : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from four masses and angular momenta
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and
       *  the third particle
       *  @param l  the angular momentum between the first and the second particle
       */
      PhaseSpace23L ( const double         m1 = 0.5 ,
                      const double         m2 = 0.5 ,
                      const double         m3 = 3   ,
                      const double         m  = 5   ,
                      const unsigned short L  = 1   ,
                      const unsigned short l  = 0   ) ;
      /// deststructor
      ~PhaseSpace23L () ;                                     // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the phase space
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the phase space
      double ps23L ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m1        () const { return m_m1 ; }
      double m2        () const { return m_m2 ; }
      double m3        () const { return m_m3 ; }
      double m         () const { return m_m  ; }
      unsigned short l () const { return m_l  ; }
      unsigned short L () const { return m_L  ; }
      // ======================================================================
      double lowEdge   () const { return m1 () + m2 () ; }
      double highEdge  () const { return m  () - m3 () ; }
      /// get the momentum of 1st particle in rest frame of (1,2)
      double         q ( const double x ) const ;
      /// get the momentum of 3rd particle in rest frame of mother
      double         p ( const double x ) const ;
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
      /// the first mass
      double m_m1 ;            // the first mass
      /// the second mass
      double m_m2 ;            // the second mass
      /// the third  mass
      double m_m3 ;            // the third mass
      /// the mass of mother particle
      double m_m  ;            // the mass of mother particle
      /// the orbital momentum between the first and the second particle
      unsigned short m_l ; // the orbital momentum between the 1st and 2nd
      /// the orbital momentum between the pair an dthe third particle
      unsigned short m_L ; // the orbital momentum between the (12) and 3rd
      // ======================================================================
      /// helper normalization parameter
      double m_norm ; // helper normalization parameter
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /// base class for formfactors
    class FormFactor ;
    // ========================================================================
    namespace FormFactors
    {
      // ======================================================================
      /** @typedef rho_fun
       *  the \f$\rho(\omega)\f$ function from Jackson
       *  Arguments
       *    - the        mass
       *    - the pole   mass
       *    - the first  daughter mass
       *    - the second daughter mass
       */
      typedef double (*rho_fun) ( double , double , double , double ) ;
      // ======================================================================
      /** parameterization for \f$\rho(\omega)\f$-function from (A.1)
       *  J.D.Jackson,
       *  "Remarks on the Phenomenological Analysis of Resonances",
       *  In Nuovo Cimento, Vol. XXXIV, N.6
       */
      enum JacksonRho {
        Jackson_0  = 0 ,/// \f$\rho(\omega) = 1 \f$
        Jackson_A2 ,/// \f$ 1^- \rightarrow 0^- 0^- \f$ , l = 1
        Jackson_A3 ,/// \f$          1^- \rightarrow 0^- 1^- \f$ , l = 1
        Jackson_A4 ,/// \f$ \frac{3}{2}+ \rightarrow 0^- \frac{1}{2}^+ \f$ , l = 1
        Jackson_A5 ,/// \f$ \frac{3}{2}- \rightarrow 0^- \frac{1}{2}^+ \f$ , l = 2
        Jackson_A7 /// recommended for rho0 -> pi+ pi-
      } ;
      // ======================================================================
    }
    // ========================================================================
    /** @class BreitWigner
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *
     *  http://www.springerlink.com/content/q773737260425652/
     *
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  BreitWigner
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /// constructor from all parameters
      BreitWigner ( const double         m0     = 0.770 ,
                    const double         gam0   = 0.150 ,
                    const double         m1     = 0.139 ,
                    const double         m2     = 0.139 ,
                    const unsigned short L      = 0     ) ;
      /// constructor from all parameters
      BreitWigner ( const double         m0       ,
                    const double         gam0     ,
                    const double         m1       ,
                    const double         m2       ,
                    const unsigned short L        ,
                    const FormFactors::JacksonRho     r ) ;
      /// constructor from all parameters
      BreitWigner ( const double         m0       ,
                    const double         gam0     ,
                    const double         m1       ,
                    const double         m2       ,
                    const unsigned short L        ,
                    const FormFactor&    f ) ;
      /// copy constructor
      BreitWigner ( const BreitWigner&  bw ) ;
      ///
      // ======================================================================
      /// move constructor
      BreitWigner (       BreitWigner&& bw ) ;
      // ======================================================================
      /// destructor
      virtual ~BreitWigner () ;
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the Breit-Wigner shape
       *  \f$\frac{1}{\pi}\frac{\omega\Gamma(\omega)}
       *   { (\omega_0^2-\omega^2)^2-\omega_0^2\Gammma^2(\omega)-}\f$
       */
      virtual double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get breit-wigner amplitude
      std::complex<double> amplitude ( const double x ) const ;
      // ======================================================================
      double breit_wigner  ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double         m0     () const { return m_m0      ; }
      double         mass   () const { return   m0   () ; }
      double         peak   () const { return   m0   () ; }
      double         gam0   () const { return m_gam0    ; }
      double         gamma0 () const { return   gam0 () ; }
      double         gamma  () const { return   gam0 () ; }
      double         width  () const { return   gam0 () ; }
      // ======================================================================
    public:
      // ======================================================================
      double         m1     () const { return m_m1 ; }
      double         m2     () const { return m_m2 ; }
      unsigned short L      () const { return m_L  ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0     ( const double x ) ;
      bool setMass   ( const double x ) { return setM0     ( x ) ; }
      bool setPeak   ( const double x ) { return setM0     ( x ) ; }
      bool setGamma0 ( const double x ) ;
      bool setGamma  ( const double x ) { return setGamma0 ( x ) ; }
      bool setWidth  ( const double x ) { return setGamma0 ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the current width
      double gamma ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of formfactor at given m
      double            formfactor ( const double m ) const ;
      /// get the formfactor itself
      const Ostap::Math::FormFactor*
        formfactor () const { return m_formfactor ; }
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
      /// the mass
      double m_m0   ; // the mass
      /// the width
      double m_gam0 ; // the width
      // ======================================================================
    private:
      // ======================================================================
      /// the mass of the first  particle
      double            m_m1         ;
      /// the mass of the second particle
      double            m_m2         ;
      /// the orbital momentum
      unsigned int      m_L          ; // the orbital momentum
      /// the formfactor
      const Ostap::Math::FormFactor* m_formfactor ; // the formfactor
      // ======================================================================
    private:
      // ======================================================================
      /// assignement operator is disabled
      BreitWigner& operator=( const BreitWigner& ) ; // no assignement
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Rho0
     *  \f$ \rho^{0} \rightarrow \pi^+ \pi^- \f$
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *  @see Ostap::Math::BreitWigner::Jackson_A7
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Rho0 : public Ostap::Math::BreitWigner
    {
    public:
      // ======================================================================
      // constructor from all parameters
      Rho0  ( const double m0       = 770   ,     // MeV
              const double gam0     = 150   ,     // MeV
              const double pi_mass  = 139.6 ) ;   // MeV
      /// destructor
      virtual ~Rho0 () ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Kstar0
     *  \f$ K^{*0} \rightarrow K^+ \pi^- \f$
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *  @see Ostap::Math::BreitWigner::Jackson_A2
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2014-04-27
     */
    class  Kstar0 : public Ostap::Math::BreitWigner
    {
    public:
      // ======================================================================
      // constructor from all parameters
      Kstar0  ( const double m0       = 770   ,     // MeV
                const double gam0     = 150   ,     // MeV
                const double k_mass   = 493.7 ,     // MeV
                const double pi_mass  = 139.6 ) ;   // MeV
      /// destructor
      virtual ~Kstar0 () ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Phi0
     *  \f$ \phi \rightarrow K^+ K^- \f$
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *  @see Ostap::Math::BreitWigner::Jackson_A2
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2014-04-27
     */
    class  Phi0 : public Ostap::Math::BreitWigner
    {
    public:
      // ======================================================================
      // constructor from all parameters
      Phi0  ( const double m0       = 1019.5 ,     // MeV
              const double gam0     =    4.3 ,     // MeV
              const double k_mass   =  493.7 ) ;   // MeV
      /// destructor
      virtual ~Phi0 () ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Rho0FromEtaPrime
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Rho0FromEtaPrime : public Ostap::Math::Rho0
    {
    public:
      // ======================================================================
      /// constructor from all parameters
      Rho0FromEtaPrime  ( const double m0        = 770   ,    // MeV
                          const double gam0      = 150   ,    // MeV
                          const double pi_mass   = 139.6 ,    // MeV
                          const double eta_prime = 957.7 ) ;  // MeV
      /// constructor from all parameters
      Rho0FromEtaPrime  ( const Ostap::Math::Rho0& rho   ,
                          const double eta_prime = 957.7 ) ;  // MeV
      /// destructor
      virtual ~Rho0FromEtaPrime () ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the function
      double operator() ( const double x ) const  override;
      // ======================================================================
    private:
      // ======================================================================
      double m_eta_prime ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Flatte
     *
     *  S.M. Flatte
     *  "Coupled-channel analysis of the \f$\pi\eta\f$ and \f$K\bar{K}\f$
     *  systems near \f$K\bar{K}\f$threshold"
     *  Physics Letters B, Volume 63, Issue 2, 19 July 1976, Pages 224-227
     *
     *  http://www.sciencedirect.com/science/article/pii/0370269376906547
     *
     *  \f$\pi\pi\f$-channel
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Flatte
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor  from all parameters
       *  \f$ f \rigtharrow A_1 + A_2\f$
       *  @param m0    the mass
       *  @param m0g1  parameter \f$ m_0\times g_1\f$
       *  @param g2og2 parameter \f$ g2/g_1       \f$
       *  @param mA1   mass of A1
       *  @param mA2   mass of A2
       *  @param mB1   mass of B1
       *  @param mB2   mass of B2
       */
      Flatte  ( const double m0    = 980   ,
                const double m0g1  = 165   ,
                const double g2og1 = 4.21  ,
                const double mA1   = 139.6 ,
                const double mA2   = 139.6 ,
                const double mB1   = 493.7 ,
                const double mB2   = 493.7 ) ;
      /// destructor
      virtual ~Flatte () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of Flatte function
      virtual double operator() ( const double x ) const ;
      /// get the value of Flatte amplitude
      std::complex<double> amplitude ( const double x ) const
      { return flatte_amp ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude for pipi-channel
      std::complex<double> flatte_amp  ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the curve for pipi-channel
      double flatte  ( const double x ) const ;
      /// get the curve for   KK-channel
      double flatte2 ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m0     () const { return m_m0      ; }
      double mass   () const { return   m0   () ; }
      double peak   () const { return   m0   () ; }
      double m0g1   () const { return m_m0g1    ; }
      double g2og1  () const { return m_g2og1   ; }
      double mA1    () const { return m_A1      ; }
      double mA2    () const { return m_A2      ; }
      double mB1    () const { return m_B1      ; }
      double mB2    () const { return m_B2      ; }
      // ======================================================================
    public:
      // ======================================================================
      double thresholdA () const { return mA1() + mA2() ; }
      double thresholdB () const { return mB1() + mB2() ; }
      double threshold  () const
      { return std::min ( thresholdA () , thresholdB () ) ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0     ( const double x ) ;
      bool setMass   ( const double x ) { return setM0 ( x ) ; }
      bool setPeak   ( const double x ) { return setM0 ( x ) ; }
      bool setM0G1   ( const double x ) ;
      bool setG2oG1  ( const double x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      virtual double integral () const ;
      /// get the integral between low and high limits
      virtual double integral  ( const double low  ,
                                 const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_m0     ;
      double m_m0g1   ;
      double m_g2og1  ;
      double m_A1     ;
      double m_A2     ;
      double m_B1     ;
      double m_B2     ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Flatte2
     *
     *  S.M. Flatte
     *  "Coupled-channel analysis of the \f$\pi\eta\f$ and \f$K\bar{K}\f$
     *  systems near \f$K\bar{K}\f$threshold"
     *  Physics Letters B, Volume 63, Issue 2, 19 July 1976, Pages 224-227
     *
     *  http://www.sciencedirect.com/science/article/pii/0370269376906547
     *
     *  KK-channel
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Flatte2 : public Ostap::Math::Flatte
    {
    public:
      // ======================================================================
      /** constructor  from all parameters
       *  \f$ f \rigtharrow B_1 + B_2\f$
       *  @param m0    the mass
       *  @param m0g1  parameter \f$ m_0\times g_1\f$
       *  @param g2og2 parameter \f$ g2/g_1       \f$
       *  @param mA1   mass of A1
       *  @param mA2   mass of A2
       *  @param mB1   mass of B1
       *  @param mB2   mass of B2
       */
      Flatte2 ( const double m0    = 980   ,
                const double m0g1  = 165   ,
                const double g2og1 = 4.21  ,
                const double mA1   = 139.6 ,
                const double mA2   = 139.6 ,
                const double mB1   = 493.7 ,
                const double mB2   = 493.7 ) ;
      /// constructor  from Flatte
      Flatte2 ( const Flatte& flatte ) ;
      /// destructor
      virtual ~Flatte2 () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of Flatte function (KK-channel)
      double operator() ( const double x ) const  override;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Voigt
     *  simple Voigtian function:
     *  convolution of Lorenzian (non-relativistic Breit-Wigner function)
     *  with Gaussian resoltuion
     *  @see http://en.wikipedia.org/wiki/Voigt_profile
     *  The implementation relied on Faddeeva function
     *  @see http://en.wikipedia.org/wiki/Faddeeva_function
     *  @see http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Voigt
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      ///  constructor  from the three parameters
      Voigt  ( const double m0     = 1      ,
               const double gamma  = 0.004  ,
               const double sigma  = 0.001  ) ;
      /// destructor
      virtual ~Voigt () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of Voigt function
      // ======================================================================
      virtual double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m0     () const { return m_m0      ; }
      double mass   () const { return   m0   () ; }
      double peak   () const { return   m0   () ; }
      double gamma  () const { return m_gamma   ; }
      double sigma  () const { return m_sigma   ; }
      // ======================================================================
      /** full width at half maximum
       *  @see http://en.wikipedia.org/wiki/Voigt_profile
       */
      double fwhm   () const ;
      // ======================================================================
    public:
      // ======================================================================
      bool setM0     ( const double x ) ;
      bool setMass   ( const double x ) { return setM0 ( x ) ; }
      bool setPeak   ( const double x ) { return setM0 ( x ) ; }
      bool setGamma  ( const double x ) ;
      bool setSigma  ( const double x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      virtual double integral () const ;
      /// get the integral between low and high limits
      virtual double integral ( const double low  ,
                                const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_m0     ;
      double m_gamma  ;
      double m_sigma  ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PseudoVoigt
     *  Simplified verison of Voigt profile
     *  @see T. Ida, M. Ando and H. Toraya,
     *       "Extended pseudo-Voigt function for approximating the Voigt profile"
     *       J. Appl. Cryst. (2000). 33, 1311-1316
     *  @see doi:10.1107/S0021889800010219
     *  @see http://dx.doi.org/10.1107/S0021889800010219
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-13
     */
    class  PseudoVoigt
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      ///  constructor  from the three parameters
      PseudoVoigt  ( const double m0     = 1      ,
                     const double gamma  = 0.004  ,
                     const double sigma  = 0.001  ) ;
      /// destructor
      virtual ~PseudoVoigt () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of Voigt function
      // ======================================================================
      virtual double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m0     () const { return m_m0      ; }
      double mass   () const { return   m0   () ; }
      double peak   () const { return   m0   () ; }
      double gamma  () const { return m_gamma   ; }
      double sigma  () const { return m_sigma   ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0     ( const double x ) ;
      bool setMass   ( const double x ) { return setM0 ( x ) ; }
      bool setPeak   ( const double x ) { return setM0 ( x ) ; }
      bool setGamma  ( const double x ) ;
      bool setSigma  ( const double x ) ;
      // ======================================================================
    public: // helper constants
      // ======================================================================
      double fwhm_gauss      () const ;
      double fwhm_lorentzian () const { return 2 * m_gamma ; }
      double rho             () const
      { return fwhm_lorentzian() / ( fwhm_lorentzian() + fwhm_gauss() ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      virtual double integral () const ;
      /// get the integral between low and high limits
      virtual double integral ( const double low  ,
                                const double high ) const ;
      // ======================================================================
    public: // get parameters of the four components
      // ======================================================================
      /// get widths of the components
      double w   ( const unsigned short i ) const { return i < 4 ? m_w  [i] : 0.0 ; }
      /// get stength of the components
      double eta ( const unsigned short i ) const { return i < 4 ? m_eta[i] : 0.0 ; }
      // ======================================================================
    public: // get the separate components
      // ======================================================================
      /// get the Gaussian component
      double gaussian   ( const double x ) const ;
      /// get the Lorentzian component
      double lorentzian ( const double x ) const ;
      /// get the Irrational  component
      double irrational ( const double x ) const ;
      /// get the squared hyperbolic secant component
      double sech2      ( const double x ) const ;
      // ======================================================================
    private:  // calculate internal data
      // ======================================================================
      void update () ;
      // ======================================================================
    private:
      // ======================================================================
      double m_m0     ;
      double m_gamma  ;
      double m_sigma  ;
      // ======================================================================
    private: // some data
      // ======================================================================
      /// the widths/gammas of four components: Gaussian,Lorentzian,rIrational and Sech2
      std::vector<double> m_w   ;
      //// the strengths of four components
      std::vector<double>  m_eta ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    namespace Jackson
    {
      // ======================================================================
      /** the simplest function: constant
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      
      double jackson_0 ( double /* m  */ ,
                         double /* m0 */ ,
                         double /* m1 */ ,
                         double /* m2 */ ) ;
      // ======================================================================
      /** the simple function for \f$ 1^- \rightarrow 0^- 0^- \f$, l = 1
       *  \f$\rho(\omega)= \omega^{-1} \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      
      double jackson_A2 ( double    m     ,
                          double /* m0 */ ,
                          double /* m1 */ ,
                          double /* m2 */ ) ;
      // ======================================================================
      /** the simple function for \f$ 1^- \rightarrow 0^- 1^- \f$, l = 1
       *  \f$\rho(\omega)= \omega \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      
      double jackson_A3 ( double    m     ,
                          double /* m0 */ ,
                          double /* m1 */ ,
                          double /* m2 */ ) ;
      // ======================================================================
      /** the simple function for
       *  \f$ \frac{3}{2}^+ \rightarrow \frac{1}{2}^+ 0^- \f$, l = 1
       *  $\rho(\omega)= \frac{ ( \omega + M )^2 - m^2 }{ \omega^2} \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @param m1 the invariant mass of the first  (spinor) particle
       *  @param m2 the invariant mass of the secodn (scalar) particle
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      
      double jackson_A4 ( double    m     ,
                          double /* m0 */ ,
                          double    m1    ,
                          double    m2     ) ;
      // ======================================================================
      /** the simple function for
       *  \f$ \frac{3}{2}^- \rightarrow \frac{1}{2}^+ 0^- \f$, l = 2
       *  $\rho(\omega)= \left[ ( \omega + M )^2 - m^2 \right]^{-1} \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @param m1 the invariant mass of the first  (spinor) particle
       *  @param m2 the invariant mass of the secodn (scalar) particle
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      
      double jackson_A5 ( double    m     ,
                          double /* m0 */ ,
                          double    m1    ,
                          double    m2     ) ;
      // ======================================================================
      /** the simple function for \f$\rho^- \rightarrow \pi^+ \pi^-\f$
       *  \f$ 1- \rightarrow 0^- 0^- \f$, l = 1
       *  $\rho(\omega)= \left[ q_0^2 + q^2 \right]^{-1}f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @param m the nominam   mass
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      
      double jackson_A7 ( double    m     ,
                          double    m0    ,
                          double    m1    ,
                          double    m2    ) ;
      // ======================================================================
    } //                                               end of namespace Jackson
    // ========================================================================
    /** @class FormFactor
     *  abstract class to implement various formfactors
     */
    class FormFactor
    {
    public :
      // ======================================================================
      /// the only important method
      virtual double operator()
      ( const double m  , const double m0 ,
        const double m1 , const double m2 ) const  = 0 ;
      /// virtual destructor
      virtual ~FormFactor () ;
      /// clone method ("virtual constructor" )
      virtual  FormFactor* clone() const = 0 ;
      // ======================================================================
    } ;
    // ========================================================================
    namespace FormFactors
    {
      // ======================================================================
      /** Formfactor for Breit-Wigner amplitude
       *  parameterization for \f$\rho(\omega)\f$-function from (A.1)
       *  J.D.Jackson,
       *  "Remarks on the Phenomenological Analysis of Resonances",
       *  In Nuovo Cimento, Vol. XXXIV, N.6
       */
      class Jackson : public Ostap::Math::FormFactor
      {
      public:
        // ====================================================================
        /// default constructor
        Jackson () ;
        /// constructor from enum
        Jackson ( const Ostap::Math::FormFactors::JacksonRho rho ) ;
        /// constructor from rho-function
        Jackson (       Ostap::Math::FormFactors::rho_fun    rho ) ;
        /// virtual destructor
        virtual ~Jackson  () ;
        /// clone method ("virtual constructor")
        Jackson* clone() const   override;
        /// the only important method
        double operator() ( const double m  , const double m0 ,
                            const double m1 , const double m2 ) const override;
        // ====================================================================
      private:
        // ====================================================================
        /// the finction itself
        Ostap::Math::FormFactors::rho_fun m_rho ; // the finction itself
        // ====================================================================
      } ;
      // ======================================================================
      /** Blatt-Weisskopf formfactor/barrier factor
       *  actually it is "traslation" of
       *  Blatt-Weiskopf barrier factor into in Jackson's "rho"-function
       */
      class BlattWeisskopf : public Ostap::Math::FormFactor
      {
      public:
        // ====================================================================
        /// orbital momentum
        enum Case {
          Zero  = 0 ,
          One   = 1 ,
          Two   = 2 ,
          Three = 3 ,
          Four  = 4 ,
          Five  = 5
        } ;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor from enum and barrier factor
        BlattWeisskopf ( const Case   L , const double b ) ;
        /// default constructor (needed for  serialization)
        BlattWeisskopf () ;
        /// virtual destructor
        virtual ~BlattWeisskopf () ;
        /// clone method ("virtual constructor")
        BlattWeisskopf* clone() const   override;
        /// the only important method
        double operator() ( const double m  , const double m0 ,
                            const double m1 , const double m2 ) const override;
        // ====================================================================
      protected:
        // ====================================================================
        /// get the barrier factor
        double   b ( const double z , const double z0 ) const ;
        // ====================================================================
      private:
        // ====================================================================
        Case   m_L ;
        double m_b ;
        // ====================================================================
      } ;
      // ======================================================================
    } // end of namespace Ostap:Math::FormFactors
    // ========================================================================
    /** @class Swanson
     *  Swanson's parameterization of S-wave cusp
     *  @see LHCb-PAPER-2016-019 appendix D
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-11
     */
    class  Swanson : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /// constructor from all parameters (numbers are arbitrary...)
      Swanson ( const double         m1     = 0.139 ,   // the first  real particle
                const double         m2     = 0.139 ,   // the second real particle
                const double         m1_0   = 0.135 ,   // the first  particle for cusp
                const double         m2_0   = 0.135 ,   // the second particle for cusp
                const double         beta_0 = 0.300 ,   // beta_0 parameter
                const unsigned short L      = 0     ) ; // orbital momentum for real particles
      /// constructor from all parameters
      Swanson ( const double         m1             ,   // the first  real particle
                const double         m2             ,   // the second real particle
                const double         m1_0           ,   // the first  particle for cusp
                const double         m2_0           ,   // the second particle for cusp
                const double         beta_0         ,   // beta_0 parameter
                const unsigned short L              ,   // orbital momentum for real particles
                const Ostap::Math::FormFactors::JacksonRho  r ) ; //  formfactor
      /// constructor from all parameters
      Swanson ( const double         m1             ,   // the first  real particle
                const double         m2             ,   // the second real particle
                const double         m1_0           ,   // the first particle for cusp
                const double         m2_0           ,   // the second particle for cusp
                const double         beta_0         ,   // beta_0 parameter
                const unsigned short L              ,   // orbital momentum for real particles
                const Ostap::Math::FormFactor&    f ) ; // formfactor
      /// constructor from all parameters
      Swanson ( const BreitWigner&   bw             ,   // breit-wigner
                const double         m1_0           ,   // the first  particle for cusp
                const double         m2_0           ,   // the second particle for cusp
                const double         beta_0         ) ; // beta_0 parameter
      /// copy constructor
      Swanson ( const Swanson&  sw ) ;
      /// destructor
      virtual ~Swanson() ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the Swanson shape
      double operator () ( const double x ) const { return swanson ( x ) ; }
      /// calculate the Swanson shape
      double swanson     ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate complex amplitude
      std::complex<double> amplitude ( const double x ) const ;
      // ======================================================================
    public: // getters
      // ======================================================================
      /// get beta_0 parameter
      double  beta0 () const { return m_beta0 ; }  // get beta_0 parameter
      /// mass of the first particle
      double  m1    () const { return m_m1    ; }  // mass of first particle
      /// mass of the second particle
      double  m2    () const { return m_m2    ; }  // mass of the second particle
      // ======================================================================
    public: // derived getters
      // ======================================================================
      double mmin () const { return m_bw.m1() + m_bw.m2() ; }
      double cusp () const { return    m_m1   +    m_m2   ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      /// set new value for beta_0
      bool setBeta0  ( const double value ) ;
      /// set new value for beta_0
      bool setBeta_0 ( const double value ) { return setBeta0 ( value ) ; }
      /// set new valeu for m1
      bool setM1_0   ( const double value ) ;
      /// set new valeu for m2
      bool setM2_0   ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high limits
      virtual double integral  ( const double low  ,
                                 const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// assignement operator is disabled
      Swanson& operator=( const Swanson& ) ; // no assignement
      // ======================================================================
    private:
      // ======================================================================
      /// use Breit-Wigner to keep parameters of real particles
      Ostap::Math::BreitWigner   m_bw ;
      /// the mass of the first  particle
      double            m_m1         ; // the mass of the first  particle
      /// the mass of the second particle
      double            m_m2         ; // the mass of the second particle
      /// beta0 parameter
      double            m_beta0      ; // beta0 parameter
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LASS
     *  The LASS parameterization (Nucl. Phys. B296, 493 (1988))
     *  describes the 0+ component of the Kpi spectrum.
     *  It consists of the K*(1430) resonance together with an
     *  effective range non-resonant component
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-10-05
     */
    class  LASS : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param a  the LASS parameter
       *  @param r  the LASS parameter
       *  @param e  the LASS parameter
       */
      LASS ( const double         m1 =  493.7  ,
             const double         m2 =  139.6  ,
             const double         m0 = 1435    , // K*(1450) mass
             const double         g0 =  279    , // K*(1430) width
             const double         a  = 1.94e-3 ,
             const double         r  = 1.76e-3 ,
             const double         e  = 1.0     ) ;
      /// destructor
      virtual ~LASS () ;                                          // destructor
      // ======================================================================
    public:
      // ======================================================================
      /// get the (complex) LASS amplitude
      std::complex<double> amplitude  ( const double x ) const ;
      /// get the phase space factor
      double               phaseSpace ( const double x ) const ;
      /// evaluate LASS
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m0  ( ) const { return m_m0 ; }
      double g0  ( ) const { return m_g0 ; }
      double a   ( ) const { return m_a  ; }
      double r   ( ) const { return m_r  ; }
      double e   ( ) const { return m_e  ; }
      // ======================================================================
      double m1  ( ) const { return m_ps2.m1 () ; }
      double m2  ( ) const { return m_ps2.m2 () ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0 ( const double value ) ;
      bool setG0 ( const double value ) ;
      bool setA  ( const double value ) ;
      bool setR  ( const double value ) ;
      bool setE  ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the pole position for scalar meson
      double   m_m0 ;
      double   m_g0 ;
      /// LASS-parameters
      double   m_a  ;
      double   m_r  ;
      double   m_e  ;
      // ======================================================================
    private:
      // ======================================================================
      /// phase space
      Ostap::Math::PhaseSpace2 m_ps2     ;    // phase space
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LASS23L
     *
     *  The LASS parameterization (Nucl. Phys. B296, 493 (1988))
     *  describes the 0+ component of the Kpi spectrum.
     *  It consists of the K*(1430) resonance together with an
     *  effective range non-resonant component.
     *
     *  This function is suitable to describe the S-wave Kpi distribtion
     *  from X -> K pi Y decay
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-01
     */
    class  LASS23L : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and the third
       *  @param a  the LASS parameter
       *  @param r  the LASS parameter
       */
      LASS23L ( const double         m1 =  493.7  ,
                const double         m2 =  139.6  ,
                const double         m3 = 3097    ,
                const double         m  = 5278    ,
                const double         m0 = 1435    ,
                const double         g0 =  279    ,
                const unsigned short L  =    1    ,
                const double         a  = 1.94e-3 ,
                const double         r  = 1.76e-3 ,
                const double         e  = 1.0     ) ;
      /** constructor from LASS and 3-rd particle
       *  @param lass the actual lass shape
       *  @param m3   the mass of third particle (Y)
       *  @param m    the mass of mother particle (X)
       *  @param L    the orbital momentum between Y and (Kpi)
       */
      LASS23L ( const LASS&          lass   ,
                const double         m3     , // the third particle, e.g. J/psi
                const double         m      , // mother particle, e.g. B
                const unsigned short L  = 1 ) ;
      /// destructor
      virtual ~LASS23L () ;                                     // destructor
      // ======================================================================
    public:
      // ======================================================================
      /// get the (complex) LASS amplitude
      std::complex<double> amplitude  ( const double x ) const ;
      /// get the phase space factor
      double               phaseSpace ( const double x ) const ;
      /// evaluate LASS
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m0  () const { return m_lass . m0 () ; } // K*(1430) mass
      double g0  () const { return m_lass . g0 () ; } // K*(1430) width
      double a   () const { return m_lass . a  () ; }
      double r   () const { return m_lass . r  () ; }
      double e   () const { return m_lass . e  () ; }
      // ======================================================================
      double m1  () const { return m_ps   . m1 () ; }
      double m2  () const { return m_ps   . m2 () ; }
      double m3  () const { return m_ps   . m3 () ; }
      double m   () const { return m_ps   . m  () ; }
      double l   () const { return m_ps   . l  () ; }
      double L   () const { return m_ps   . L  () ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0 ( const double value ) { return m_lass . setM0 ( value ) ; }
      bool setG0 ( const double value ) { return m_lass . setG0 ( value ) ; }
      bool setA  ( const double value ) { return m_lass . setA  ( value ) ; }
      bool setR  ( const double value ) { return m_lass . setR  ( value ) ; }
      bool setE  ( const double value ) { return m_lass . setE  ( value ) ; }
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
      /// Lass itself
      Ostap::Math::LASS          m_lass  ;    // lass-function itself
      /// phase space
      Ostap::Math::PhaseSpace23L m_ps    ;    // phase space
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bugg
     *  parametrisation of sigma-pole for
     *  two pion mass distribution
     *
     *  The parameterization of sigma pole by
     *  B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-01
     */
    class  Bugg
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param M  mass of sigma (very different from the pole positon!)
       *  @param g2 width parameter g2 (4pi width)
       *  @param b1 width parameter b1  (2pi coupling)
       *  @param b2 width parameter b2  (2pi coupling)
       *  @param s1 width parameter s1  (cut-off for 4pi coupling)
       *  @param s2 width parameter s2  (cut-off for 4pi coupling)
       *  @param a  parameter a (the exponential cut-off)
       *  @param m1 the mass of the first  particle
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and the third
       */
      Bugg    ( const double         M  = 0.9264        ,  // GeV
                const double         g2 = 0.0024        ,  // GeV
                const double         b1 = 0.5848        ,  // GeV
                const double         b2 = 1.6663        ,  // GeV-1
                const double         a  = 1.082         ,  // GeV^2
                const double         s1 = 2.8           ,  // GeV^2
                const double         s2 = 3.5           ,
                const double         m1 =  139.6 / 1000 ) ; // GeV
      /// destructor
      ~Bugg () ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the Bugg shape
      double operator() ( const double x ) const { return pdf ( x ) ; }
      double pdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the amlitude  (not normalized!)
      std::complex<double> amplitude (  const double x ) const ;
      /// get the phase space factor (taking into account L)
      double phaseSpace ( const double x ) const { return m_ps  ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // phase space variables
      // ======================================================================
      double m1        () const { return m_ps.m1 () ; }
      double m2        () const { return m_ps.m2 () ; }
      // ======================================================================
      double lowEdge   () const { return m_ps. lowEdge() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the running width by Bugg
      std::complex<double>
      gamma ( const double x ) const ; // get the running width by Bugg
      // ======================================================================
    public:
      // ======================================================================
      /// adler factor
      double               adler       ( const double x ) const ; // adler factor
      /// ratio of 2pi-phase spaces
      double               rho2_ratio  ( const double x ) const ;
      /// ratio of 4pi-phase spaces
      std::complex<double> rho4_ratio  ( const double x ) const ;
      /// b-factor for 2-pi coupling
      double b ( const double x ) const { return  b1 () + x * x * b2 () ; }
      // ======================================================================
    private:
      // ======================================================================
      /// approximation for  4pi-phase space
      std::complex<double> rho4        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // sigma & Bugg variables
      // ======================================================================
      double M     () const  { return m_M       ; }
      double M2    () const  { return m_M * m_M ; }
      double m0    () const  { return   M ()    ; }
      double mass  () const  { return   M ()    ; }
      double peak  () const  { return   M ()    ; }
      // ======================================================================
      double g2    () const  { return m_g2   ; }
      double b1    () const  { return m_b1   ; }
      double b2    () const  { return m_b2   ; }
      double s1    () const  { return m_s1   ; }
      double s2    () const  { return m_s2   ; }
      double a     () const  { return m_a    ; }
      // ======================================================================
      bool setM    ( const double value  ) ;
      bool setM0   ( const double value  ) { return setM ( value )  ; }
      bool setMass ( const double value  ) { return setM ( value )  ; }
      bool setPeak ( const double value  ) { return setM ( value )  ; }
      // ======================================================================
      bool setG2   ( const double value  ) ;
      bool setB1   ( const double value  ) ;
      bool setB2   ( const double value  ) ;
      bool setS1   ( const double value  ) ;
      bool setS2   ( const double value  ) ;
      bool setA    ( const double value  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      /// double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      // sigma & Bugg varibales
      // ======================================================================
      /// mass of sigma (very different from the pole positon!)
      double m_M  ; // mass of sigma (very different from the pole positon!)
      /// width parameter g2 (4pi width)
      double m_g2 ; // width parameter g2 (4-p width)
      /// width parameter b1  (2pi coupling)
      double m_b1 ; // width parameter b1  (2pi coupling)
      /// width parameter b2  (2pi coupling)
      double m_b2 ; // width parameter b2  (2pi coupling)
      /// width parameter s1  (cut-off for 4pi coupling)
      double m_s1 ; // width parameter b1  (cut-off for 4pi coupling)
      /// width parameter s2  (cut-off for 4pi coupling)
      double m_s2 ; // width parameter b2  (cut-off for 4pi coupling)
      /// parameter a (the exponential cut-off)
      double m_a  ; // parameter a (the exponential cut-off)
      // ======================================================================
    private:
      // ======================================================================
      /// phase space
      Ostap::Math::PhaseSpace2   m_ps         ; // phase space
      // ======================================================================
    private:
      /// integration workspace
      Ostap::Math::WorkSpace     m_workspace  ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bugg23L
     *  parametrisation of sigma-pole for
     *  two pion mass distribution from three body decays
     *
     *  The parameterization of sigma pole by
     *  B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-01
     */
    class  Bugg23L
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param M  mass of sigma (very different from the pole positon!)
       *  @param g2 width parameter g2 (4pi width)
       *  @param b1 width parameter b1  (2pi coupling)
       *  @param b2 width parameter b2  (2pi coupling)
       *  @param s1 width parameter s1  (cut-off for 4pi coupling)
       *  @param s2 width parameter s2  (cut-off for 4pi coupling)
       *  @param a  parameter a (the exponential cut-off)
       *  @param m1 the mass of the first  particle
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and the third
       */
      Bugg23L ( const double         M  = 0.9264        ,  // GeV
                const double         g2 = 0.0024        ,  // GeV
                const double         b1 = 0.5848        ,  // GeV
                const double         b2 = 1.6663        ,  // GeV-1
                const double         a  = 1.082         ,  // GeV^2
                const double         s1 = 2.8           ,  // GeV^2
                const double         s2 = 3.5           ,
                const double         m1 =  139.6 / 1000 ,  // MeV
                const double         m3 = 3097.0 / 1000 ,  // MeV
                const double         m  = 5278.0 / 1000 ,  // MeV
                const unsigned short L  =    1          ) ;
      /** constructor from bugg & phase space parameters
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and the third
       */
      Bugg23L ( const Ostap::Math::Bugg& bugg               ,
                const double             m3 = 3097.0 / 1000 ,  // MeV
                const double             m  = 5278.0 / 1000 ,  // MeV
                const unsigned short     L  =    1          ) ;


      /// destructor
      ~Bugg23L () ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the Bugg shape
      double operator() ( const double x ) const { return pdf ( x ) ; }
      double pdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the amlitude  (not normalized!)
      std::complex<double> amplitude (  const double x ) const
      { return m_bugg.amplitude ( x ) ; }
      /// get the phase space factor (taking into account L)
      double phaseSpace ( const double x ) const { return m_ps  ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // phase space variables
      // ======================================================================
      double m1        () const { return m_ps.m1 () ; }
      double m2        () const { return m_ps.m2 () ; }
      double m3        () const { return m_ps.m3 () ; }
      double m         () const { return m_ps.m  () ; }
      // ======================================================================
      double lowEdge   () const { return m_ps. lowEdge() ; }
      double highEdge  () const { return m_ps.highEdge() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the running width by Bugg
      std::complex<double>
      gamma ( const double x ) const { return m_bugg.gamma ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// adler factor
      double               adler       ( const double x ) const
      { return m_bugg.adler      ( x ) ; } // adler factor
      /// ratio of 2pi-phase spaces
      double               rho2_ratio  ( const double x ) const
      { return m_bugg.rho2_ratio ( x ) ; }
      /// ratio of 4pi-phase spaces
      std::complex<double> rho4_ratio  ( const double x ) const
      { return m_bugg.rho4_ratio ( x ) ; }
      /// b-factor for 2-pi coupling
      double b ( const double x ) const { return m_bugg. b ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // sigma & Bugg variables
      // ======================================================================
      double M     () const  { return m_bugg. M    () ; }
      double M2    () const  { return m_bugg. M2   () ; }
      double m0    () const  { return m_bugg. m0   () ; }
      double mass  () const  { return m_bugg. mass () ; }
      double peak  () const  { return m_bugg. peak () ; }
      // ======================================================================
      double g2    () const  { return m_bugg. g2   () ; }
      double b1    () const  { return m_bugg. b1   () ; }
      double b2    () const  { return m_bugg. b2   () ; }
      double s1    () const  { return m_bugg. s1   () ; }
      double s2    () const  { return m_bugg. s2   () ; }
      double a     () const  { return m_bugg. a    () ; }
      // ======================================================================
      bool setM    ( const double value  ) { return m_bugg.setM    ( value ) ; }
      bool setM0   ( const double value  ) { return m_bugg.setM0   ( value ) ; }
      bool setMass ( const double value  ) { return m_bugg.setMass ( value ) ; }
      bool setPeak ( const double value  ) { return m_bugg.setPeak ( value ) ; }
      // ======================================================================
      bool setG2   ( const double value  ) { return m_bugg.setG2   ( value ) ; }
      bool setB1   ( const double value  ) { return m_bugg.setB1   ( value ) ; }
      bool setB2   ( const double value  ) { return m_bugg.setB2   ( value ) ; }
      bool setS1   ( const double value  ) { return m_bugg.setS1   ( value ) ; }
      bool setS2   ( const double value  ) { return m_bugg.setS2   ( value ) ; }
      bool setA    ( const double value  ) { return m_bugg.setA    ( value ) ; }
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
      /// bugg function
      Ostap::Math::Bugg          m_bugg       ; // bugg function
      /// phase space
      Ostap::Math::PhaseSpace23L m_ps         ; // phase space
      // ======================================================================
    private:
      /// integration workspace
      Ostap::Math::WorkSpace     m_workspace  ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BW23L
     *  @see Ostap::Math::BreitWigner
     *  @see Ostap::Math::PhaseSpace23L
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2012-05-23
     */
    class  BW23L
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      // constructor from all parameters
      BW23L ( const double         m0       = 0.770 ,
              const double         gam0     = 0.150 ,
              const double         m1       = 0.139 ,
              const double         m2       = 0.139 ,
              const double         m3       = 3.096 ,
              const double         m        = 5.278 ,
              const unsigned short L1       = 0     ,
              const unsigned short L2       = 0     ) ;
      // constructor from all parameters
      BW23L ( const double         m0       ,
              const double         gam0     ,
              const double         m1       ,
              const double         m2       ,
              const double         m3       ,
              const double         m        ,
              const unsigned short L1       ,
              const unsigned short L2       ,
              const Ostap::Math::FormFactors::JacksonRho r ) ;
      /// constructor from BreitWigner
      BW23L ( const Ostap::Math::BreitWigner& bw ,
              const double                    m3 ,
              const double                    m  ,
              const unsigned short            L2 ) ;
      /// destructor
      virtual ~BW23L () ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the shape
      double operator() ( const double x ) const ;
      /// get the amplitude
      std::complex<double>
      amplitude ( const double x ) const { return m_bw.amplitude ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double m0     () const { return m_bw . m0   () ; }
      double mass   () const { return        m0   () ; }
      double peak   () const { return        m0   () ; }
      double gam0   () const { return m_bw . gam0 () ; }
      double gamma0 () const { return        gam0 () ; }
      double gamma  () const { return        gam0 () ; }
      double width  () const { return        gam0 () ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0     ( const double x ) { return m_bw.setM0     ( x ) ; }
      bool setMass   ( const double x ) { return setM0          ( x ) ; }
      bool setPeak   ( const double x ) { return setM0          ( x ) ; }
      bool setGamma0 ( const double x ) { return m_bw.setGamma0 ( x ) ; }
      bool setGamma  ( const double x ) { return setGamma0      ( x ) ; }
      bool setWidth  ( const double x ) { return setGamma0      ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double lowEdge   () const { return m_ps. lowEdge() ; }
      double highEdge  () const { return m_ps.highEdge() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the current width
      double gamma ( const double x ) const { return m_bw.gamma ( x ) ; }
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
      /// the breit wigner
      Ostap::Math::BreitWigner   m_bw        ;    // the breit wigner
      /// the phase space
      Ostap::Math::PhaseSpace23L m_ps        ;    // the phase space
      /// integration workspace
      Ostap::Math::WorkSpace     m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Flatte23L
     *  \f$\pi\pi\f$-channel
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Flatte23L
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor  from all parameters
       *  @param m0    the mass
       *  @param m0g1  parameter \f$ m_0\times g_1\f$
       *  @param g2og2 parameter \f$ g2/g_1       \f$
       *  @param mA    A mass
       *  @param mB    B mass
       *  @param m3    the mass of the third particle
       *  @param m     the mass of mother particle
       *  @param L     the orbital momentum between the pair and the third particle
       */
      Flatte23L  ( const double         m0    =  980.0   ,     // MeV
                   const double         m0g1  =  165     ,     // MeV^2
                   const double         g2og1 =    4.21  ,     // dimensionless
                   const double         mK    =  493.7   ,     // MeV
                   const double         mPi   =  139.6   ,     // MeV
                   const double         m3    = 3096.9   ,     // MeV
                   const double         m     = 5366.0   ,     // MeV
                   const unsigned short L     = 1        ) ;
      // ======================================================================
      /** constructor  from flatte function
       *  @param m3    the mass of the third particle
       *  @param m     the mass of mother particle
       *  @param L     the orbital momentum between the pair and the third particle
       */
      Flatte23L  ( const Flatte&        fun              ,     // MeV
                   const double         m3    = 3096.9   ,     // MeV
                   const double         m     = 5366.0   ,     // MeV
                   const unsigned short L     = 1        ) ;
      /// destructor
      virtual ~Flatte23L () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of Flatte function
      // ======================================================================
      double operator() ( const double x ) const ;
      // ======================================================================
      /// get the value of complex Flatte amplitude (pipi-channel)
      std::complex<double> amplitude ( const double x ) const
      { return m_flatte.flatte_amp ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double m0     () const { return m_flatte . m0    () ; }
      double mass   () const { return            m0    () ; }
      double peak   () const { return            m0    () ; }
      double m0g1   () const { return m_flatte . m0g1  () ; }
      double g2og1  () const { return m_flatte . g2og1 () ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setM0     ( const double x ) { return m_flatte . setM0    ( x ) ; }
      bool setMass   ( const double x ) { return            setM0    ( x ) ; }
      bool setPeak   ( const double x ) { return            setM0    ( x ) ; }
      bool setM0G1   ( const double x ) { return m_flatte . setM0G1  ( x ) ; }
      bool setG2oG1  ( const double x ) { return m_flatte . setG2oG1 ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double lowEdge   () const { return m_ps .  lowEdge () ; }
      double highEdge  () const { return m_ps . highEdge () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      virtual double integral () const ;
      /// get the integral between low and high limits
      virtual double integral ( const double low  ,
                                const double high ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual Flatte function
      Ostap::Math::Flatte        m_flatte ; // the actual Flatte function
      /// phase space factor
      Ostap::Math::PhaseSpace23L m_ps     ; // phase space factor
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Gounaris23L
     *  parametrisation of rho0 for
     *  two pion mass distribution
     *
     *  G.J.Gounaris and J.J.Sakurai,
     *  "Finite width corrections to the vector meson dominance
     *  predictions for \f$\rho\rightarrow e^+e^-\f$",
     *  Phys.Rev.Lett. 21 (1968) 244
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-01
     */
    class  Gounaris23L
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param M  mass of rho
       *  @param g0 width parameter
       *  @param m1 the mass of the first  particle (the same as the second)
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and the third
       */
      Gounaris23L ( const double         M  = 0.770         ,  // GeV
                    const double         g0 = 0.150         ,  // GeV
                    const double         m1 =  139.6 / 1000 ,  // MeV
                    const double         m3 = 3097.0 / 1000 ,  // MeV
                    const double         m  = 5278.0 / 1000 ,  // MeV
                    const unsigned short L  =    1          ) ;
      /// destructor
      ~Gounaris23L () ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the Gounaris-Sakurai shape
      double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the amlitude  (not normalized!)
      std::complex<double> amplitude (  const double x ) const ;
      /// get the phase space factor (taking into account L)
      double phaseSpace ( const double x ) const { return m_ps  ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // phase space variables
      // ======================================================================
      double m1        () const { return m_ps.m1 () ; }
      double m2        () const { return m_ps.m2 () ; }
      double m3        () const { return m_ps.m3 () ; }
      double m         () const { return m_ps.m  () ; }
      // ======================================================================
      double lowEdge   () const { return m_ps. lowEdge() ; }
      double highEdge  () const { return m_ps.highEdge() ; }
      // ======================================================================
    private:
      // ======================================================================
      /// get h-factor
      double h       ( const double x ) const ;
      /// get h-factor
      double h       ( const double x , const double k ) const ;
      /// get h'-factor
      double h_prime ( const double x ) const ;
      /// get h'-factor
      double h_prime ( const double x , const double k ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // Gounaris & Sakurai variables
      // ======================================================================
      double M      () const  { return m_M     ; }
      double m0     () const  { return   M  () ; }
      double mass   () const  { return   M  () ; }
      double peak   () const  { return   M  () ; }
      // ======================================================================
      double g0     () const  { return m_g0    ; }
      double gamma  () const  { return   g0 () ; }
      double width  () const  { return   g0 () ; }
      // ======================================================================
      bool setM     ( const double value  ) ;
      bool setM0    ( const double value  ) { return setM  ( value ) ; }
      bool setMass  ( const double value  ) { return setM  ( value ) ; }
      bool setPeak  ( const double value  ) { return setM  ( value ) ; }
      // ======================================================================
      bool setG0    ( const double value  ) ;
      bool setGamma ( const double value  ) { return setG0 ( value ) ; }
      bool setWidth ( const double value  ) { return setG0 ( value ) ; }
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
      //  Gounaris and Sakurai variables
      // ======================================================================
      /// mass of rho
      double m_M  ; // mass of sigma (very different from the pole positon!)
      /// width parameter
      double m_g0 ; // width parameter
      // ======================================================================
    private:
      // ======================================================================
      /// phase space
      Ostap::Math::PhaseSpace23L m_ps         ; // phase space
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace     m_workspace  ;    // integration workspace
      // ======================================================================
    } ;
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
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor from mass, resolution and "n"-parameter
       *  @param M     mass
       *  @param sigma width parameter
       *  @param N     n-parameter  ( actually  n=1+|N| )
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
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor from mass, resolution and "n"-parameter
       *  @param M     mass
       *  @param sigma width parameter
       *  @param N     n-parameter  ( actually  n=1+|N| )
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
    /** @class GammaDist
     *  Gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  GammaDist
      : public std::unary_function<double,double>
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
      : public std::unary_function<double,double>
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
      : public std::unary_function<double,double>
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
      : public std::unary_function<double,double>
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
      : public std::unary_function<double,double>
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
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param alpha \f$\alpha\f$-parameter
       *  @param beta  \f$\beta\f$-parameter
       *  @param scale scale-parameter
       *  @param low   shift-parameter
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
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param alpha \f$\alpha\f$-parameter
       *  @param beta  \f$\beta\f$-parameter
       *  @param scale scale-parameter
       *  @param low   shift-parameter
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
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param location \f$\mu\f$-parameter       \f$-\inf<\mu<+\inf\f$
       *  @param scale    \f$\sigma\f$-parameter    \f$0<\sigma\f$
       *  @param epsilon  \f$\epsilon\f$-parameter  \f$-\inf<\epsilon<+\inf\f$
       *  @param delta    \f$\delta\f$-parameter    \f$0<\epsilon<+\inf\f$
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
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param xi     \f$\xi\f$-parameter       \f$-\inf<\xi<+\inf\f$
       *  @param lambda \f$\lambda\f$-parameter   \f$   0<\lambda<+\inf\f$
       *  @param delta  \f$\delta\f$-parameter    \f$   0<\delta<+\inf\f$
       *  @param gamma  \f$\gamma\f$-parameter    \f$-\inf<\epsilon<+\inf\f$
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
     *  \f$  f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\deltax/2}}}{2})\f$,
     *  where \f$\delta x = \left| x - \mu \right|/\sigma\f$
     *  Function is taken from http://arxiv.org/abs/arXiv:1507.07099
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-08-21
     */
    class  Atlas : public std::unary_function<double,double>
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
     *  \f$ f(x,\mu,\sigma) \propto \frac{1}{2} \sech ( \frac{\pi}{2}\frac{x-\mu}{\sigma} )\f$
     *  @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-04-25
     */
    class  Sech : public std::unary_function<double,double>
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
    class  Logistic : public std::unary_function<double,double>
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
    /** @class Argus
     *  http://en.wikipedia.org/wiki/ARGUS_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     */
    class  Argus
      : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor with all parameters
       *  @param alpha \f$\alpha\f$-parameter
       *  @param beta  \f$\beta\f$-parameter
       *  @param scale scale-parameter
       *  @param low   shift-parameter
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
    class  ExpoPositive :  public std::unary_function<double,double>
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
    // 2D-models
    // ========================================================================
    /** @class PS2DPol
     *  The 2D-function:
     *  \f$ f(x,y) = Ps(x)*Ps(y)*P_{pos}(x,y) \f$, where
     *  \f$Ps\f$ denotes phase-space function and
     * \f$P_{pos}\f$ denotes the positive polynomial
     */
    class  PS2DPol
      : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPol ( const PhaseSpaceNL&   psx = PhaseSpaceNL () ,
                const PhaseSpaceNL&   psy = PhaseSpaceNL () ,
                const unsigned short  Nx  =  1 ,
                const unsigned short  Ny  =  1 ) ;
      /// constructor from the order
      PS2DPol ( const PhaseSpaceNL&   psx     ,
                const PhaseSpaceNL&   psy     ,
                const unsigned short  Nx      ,
                const unsigned short  Ny      ,
                const double          xmin    ,
                const double          xmax    ,
                const double          ymin    ,
                const double          ymax    ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_positive.npars () ; }
      /// set k-parameter
      bool setPar       ( const unsigned int k , const double value )
      { return m_positive.setPar ( k , value ) ;}
      /// set k-parameter
      bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned int k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get nX & nY
      unsigned short nX () const { return m_positive.nX () ; }
      unsigned short nY () const { return m_positive.nY () ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral ( const double xlow , const double xhigh ,
                        const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL& psX         () const { return m_psx      ; }
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_psy      ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceX () const { return psX ()     ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY ()     ; }
      const Ostap::Math::Positive2D&   positive    () const { return m_positive ; }
      const Ostap::Math::Positive2D&   polynom     () const { return m_positive ; }
      // ====================================== ===============================
    private:
      // ======================================================================
      /// the actual (positive) bernstein polynomial in 2D
      Ostap::Math::Positive2D   m_positive ; // the actual bernstein polynomial
      /// Phase space
      Ostap::Math::PhaseSpaceNL m_psx      ;
      Ostap::Math::PhaseSpaceNL m_psy      ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace   ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PS2DPolSym
     *  The symmetric 2D-function:
     *  \f$ f(x,y) = Ps(x)*Ps(y)*P_{sym}(x,y) \f$, where
     *  \f$Ps\f$ denotes phase-space function and
     * \f$P_{sym}\f$ denotes the symmetric positive polynomial
     */
    class  PS2DPolSym
      : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPolSym ( const PhaseSpaceNL&   ps = PhaseSpaceNL() ,
                   const unsigned short  N  =  1             ) ;
      /// constructor from the order
      PS2DPolSym ( const PhaseSpaceNL&   ps      ,
                   const unsigned short  N       ,
                   const double          xmin    ,
                   const double          xmax    ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_positive.npars () ; }
      /// set k-parameter
      bool setPar       ( const unsigned int k , const double value )
      { return m_positive.setPar ( k , value ) ;}
      /// set k-parameter
      bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned int k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral ( const double xlow , const double xhigh ,
                        const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL&  psX         () const { return m_ps       ; }
      const Ostap::Math::PhaseSpaceNL&  psY         () const { return m_ps       ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceX () const { return psX()      ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceY () const { return psY()      ; }
      const Ostap::Math::Positive2DSym& positive    () const { return m_positive ; }
      const Ostap::Math::Positive2DSym& polynom     () const { return m_positive ; }
      // ====================================== ===============================
    private:
      // ======================================================================
      /// the actual (positive) bernstein polynomial in 2D
      Ostap::Math::Positive2DSym m_positive ; // the actual bernstein polynomial
      /// Phase space
      Ostap::Math::PhaseSpaceNL m_ps        ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace for numerical integration
      Ostap::Math::WorkSpace m_workspace    ;
      // ======================================================================
    };
    // ========================================================================
    /** @class ExpoPS2DPol
     *  The 2D-function:
     *  \f$ f(x,y) = exp(tau*x)*Ps(y)*P_{pos}(x,y) \f$, where
     *  \f$Ps\f$ denotes phase-space function and
     * \f$P_{pos}\f$ denotes the positive polynomial
     */
    class  ExpoPS2DPol
      : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      ExpoPS2DPol ( const PhaseSpaceNL&   psy  = PhaseSpaceNL() ,
                    const double          xmin = 0 ,
                    const double          xmax = 1 ,
                    const unsigned short  Nx   = 1 ,
                    const unsigned short  Ny   = 1 ) ;
      /// constructor from the order
      ExpoPS2DPol ( const PhaseSpaceNL&   psy     ,
                    const double          xmin    ,
                    const double          xmax    ,
                    const unsigned short  Nx      ,
                    const unsigned short  Ny      ,
                    const double          ymin    ,
                    const double          ymax    ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_positive.npars () ; }
      /// set k-parameter
      bool setPar       ( const unsigned int k , const double value )
      { return m_positive.setPar ( k , value ) ;}
      /// set k-parameter
      bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned int k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get nX & nY
      unsigned short nX () const { return m_positive.nX () ; }
      unsigned short nY () const { return m_positive.nY () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get tau
      double         tau () const { return m_tau ;}
      /// set tau
      bool           setTau ( const double val ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral ( const double xlow , const double xhigh ,
                        const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_psy      ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY ()     ; }
      const Ostap::Math::Positive2D&   positive    () const { return m_positive ; }
      const Ostap::Math::Positive2D&   polynom     () const { return m_positive ; }
      // ====================================== ===============================
    private:
      // ======================================================================
      /// the actual (positive) bernstein polynomial in 2D
      Ostap::Math::Positive2D   m_positive ; // the actual bernstein polynomial
      /// Phase space
      Ostap::Math::PhaseSpaceNL m_psy      ;
      /// exponential
      double                    m_tau      ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace   ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Expo2DPol
     *  The 2D-function:
     *  \f$ f(x,y) = exp(x)*expo(y)*P_{pos}(x,y) \f$, where
     * \f$P_{pos}\f$ denotes the positive polynomial
     */
    class  Expo2DPol
      : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Expo2DPol ( const double          xmin = 0  ,
                  const double          xmax = 1  ,
                  const double          ymin = 0  ,
                  const double          ymax = 1  ,
                  const unsigned short  Nx   =  1 ,
                  const unsigned short  Ny   =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_positive.npars () ; }
      /// set k-parameter
      bool setPar       ( const unsigned int k , const double value )
      { return m_positive.setPar ( k , value ) ;}
      /// set k-parameter
      bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned int k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get nX & nY
      unsigned short nX () const { return m_positive.nX () ; }
      unsigned short nY () const { return m_positive.nY () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get tau
      double         tauX    () const { return m_tauX ;}
      double         tauY    () const { return m_tauY ;}
      /// set tau
      bool           setTauX ( const double val ) ;
      bool           setTauY ( const double val ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral ( const double xlow , const double xhigh ,
                        const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Positive2D&   positive    () const { return m_positive ; }
      const Ostap::Math::Positive2D&   polynom     () const { return m_positive ; }
      // ====================================== ===============================
    private:
      // ======================================================================
      /// the actual (positive) bernstein polynomial in 2D
      Ostap::Math::Positive2D   m_positive ; // the actual bernstein polynomial
      /// exponential
      double                    m_tauX     ;
      double                    m_tauY     ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Expo2DPolSym
     *  The 2D-function:
     *  \f$ f(x,y) = exp(x)*expo(y)*P_{sym}(x,y) \f$, where
     * \f$P_{pos}\f$ denotes the symmetric positive polynomial
     */
    class  Expo2DPolSym
      : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Expo2DPolSym ( const double          xmin = 0 ,
                     const double          xmax = 1 ,
                     const unsigned short  N    = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_positive.npars () ; }
      /// set k-parameter
      bool setPar       ( const unsigned int k , const double value )
      { return m_positive.setPar ( k , value ) ;}
      /// set k-parameter
      bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned int k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get nX & nY
      unsigned short n  () const { return m_positive.nX () ; }
      unsigned short nX () const { return m_positive.nX () ; }
      unsigned short nY () const { return m_positive.nY () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get tau
      double         tau     () const { return m_tau  ;}
      /// set tau
      bool           setTau  ( const double val ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral ( const double xlow , const double xhigh ,
                        const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public:  // expose some internmals
      // ======================================================================
      const Ostap::Math::Positive2DSym& positive () const { return m_positive ; }
      const Ostap::Math::Positive2DSym& polynom  () const { return m_positive ; }
      // ====================================== ===============================
    private:
      // ======================================================================
      /// the actual (positive) bernstein polynomial in 2D
      Ostap::Math::Positive2DSym m_positive ; // the actual bernstein polynomial
      /// exponential
      double                     m_tau      ;
      // ======================================================================
    };
    // ========================================================================
    /** @class GenSigmoid
     *  Sigmoid function, modulated by the positive polynomial
     *  \f$ f(x) = ( 1 + tanh( \alpha ( x  - x_0) ) \times P_{pos} (x) \f$
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  Sigmoid : public std::unary_function<double,double>
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
    class  TwoExpos : public std::unary_function<double,double>
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
    class  TwoExpoPositive : public std::unary_function<double,double>
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
    class  Tsallis : public std::unary_function<double,double>
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
    class  QGSM: public std::unary_function<double,double>
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
    // forward declarations
    class FourierSum ;
    class CosineSum  ;
    // ========================================================================
    /** @class FourierSum
     *  Fourier sum
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-07-26
     */
    class  FourierSum : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** @param degree  degree
       *  @param xmin    low  edge
       *  @param xmax    high edge
       *  @param fejer   use fejer summation
       */
      FourierSum ( const unsigned short degree = 0     ,   // degree
                   const double         xmin   = 0     ,   // low edge
                   const double         xmax   = 1     ,   // high edge
                   const bool           fejer  = false );  // use Fejer summation
      /// constructor from cosine serie
      FourierSum ( const CosineSum&  sum ) ;
      /// constructor from Fourier series and fejer flag
      FourierSum ( const FourierSum& sum  , const bool fejer ) ;
      // ======================================================================
    protected:  // protected constructor from parameters
      // ======================================================================
      /// protected constructor from parameters
      FourierSum ( const std::vector<double>& pars  ,
                   const double               xmin  ,
                   const double               xmax  ,
                   const double               fejer );
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const
      { return m_fejer ? fejer_sum ( x ) : fourier_sum ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate Fourier sum
      double fourier_sum ( const double x ) const ;
      /// calculate Fejer sum
      double fejer_sum   ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin  ; }
      /// get upper edge
      double xmax  () const { return m_xmax  ; }
      /// use Fejer summation?
      bool   fejer () const { return m_fejer ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const { return    t / m_scale   + m_delta ; }
      double t ( const double x ) const { return  ( x - m_delta ) * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// degree  of polynomial
      unsigned short degree () const { return ( m_pars.size() - 1 ) / 2 ; }
      /// number of parameters
      unsigned short npars  () const { return m_pars.size()     ; }
      /// all zero ?
      bool           zero   () const ;
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setPar          ( const unsigned short k , const double value ) ;
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setParameter    ( const unsigned short k , const double value )
      { return setPar      ( k , value ) ; }
      /// get the parameter value
      double  par          ( const unsigned short k ) const
      { return ( k < m_pars.size() ) ? m_pars[k] : 0.0 ; }
      /// get the parameter value
      double  parameter    ( const unsigned short k ) const { return par ( k ) ; }
      /// get all parameters:
      const std::vector<double>& pars () const { return m_pars ; }
      /// get k-th cos-parameter
      double a ( const unsigned short k ) const { return par ( 2 * k     ) ; }
      /// get k-th sin-parameter
      double b ( const unsigned short k ) const
      { return 1 <= k ? par   ( 2 * k - 1 ) : 0 ; }
      // set cosine terms
      bool setA ( const unsigned short k , const double value )
      { return setPar ( 2 * k , value ) ; }
      // set cosine terms
      bool setB ( const unsigned short k , const double value )
      { return 1<= k ? setPar ( 2 * k - 1 , value ) : false ; }
      /** get the magnitude of nth-harmonic
       *  \f$m_k = \sqrt( a^2_k + b^2_k) \f$
       */
      double mag    ( const unsigned short k ) const ;
      /// get the phase for n-th harmonic
      double phase  ( const unsigned short k ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get Fejer sum
      FourierSum fejer_sum   () const ;                       // get Fejer sum
      // ======================================================================
    public:
      // ======================================================================
      /// get the derivative at point x
      double     derivative   ( const double x ) const ;
      /// get the derivative as function
      FourierSum derivative   ( ) const ;
      /// get nth derivative as function
      FourierSum derivative_n ( const unsigned short n ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get integral between low and high
      double     integral   ( const double low , const double high ) const ;
      /** get integral as function
       *  @param c0  integration constant
       */
      FourierSum integral   ( const double c0 = 0 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** convolute with gaussian
       *  @param sigma resoltuion parameter for gaussian
       *  @return convolution witgh gaussian
       */
      FourierSum   convolve     ( const double sigma     ) const ;
      /** deconvolute with optional regularization
       *  @param sigma sigma of gaussian
       *  @param delta parameter of Tikhonov's regularization
       *  for delta<=0, no regularization
       *  @return regularised deconvolution
       */
      FourierSum deconvolve     ( const double sigma     ,
                                  const double delta = 0 ) const ;
      /**  get the effective cut-off (==number of effective harmonics)
       *   of Tikhonov's regularization
       *   \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
       *   @param sigma  gaussian resoltuion
       *   @param delta  regularization parameter
       *   @return number of effective harmonic
       */
      double     regularization ( const double sigma     ,
                                  const double delta     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: scale it!
      FourierSum& operator *= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: scale it!
      FourierSum& operator /= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: add constant
      FourierSum& operator += ( const double a ) ;     // add constant
      /// simple  manipulations with polynoms: subtract constant
      FourierSum& operator -= ( const double a ) ;     // subtract constant
      // ======================================================================
    public:
      // ======================================================================
      /** sum of two Fourier series (with the same interval!)
       *  @param other the first fourier sum
       *  @return the sum of two Fourier series
       */
      FourierSum sum ( const FourierSum& other ) const ;
      // ======================================================================
      /** get "shifted" fourier sum
       *  \f$ g(x) \equiv f ( x - a ) \f$
       *  @param a the bias aprameter
       *  @return the shifted fourier sum
       */
      FourierSum shift ( const double a ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// actual vector of coefficients
      std::vector<double> m_pars ; // actual vector of coefficients
      /// low edge
      double m_xmin  ;             // the low edge
      /// high edge
      double m_xmax  ;             // the high edge
      /// scale factor
      double m_scale ;             // scale factor
      /// delta
      double m_delta ;             // delta
      /// summation algorithm
      bool m_fejer   ;             // summation algorithm
       // ======================================================================
    } ;
    // ========================================================================
    /** @class CosineSum
     *  Fourier sum over cosines
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-07-26
     */
    class  CosineSum : public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** @param degree  degree
       *  @param xmin    low  edge
       *  @param xmax    high edge
       *  @param fejer   use fejer summation
       */
      CosineSum ( const unsigned short degree = 0     ,    // degree
                  const double         xmin   = 0     ,    // low edge
                  const double         xmax   = 1     ,    // high edge
                  const bool           fejer  = false ) ;  // use Fejer summation
      /// constructor from Fourier sum
      CosineSum ( const FourierSum&    sum            ) ;
      /// constructor from Fourier series and fejer flag
      CosineSum ( const CosineSum&     sum  , const bool fejer ) ;
      // ======================================================================
    protected:  // protected constructor from parameters
      // ======================================================================
      /// protected constructor from parameters
      CosineSum ( const std::vector<double>& pars  ,
                  const double               xmin  ,
                  const double               xmax  ,
                  const double               fejer );
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const
      { return m_fejer ? fejer_sum ( x ) : fourier_sum ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate Fourier sum
      double fourier_sum ( const double x ) const ;
      /// calculate Fejer sum
      double fejer_sum   ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin  ; }
      /// get upper edge
      double xmax  () const { return m_xmax  ; }
      /// use Fejer summation?
      bool   fejer () const { return m_fejer ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const { return    t / m_scale  + m_xmin  ; }
      double t ( const double x ) const { return  ( x - m_xmin ) * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// degree  of polynomial
      unsigned short degree () const { return m_pars.size() - 1 ; }
      /// number of parameters
      unsigned short npars  () const { return m_pars.size()     ; }
      /// all zero ?
      bool           zero   () const ;
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setPar          ( const unsigned short k , const double value ) ;
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setParameter    ( const unsigned short k , const double value )
      { return setPar      ( k , value ) ; }
      /// get the parameter value
      double  par          ( const unsigned short k ) const
      { return ( k < m_pars.size() ) ? m_pars[k] : 0.0 ; }
      /// get the parameter value
      double  parameter    ( const unsigned short k ) const { return par ( k ) ; }
      /// get all parameters:
      const std::vector<double>& pars () const { return m_pars ; }
      /// get k-th cos-parameter
      double a    ( const unsigned short k ) const { return par ( k     ) ; }
      // set cosine terms
      bool   setA ( const unsigned short k , const double value )
      { return setPar ( k , value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get Fejer sum
      CosineSum fejer_sum   () const ;                         // get Fejer sum
      // ======================================================================
    public:
      // ======================================================================
      /// get the derivative at point x
      double     derivative ( const double x ) const ;
      /// get the derivative as function
      FourierSum derivative ( ) const ;
      /// get nth derivative as function
      FourierSum derivative_n ( const unsigned short n ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get integral between low and high
      double     integral   ( const double low , const double high ) const ;
      /** get integral as function
       *  @param c0  integration constant
       */
      FourierSum integral   ( const double c0 = 0 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** convolute with gaussian
       *  @param sigma resoltuion parameter for gaussian
       *  @return convolution witgh gaussian
       */
      CosineSum   convolve     ( const double sigma     ) const ;
      /** deconvolute with optional regularization
       *  @param sigma sigma of gaussian
       *  @param delta parameter of Tikhonov's regularization
       *  for delta<=0, no regularization
       *  @return regularised deconvolution
       */
      CosineSum deconvolve     ( const double sigma     ,
                                 const double delta = 0 ) const ;
      /** get the effective cut-off (==number of terms/harmonics)
       *  of Tikhonov's regularization
       *  \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
       *  @param sigma  gaussian resoltuion
       *  @param delta  regularization parameter
       *  @return number of effective harmonic
       */
      double    regularization ( const double sigma     ,
                                 const double delta     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: scale it!
      CosineSum& operator *= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: scale it!
      CosineSum& operator /= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: add constant
      CosineSum& operator += ( const double a ) ;     // add constant
      /// simple  manipulations with polynoms: subtract constant
      CosineSum& operator -= ( const double a ) ;     // subtract constant
      // ======================================================================
    public:
      // ======================================================================
      /** sum of two Fourier series (with the same interval!)
       *  @param other the first fourier sum
       *  @return the sum of two Fourier series
       */
      CosineSum sum ( const CosineSum& other ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// actual vector of coefficients
      std::vector<double> m_pars ; // actual vector of coefficients
      /// low edge
      double m_xmin  ;             // the low edge
      /// high edge
      double m_xmax  ;             // the high edge
      /// scale factor
      double m_scale ;             // scale factor
      /// summation algorithm
      bool m_fejer   ;             // summation algorithm
      // ======================================================================
    } ;
    // ========================================================================
    /** make a sum of two fourier series (with the same interval!)
     *  @param s1 the first fourier sum
     *  @param s2 the first fourier sum
     *  @return s1+s2
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-26
     */
    
    FourierSum sum ( const FourierSum& s1 , const FourierSum& s2 ) ;
    // ========================================================================
    /** make a sum of two fourier cosine series (with the same interval!)
     *  @param s1 the first fourier cosine sum
     *  @param s2 the first fourier cosine sum
     *  @return s1+s2
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-26
     */
    CosineSum sum ( const CosineSum& s1 , const CosineSum& s2 ) ;
    // ========================================================================
    /// sum of two fourier series
    inline FourierSum operator+( const FourierSum& a ,
                                 const FourierSum& b ) { return a.sum ( b ) ; }
    /// sum of two cosine series
    inline CosineSum  operator+( const CosineSum&  a ,
                                 const CosineSum&  b ) { return a.sum ( b ) ; }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_FUNCTIONS_H
// ============================================================================
