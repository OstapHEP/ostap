// ============================================================================
#ifndef OSTAP_BREITWIGNER_H 
#define OSTAP_BREITWIGNER_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
#include "Ostap/PhaseSpace.h"
// ============================================================================
/** @file Ostap/BreitWigner.h
 *
 *  set of useful models for describing signal peaks with the natural width \
 *  - Breit-Wigner
 *  - Flatte 
 *  - LASS  (kappa) 
 *  - Bugg  (sigma-pole)
 *  - Gounaris-Sakurai
 *  - Voight &Co 
 *  - ...
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
    /// base class for formfactors
    class FormFactor ;
    // ========================================================================
    /** @namespace Ostap::Math::FormFactors
     *  form-factor functions for Breit-Wigner
     *   - Blatt-Weisskopf form-factors 
     *   - Jackson's form-factors 
     */
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
    } //                          the end of namespace Ostap::Math::FormFactors
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
      /// clone it 
      virtual BreitWigner* clone() const ;
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the Breit-Wigner shape
       *  \f$\frac{1}{\pi}\frac{\omega\Gamma(\omega)}
       *   { (\omega_0^2-\omega^2)^2-\omega_0^2\Gamma^2(\omega)-}\f$
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
      /// pole position 
      double         m0     () const { return m_m0      ; }
      /// pole position 
      double         mass   () const { return   m0   () ; }
      /// pole position 
      double         peak   () const { return   m0   () ; }
      /// the width at the pole 
      double         gam0   () const { return m_gam0    ; }
      /// the width at the pole 
      double         gamma0 () const { return   gam0 () ; }
      /// the width at the pole 
      double         gamma  () const { return   gam0 () ; }
      /// the width at the pole 
      double         width  () const { return   gam0 () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the mass of the first daughter 
      double         m1     () const { return m_m1 ; }
      /// the mass of the second daughter 
      double         m2     () const { return m_m2 ; }
      /// relative orbital momentum 
      unsigned short L      () const { return m_L  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set pole position 
      bool setM0     ( const double x ) ;
      /// set pole position 
      bool setMass   ( const double x ) { return setM0     ( x ) ; }
      /// set pole position 
      bool setPeak   ( const double x ) { return setM0     ( x ) ; }
      // set width at pole 
      bool setGamma0 ( const double x ) ;
      // set width at pole 
      bool setGamma  ( const double x ) { return setGamma0 ( x ) ; }
      // set width at pole 
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
    class Rho0FromEtaPrime : public Ostap::Math::Rho0
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
      /// clone it! 
      Rho0FromEtaPrime* clone() const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the function
      double operator() ( const double x ) const  override;
      // ======================================================================
    private:
      // ======================================================================
      /// the mass of eta'
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
    {
    public:
      // ======================================================================
      /** constructor  from all parameters
       *  \f$ f \rightarrow A_1 + A_2\f$
       *  @param m0    the mass
       *  @param m0g1  parameter \f$ m_0\times g_1\f$
       *  @param g2og1 parameter \f$ g2/g_1       \f$
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
      /// copy constructor 
      Flatte ( const Flatte&  right ) = default ;
      /// destructor
      virtual ~Flatte () ;
      // ======================================================================
      /// clone it! 
      virtual Flatte* clone() const ;
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
      /// pole position
      double m0     () const { return m_m0      ; }
      /// pole position
      double mass   () const { return   m0   () ; }
      /// pole position
      double peak   () const { return   m0   () ; }
      /// m*g1 
      double m0g1   () const { return m_m0g1    ; }
      /// g2/g1 
      double g2og1  () const { return m_g2og1   ; }
      /// mass of the first  daughter 
      double mA1    () const { return m_A1      ; }
      /// mass of the second daughter 
      double mA2    () const { return m_A2      ; }
      /// mass of the first  daughter in for the coupled channel
      double mB1    () const { return m_B1      ; }
      /// mass of the second daughter in for the coupled channel
      double mB2    () const { return m_B2      ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the thereshold 
      double thresholdA () const { return mA1() + mA2() ; }
      /// the threshold for the coupled channel 
      double thresholdB () const { return mB1() + mB2() ; }
      /// minimal threshold 
      double threshold  () const
      { return std::min ( thresholdA () , thresholdB () ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set pole position
      bool setM0     ( const double x ) ;
      /// set pole position
      bool setMass   ( const double x ) { return setM0 ( x ) ; }
      /// set pole position
      bool setPeak   ( const double x ) { return setM0 ( x ) ; }
      /// set m*g1 
      bool setM0G1   ( const double x ) ;
      /// set g2/g1 
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
       *  \f$ f \rightarrow B_1 + B_2\f$
       *  @param m0    the mass
       *  @param m0g1  parameter \f$ m_0\times g_1\f$
       *  @param g2og1 parameter \f$ g2/g_1       \f$
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
      Flatte2 ( const Flatte&  flatte ) ;
      /// copy constructor 
      Flatte2 ( const Flatte2& flatte ) = default ;
      /// destructor
      virtual ~Flatte2 () ;
      /// clone it!
      Flatte2* clone() const override ;
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
      /// pole position 
      bool setM0     ( const double x ) ;
      /// pole position 
      bool setMass   ( const double x ) { return setM0 ( x ) ; }
      /// pole position 
      bool setPeak   ( const double x ) { return setM0 ( x ) ; }
      /// width at the pole 
      bool setGamma  ( const double x ) ;
      /// width at the pole 
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
     *  Simplified version of Voigt profile
     *  @see T. Ida, M. Ando and H. Toraya,
     *       "Extended pseudo-Voigt function for approximating the Voigt profile"
     *       J. Appl. Cryst. (2000). 33, 1311-1316
     *  @see doi:10.1107/S0021889800010219
     *  @see http://dx.doi.org/10.1107/S0021889800010219
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-13
     */
    class  PseudoVoigt
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
      /// pole position 
      double m0     () const { return m_m0      ; }
      /// pole position 
      double mass   () const { return   m0   () ; }
      /// pole position 
      double peak   () const { return   m0   () ; }
      /// width at the pole 
      double gamma  () const { return m_gamma   ; }
      /// resolution  
      double sigma  () const { return m_sigma   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set pole position
      bool setM0     ( const double x ) ;
      /// set pole position
      bool setMass   ( const double x ) { return setM0 ( x ) ; }
      /// set pole position
      bool setPeak   ( const double x ) { return setM0 ( x ) ; }
      /// set width at the pole 
      bool setGamma  ( const double x ) ; 
      /// set resolution 
      bool setSigma  ( const double x ) ;
      // ======================================================================
    public: // helper constants
      // ======================================================================
      /// FWHM-gaussian
      double fwhm_gauss      () const ;
      /// FWHM-lorenzian
      double fwhm_lorentzian () const { return 2 * m_gamma ; }
      /// rho 
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
    /** @namespace Ostap::Math::Jackson
     *   - Jackson's form-factors 
     */
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
       *  \f$\rho(\omega)= \frac{ ( \omega + M )^2 - m^2 }{ \omega^2} \f$
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
       *  \f$\rho(\omega)= \left[ ( \omega + M )^2 - m^2 \right]^{-1} \f$
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
      /** the simple function for \f$\rho^0 \rightarrow \pi^+ \pi^-\f$ and           
       *  \f$ 1- \rightarrow 0^- 0^- \f$, l = 1
       *  \f$ \rho(\omega)= \left[ q_0^2 + q^2 \right]^{-1} \f$
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
        /// orbital momentum
        Case   m_L ; // orbital momentum
        /// Break-up 
        double m_b ; // Break-up 
        // ====================================================================
      } ;
      // ======================================================================
    } // end of namespace Ostap:Math::FormFactors
    // ========================================================================
    // VARIOUS BEASTS 
    // ========================================================================
    /** @class LASS
     *  The LASS parameterization (Nucl. Phys. B296, 493 (1988))
     *  describes the 0+ component of the Kpi spectrum.
     *  It consists of the K*(1430) resonance together with an
     *  effective range non-resonant component
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-10-05
     */
    class  LASS
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m0 the mass of K* 
       *  @param g0 the width of  K*
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
      /// pole position
      double m0  ( ) const { return m_m0 ; }
      /// width at the pole 
      double g0  ( ) const { return m_g0 ; }
      /// a 
      double a   ( ) const { return m_a  ; }
      /// r 
      double r   ( ) const { return m_r  ; }
      /// elasticity 
      double e   ( ) const { return m_e  ; }
      // ======================================================================
      /// the mass of the first daughter 
      double m1  ( ) const { return m_ps2.m1 () ; }
      /// the mass of the second daughter 
      double m2  ( ) const { return m_ps2.m2 () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set pole posiiton 
      bool setM0 ( const double value ) ;
      /// set with 
      bool setG0 ( const double value ) ;
      /// a 
      bool setA  ( const double value ) ;
      /// r
      bool setR  ( const double value ) ;
      /// elatisicity 
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
    {
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param M  mass of sigma (very different from the pole positon!)
       *  @param g2 width parameter g2 (4pi width)
       *  @param b1 width parameter b1  (2pi coupling)
       *  @param b2 width parameter b2  (2pi coupling)
       *  @param a  parameter a (the exponential cut-off)
       *  @param s1 width parameter s1  (cut-off for 4pi coupling)
       *  @param s2 width parameter s2  (cut-off for 4pi coupling)
       *  @param m1 the mass of the first  particle
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
      /// calculate the Bugg shape      
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
      /// pole position  
      double M     () const  { return m_M       ; }
      /// m^2
      double M2    () const  { return m_M * m_M ; }
      /// pole positon  
      double m0    () const  { return   M ()    ; }
      /// pole positon      
      double mass  () const  { return   M ()    ; }
      /// pole positon
      double peak  () const  { return   M ()    ; }
      // ======================================================================
      /// g2 
      double g2    () const  { return m_g2   ; }
      /// b1 
      double b1    () const  { return m_b1   ; }
      /// b2 
      double b2    () const  { return m_b2   ; }
      /// s1 
      double s1    () const  { return m_s1   ; }
      /// s2 
      double s2    () const  { return m_s2   ; }
      /// a 
      double a     () const  { return m_a    ; }
      // ======================================================================
      /// set pole position
      bool setM    ( const double value  ) ;
      /// set pole position
      bool setM0   ( const double value  ) { return setM ( value )  ; }
      /// set pole position
      bool setMass ( const double value  ) { return setM ( value )  ; }
      /// set pole position
      bool setPeak ( const double value  ) { return setM ( value )  ; }
      // ======================================================================
      /// g2
      bool setG2   ( const double value  ) ;
      /// b1 
      bool setB1   ( const double value  ) ;
      /// b2 
      bool setB2   ( const double value  ) ;
      /// s1 
      bool setS1   ( const double value  ) ;
      /// s2 
      bool setS2   ( const double value  ) ;
      /// a
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
    /** @class Swanson
     *  Swanson's parameterization of S-wave cusp
     *  @see LHCb-PAPER-2016-019 appendix D
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-11
     */
    class  Swanson
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
    

    // ========================================================================
    // "2-from-3" variants 
    // ========================================================================

    // ========================================================================
    /** @class BW23L
     *  @see Ostap::Math::BreitWigner
     *  @see Ostap::Math::PhaseSpace23L
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2012-05-23
     */
    class  BW23L
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
      /// copy 
      BW23L ( const Ostap::Math::BW23L& bw  ) ;
      /// destructor
      virtual ~BW23L () ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the shape
      double operator() ( const double x ) const ;
      /// get the amplitude
      std::complex<double>
      amplitude ( const double x ) const { return m_bw->amplitude ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// pole position 
      double m0     () const { return m_bw-> m0   () ; }
      /// pole position 
      double mass   () const { return        m0   () ; }
      /// pole position 
      double peak   () const { return        m0   () ; }
      /// width  at the pole 
      double gam0   () const { return m_bw-> gam0 () ; }
      /// width  at the pole 
      double gamma0 () const { return        gam0 () ; }
      /// width  at the pole 
      double gamma  () const { return        gam0 () ; }
      /// width  at the pole 
      double width  () const { return        gam0 () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set pole position 
      bool setM0     ( const double x ) { return m_bw->setM0     ( x ) ; }
      /// set pole position 
      bool setMass   ( const double x ) { return setM0           ( x ) ; }
      /// set width
      bool setPeak   ( const double x ) { return setM0           ( x ) ; }
      /// set width
      bool setGamma0 ( const double x ) { return m_bw->setGamma0 ( x ) ; }
      /// set width
      bool setGamma  ( const double x ) { return setGamma0       ( x ) ; }
      /// set width
      bool setWidth  ( const double x ) { return setGamma0       ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// low  edge of phasespace 
      double lowEdge   () const { return m_ps. lowEdge() ; }
      /// high edge of phasespace 
      double highEdge  () const { return m_ps.highEdge() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the current width
      double gamma ( const double x ) const { return m_bw->gamma ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the Breit-Wigner
      const Ostap::Math::BreitWigner& breitwigner() const { return *m_bw ; }
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
      std::unique_ptr<Ostap::Math::BreitWigner> m_bw;    // the breit wigner
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
    {
    public:
      // ======================================================================
      /** constructor  from all parameters
       *  @param m0    the mass
       *  @param m0g1  parameter \f$ m_0\times g_1\f$
       *  @param g2og1 parameter \f$ g2/g_1       \f$
       *  @param mK    A mass
       *  @param mPi    B mass
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
       *  @param fun preconfigured flatte function 
       *  @param m3    the mass of the third particle
       *  @param m     the mass of mother particle
       *  @param L     the orbital momentum between the pair and the third particle
       */
      Flatte23L  ( const Flatte&        fun              ,     // MeV
                   const double         m3    = 3096.9   ,     // MeV
                   const double         m     = 5366.0   ,     // MeV
                   const unsigned short L     = 1        ) ;
      /// copy 
      Flatte23L  ( const Flatte23L& right ) ;
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
      { return m_flatte->flatte_amp ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// pole position
      double m0     () const { return m_flatte-> m0    () ; }
      /// pole position
      double mass   () const { return            m0    () ; }
      /// pole position
      double peak   () const { return            m0    () ; }
      /// m*g1 
      double m0g1   () const { return m_flatte-> m0g1  () ; }
      /// g2/g1 
      double g2og1  () const { return m_flatte-> g2og1 () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// ste pole   position 
      bool setM0     ( const double x ) { return m_flatte-> setM0    ( x ) ; }
      /// ste pole   position 
      bool setMass   ( const double x ) { return            setM0    ( x ) ; }
      /// ste pole   position 
      bool setPeak   ( const double x ) { return            setM0    ( x ) ; }
      /// m*g1 
      bool setM0G1   ( const double x ) { return m_flatte-> setM0G1  ( x ) ; }
      /// g2/g1
      bool setG2oG1  ( const double x ) { return m_flatte-> setG2oG1 ( x ) ; }
      // ======================================================================      
    public:
      // ======================================================================
      /// low  edge of phase space 
      double lowEdge   () const { return m_ps .  lowEdge () ; }
      /// high edge of phase space 
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
    public:
      // ======================================================================
      /// get the Flatte
      const Ostap::Math::Flatte& flatte () const { return *m_flatte ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the actual Flatte function
      std::unique_ptr<Ostap::Math::Flatte> m_flatte ; // the actual Flatte function
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
    class  LASS23L 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param m0 the mass of K*
       *  @param g0 the width of K* 
       *  @param L  the angular momentum between the first pair and the third
       *  @param a  the LASS parameter
       *  @param r  the LASS parameter
       *  @param e  the LASS parameter
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
      /// pole position 
      double m0  () const { return m_lass . m0 () ; } // K*(1430) mass
      /// width at the pole 
      double g0  () const { return m_lass . g0 () ; } // K*(1430) width
      /// a 
      double a   () const { return m_lass . a  () ; }
      /// r 
      double r   () const { return m_lass . r  () ; }
      ///  e 
      double e   () const { return m_lass . e  () ; }
      // ======================================================================
      /// mass of the first daughter 
      double m1  () const { return m_ps   . m1 () ; }
      /// mass of the second daughter 
      double m2  () const { return m_ps   . m2 () ; }
      /// the tthird mass 
      double m3  () const { return m_ps   . m3 () ; }
      /// m 
      double m   () const { return m_ps   . m  () ; }
      /// l 
      double l   () const { return m_ps   . l  () ; }
      /// L 
      double L   () const { return m_ps   . L  () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set pole position 
      bool setM0 ( const double value ) { return m_lass . setM0 ( value ) ; }
      /// set g0 
      bool setG0 ( const double value ) { return m_lass . setG0 ( value ) ; }
      /// a 
      bool setA  ( const double value ) { return m_lass . setA  ( value ) ; }
      /// r 
      bool setR  ( const double value ) { return m_lass . setR  ( value ) ; }
      ///  e 
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
    {
    public:
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param M  mass of sigma (very different from the pole positon!)
       *  @param g2 width parameter g2 (4pi width)
       *  @param b1 width parameter b1  (2pi coupling)
       *  @param b2 width parameter b2  (2pi coupling)
       *  @param a  parameter a (the exponential cut-off)
       *  @param s1 width parameter s1  (cut-off for 4pi coupling)
       *  @param s2 width parameter s2  (cut-off for 4pi coupling)
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
       *  @param bugg precofoigured bugg-function
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
      /// calculate the Bugg shape
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
      /// the first mass 
      double m1        () const { return m_ps.m1 () ; }
      /// the second mass 
      double m2        () const { return m_ps.m2 () ; }
      /// the third mass 
      double m3        () const { return m_ps.m3 () ; }
      /// mass 
      double m         () const { return m_ps.m  () ; }
      // ======================================================================
      /// low   edge of phase space 
      double lowEdge   () const { return m_ps. lowEdge() ; }
      /// high edge of phase space 
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
      /// M
      double M     () const  { return m_bugg. M    () ; }
      /// M2
      double M2    () const  { return m_bugg. M2   () ; }
      /// pole postition
      double m0    () const  { return m_bugg. m0   () ; }
      /// pole position 
      double mass  () const  { return m_bugg. mass () ; }
      /// pole position 
      double peak  () const  { return m_bugg. peak () ; }
      // ======================================================================
      /// g2 
      double g2    () const  { return m_bugg. g2   () ; }
      /// b1 
      double b1    () const  { return m_bugg. b1   () ; }
      ///  b2 
      double b2    () const  { return m_bugg. b2   () ; }
      ///  s1 
      double s1    () const  { return m_bugg. s1   () ; }
      /// s2 
      double s2    () const  { return m_bugg. s2   () ; }
      /// a 
      double a     () const  { return m_bugg. a    () ; }
      // ======================================================================
      /// set pole position
      bool setM    ( const double value  ) { return m_bugg.setM    ( value ) ; }
      /// set pole position
      bool setM0   ( const double value  ) { return m_bugg.setM0   ( value ) ; }
      /// set pole position
      bool setMass ( const double value  ) { return m_bugg.setMass ( value ) ; }
      /// set pole position
      bool setPeak ( const double value  ) { return m_bugg.setPeak ( value ) ; }
      // ======================================================================
      /// g2 
      bool setG2   ( const double value  ) { return m_bugg.setG2   ( value ) ; }
      /// b1 
      bool setB1   ( const double value  ) { return m_bugg.setB1   ( value ) ; }
      ///  b2 
      bool setB2   ( const double value  ) { return m_bugg.setB2   ( value ) ; }
      /// s1 
      bool setS1   ( const double value  ) { return m_bugg.setS1   ( value ) ; }
      /// s2 
      bool setS2   ( const double value  ) { return m_bugg.setS2   ( value ) ; }
      /// a 
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
      /// m1 
      double m1        () const { return m_ps.m1 () ; }
      /// m2 
      double m2        () const { return m_ps.m2 () ; }
      /// m3 
      double m3        () const { return m_ps.m3 () ; }
      /// mass 
      double m         () const { return m_ps.m  () ; }
      // ======================================================================
      /// low edge of the phasespace 
      double lowEdge   () const { return m_ps. lowEdge() ; }
      /// high dge of the phasespace 
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
      /// pole position
      double M      () const  { return m_M     ; }
      /// pole position
      double m0     () const  { return   M  () ; }
      /// pole position
      double mass   () const  { return   M  () ; }
      /// pole position
      double peak   () const  { return   M  () ; }
      // ======================================================================
      /// width at the pole
      double g0     () const  { return m_g0    ; }
      /// width at the pole
      double gamma  () const  { return   g0 () ; }
      /// width at the pole
      double width  () const  { return   g0 () ; }
      // ======================================================================
      /// set pole position
      bool setM     ( const double value  ) ;
      /// set pole position
      bool setM0    ( const double value  ) { return setM  ( value ) ; }
      /// set pole position
      bool setMass  ( const double value  ) { return setM  ( value ) ; }
      /// set pole position
      bool setPeak  ( const double value  ) { return setM  ( value ) ; }
      // ======================================================================
      /// set width
      bool setG0    ( const double value  ) ;
      /// set width
      bool setGamma ( const double value  ) { return setG0 ( value ) ; }
      /// set width
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
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_BREITWIGNER_H
// ============================================================================
