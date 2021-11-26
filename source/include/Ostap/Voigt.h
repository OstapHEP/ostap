// ============================================================================
#ifndef OSTAP_VOIGT_H 
#define OSTAP_VOIGT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
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
     *  @see https://doi.org/10.1107/S0021889800010219
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
  } //                                        The end oif namespace Ostap::Math
  // ==========================================================================
} //                                                The end of namespoace Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_VOIGT_H
// ============================================================================
