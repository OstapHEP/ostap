// ============================================================================
#ifndef OSTAP_SPECTRA_H
#define OSTAP_SPECTRA_H 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
// ============================================================================
/** @file Ostap/Spectra.h
 *
 *  Set of simple phenomenology models to describe the  pt-spectra
 *  - Tsallis
 *  - QGSM 
 *  - Hagedorn 
 *
 *  @see Ostap::Math::Tsallis
 *  @see Ostap::Math::QGSM
 *  @see Ostap::Math::Hagedorn
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
      inline double operator() ( const double pt ) const { return evaluate ( pt ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get mass-parameter
      inline double mass () const { return m_mass  ; } // get mass-parameter
      /// get n-parameter
      inline double n    () const { return m_n     ; } // get n-parameter
      /// get T-parameter
      inline double T    () const { return m_T     ; } // get T-parameter
      // ======================================================================
      // aliases
      // ======================================================================
      /// get mass-parameter
      inline double m    () const { return  mass () ; } // get mass-parameter
      /// get mass-parameter
      inline double M    () const { return  mass () ; } // get mass-parameter
      /// get n-parameter
      inline double N    () const { return  n    () ; } // get n-parameter
      // ======================================================================
    public:
      // ======================================================================
      /// q-parameter for Tsallis enthropy 
      inline double q    () const { return m_n / ( m_n - 1 ) ; }
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
      inline double xmin() const { return 0 ; }
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
      inline double mass () const { return m_mass  ; } // get mass-parameter
      /// get n-parameter
      inline double b    () const { return m_b     ; } // get b-parameter
      // ======================================================================
      // aliases
      // ======================================================================
      /// get mass-parameter
      inline double m    () const { return  mass () ; } // get mass-parameter
      /// get mass-parameter
      inline double M    () const { return  mass () ; } // get mass-parameter
      /// get b-parameter
      inline double B    () const { return  b    () ; } // get n-parameter
      /// get b-parameter
      inline double B_0   () const { return  b    () ; } // get n-parameter
      /// get b-parameter
      inline double b0   () const { return  b    () ; } // get n-parameter
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
      inline double xmin() const { return 0 ; }
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
     *  simple function to describe pT spectra of particles 
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
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namesapce Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_SPECTRA_H
// ============================================================================
