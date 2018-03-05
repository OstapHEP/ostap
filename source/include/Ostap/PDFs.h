// ============================================================================
#ifndef OSTAP_PDFS_H
#define OSTAP_PDFS_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Peaks.h"
#include "Ostap/BreitWigner.h"
#include "Ostap/Models.h"
#include "Ostap/BSpline.h"
// ============================================================================
// ROOT
// ============================================================================
using std::size_t ;
// ============================================================================
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsReal.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @namespace  Models PDFs.h Ostap/PDFs.h
   *
   *  Naturally "wide" models:
   *
   *  - BreitWigner, Rho0, Kstar, Phi, ...
   *  - BreitWigner from 3-body decay of mother particle: BW23L
   *  - LASS (kappa pole)
   *  - LASS from 3-body decay of mother particle: LASS23L
   *  - Bugg (sigma pole)
   *  - Bugg from 3-body decay of mother particle: Bugg23L
   *  - Voigt
   *  - Swanson's S-wave cusp
   *
   *  Empirical resolution models:
   *
   *  - Crystal Ball
   *  - right side Crystal Ball
   *  - double-sided Crystal Ball
   *  - Needham: Crystal Ball with \f$\alpha(\sigma)\f$
   *  - Apolonios
   *  - Apolonios2 (bifurcated apolonious)
   *  - Bifurcated Gauissian
   *  - Generalized Gaussian v1
   *  - Generalized Gaussian v2
   *  - Skew Gaussian
   *  - Bukin
   *  - Student-T
   *  - bifurcated Student-T
   *  - Gram-Charlier-A
   *
   *  Smooth phase-space induced models for background
   *
   *  - 2-body phase space
   *  - L-body phase space at low-edge
   *  - L-body phase space in N-body decays at high-edge
   *  - L-body phase space from N-body decay
   *  - L-body phase space from N-body decay times the positive polynomial
   *  - 2-body phase space from 3-body decays taking into accout orbital momenta
   *
   *  Various smooth empirical models for background
   *
   *  - positive polynomial
   *  - exponential times positive polynomial
   *  - gamma distribution
   *  - generalized gamma distribution
   *  - Amoroso function
   *  - log  (Gamma-distribution)
   *  - log10(Gamma-distribution)
   *  - Log-Gamma distribution
   *  - Beta-prime distribution
   *
   *  Non-factorazeable smooth 2D-models
   *
   *  - generic   positive non-factorizable polynomial in 2D
   *   \f$ P^+(x,y) = \sum_i \sum_j \alpha^2_{i,j} B^n_i(x) B^k_j(y) \f$
   *  - symmetric positive non-factorizable polynomial in 2D \f$ P^+_{sym}(x,y) \f$
   *  - \f$ f(x,y)       = \Phi_1(x)\times\Phi_2(y)\timesP^+(x,y)       \f$
   *  - \f$ f_{sym}(x,y) = \Phi  (x)\times\Phi  (y)\timesP^+_{sym}(x,y) \f$
   *  - \f$ f(x,y)       = exp   (x)\times\Phi  (y)\timesP^+(x,y)       \f$
   *  - \f$ f(x,y)       = exp   (x)\times exp  (y)\timesP^+(x,y)       \f$
   *  - \f$ f_{sym}(x,y) = exp   (x)\times exp  (y)\timesP^+_{sym}(x,y) \f$
   *
   *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
   *  @date   2011-11-30
   */
  namespace Models
  {
    // ========================================================================
    // Naturally "wide" models
    // ========================================================================

    // ========================================================================
    /** @class BreitWigner
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *
     *  http://www.springerlink.com/content/q773737260425652/
     *
     *  @see Ostap::Math::BreitWigner
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  BreitWigner : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::BreitWigner, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BreitWigner ( const char*          name      ,
                    const char*          title     ,
                    RooAbsReal&          x         ,
                    RooAbsReal&          mass      ,
                    RooAbsReal&          width     ,
                    const double         m1        ,
                    const double         m2        ,
                    const unsigned short L     = 0 ) ;
      /// constructor from all parameters
      BreitWigner ( const char*          name      ,
                    const char*          title     ,
                    RooAbsReal&          x         ,
                    RooAbsReal&          mass      ,
                    RooAbsReal&          width     ,
                    const double         m1        ,
                    const double         m2        ,
                    const unsigned short L                         ,
                    const Ostap::Math::FormFactors::JacksonRho rho ) ;
      /// constructor from main parameters and "shape"
      BreitWigner ( const char*          name      ,
                    const char*          title     ,
                    RooAbsReal&          x         ,
                    RooAbsReal&          mass      ,
                    RooAbsReal&          width     ,
                    const Ostap::Math::BreitWigner& bw ) ;
      /// "copy" constructor
      BreitWigner ( const BreitWigner& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~BreitWigner() ;
      /// clone
      BreitWigner* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      BreitWigner() {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      std::complex<double>            amplitude () const  ;
      /// access to underlying function
      const Ostap::Math::BreitWigner& function  () const { return m_bw ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mass  ;
      RooRealProxy m_width ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::BreitWigner m_bw ;            // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Rho0
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *
     *  http://www.springerlink.com/content/q773737260425652/
     *
     *  @see Ostap::Models::BreitWigner
     *  @see Ostap::Math::BreitWigner
     *  @see Ostap::Math::Rho0
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Rho0 : public Ostap::Models::BreitWigner
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Rho0, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Rho0 ( const char*          name      ,
             const char*          title     ,
             RooAbsReal&          x         ,
             RooAbsReal&          mass      ,
             RooAbsReal&          width     ,
             const double         pi_mass   ) ;
      /// "copy" constructor
      Rho0 ( const Rho0& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~Rho0 () ;
      /// clone
      Rho0* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Rho0 () {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Kstar
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *
     *  http://www.springerlink.com/content/q773737260425652/
     *
     *  @see Ostap::Models::BreitWigner
     *  @see Ostap::Math::BreitWigner
     *  @see Ostap::Math::Kstar
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Kstar : public Ostap::Models::BreitWigner
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Kstar, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Kstar ( const char*          name      ,
              const char*          title     ,
              RooAbsReal&          x         ,
              RooAbsReal&          mass      ,
              RooAbsReal&          width     ,
              const double         k_mass    ,
              const double         pi_mass   ) ;
      /// "copy" constructor
      Kstar ( const Kstar& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~Kstar () ;
      /// clone
      Kstar* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Kstar () {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Phi
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *
     *  http://www.springerlink.com/content/q773737260425652/
     *
     *  @see Ostap::Models::BreitWigner
     *  @see Ostap::Math::BreitWigner
     *  @see Ostap::Math::Rho0
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Phi : public Ostap::Models::BreitWigner
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Phi, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Phi  ( const char*          name      ,
             const char*          title     ,
             RooAbsReal&          x         ,
             RooAbsReal&          mass      ,
             RooAbsReal&          width     ,
             const double         k_mass    ) ;
      /// "copy" constructor
      Phi ( const Phi& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~Phi () ;
      /// clone
      Phi* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Phi () {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BW23L
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *
     *  http://www.springerlink.com/content/q773737260425652/
     *
     *  @see Ostap::Math::BreitWigner
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  BW23L : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::BW23L, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BW23L ( const char*          name      ,
              const char*          title     ,
              RooAbsReal&          x         ,
              RooAbsReal&          mass      ,
              RooAbsReal&          width     ,
              const double         m1        ,
              const double         m2        ,
              const unsigned short l         ,
              //
              const double         m3        ,
              const double         m         ,
              const double         L         ) ;
      /// constructor from all parameters
      BW23L ( const char*          name      ,
              const char*          title     ,
              RooAbsReal&          x         ,
              RooAbsReal&          mass      ,
              RooAbsReal&          width     ,
              const double         m1        ,
              const double         m2        ,
              const unsigned short l                         ,
              const Ostap::Math::FormFactors::JacksonRho rho ,
              //
              const double         m3        ,
              const double         m         ,
              const double         L         ) ;
      /// constructor from main parameters and "shape"
      BW23L ( const char*          name      ,
              const char*          title     ,
              RooAbsReal&          x         ,
              RooAbsReal&          mass      ,
              RooAbsReal&          width     ,
              const Ostap::Math::BW23L& bw   ) ; // shape
      /// "copy" constructor
      BW23L ( const BW23L& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~BW23L() ;
      /// clone
      BW23L* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      BW23L () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t     evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      std::complex<double>      amplitude () const  ;
      /// access to underlying function
      const Ostap::Math::BW23L& function  () const { return m_bw ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mass  ;
      RooRealProxy m_width ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::BW23L m_bw ;            // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Flatte
     *
     *  S.M.Flatte,
     *  "Coupled-channel analysis of the \f$\pi\eta\f$
     *  and \f$K\bar{K}\f$ systems near \f$K\bar{K}\f$ threshold
     *  Phys. Lett. B63, 224 (1976
     *
     *  http://www.sciencedirect.com/science/article/pii/0370269376906547
     *
     *  \f$\pi\pi\f$-channel
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Flatte : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Flatte, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Flatte ( const char*                name      ,
               const char*                title     ,
               RooAbsReal&                x         ,
               RooAbsReal&                m0        ,
               RooAbsReal&                m0g1      ,
               RooAbsReal&                g2og1     ,
               const Ostap::Math::Flatte& flatte    ) ;
      /// "copy" constructor
      Flatte ( const Flatte& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~Flatte () ;
      /// clone
      Flatte* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Flatte () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      virtual std::complex<double> amplitude () const  ;
      /// access to underlying function
      const Ostap::Math::Flatte&   function  () const { return m_flatte ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_m0    ;
      RooRealProxy m_m0g1  ;
      RooRealProxy m_g2og1 ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Flatte m_flatte ;             // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Flatte2
     *
     *  S.M.Flatte,
     *  "Coupled-channel analysis of the \f$\pi\eta\f$
     *  and \f$K\bar{K}\f$ systems near \f$K\bar{K}\f$ threshold
     *  Phys. Lett. B63, 224 (1976
     *
     *  http://www.sciencedirect.com/science/article/pii/0370269376906547
     *
     *  \f$\pi\pi\f$-channel
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Flatte2 : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Flatte2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Flatte2 ( const char*                name      ,
                const char*                title     ,
                RooAbsReal&                x         ,
                RooAbsReal&                m0        ,
                RooAbsReal&                m0g1      ,
                RooAbsReal&                g2og1     ,
                const Ostap::Math::Flatte& flatte    ) ;
      /// "copy" constructor
      Flatte2 ( const Flatte2& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~Flatte2 () ;
      /// clone
      Flatte2* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Flatte2 () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      std::complex<double>          amplitude () const  ;
      /// access to underlying function
      const Ostap::Math::Flatte2&   function  () const { return m_flatte2 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_m0    ;
      RooRealProxy m_m0g1  ;
      RooRealProxy m_g2og1 ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Flatte2 m_flatte2 ;          // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LASS
     *  S-wave Kpi amplitude for S-wave Kpi distribtion
     *  @see Ostap::Math::LASS
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-10-05
     */
    class  LASS : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::LASS, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      LASS  ( const char*          name          ,
              const char*          title         ,
              RooAbsReal&          x             ,
              RooAbsReal&          m1430         , // mass  of K*(1430)
              RooAbsReal&          g1430         , // width of K*(1430)
              RooAbsReal&          a             ,
              RooAbsReal&          r             ,
              RooAbsReal&          e             ,
              const double         m1    = 493.7 ,   // mass of K
              const double         m2    = 139.6 ) ; // mass of pi
      /// "copy constructor"
      LASS  ( const LASS& right , const char* name = 0 )  ;
      /// destructor
      virtual ~LASS () ;
      /// clone
      LASS * clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      LASS () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      std::complex<double>     amplitude () const ;
      /// access to underlying function
      const Ostap::Math::LASS& function  () const { return m_lass ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// the mass
      RooRealProxy m_x     ;
      /// K*(1430) parameters
      RooRealProxy m_m0    ;
      RooRealProxy m_g0    ;
      /// LASS parameters
      RooRealProxy m_a     ;
      RooRealProxy m_r     ;
      RooRealProxy m_e     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::LASS m_lass ;              // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LASS23L
     *  S-wave Kpi amplitude for Kpi from B-> Kpi X decays
     *  @see Ostap::Math::LASS23L
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-02
     */
    class  LASS23L : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::LASS23L, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      LASS23L ( const char*          name          ,
                const char*          title         ,
                RooAbsReal&          x             ,
                RooAbsReal&          m1430         ,
                RooAbsReal&          g1430         ,
                RooAbsReal&          a             ,
                RooAbsReal&          r             ,
                RooAbsReal&          e             ,
                const double         m1    = 493.7 ,
                const double         m2    = 139.6 ,
                const double         m3    = 3097  ,
                const double         m     = 5278  ,
                const unsigned short L     = 1     ) ;
      /// "copy constructor"
      LASS23L ( const LASS23L& right , const char* name = 0 )  ;
      /// destructor
      virtual ~LASS23L() ;
      /// clone
      LASS23L* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      LASS23L () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// get the complex amplitude
      std::complex<double>        amplitude () const ; // get the complex amplitude
      /// access to underlying function
      const Ostap::Math::LASS23L& function  () const { return m_lass ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// the mass
      RooRealProxy m_x     ;
      /// K*(1430) parameters:
      RooRealProxy m_m0    ;
      RooRealProxy m_g0    ;
      /// LASS parameters
      RooRealProxy m_a     ;
      RooRealProxy m_r     ;
      RooRealProxy m_e     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::LASS23L m_lass ;              // the actual function
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
     *  @see Ostap::Math::Bugg
     */
    class  Bugg : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Bugg, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Bugg  ( const char*          name               ,
              const char*          title              ,
              RooAbsReal&          x                  ,
              RooAbsReal&          M                  ,   // sigma M
              RooAbsReal&          g2                 ,   // sigma G2
              RooAbsReal&          b1                 ,   // sigma B1
              RooAbsReal&          b2                 ,   // sigma B2
              RooAbsReal&          a                  ,   // sigma a
              RooAbsReal&          s1                 ,   // sigma s1
              RooAbsReal&          s2                 ,   // sigma s2
              const double         m1    = 139.6/1000 ) ; // mass of pi GeV
      /// "copy constructor"
      Bugg  ( const Bugg& right , const char* name = 0 )  ;
      /// destructor
      virtual ~Bugg () ;
      /// clone
      Bugg* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Bugg () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      std::complex<double>     amplitude () const ;
      /// access to underlying function
      const Ostap::Math::Bugg& function  () const { return m_bugg ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// the mass
      RooRealProxy m_x     ;
      /// sigma/bugg parameters
      RooRealProxy m_M     ;
      RooRealProxy m_g2    ;
      RooRealProxy m_b1    ;
      RooRealProxy m_b2    ;
      RooRealProxy m_a     ;
      RooRealProxy m_s1    ;
      RooRealProxy m_s2    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Bugg m_bugg ;              // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bugg23L
     *  parametrisation of sigma-pole for
     *  two pion mass distribution form three body decays
     *
     *  The parameterization of sigma pole by
     *  B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-01
     *  @see Ostap::Math::Bugg23L
     */
    class  Bugg23L : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Bugg23L, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Bugg23L ( const char*          name               ,
                const char*          title              ,
                RooAbsReal&          x                  ,
                RooAbsReal&          M                  ,   // sigma M
                RooAbsReal&          g2                 ,   // sigma G2
                RooAbsReal&          b1                 ,   // sigma B1
                RooAbsReal&          b2                 ,   // sigma B2
                RooAbsReal&          a                  ,   // sigma a
                RooAbsReal&          s1                 ,   // sigma s1
                RooAbsReal&          s2                 ,   // sigma s2
                const double         m1    = 139.6/1000 ,   // mass of pi GeV
                const double         m3 = 3097.0 / 1000 ,   //  GeV
                const double         m  = 5278.0 / 1000 ,   // GeV
                const unsigned short L  =    1          ) ;
      /// "copy constructor"
      Bugg23L ( const Bugg23L& right , const char* name = 0 )  ;
      /// destructor
      virtual ~Bugg23L () ;
      /// clone
      Bugg23L* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Bugg23L () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      std::complex<double>        amplitude () const ;
      /// access to underlying function
      const Ostap::Math::Bugg23L& function  () const { return m_bugg ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// the mass
      RooRealProxy m_x     ;
      /// sigma/bugg parameters
      RooRealProxy m_M     ;
      RooRealProxy m_g2    ;
      RooRealProxy m_b1    ;
      RooRealProxy m_b2    ;
      RooRealProxy m_a     ;
      RooRealProxy m_s1    ;
      RooRealProxy m_s2    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Bugg23L m_bugg ;              // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Voigt
     *  "Voigt"-function
     *  @see Ostap::Math::Voigt
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-12-05
     */
    class  Voigt : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Voigt, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Voigt
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          gamma     ,
        RooAbsReal&          sigma     ) ;
      /// "copy" constructor
      Voigt ( const Voigt& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Voigt () ;
      /// clone
      Voigt* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Voigt () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Voigt& function() const { return m_voigt ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_gamma  ;
      RooRealProxy m_sigma  ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Voigt m_voigt ;                      // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class PseudoVoigt
     *  "PseudoVoigt"-function
     *  @see Ostap::Math::PseudoVoigt
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2016-06-14
     */
    class  PseudoVoigt : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PseudoVoigt, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PseudoVoigt
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          gamma     ,
        RooAbsReal&          sigma     ) ;
      /// "copy" constructor
      PseudoVoigt ( const PseudoVoigt& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~PseudoVoigt () ;
      /// clone
      PseudoVoigt* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PseudoVoigt () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PseudoVoigt& function() const { return m_voigt ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_gamma  ;
      RooRealProxy m_sigma  ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PseudoVoigt m_voigt ;                // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class Swanson
     *  Swanson's S-wau cusp
     *  @see Ostap::Math::Swanson
     *  @see LHCb-PAPER-2016-019, Appendix D
     *  @see E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952,
     *  @see http://arxiv.org/abs/1504.07952
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-12
     */
    class  Swanson : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Swanson, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Swanson
        ( const char*          name          ,
          const char*          title         ,
          RooAbsReal&          x             ,
          RooAbsReal&          beta0         ,
          const Ostap::Math::Swanson& sw     ) ;
      /// constructor from all parameters
      Swanson
        ( const char*          name          ,
          const char*          title         ,
          RooAbsReal&          x             ,
          RooAbsReal&          beta0         ,
          const double         m1_0          ,
          const double         m2_0          ,
          const Ostap::Math::BreitWigner& bw ) ;
      /// "copy" constructor
      Swanson ( const Swanson& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Swanson () ;
      /// clone
      Swanson* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Swanson () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Swanson& function() const { return m_swanson ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_beta0  ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Swanson m_swanson ;                 // the function
      // ======================================================================
    };
    // ========================================================================


    // ========================================================================
    // Resolution models
    // ========================================================================

    // ========================================================================
    /** @class CrystalBall
     *  The special parametrization of ``Crystal Ball-function''
     *  @see Ostap::Math::CrystalBall
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-05-13
     */
    class  CrystalBall : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::CrystalBall, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      CrystalBall
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigma     ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          n         ) ;  // n-1
      /// "copy" constructor
      CrystalBall ( const CrystalBall& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~CrystalBall () ;
      /// clone
      CrystalBall* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      CrystalBall () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::CrystalBall& function() const { return m_cb ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_n      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CrystalBall m_cb ;                  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallRS
     *  The special parametrization of ``Crystal Ball-function''
     * rigth-side crystal ball
     *  @see Ostap::Math::CrystalBallRightSide
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-05-13
     */
    class  CrystalBallRS : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::CrystalBallRS, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      CrystalBallRS
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigma     ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          n         ) ;  // n-1
      /// "copy" constructor
      CrystalBallRS ( const CrystalBallRS& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~CrystalBallRS () ;
      /// clone
      CrystalBallRS* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      CrystalBallRS () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::CrystalBallRightSide& function() const { return m_cb ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_n      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CrystalBallRightSide m_cb ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallDS
     *  double-sided ``Crystal Ball-function''
     *  for description of gaussian with the tail
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @see Ostap::Math::CrystalBallDoubleSided
     *  @date 2011-05-25
     */
    class  CrystalBallDS : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::CrystalBallDS, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      CrystalBallDS
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigma     ,  //
        RooAbsReal&          alphaL    ,  // alpha_L
        RooAbsReal&          nL        ,  //     n_L - 1
        RooAbsReal&          alphaR    ,  // alpha_R - 1
        RooAbsReal&          nR        ); //     n_R
      /// "copy" constructor
      CrystalBallDS ( const CrystalBallDS& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~CrystalBallDS() ;
      /// clone
      CrystalBallDS* clone ( const char* name ) const override;
      // ======================================================================
    public : // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      CrystalBallDS() {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::CrystalBallDoubleSided& function() const
      { return m_cb2 ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_alphaL ;
      RooRealProxy m_nL     ;
      RooRealProxy m_alphaR ;
      RooRealProxy m_nR     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CrystalBallDoubleSided m_cb2 ;       // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Needham
     *  The special parametrization by Matthew NEEDHAM of
     *  ``Crystal Ball-function'' nicely suitable for \f$J/\psi\f$-peak
     *  @thank Matthew Needham
     *  @see Ostap::Math::Needham
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-05-13
     */
    class  Needham : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Needham, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Needham
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigma     ,
        RooAbsReal&          a0        ,
        RooAbsReal&          a1        ,
        RooAbsReal&          a2        ) ;
      /// "copy" constructor
      Needham ( const Needham& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~Needham () ;
      /// clone
      Needham* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Needham () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Needham& function() const { return m_needham ; }
      /// get current alpha
      double                      alpha   () const ;
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_a0     ;
      RooRealProxy m_a1     ;
      RooRealProxy m_a2     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Needham m_needham ;                  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Apolonios
     *  A modified gaussian with power-law tail on rigth ride and exponential
     *  tail on low-side
     *  The function is proposed by Diego Martinez Santos
     *  @see http://arxiv.org/abs/1312.5000
     *  Here a bit modified version is used with redefined parameter <code>n</code>
     *  to be coherent with local definitions of Crystal Ball
     *
     *  @see Ostap::Math::Apolonios
     *  @author Vanya BELYAEV Ivane.BElyaev@itep.ru
     *  @date 2013-12-01
     */
    // ========================================================================
    class  Apolonios : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Apolonios, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Apolonios
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mean      ,
        RooAbsReal&          sigma     ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          n         ,
        RooAbsReal&          b         ) ;
      /// "copy" constructor
      Apolonios  ( const Apolonios& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Apolonios () ;
      /// clone
      Apolonios* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Apolonios () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Apolonios& function() const { return m_apo ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_n      ;
      RooRealProxy m_b      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Apolonios m_apo ;                // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Apolonios2
     *  A modified gaussian with exponential
     *  tails on low-side
     *
     *  @see Ostap::Math::Apolonios2
     *  @author Vanya BELYAEV Ivane.BElyaev@itep.ru
     *  @date 2013-12-01
     */
    // ========================================================================
    class  Apolonios2 : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Apolonios2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Apolonios2
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mean      ,
        RooAbsReal&          sigmaL    ,
        RooAbsReal&          sigmaR    ,
        RooAbsReal&          beta      ) ;
      /// "copy" constructor
      Apolonios2  ( const Apolonios2& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Apolonios2 () ;
      /// clone
      Apolonios2* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Apolonios2 () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Apolonios2& function() const { return m_apo2 ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigmaL ;
      RooRealProxy m_sigmaR ;
      RooRealProxy m_beta   ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Apolonios2 m_apo2 ;                // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BifurcatedGauss
     *  @see Ostap::Math::BifurkatedGauss
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2013-08-27
     */
    // ========================================================================
    class  BifurcatedGauss : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::BifurcatedGauss, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BifurcatedGauss
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          peak      ,
        RooAbsReal&          sigmaL    ,
        RooAbsReal&          sigmaR    ) ;
      /// "copy" constructor
      BifurcatedGauss ( const BifurcatedGauss& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~BifurcatedGauss () ;
      /// clone
      BifurcatedGauss* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      BifurcatedGauss () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::BifurcatedGauss& function() const { return m_bg ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_peak   ;
      RooRealProxy m_sigmaL ;
      RooRealProxy m_sigmaR ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::BifurcatedGauss m_bg ;               // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenGaussV1
     *  Simple class that implements the generalized normal distribution v1
     *  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1
     *  @see Ostap::Math::GenGaussV1
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2013-08-27
     */
    // ========================================================================
    class  GenGaussV1 : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::GenGaussV1, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GenGaussV1
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          beta      ) ;
      /// "copy" constructor
      GenGaussV1 ( const GenGaussV1& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~GenGaussV1 () ;
      /// clone
      GenGaussV1* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GenGaussV1 () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::GenGaussV1& function() const { return m_ggv1 ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_mu     ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_beta   ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenGaussV1 m_ggv1 ;                 // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenGaussV2
     *  Simple class that implements the generalized normal distribution v2
     *  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_2
     *  @see Ostap::Math::GenGaussV2
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2013-08-27
     */
    // ========================================================================
    class  GenGaussV2 : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::GenGaussV2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GenGaussV2
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          xi        ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          kappa     ) ;
      /// "copy" constructor
      GenGaussV2 ( const GenGaussV2& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~GenGaussV2 () ;
      /// clone
      GenGaussV2* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GenGaussV2 () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::GenGaussV2& function() const { return m_ggv2 ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_xi     ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_kappa  ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenGaussV2 m_ggv2 ;                 // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class SkewGauss
     *  Simple class that implements the skew normal distribution
     *  @see http://en.wikipedia.org/wiki/Skew_normal_distribution
     *  @see Ostap::Math::SkewGauss
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2013-08-27
     */
    // ========================================================================
    class  SkewGauss : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::SkewGauss, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      SkewGauss
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          xi        ,
        RooAbsReal&          omega     ,
        RooAbsReal&          alpha     ) ;
      /// "copy" constructor
      SkewGauss ( const SkewGauss& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~SkewGauss () ;
      /// clone
      SkewGauss* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      SkewGauss () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::SkewGauss& function() const { return m_sg ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_xi     ;
      RooRealProxy m_omega  ;
      RooRealProxy m_alpha  ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::SkewGauss m_sg ;                     // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bukin
     *  "Bukin"-function, aka "Modified Novosibirsk function"
     *  @see http://arxiv.org/abs/1107.5751
     *  @see http://dx.doi.org/10.1007/JHEP06(2012)141
     *  @see Ostap::Math::Bukin
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-12-05
     */
    class  Bukin : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Bukin, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Bukin
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          peak      ,   // peak position
        RooAbsReal&          sigma     ,   // "width"
        RooAbsReal&          xi        ,   // asymmetry pareamneter
        RooAbsReal&          rhoL      ,   // left tail
        RooAbsReal&          rhoR      ) ; // right tail
      /// "copy" constructor
      Bukin ( const Bukin& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Bukin () ;
      /// clone
      Bukin* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Bukin () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Bukin& function() const { return m_bukin ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_peak   ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_xi     ;
      RooRealProxy m_rhoL   ;
      RooRealProxy m_rhoR   ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Bukin m_bukin ;                      // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class StudentT
     *  Student-T distribution
     *
     *  \f[  f(y) = \frac{1}{\sqrt{\pi n}} \frac { \Gamma( \frac{n+1}{2}) } { \Gamma( \frac{n}{2}  ) }
     *  \left( 1 + \frac{y^2}{n} \right)^{ -\frac{n+1}{2}} \f],
     *  where \f$ y = \frac{x - \mu}{\sigma} \f$
     *
     *  @see Ostap::Math::StudentT
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2013-01-05
     */
    class  StudentT: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::StudentT, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      StudentT ( const char*          name      ,
                 const char*          title     ,
                 RooAbsReal&          x         ,
                 RooAbsReal&          mu        ,
                 RooAbsReal&          sigma     ,
                 RooAbsReal&          n         ) ;
      /// "copy constructor"
      StudentT ( const StudentT&      right     ,
                 const char*          name  = 0 )  ;
      /// destructor
      virtual ~StudentT() ;
      /// clone
      StudentT* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      StudentT () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::StudentT& function() const { return m_stt ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_sigma    ;
      RooRealProxy m_n        ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::StudentT m_stt ;           // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BifurcatedStudentT
     *  @see Ostap::Math::BifurcatedStudentT
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2013-01-05
     */
    class  BifurcatedStudentT: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::BifurcatedStudentT, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BifurcatedStudentT ( const char*          name      ,
                           const char*          title     ,
                           RooAbsReal&          x         ,
                           RooAbsReal&          mu        ,
                           RooAbsReal&          sigmaL    ,
                           RooAbsReal&          sigmaR    ,
                           RooAbsReal&          nL        ,
                           RooAbsReal&          nR        ) ;
      /// "copy constructor"
      BifurcatedStudentT ( const BifurcatedStudentT& right     ,
                           const char*               name  = 0 )  ;
      /// destructor
      virtual ~BifurcatedStudentT() ;
      /// clone
      BifurcatedStudentT* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      BifurcatedStudentT () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::BifurcatedStudentT& function() const { return m_stt ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_sigmaL   ;
      RooRealProxy m_sigmaR   ;
      RooRealProxy m_nL       ;
      RooRealProxy m_nR       ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::BifurcatedStudentT m_stt ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GramCharlierA
     *  The peak with Gram-Charlier type A parameterization
     *  @see Ostap::Math::GramCharlierA
     *  http://en.wikipedia.org/wiki/Edgeworth_series
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-12-05
     */
    class  GramCharlierA : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::GramCharlierA, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GramCharlierA
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mean      ,
        RooAbsReal&          sigma     ,
        RooAbsReal&          kappa3    ,
        RooAbsReal&          kappa4    );
      /// "copy" constructor
      GramCharlierA ( const GramCharlierA& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~GramCharlierA () ;
      /// clone
      GramCharlierA* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GramCharlierA () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::GramCharlierA& function() const { return m_gca ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_kappa3 ;
      RooRealProxy m_kappa4 ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GramCharlierA m_gca ;                // the function
      // ======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    // Smooth functions for background
    // ========================================================================

    // ========================================================================
    /** @class PhaseSpace2
     *  simple model for 2-body phase space
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpace2 : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::PhaseSpace2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpace2 ( const char*          name      ,
                    const char*          title     ,
                    RooAbsReal&          x         ,
                    const double         m1        ,
                    const double         m2        ) ;
      /// "copy constructor"
      PhaseSpace2 ( const PhaseSpace2& right     ,
                    const char*        name  = 0 )  ;
      /// destructor
      virtual ~PhaseSpace2() ;
      /// clone
      PhaseSpace2* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PhaseSpace2 () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function       /// access to underlying function
      const Ostap::Math::PhaseSpace2& function() const { return m_ps2 ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PhaseSpace2 m_ps2 ;           // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceLeft
     *  simple model for left-edge of N-body phase-space
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpaceLeft : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::PhaseSpaceLeft, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpaceLeft ( const char*          name      ,
                       const char*          title     ,
                       RooAbsReal&          x         ,
                       RooAbsReal&          threshold ,
                       const unsigned short N         ) ;
      /// "copy constructor"
      PhaseSpaceLeft ( const PhaseSpaceLeft& right     ,
                       const char*           name  = 0 )  ;
      /// destructor
      virtual ~PhaseSpaceLeft() ;
      /// clone
      PhaseSpaceLeft* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PhaseSpaceLeft () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PhaseSpaceLeft& function() const { return m_left ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x         ;
      RooRealProxy m_threshold ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PhaseSpaceLeft m_left ;        // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceRight
     *  simple model for right-edge of L-body phase-space in N-body decays
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpaceRight : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::PhaseSpaceRight, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpaceRight ( const char*          name      ,
                        const char*          title     ,
                        RooAbsReal&          x         ,
                        RooAbsReal&          threshold ,
                        const unsigned short L         ,
                        const unsigned short N         ) ;
      /// "copy constructor"
      PhaseSpaceRight ( const PhaseSpaceRight& right     ,
                        const char*            name  = 0 )  ;
      /// destructor
      virtual ~PhaseSpaceRight () ;
      /// clone
      PhaseSpaceRight* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PhaseSpaceRight () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PhaseSpaceRight& function() const { return m_right ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x         ;
      RooRealProxy m_threshold ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PhaseSpaceRight m_right ;     // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceNL
     *
     *  The phase space function for L-body systema from N-body decay
     *
     *  @see Ostap::Math::PhaseSpaceNL
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpaceNL : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::PhaseSpaceNL, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpaceNL ( const char*          name      ,
                     const char*          title     ,
                     RooAbsReal&          x         ,
                     RooAbsReal&          low       ,
                     RooAbsReal&          high      ,
                     const unsigned short N         ,
                     const unsigned short L         ) ;
      /// "copy" constructor
      PhaseSpaceNL ( const PhaseSpaceNL& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~PhaseSpaceNL () ;
      /// clone
      PhaseSpaceNL* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PhaseSpaceNL () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PhaseSpaceNL& function() const { return m_ps ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_low   ;
      RooRealProxy m_high  ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PhaseSpaceNL m_ps ;           // the actual function
      // ======================================================================
    };
    // ========================================================================
    /** @class PhaseSpacePol
     *  The mass-ditribtion of L-particles from N-body phase space decays,
     *  modulate with non-negative polynomial
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-10-05
     */
    class PhaseSpacePol : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::PhaseSpacePol, 1) ;
      // ======================================================================
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const double                     low    ,
        const double                     high   ,
        const unsigned short             N      ,
        const unsigned short             L      ,
        RooAbsReal&                      phi1   ) ;
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const double                     low    ,
        const double                     high   ,
        const unsigned short             N      ,
        const unsigned short             L      ,
        RooAbsReal&                      phi1   ,
        RooAbsReal&                      phi2   ) ;
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const double                     low    ,
        const double                     high   ,
        const unsigned short             N      ,
        const unsigned short             L      ,
        RooAbsReal&                      phi1   ,
        RooAbsReal&                      phi2   ,
        RooAbsReal&                      phi3   ) ;
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const double                     low    ,
        const double                     high   ,
        const unsigned short             N      ,
        const unsigned short             L      ,
        RooArgList&                      phis   ) ;
      // ======================================================================
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const Ostap::Math::PhaseSpaceNL& ps     ,
        RooAbsReal&                      phi1   ) ;
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const Ostap::Math::PhaseSpaceNL& ps     ,
        RooAbsReal&                      phi1   ,
        RooAbsReal&                      phi2   ) ;
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const Ostap::Math::PhaseSpaceNL& ps     ,
        RooAbsReal&                      phi1   ,
        RooAbsReal&                      phi2   ,
        RooAbsReal&                      phi3   ) ;
      /// constructor from all parameters
      PhaseSpacePol
      ( const char*                      name   ,
        const char*                      title  ,
        RooAbsReal&                      x      ,
        const Ostap::Math::PhaseSpaceNL& ps     ,
        RooArgList&                      phis   ) ;
      // ======================================================================
      // "copy" constructor
      // ======================================================================
      PhaseSpacePol ( const PhaseSpacePol& right    ,
                      const char*          name = 0 ) ;
      /// destructor
      virtual ~PhaseSpacePol () ;
      /// clone
      PhaseSpacePol* clone( const char* name )  const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PhaseSpacePol () {} ;
      // ======================================================================
    public:
      // ======================================================================
      Double_t evaluate () const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpacePol& function() const { return m_ps ; }
      // ======================================================================
   private:
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
      TIterator* m_iterator;  //! do not persist
      // ======================================================================
    private:
      // ======================================================================
      /// the actual phase space function
      mutable Ostap::Math::PhaseSpacePol m_ps ;  // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpace23L
     *  simple model for 2-body phase space from 3-body decays with
     *  the orbital momenta:
     *
     *   \f$ f \propto q^{2\ell+1}p^{2L+1}\f$, where
     *     \f$\ell\f$ is the orbital momentum of the pair of particles,
     *    and \f$L\f$ is the orbital momentum between the pair and
     *    the third particle.
     *   E.g. taking \f$\ell=0, L=1\f$, one can get the S-wave contribution for
     *   \f$\pi^+\pi^-\f$-mass from \f$B^0\rightarrowJ/\psi\pi^+\pi^-\f$ decay.
     *
     *  @see Ostap::Math::PhaseSpace23L
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-01
     */
    class PhaseSpace23L : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::PhaseSpace23L, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and the third particle
       *  @param l  the angular momentum between the first and the second particle
       */
      PhaseSpace23L ( const char*          name      ,
                      const char*          title     ,
                      RooAbsReal&          x         ,
                      const double         m1        ,
                      const double         m2        ,
                      const double         m3        ,
                      const double         m         ,
                      const unsigned short L         ,
                      const unsigned short l     = 0 ) ;
      /// "copy constructor"
      PhaseSpace23L ( const PhaseSpace23L& right , const char* name = 0 )  ;
      /// destructor
      virtual ~PhaseSpace23L() ;
      /// clone
      PhaseSpace23L* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PhaseSpace23L () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PhaseSpace23L& function() const { return m_ps23L ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      Ostap::Math::PhaseSpace23L m_ps23L ;               // the actual function
      // ======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    // Smooth empirical models for background
    // ========================================================================

    // ========================================================================
    /** @class PolyPositive
     *  PolyPositive polynomial
     *  @see Ostap::Math::Positive
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  PolyPositive: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PolyPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      PolyPositive
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        const RooArgList&    coeffs    ,
        const double         xmin      ,
        const double         xmax      ) ;
      /// copy
      PolyPositive
      ( const PolyPositive&     right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~PolyPositive() ;
      /// clone
      PolyPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PolyPositive () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Positive& function() const { return m_positive ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Positive m_positive ;               // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PolyPositiveEven
     *  positive even polynomial
     *  @see Ostap::Math::PositiveEven
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2016-10-03
     */
    class  PolyPositiveEven: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PolyPositiveEven, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      PolyPositiveEven
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        const RooArgList&    coeffs    ,
        const double         xmin      ,
        const double         xmax      ) ;
      /// copy
      PolyPositiveEven
      ( const PolyPositiveEven& right     ,
        const char*             name = 0  ) ;
      /// destructor
      virtual ~PolyPositiveEven() ;
      /// clone
      PolyPositiveEven* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PolyPositiveEven () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PositiveEven& function() const { return m_even ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PositiveEven m_even ;               // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PolyMonothonic
     *  positive monothonic  polynomial
     *  @see Ostap::Math::Monothonic
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  PolyMonothonic: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PolyMonothonic, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      PolyMonothonic
        ( const char*          name       ,
          const char*          title      ,
          RooAbsReal&          x          ,
          const RooArgList&    coeffs     ,
          const double         xmin       ,
          const double         xmax       ,
          const bool           increasing ) ;
      /// copy
      PolyMonothonic
        ( const PolyMonothonic&     right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~PolyMonothonic() ;
      /// clone
      PolyMonothonic* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PolyMonothonic () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Monothonic& function() const { return m_monothonic ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Monothonic m_monothonic ;            // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PolyConvex
     *  Positive polynomial with fixes sign first and second derivatives
     *  @see Ostap::Math::Convex
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  PolyConvex: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PolyConvex, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      PolyConvex
        ( const char*          name       ,
          const char*          title      ,
          RooAbsReal&          x          ,
          const RooArgList&    coeffs     ,
          const double         xmin       ,
          const double         xmax       ,
          const bool           increasing ,
          const bool           convex     ) ;
      /// copy
      PolyConvex
        ( const PolyConvex&    right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~PolyConvex () ;
      /// clone
      PolyConvex* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PolyConvex () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Convex& function() const { return m_convex ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Convex m_convex ;                    // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PolyConvexOnly
     *  Positive polynomial with fixes sign second derivatives
     *  @see Ostap::Math::ConvexOnly
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  PolyConvexOnly: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PolyConvexOnly, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      PolyConvexOnly
        ( const char*          name       ,
          const char*          title      ,
          RooAbsReal&          x          ,
          const RooArgList&    coeffs     ,
          const double         xmin       ,
          const double         xmax       ,
          const bool           convex     ) ;
      /// copy
      PolyConvexOnly
        ( const PolyConvexOnly& right     ,
          const char*           name = 0  ) ;
      /// destructor
      virtual ~PolyConvexOnly () ;
      /// clone
      PolyConvexOnly* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PolyConvexOnly () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::ConvexOnly& function() const { return m_convex ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ConvexOnly m_convex ;                // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ExpoPositive
     *  exponential multiplied on positive polynomial
     *  @see Ostap::Math::Positive
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  ExpoPositive: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::ExpoPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      ExpoPositive
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          tau       ,
        const RooArgList&    coeffs    ,
        const double         xmin      ,
        const double         xmax      ) ;
      /// copy
      ExpoPositive
      ( const ExpoPositive&     right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~ExpoPositive() ;
      /// clone
      ExpoPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      ExpoPositive () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::ExpoPositive& function() const { return m_positive ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_tau  ;
      RooListProxy m_phis ;
      // ======================================================================
      TIterator* m_iterator;  //! do not persist
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ExpoPositive m_positive ;           // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Sigmoid
     *  The product of sigmoid function and positive polynomial
     *  @see Ostap::Math::Sigmoid
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  PolySigmoid: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PolySigmoid, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      PolySigmoid
        ( const char*          name       ,
          const char*          title      ,
          RooAbsReal&          x          ,
          const RooArgList&    coeffs     ,
          const double         xmin       ,
          const double         xmax       ,
          RooAbsReal&          alpha      ,
          RooAbsReal&          x0         ) ;
      /// copy
      PolySigmoid
        ( const PolySigmoid&   right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~PolySigmoid() ;
      /// clone
      PolySigmoid* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PolySigmoid () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Sigmoid& sigmoid () const { return m_sigmoid    ; }
      /// access to underlying function
      const Ostap::Math::Sigmoid& function() const { return   sigmoid () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     ;
      RooListProxy m_phis  ;
      RooRealProxy m_alpha ;
      RooRealProxy m_x0    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Sigmoid m_sigmoid ;                 // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class TwoExpoPositive
     *  simple difference of two exponents modulated with positive polynomials
     *  \f$ f(x) = e_2(x) * p_n(x) \f$, where
     *  \f$ e_2(x) \propto
     *        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} =
     *        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
     *  @see Ostap::Math::TwoExpoPositive
     *  @see Ostap::Math::ExpoPositive
     *  @see Ostap::Math::Positive
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  TwoExpoPositive: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::TwoExpoPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      TwoExpoPositive
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          delta     ,
        RooAbsReal&          x0        ,
        const RooArgList&    coeffs    ,
        const double         xmin      ,
        const double         xmax      ) ;
      /// copy
      TwoExpoPositive
        ( const TwoExpoPositive& right     ,
          const char*           name = 0  ) ;
      /// destructor
      virtual ~TwoExpoPositive() ;
      /// clone
      TwoExpoPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      TwoExpoPositive() {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::TwoExpoPositive& function() const { return m_2expopos ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_alpha ;
      RooRealProxy m_delta ;
      RooRealProxy m_x0    ;
      RooListProxy m_phis  ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::TwoExpoPositive m_2expopos;         // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GammaDist
     *  Gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     *  @see Ostap::Math::GammaDist
     */
    class  GammaDist : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::GammaDist, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GammaDist ( const char*          name      ,
                  const char*          title     ,
                  RooAbsReal&          x         ,
                  RooAbsReal&          k         ,
                  RooAbsReal&          theta     ) ;
      /// "copy constructor"
      GammaDist ( const GammaDist&     right     ,
                  const char*          name  = 0 )  ;
      /// destructor
      virtual ~GammaDist () ;
      /// clone
      GammaDist* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GammaDist () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::GammaDist& function() const { return m_gamma ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_k        ;
      RooRealProxy m_theta    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GammaDist m_gamma ; // the actual function
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
     *  @see Ostap::Math::GenGammaDist
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     *  @see Ostap::Math::GenGammaDist
     */
    class  GenGammaDist : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::GenGammaDist, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GenGammaDist
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          k         ,
        RooAbsReal&          theta     ,
        RooAbsReal&          p         ,
        RooAbsReal&          low       ) ;
      /// "copy constructor"
      GenGammaDist ( const GenGammaDist&  right     ,
                     const char*          name  = 0 )  ;
      /// destructor
      virtual ~GenGammaDist () ;
      /// clone
      GenGammaDist* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GenGammaDist () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::GenGammaDist& function() const { return m_ggamma ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_k        ;
      RooRealProxy m_theta    ;
      RooRealProxy m_p        ;
      RooRealProxy m_low      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenGammaDist m_ggamma ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Amoroso
     *  Another view on generalized gamma distribtion
     *  http://arxiv.org/pdf/1005.3274
     *  @see Ostap::Math::Amoroso
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-12-05
     */
    class  Amoroso : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Amoroso, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Amoroso
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          theta     ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          beta      ,
        RooAbsReal&          a         ) ;
      /// "copy" constructor
      Amoroso ( const Amoroso& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Amoroso () ;
      /// clone
      Amoroso* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Amoroso () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Amoroso& function() const { return m_amoroso ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_theta  ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_beta   ;
      RooRealProxy m_a      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Amoroso m_amoroso ;                  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LogGammaDist
     *  Distribution for log(x), where x follows
     *  gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     *  @see Ostap::Math::LogGammaDist
     */
    class  LogGammaDist : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::LogGammaDist, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      LogGammaDist ( const char*          name      ,
                     const char*          title     ,
                     RooAbsReal&          x         ,
                     RooAbsReal&          k         ,
                     RooAbsReal&          theta     ) ;
      /// "copy constructor"
      LogGammaDist ( const LogGammaDist&  right     ,
                     const char*          name  = 0 )  ;
      /// destructor
      virtual ~LogGammaDist () ;
      /// clone
      LogGammaDist* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      LogGammaDist () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::LogGammaDist& function() const { return m_gamma ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_k        ;
      RooRealProxy m_theta    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::LogGammaDist m_gamma ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Log10GammaDist
     *  Distribution for log10(x), where x follows
     *  gamma-distribution shape/scale parameters
     *  http://en.wikipedia.org/wiki/Gamma_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     *  @see Ostap::Math::Log10GammaDist
     */
    class  Log10GammaDist : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Log10GammaDist, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Log10GammaDist ( const char*           name      ,
                       const char*           title     ,
                       RooAbsReal&           x         ,
                       RooAbsReal&           k         ,
                       RooAbsReal&           theta     ) ;
      /// "copy constructor"
      Log10GammaDist ( const Log10GammaDist& right     ,
                       const char*           name  = 0 )  ;
      /// destructor
      virtual ~Log10GammaDist () ;
      /// clone
      Log10GammaDist* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Log10GammaDist () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Log10GammaDist& function() const { return m_gamma ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_k        ;
      RooRealProxy m_theta    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Log10GammaDist m_gamma ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LogGamma
     *  - http://arxiv.org/pdf/1005.3274
     *  - Prentice, R. L. (1974). A log gamma model and its maximum likelihood
     *                            estimation. Biometrika 61, 539
     *  - Johnson, N. L., Kotz, S., and Balakrishnan, N. (1995). Continuous
     *            univariate distributions, 2nd ed. Vol. 2. Wiley, New York.
     *  - Bartlett, M. S. and G., K. M. (1946). The statistical analysis of
     *                  variance-heterogeneity and the logarithmic transformation.
     *                 J. Roy. Statist. Soc. Suppl. 8, 1, 128.
     *
     *  dot not mix with Ostap::Models::LogGammaDist
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     *  @see Ostap::Math::LogGamma
     */
    class  LogGamma : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::LogGamma, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      LogGamma ( const char*          name      ,
                 const char*          title     ,
                 RooAbsReal&          x         ,
                 RooAbsReal&          nu        ,
                 RooAbsReal&          lambda    ,
                 RooAbsReal&          alpha     ) ;
      /// "copy constructor"
      LogGamma ( const LogGamma&      right     ,
                 const char*          name  = 0 )  ;
      /// destructor
      virtual ~LogGamma () ;
      /// clone
      LogGamma* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      LogGamma () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::LogGamma& function() const { return m_lgamma ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_nu       ;
      RooRealProxy m_lambda   ;
      RooRealProxy m_alpha    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::LogGamma m_lgamma ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BetaPrime
     *  http://en.wikipedia.org/wiki/Beta_prime_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-05-11
     *  @see Ostap::Math::BetaPrime
     */
    class  BetaPrime : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::BetaPrime, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BetaPrime
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           alpha     ,
        RooAbsReal&           beta      ,
        RooAbsReal&           scale     ,
        RooAbsReal&           shift     ) ;
      /// "copy constructor"
      BetaPrime
      ( const BetaPrime&      right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~BetaPrime () ;
      /// clone
      BetaPrime* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      BetaPrime () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::BetaPrime& function() const { return m_betap ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_alpha    ;
      RooRealProxy m_beta     ;
      RooRealProxy m_scale    ;
      RooRealProxy m_shift    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::BetaPrime m_betap ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Landau
     *  http://en.wikipedia.org/wiki/Landau_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-06-07
     *  @see Ostap::Math::Landau
     */
    class  Landau : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Landau, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Landau
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           scale     ,
        RooAbsReal&           shift     ) ;
      /// "copy constructor"
      Landau
      ( const Landau&         right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~Landau  () ;
      /// clone
      Landau* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Landau () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Landau& function() const { return m_landau ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_scale    ;
      RooRealProxy m_shift    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Landau m_landau ; // the actual function
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

     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-08-02
     *  @see Ostap::Math::SinhAsinh
     */
    class  SinhAsinh : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::SinhAsinh, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      SinhAsinh
        ( const char*           name      ,
          const char*           title     ,
          RooAbsReal&           x         ,
          RooAbsReal&           mu        ,
          RooAbsReal&           sigma     ,
          RooAbsReal&           epsilon   ,
          RooAbsReal&           delta     ) ;
      /// "copy constructor"
      SinhAsinh
        ( const SinhAsinh&      right     ,
          const char*           name  = 0 )  ;
      /// destructor
      virtual ~SinhAsinh  () ;
      /// clone
      SinhAsinh* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      SinhAsinh () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::SinhAsinh& function() const { return m_sinhasinh ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_sigma    ;
      RooRealProxy m_epsilon  ;
      RooRealProxy m_delta    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::SinhAsinh m_sinhasinh ; // the actual function
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
     *  @see Ostap::Math::JohnsonSU
     */
    class  JohnsonSU : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::JohnsonSU, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      JohnsonSU
        ( const char*           name      ,
          const char*           title     ,
          RooAbsReal&           x         ,
          RooAbsReal&           xi        ,
          RooAbsReal&           lam       ,
          RooAbsReal&           delta     ,
          RooAbsReal&           gamma     ) ;
      /// "copy constructor"
      JohnsonSU
        ( const JohnsonSU&      right     ,
          const char*           name  = 0 )  ;
      /// destructor
      virtual ~JohnsonSU  () ;
      /// clone
      JohnsonSU* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      JohnsonSU () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::JohnsonSU& function() const { return m_johnsonSU ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_xi       ;
      RooRealProxy m_lambda   ;
      RooRealProxy m_delta    ;
      RooRealProxy m_gamma    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::JohnsonSU m_johnsonSU ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Atlas
     *  Modified gaussian function
     *  \f$  f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\deltax/2}}}{2})\f$,
     *  where \f$\delta x = \left| x - \mu \right|/\sigma\f$
     *  Fuction is taken from http://arxiv.org/abs/arXiv:1507.07099
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-08-21
     *  @see Ostap::Math::Atlas
     */
    class  Atlas : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Atlas, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Atlas
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           mu        ,
        RooAbsReal&           sigma     ) ;
      /// "copy constructor"
      Atlas
      ( const Atlas&          right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~Atlas  () ;
      /// clone
      Atlas* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Atlas () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Atlas& function() const { return m_atlas ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_sigma    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Atlas m_atlas ; // the actual function
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
     *  @see Ostap::Math::Sech
     */
    class  Sech : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Sech, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Sech
        ( const char*           name      ,
          const char*           title     ,
          RooAbsReal&           x         ,
          RooAbsReal&           mu        ,
          RooAbsReal&           sigma     ) ;
      /// "copy constructor"
      Sech
        ( const Sech&          right     ,
          const char*           name  = 0 )  ;
      /// destructor
      virtual ~Sech  () ;
      /// clone
      Sech* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Sech () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Sech& function() const { return m_sech ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_sigma    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Sech m_sech ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Logistic
     *  aka "Sech-square"
     *  \f$ f(x;\mu;s) = \frac{1}{4s}sech^2\left(\frac{x-\mu}{2s}\right)\f$
     *  where
     *  \f$  s = \sigma \frac{\sqrt{3}}{\pi}\f$
     *  @see https://en.wikipedia.org/wiki/Logistic_distribution
     *  @see Ostap::Math::Logistic
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-14
     */
    class  Logistic: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Logistic, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Logistic
        ( const char*           name      ,
          const char*           title     ,
          RooAbsReal&           x         ,
          RooAbsReal&           mu        ,
          RooAbsReal&           sigma     ) ;
      /// "copy constructor"
      Logistic
        ( const Logistic&       right     ,
          const char*           name  = 0 )  ;
      /// destructor
      virtual ~Logistic  () ;
      /// clone
      Logistic* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Logistic () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Logistic& function() const { return m_logistic ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_sigma    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Logistic m_logistic ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Argus
     *  http://en.wikipedia.org/wiki/ARGUS_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-06-07
     *  @see Ostap::Math::Argus
     */
    class  Argus : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Argus, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Argus
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           shape     ,
        RooAbsReal&           high      ,
        RooAbsReal&           low       ) ;
      /// "copy constructor"
      Argus
      ( const Argus&          right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~Argus  () ;
      /// clone
      Argus* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Argus () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Argus& function() const { return m_argus ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_shape    ;
      RooRealProxy m_high     ;
      RooRealProxy m_low      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Argus m_argus ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  Slash 
     *  ``Slash''-distribution -  symmetric peak with veyr heavy tail
     *  @see https://en.wikipedia.org/wiki/Slash_distribution
     *  Tails arew so heavy that moments (e.g. variance) do not exist 
     *  @see Ostap::Math::Slash
     */
    class Slash : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Slash, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Slash
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           mu        ,
        RooAbsReal&           scale     ) ;
      /// "copy constructor"
      Slash
      ( const Slash&          right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~Slash  () ;
      /// clone
      Slash* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Slash () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Slash& function () const { return m_slash ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mu    ;
      RooRealProxy m_scale ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Slash m_slash ; // the actual function
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class AsymmetricLaplace
     *  @see https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution
     *  @see  Ostap::Math::AsymmetricLaplace
     */
    class AsymmetricLaplace: public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDef(Ostap::Models::AsymmetricLaplace, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      AsymmetricLaplace
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           mu        ,
        RooAbsReal&           lambdaL   ,
        RooAbsReal&           lambdaR   ) ;
      /// "copy constructor"
      AsymmetricLaplace
      ( const AsymmetricLaplace& right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~AsymmetricLaplace  () ;
      /// clone
      AsymmetricLaplace* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      AsymmetricLaplace () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::AsymmetricLaplace& function () const { return m_laplace ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x       ;
      RooRealProxy m_mu      ;
      RooRealProxy m_lambdaL ;
      RooRealProxy m_lambdaR ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::AsymmetricLaplace m_laplace ; // the actual function
      // ======================================================================
    } ;   
    // ========================================================================
    /** @class Tsallis
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
     *
     *  @see Ostap::Math::Tsallis
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-07-11
     */
    class  Tsallis : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Tsallis, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Tsallis
        ( const char*           name      ,
          const char*           title     ,
          RooAbsReal&           x         ,
          RooAbsReal&           n         ,   // parameter N
          RooAbsReal&           T         ,   // parameter T
          RooAbsReal&           mass      ) ; // particle mass (fixed)
      /// "copy constructor"
      Tsallis
        ( const Tsallis&        right     ,
          const char*           name  = 0 )  ;
      /// destructor
      virtual ~Tsallis  () ;
      /// clone
      Tsallis* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Tsallis () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Tsallis& function () const { return m_tsallis ; }
      const Ostap::Math::Tsallis& tsallis  () const { return m_tsallis ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_n        ;
      RooRealProxy m_T        ;
      RooRealProxy m_mass     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Tsallis m_tsallis ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class QGSM
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
     *  @see Ostap::Math::QGSM
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-07-11
     */
    class  QGSM : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::QGSM, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      QGSM
        ( const char*           name      ,
          const char*           title     ,
          RooAbsReal&           x         ,
          RooAbsReal&           b         ,   // parameter b
          RooAbsReal&           mass      ) ; // particle mass (fixed)
      /// "copy constructor"
      QGSM
        ( const QGSM&           right     ,
          const char*           name  = 0 )  ;
      /// destructor
      virtual ~QGSM  () ;
      /// clone
      QGSM* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      QGSM  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::QGSM& function() const { return m_qgsm ; }
      const Ostap::Math::QGSM& qgsm    () const { return m_qgsm ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_b        ;
      RooRealProxy m_mass     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::QGSM m_qgsm ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class TwoExpos
     *  simple difference of two exponents
     *  \f$ f \propto
     *        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} =
     *        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  TwoExpos : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::TwoExpos, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      TwoExpos
        ( const char*           name  ,
          const char*           title ,
          RooAbsReal&           x     ,
          RooAbsReal&           alpha ,
          RooAbsReal&           delta ,
          RooAbsReal&           x0    ) ;
      /// "copy constructor"
      TwoExpos
        ( const TwoExpos&          right     ,
          const char*           name  = 0 )  ;
      /// destructor
      virtual ~TwoExpos () ;
      /// clone
      TwoExpos* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      TwoExpos  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::TwoExpos& function() const { return m_2expos ; }
      const Ostap::Math::TwoExpos& twoexpos() const { return m_2expos ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_alpha ;
      RooRealProxy m_delta ;
      RooRealProxy m_x0    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::TwoExpos m_2expos ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  DoubleGauss
     *  Simple double Gaussian PDF
     *  suitable as resolution model
     */
    class DoubleGauss: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::DoubleGauss, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      DoubleGauss ( const char*          name      , 
                    const char*          title     ,
                    RooAbsReal&          x         ,
                    RooAbsReal&          sigma     ,   // narrow sigma 
                    RooAbsReal&          fraction  ,   // fraction of narrow sigma 
                    RooAbsReal&          scale     ,   // wide/narrow sigma ratio    
                    RooAbsReal&          mean      ) ; // mean, presumably fixed at 0
      /// "copy" constructor 
      DoubleGauss ( const DoubleGauss& , const char* name = 0 ) ;
      /// virtual destructor 
      virtual ~DoubleGauss(){} ;
      /// clone 
      DoubleGauss* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      DoubleGauss () {} ;
      // =====================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // =====================================================================
      Int_t    getAnalyticalIntegral 
      ( RooArgSet&  allVars       , 
        RooArgSet&  analVars      , 
        const char* rangeName = 0 ) const override ;
      Double_t analyticalIntegral
      ( Int_t       code          , 
        const char* rangeName = 0 ) const override ;
      // =====================================================================
    public:
      // =====================================================================
      // the actual evaluation of function 
      Double_t evaluate() const override ;
      // =====================================================================
    public:
      // =====================================================================
      /// access to underlying function
      const Ostap::Math::DoubleGauss& function () const { return m_2gauss ; }
      // =====================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x         ;
      RooRealProxy m_sigma     ;
      RooRealProxy m_fraction  ;
      RooRealProxy m_scale     ;
      RooRealProxy m_mean      ;
      // ======================================================================
    protected: // the function
      // ======================================================================
      mutable Ostap::Math::DoubleGauss  m_2gauss ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Gumbel
     *  Gumbel  distribution:
     *  https://en.wikipedia.org/wiki/Gumbel_distribution
     *  \f$  G(x;\mu,\beta) = \frac{1}{\left|\beta\right| e^{-(z+e^{-z}}\f$, 
     *  where \f$ z = \frac{x-\mu}{\beta}\f$.
     *  Importnat  cases if \f$ E(x) = e^{-\tau x}\f$, and:
     *  - \f$ z \equiv  \log(x)\f$, then \f$ F(z) = E(x) = G(z, -log(\tau) , 1 ) \f$, 
     *  - \f$ z \equiv -\log(x)\f$, then \f$ F(z) = E(x) = G(z, -log(\tau) , 1 ) \f$.
     *  As a direct sequence,  a sum of exponential componets is transformed to 
     *  a sum of ``peak-like'' Gumbel  stuctures
     *  @seee Ostap::Math::Gumbel
     */
    class  Gumbel : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Gumbel, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  z    the  variable 
       *  @param  eta  the shape eta-parameter 
       *  @param  b    the scale b-parameter 
       *  @param  xmin the bias  parameter
       */
      Gumbel ( const char*          name      , 
               const char*          title     ,
               RooAbsReal&          x         ,   // observable 
               RooAbsReal&          mu        ,   // mode/shift
               RooAbsReal&          beta      ) ; // beta/scale 
      /// "copy" constructor 
      Gumbel ( const Gumbel& , const char* name = 0 ) ;
      /// clone 
      Gumbel* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      Gumbel () {} ;
      // =====================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // =====================================================================
      Int_t    getAnalyticalIntegral 
      ( RooArgSet&  allVars       , 
        RooArgSet&  analVars      , 
        const char* rangeName = 0 ) const override ;
      Double_t analyticalIntegral
      ( Int_t       code          , 
        const char* rangeName = 0 ) const override ;
      // =====================================================================
    public:
      // =====================================================================
      // the actual evaluation of function 
      Double_t evaluate() const override ;
      // =====================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Gumbel& function () const { return m_gumbel ; }
      const Ostap::Math::Gumbel& gumbel   () const { return m_gumbel ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_mu   ;
      RooRealProxy m_beta ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::Gumbel m_gumbel ;
      // =====================================================================
    } ;
    // ========================================================================
    /** @class Weibull
     *  3-parameter  Weibull distribution 
     *  \f$ f(x,\lambda,k,x_0) = \frac{k}{\lambda}  y^{k-1} e^{-y^k}\f$, where 
     *  \f$ y \equiv \frac{x-x_0}{\lambda}\f$
     *  @see https://en.wikipedia.org/wiki/Weibull_distribution
     *  @seee Ostap::Math::Weibull
     */
    class Weibull : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::Weibull, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  x    the  variable 
       *  @param  scale the scale parameter
       *  @param  shape the shape parameter 
       *  @param  shift the shift parameter 
       */
      Weibull ( const char*          name      , 
                const char*          title     ,
                RooAbsReal&          x         ,   // observable 
                RooAbsReal&          scale     ,   // scale/lambda 
                RooAbsReal&          shape     ,   // shape/k 
                RooAbsReal&          shift     ) ; // shift/x0 
      /// "copy" constructor 
      Weibull ( const Weibull& , const char* name = 0 ) ;
      /// clone 
      Weibull* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      Weibull () {} ;
      // =====================================================================
    public:
      // =====================================================================
      Int_t    getAnalyticalIntegral 
        ( RooArgSet&  allVars       , 
          RooArgSet&  analVars      , 
          const char* rangeName = 0 ) const override ;
      Double_t analyticalIntegral
        ( Int_t       code          , 
          const char* rangeName = 0 ) const override ;
      // =====================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // =====================================================================
      // the actual evaluation of function 
      Double_t evaluate() const override ;
      // =====================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Weibull& function () const { return m_weibull ; }
      const Ostap::Math::Weibull& weibull  () const { return m_weibull ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_scale ;
      RooRealProxy m_shape ;
      RooRealProxy m_shift ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::Weibull m_weibull ;
      // =====================================================================
    } ;
    // ========================================================================
    /** @class RaisingCosine
     *  "Raising cosine" distribution
     *  \f$ f(x,\mu,s) = \frac{1}{2s}   \left( 1   +\cos \pi y \right)  \f$, 
     *  where \f$  y  \equiv = \frac{x-\mu}{s}\f$ 
     *  @see https://en.wikipedia.org/wiki/Raised_cosine_distribution
     *  @see Ostap::Math::RaisngCosine 
     */
    class RaisingCosine: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::RaisingCosine, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  x      the variable 
       *  @param  mean   the mean/mode/median/location 
       *  @param  scale  the scale parameter 
       */
      RaisingCosine ( const char*          name      , 
                      const char*          title     ,
                      RooAbsReal&          x         ,   // observable 
                      RooAbsReal&          mean      ,   // mean/mode/location
                      RooAbsReal&          scale     ) ; // scale parameter 
      /// "copy" constructor 
      RaisingCosine ( const RaisingCosine& , const char* name = 0 ) ;
      /// clone 
      RaisingCosine* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      RaisingCosine () {} ;
      // =====================================================================
    public:
      // =====================================================================
      Int_t    getAnalyticalIntegral 
        ( RooArgSet&  allVars       , 
          RooArgSet&  analVars      , 
          const char* rangeName = 0 ) const override ;
      Double_t analyticalIntegral
        ( Int_t       code          , 
          const char* rangeName = 0 ) const override ;
      // =====================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // =====================================================================
      // the actual evaluation of function 
      Double_t evaluate() const override ;
      // =====================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::RaisingCosine& function () const { return m_rcos ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mean  ;
      RooRealProxy m_scale ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::RaisingCosine m_rcos ;
      // =====================================================================
    } ;
    // ========================================================================
    /** @class QGaussian
     *  q-Gaussian distribution:
     *  \f$ f(x) = \frac{ \sqrt{\beta}}{C_q} e_q (-\beta (x-\mu)^2)$, 
     *  where  \f$ e_q (x) = \left( 1 + (1-q)x\right)^{\frac{1}{1-q}}\f$ 
     *  @see https://en.wikipedia.org/wiki/Q-Gaussian_distribution
     *  If is equal to 
     *  - scaled version of Student' t-distribution for 1<q<3
     *  - Gaussian distribution for q = 1 
     *  - has finite  support for q<1 
     *  @see Ostap::Math::QGaussian
     *  Here we use \f$ \beta = \frac{1}{2\sigma^2}\f$
     */
    class QGaussian: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Models::QGaussian, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  x      the variable 
       *  @param  mean   the mean/mode/median/location 
       *  @param  q      the q-value 
       *  @param  scale  the scale parameter 
       */
      QGaussian ( const char*          name      , 
                  const char*          title     ,
                  RooAbsReal&          x         ,   // observable 
                  RooAbsReal&          mean      ,   // mean/mode/location
                  RooAbsReal&          q         ,   // q-value
                  RooAbsReal&          scale     ) ; // scale parameter 
      /// "copy" constructor 
      QGaussian ( const QGaussian& , const char* name = 0 ) ;
      /// clone 
      QGaussian* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      QGaussian () {} ;
      // =====================================================================
    public:
      // =====================================================================
      Int_t    getAnalyticalIntegral 
        ( RooArgSet&  allVars       , 
          RooArgSet&  analVars      , 
          const char* rangeName = 0 ) const override ;
      Double_t analyticalIntegral
        ( Int_t       code          , 
          const char* rangeName = 0 ) const override ;
      // =====================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // =====================================================================
      // the actual evaluation of function 
      Double_t evaluate() const override ;
      // =====================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::QGaussian& function () const { return m_qgauss ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mean  ;
      RooRealProxy m_q     ;
      RooRealProxy m_scale ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::QGaussian m_qgauss ;
      // =====================================================================
    } ;
    // ========================================================================

    // ========================================================================
    // 1D-splines
    // ========================================================================

    // ========================================================================
    /** @class PositiveSpline
     *  The special spline for non-negative function
     *  Actually it is a sum of M-splines with non-negative coefficients
     *  \f$ f(x) = \sum_i \alpha_i * M_i^k(x) \f$,
     *  with constraints  \f$  \sum_i \alpha_i=1\f$
     *  and \f$ 0 \le \alpha_i\f$.
     *  @see http://en.wikipedia.org/wiki/M-spline
     *  @see http://en.wikipedia.org/wiki/B-spline
     *  @see Ostap::Math::PositiveSpline
     */
    class   PositiveSpline : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PositiveSpline, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor with the spline
       *  @param name  the name
       *  @param title the  title
       *  @param x     the  variable
       *  @param spine the spline
       *  @param phis  vector of parameters
       */
      PositiveSpline
        ( const char*                        name,
          const char*                        title     ,
          RooAbsReal&                        x         ,
          const Ostap::Math::PositiveSpline& spline    ,   // the spline
          RooArgList&                        phis      ) ; // parameters
      /// copy
      PositiveSpline
        ( const PositiveSpline& right     ,
          const char*           name = 0  ) ;
      /// destructor
      virtual ~PositiveSpline() ;
      /// clone
      PositiveSpline* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PositiveSpline  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PositiveSpline& function() const { return m_spline ; }
      const Ostap::Math::PositiveSpline& spline  () const { return m_spline ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PositiveSpline m_spline ;            // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class MonothonicSpline
     *  The special spline for non-negative monothonic function
     *  @see http://en.wikipedia.org/wiki/I-spline
     *  @see http://en.wikipedia.org/wiki/M-spline
     *  @see http://en.wikipedia.org/wiki/B-spline
     *  @see Ostap::Math::PositiveSpline
     */
    class   MonothonicSpline : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::MonothonicSpline, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor with the spline
       *  @param name  the name
       *  @param title the  title
       *  @param x     the  variable
       *  @param spine the spline
       *  @param phis  vector of parameters
       */
      MonothonicSpline
        ( const char*                          name,
          const char*                          title     ,
          RooAbsReal&                          x         ,
          const Ostap::Math::MonothonicSpline& spline    ,   // the spline
          RooArgList&                          phis      ) ; // parameters
      /// copy
      MonothonicSpline
        ( const MonothonicSpline& right     ,
          const char*             name = 0  ) ;
      /// destructor
      virtual ~MonothonicSpline() ;
      /// clone
      MonothonicSpline* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      MonothonicSpline  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::MonothonicSpline& function() const { return m_spline ; }
      const Ostap::Math::MonothonicSpline& spline  () const { return m_spline ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::MonothonicSpline m_spline ;          // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class ConvexPnlySpline
     *  The special spline for non-negative
     *  convex or concave function
     *  @see Ostap::Math::ConvexSpline
     */
    class   ConvexOnlySpline : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::ConvexOnlySpline, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor with the spline
       *  @param name  the name
       *  @param title the  title
       *  @param x     the  variable
       *  @param spine the spline
       *  @param phis  vector of parameters
       */
      ConvexOnlySpline
        ( const char*                          name,
          const char*                          title     ,
          RooAbsReal&                          x         ,
          const Ostap::Math::ConvexOnlySpline& spline    ,   // the spline
          RooArgList&                          phis      ) ; // parameters
      /// copy
      ConvexOnlySpline
        ( const ConvexOnlySpline& right     ,
          const char*             name = 0  ) ;
      /// destructor
      virtual ~ConvexOnlySpline() ;
      /// clone
      ConvexOnlySpline* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      ConvexOnlySpline  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::ConvexOnlySpline& function() const { return m_spline ; }
      const Ostap::Math::ConvexOnlySpline& spline  () const { return m_spline ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ConvexOnlySpline m_spline ;          // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class ConvexSpline
     *  The special spline for non-negative monothonic
     *  convex or concave function
     *  @see Ostap::Math::ConvexSpline
     */
    class   ConvexSpline : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::ConvexSpline, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor with the spline
       *  @param name  the name
       *  @param title the  title
       *  @param x     the  variable
       *  @param spine the spline
       *  @param phis  vector of parameters
       */
      ConvexSpline
        ( const char*                      name,
          const char*                      title     ,
          RooAbsReal&                      x         ,
          const Ostap::Math::ConvexSpline& spline    ,   // the spline
          RooArgList&                      phis      ) ; // parameters
      /// copy
      ConvexSpline
        ( const ConvexSpline& right     ,
          const char*         name = 0  ) ;
      /// destructor
      virtual ~ConvexSpline() ;
      /// clone
      ConvexSpline* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      ConvexSpline  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public:  // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      Double_t analyticalIntegral
        ( Int_t          code         ,
          const char*    rangeName    ) const override;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::ConvexSpline& function() const { return m_spline ; }
      const Ostap::Math::ConvexSpline& spline  () const { return m_spline ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ConvexSpline m_spline ;          // the function
      // ======================================================================
    };
    // ========================================================================
  } //                                        end of namespace Ostap::Models
  // ==========================================================================
} //                                                  end of namespace Analysis
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // ANALYSIS_MODELS_H
// ============================================================================
