// ============================================================================
#ifndef OSTAP_PDFS_H
#define OSTAP_PDFS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Peaks.h"
#include "Ostap/BreitWigner.h"
#include "Ostap/Voigt.h"
#include "Ostap/Models.h"
#include "Ostap/BSpline.h"
#include "Ostap/Positive.h"
#include "Ostap/Rational.h"
#include "Ostap/HistoInterpolators.h"
// ============================================================================
// ROOT
// ============================================================================
#include "RVersion.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
// ============================================================================
// forward declarations 
// ============================================================================
class    RooRealVar ; // ROOT,RooFit 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @namespace  Ostap::Models Ostap/PDFs.h
   *
   *  Naturally "wide" models:
   *
   *  - BreitWigner, Rho0, Kstar, Phi, ...
   *  - BreitWigner from 3-body decay of mother particle: BW3L
   *  - LASS (kappa pole)
   *  - Bugg (sigma pole)
   *  - Voigt
   *  - Swanson's S-wave cusp
   *
   *  Empirical resolution models:
   *
   *  - Crystal Ball
   *  - right side Crystal Ball
   *  - double-sided Crystal Ball
   *  - Needham: Crystal Ball with \f$\alpha(\sigma)\f$
   *  - Apollonios    (bifurcated apolonious)
   *  - ApolloniosL   (bifurcated apolonious with power-law tail)
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
   *  - Beta       distribution
   *
   *  Non-factorazeable smooth 2D-models
   *
   *  - generic   positive non-factorizable polynomial in 2D
   *   \f$ P^+(x,y) = \sum_i \sum_j \alpha^2_{i,j} B^n_i(x) B^k_j(y) \f$
   *  - symmetric positive non-factorizable polynomial in 2D \f$ P^+_{sym}(x,y) \f$
   *  - \f$ f(x,y)       = \Phi_1(x)\times\Phi_2(y)\times P^+(x,y)       \f$
   *  - \f$ f_{sym}(x,y) = \Phi  (x)\times\Phi  (y)\times P^+_{sym}(x,y) \f$
   *  - \f$ f(x,y)       = exp   (x)\times\Phi  (y)\times P^+(x,y)       \f$
   *  - \f$ f(x,y)       = exp   (x)\times exp  (y)\times P^+(x,y)       \f$
   *  - \f$ f_{sym}(x,y) = exp   (x)\times exp  (y)\times P^+_{sym}(x,y) \f$
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
      ClassDefOverride(Ostap::Models::BreitWigner, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BreitWigner
      ( const char*            name      ,
        const char*            title     ,
        RooAbsReal&            x         ,
        RooAbsReal&            mass      ,
        RooAbsReal&            width     ,
        const double           m1        ,
        const double           m2        ,
        const unsigned short   L     = 0 ) ;
      /// constructor from all parameters
      BreitWigner 
      ( const char*            name      ,
        const char*            title     ,
        RooAbsReal&            x         ,
        RooAbsReal&            mass      ,
        RooAbsReal&            width     ,
        const double           m1        ,
        const double           m2        ,
        const unsigned short   L                         ,
        const Ostap::Math::FormFactors::JacksonRho rho ) ;
      /// constructor from main parameters and "shape"
      BreitWigner
      ( const char*            name      ,
        const char*            title     ,
        RooAbsReal&            x         ,
        RooAbsReal&            mass      ,
        RooAbsReal&            width     ,
        const Ostap::Math::BW& bw        ) ;
      /// constructor from main parameters and "shape"
      BreitWigner 
      ( const char*            name      ,
        const char*            title     ,
        RooAbsReal&            x         ,
        RooAbsReal&            mass      ,
        RooArgList&            widths    ,
        const Ostap::Math::BW& bw        ) ;
      /// "copy" constructor
      BreitWigner ( const BreitWigner& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~BreitWigner() ;
      /// clone
      BreitWigner* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
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
      virtual void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x   .arg() ; }
      const RooAbsReal& mass   () const { return m_mass.arg() ; }
      const RooArgList& widths () const { return m_widths     ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the Breit-Wigner amplitude
      virtual std::complex<double> amplitude () const  ;
      /// access to underlying function
      const Ostap::Math::BW& function     () const { setPars () ; return *m_bw ; }
      /// access to underlying function
      const Ostap::Math::BW& breitwigner  () const { setPars () ; return *m_bw ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the raw Breit-Wigner amplitude
      std::complex<double> bw_amplitude () const ;
      // ======================================================================
    public: // evaluate the breit wigner 
      // ======================================================================
      /** Get Breit-Wigner \f$ a\f$ : 
       *  \f[ F_a(m) = 2m \varrho(s) N^2_a(s,m_0) 
       *    \frac{\Gamma_{tot}}{\Gamma_{0,a}} \left| \mathcal{A}  \right|^2 \f] 
       *  @param m the mass point 
       *  @param A the amplitide at this point 
       */
      double breit_wigner 
      ( const double                m , 
        const std::complex<double>& A ) const ;
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_mass   ;
      RooListProxy m_widths ;
      // ======================================================================
    protected :
      // ======================================================================
      /// the actual function
      std::unique_ptr<Ostap::Math::BW> m_bw{}  ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BreitWignerMC
     *  @see Ostap::Math::BreitWignerMC
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  BreitWignerMC : public BreitWigner 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::BreitWignerMC, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BreitWignerMC
      ( const char*                       name      ,
        const char*                       title     ,
        RooAbsReal&                       x         ,
        RooAbsReal&                       mass      , 
        RooArgList&                       widths    ,
        const Ostap::Math::BreitWignerMC& bw        ) ;
      /// "copy" constructor
      BreitWignerMC ( const BreitWignerMC& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~BreitWignerMC() ;
      /// clone
      BreitWignerMC* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      BreitWignerMC () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const override ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::BreitWignerMC& breitwigner_MC () const ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BWI 
     *  Breit-Wigner with some embedded interference: 
     *  \f[ A^{\prime}(x) = b (x) s_1(x) \mathrm{e}^{i s_2(x)\theta(x)}
     *                    + A(x)_{\mathrm{BW}} \f], 
     *  where 
     *  - \f$b(x)\f$ smooth function describing the magnitude 
     *               of the coherent background 
     *  - \f$ s_1(x)\f$ optional scale factor/function  
     *  - \f$ s_2(x)\f$ optional scale factor/function  
     *  - \f$ \theta(x)\f$ the phase of the corerent background 
     *  - \f$ A(x)_{\mathrm{BW}} \f$ is Breit-Wigner amplitude 
     */
    class BWI final : public BreitWigner
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::BWI, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the Breit-Wigner
      BWI
      ( const char*                         name      , 
        const char*                         title     ,
        const Ostap::Models::BreitWigner&   bw        ,
        RooAbsReal&                         magnitude ,  
        RooAbsReal&                         phase     ,
        RooAbsReal&                         scale1    ,
        RooAbsReal&                         scale2    ) ;
      /// constructor from the Breit-Wigner
      BWI
      ( const char*                         name      , 
        const char*                         title     , 
        const Ostap::Models::BreitWigner&   bw        ,
        RooAbsReal&                         magnitude ,  
        RooAbsReal&                         phase     ,
        RooAbsReal&                         scale1    ) ;
      /// constructor from the Breit-Wigner
      BWI
      ( const char*                         name      , 
        const char*                         title     , 
        const Ostap::Models::BreitWigner&   bw        ,
        RooAbsReal&                         magnitude ,  
        RooAbsReal&                         phase     ) ;
      /// "copy" constructor
      BWI ( const BWI& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~BWI() ;
      /// clone
      BWI* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for the proper (de)serialization
      BWI () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the amplitude
      std::complex<double> amplitude () const override ;
      // the actual evaluation of function
      Double_t             evaluate  () const override;
      // ======================================================================
    public: // integrals
      // ======================================================================
      Int_t    getAnalyticalIntegral
        ( RooArgSet&     allVars      ,
          RooArgSet&     analVars     ,
          const char* /* rangename */ ) const override;
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& magnitude () const { return m_magnitude .arg() ; }
      const RooAbsReal& phase     () const { return m_phase     .arg() ; }
      const RooAbsReal& scale1    () const { return m_scale1    .arg() ; }
      const RooAbsReal& scale2    () const { return m_scale2    .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy                m_magnitude ;  // background magnitude  
      RooRealProxy                m_phase     ;  // background phase 
      RooRealProxy                m_scale1    ;  // background factor  
      RooRealProxy                m_scale2    ;  // background factor  
      // ======================================================================
    };
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
    class  Flatte : public BreitWigner 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Flatte, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Flatte 
      ( const char*                name   ,
        const char*                title  ,
        RooAbsReal&                x      ,
        RooAbsReal&                m0     ,
        RooAbsReal&                g1     ,
        RooAbsReal&                g2     ,
        RooAbsReal&                g0     ,
        const Ostap::Math::Flatte& flatte ) ;
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
      /// access to underlying function
      const Ostap::Math::Flatte& flatte    () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const override ; // set all parameters
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FlatteBugg
     * Bugg's modification of Flatte channel 
     * @see D.V. Bugg, "Re-analysis of data on a(0)(1450) and a(0)(980)"
     *           Phys.Rev.D 78 (2008) 074023
     *  @see https://doi.org/10.1103/PhysRevD.78.074023
     *  @see https://arxiv.org/abs/0808.2706
     *  Well suitable for \f$f_0(980)\rightarrow \pi^+ \pi^-\f$
     *  \f$\pi\pi\f$-channel
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  FlatteBugg : public BreitWigner 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::FlatteBugg, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      FlatteBugg 
      ( const char*                    name   ,
        const char*                    title  ,
        RooAbsReal&                    x      ,
        RooAbsReal&                    m0     ,
        RooAbsReal&                    g1     ,
        RooAbsReal&                    g2     ,
        RooAbsReal&                    g0     ,
        const Ostap::Math::FlatteBugg& flatte ) ;
      /// "copy" constructor
      FlatteBugg ( const FlatteBugg& , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~FlatteBugg () ;
      /// clone
      FlatteBugg* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      FlatteBugg () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::FlatteBugg& flatte_bugg    () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const override ; // set all parameters
      // ======================================================================
    } ;

    // ========================================================================
    /** @class LASS 
     */
    class LASS : public BreitWigner 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::LASS, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      LASS
      ( const char*                name   ,
        const char*                title  ,
        RooAbsReal&                x      ,
        RooAbsReal&                m0     ,
        RooAbsReal&                g0     ,
        RooAbsReal&                a      ,
        RooAbsReal&                b      ,
        RooAbsReal&                e      , 
        const Ostap::Math::LASS&   lass   ) ;
      /// "copy" constructor
      LASS ( const LASS & , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~LASS () ;
      /// clone
      LASS * clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      LASS  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::LASS& lass () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const override ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& a () const { return m_a .arg() ; }
      const RooAbsReal& b () const { return m_b .arg() ; }
      const RooAbsReal& e () const { return m_e .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// a-parameter 
      RooRealProxy m_a    ; // a-parameter 
      /// b-parameter 
      RooRealProxy m_b    ; // b-parameter 
      /// eliasticity 
      RooRealProxy m_e    ; // eliasticity 
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class BWPS
     *  @see Ostap::Math::BWMC
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  BWPS : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::BWPS , 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BWPS 
      ( const char*                name   ,
        const char*                title  ,
        RooAbsReal&                x      ,
        RooAbsReal&                m0     ,
        RooAbsReal&                gamma  ,
        RooArgList&                phis   , 
        const Ostap::Math::BWPS&   bwps   ) ;
      /// constructor from all parameters
      BWPS
      ( const char*                name   ,
        const char*                title  ,
        RooAbsReal&                x      ,
        RooAbsReal&                m0     ,
        RooArgList&                gamma  ,
        RooArgList&                phis   , 
        const Ostap::Math::BWPS&   bwps   ) ;
      /// "copy" constructor 
      BWPS ( const BWPS& , const char* name = 0 ) ;
      /// virtual destructor 
      virtual ~BWPS() ;
      /// clone method
      BWPS* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// fictive public constructor, needed for  (de)serialization
      BWPS (){} ; //  = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// get the function 
      const Ostap::Math::BWPS& bwps      () const { setPars() ; return m_bwps ; }
      const Ostap::Math::BWPS& function  () const { setPars() ; return m_bwps ; }
      /// get the amplitude 
      std::complex<double>     amplitude () const ;
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
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& m0    () const { return m_m0    .arg() ; }
      const RooArgList& gamma () const { return m_gamma        ; }
      const RooArgList& phis  () const { return m_phis         ; }
      // ======================================================================
    protected: 
      // ======================================================================
      /// the function  itself 
      mutable Ostap::Math::BWPS m_bwps ; // the function  itself 
      // ======================================================================
    protected : 
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_m0    {} ;
      RooListProxy m_gamma {} ;
      RooListProxy m_phis  {} ;
      // ======================================================================
    };
    // ========================================================================
    /** @class BW3L
     *  @see Ostap::Math::BW3L
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  BW3L : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::BW3L , 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BW3L
      ( const char*                name   ,
        const char*                title  ,
        RooAbsReal&                x      ,
        RooAbsReal&                m0     ,
        RooAbsReal&                gamma  ,
        const Ostap::Math::BW3L&   bwps   ) ;
      /// constructor from all parameters
      BW3L
      ( const char*                name   ,
        const char*                title  ,
        RooAbsReal&                x      ,
        RooAbsReal&                m0     ,
        RooArgList&                gamma  ,
        const Ostap::Math::BW3L&   bwps   ) ;
      /// "copy" constructor 
      BW3L ( const BW3L& , const char* name = 0 ) ;
      /// virtual destructor 
      virtual ~BW3L() ;
      /// clone method
      BW3L* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// fictive public constructor, needed for  (de)serialization
      BW3L (){} ; //  = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// get the function 
      const Ostap::Math::BW3L& bw3l      () const { setPars() ; return m_bw3l ; }
      const Ostap::Math::BW3L& function  () const { setPars() ; return m_bw3l ; }
      // ======================================================================      
      /// get the amplitude 
      std::complex<double>     amplitude () const ;
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
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& m0    () const { return m_m0    .arg() ; }
      const RooArgList& gamma () const { return m_gamma        ; }
      // ======================================================================
    protected: 
      // ======================================================================
      /// the function  itself 
      mutable Ostap::Math::BW3L m_bw3l ; // the function  itself 
      // ======================================================================
    protected : 
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_m0    {} ;
      RooListProxy m_gamma {} ;
      // ======================================================================
    };
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
      ClassDefOverride(Ostap::Models::Voigt, 1) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& m0    () const { return m_m0    .arg() ; }
      const RooAbsReal& gamma () const { return m_gamma .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma .arg() ; }
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
      ClassDefOverride(Ostap::Models::PseudoVoigt, 1) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& m0    () const { return m_m0    .arg() ; }
      const RooAbsReal& gamma () const { return m_gamma .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma .arg() ; }
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
//     // ========================================================================
//     /** @class Swanson
//      *  Swanson's S-wau cusp
//      *  @see Ostap::Math::Swanson
//      *  @see LHCb-PAPER-2016-019, Appendix D
//      *  @see E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952,
//      *  @see http://arxiv.org/abs/1504.07952
//      *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
//      *  @date 2016-06-12
//      */
//     class  Swanson : public RooAbsPdf
//     {
//       // ======================================================================
//     public :
//       // ======================================================================
//       ClassDefOverride(Ostap::Models::Swanson, 1) ;
//       // ======================================================================
//     public:
//       // ======================================================================
//       /// constructor from all parameters
//       Swanson
//         ( const char*          name          ,
//           const char*          title         ,
//           RooAbsReal&          x             ,
//           RooAbsReal&          beta0         ,
//           const Ostap::Math::Swanson& sw     ) ;
//       /// constructor from all parameters
//       Swanson
//         ( const char*          name          ,
//           const char*          title         ,
//           RooAbsReal&          x             ,
//           RooAbsReal&          beta0         ,
//           const double         m1_0          ,
//           const double         m2_0          ,
//           const Ostap::Math::BreitWigner& bw ) ;
//       /// "copy" constructor
//       Swanson ( const Swanson& right , const char* name = 0  ) ;
//       /// virtual destructor
//       virtual ~Swanson () ;
//       /// clone
//       Swanson* clone ( const char* name ) const override;
//       // ======================================================================
//     public: // some fake functionality
//       // ======================================================================
//       // fake default contructor, needed just for proper (de)serialization
//       Swanson () {} ;
//       // ======================================================================
//     public:
//       // ======================================================================
//       // the actual evaluation of function
//       Double_t evaluate() const override;
//       // ======================================================================
//     public: // integrals
//       // ======================================================================
//       Int_t    getAnalyticalIntegral
//         ( RooArgSet&     allVars      ,
//           RooArgSet&     analVars     ,
//           const char* /* rangename */ ) const override;
//       Double_t analyticalIntegral
//         ( Int_t          code         ,
//           const char*    rangeName    ) const override;
//       // ======================================================================
//     public:
//       // ======================================================================
//       /// set all parameters
//       void setPars () const ; // set all parameters
//       // ======================================================================
//     public:
//       // ======================================================================
//       /// access to underlying function
//       const Ostap::Math::Swanson& function() const { return m_swanson ; }
//       // ======================================================================
//     protected:
//       // ======================================================================
//       RooRealProxy m_x      ;
//       RooRealProxy m_beta0  ;
//       // ======================================================================
//     private:
//       // ======================================================================
//       /// the actual function
//       mutable Ostap::Math::Swanson m_swanson ;                 // the function
//       // ======================================================================
//     };
//     // ========================================================================

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
      ClassDefOverride(Ostap::Models::CrystalBall, 1) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& m0    () const { return m_m0    .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha .arg() ; }
      const RooAbsReal& n     () const { return m_n     .arg() ; }
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
      ClassDefOverride(Ostap::Models::CrystalBallRS, 1) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& m0    () const { return m_m0    .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha .arg() ; }
      const RooAbsReal& n     () const { return m_n     .arg() ; }
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
      ClassDefOverride(Ostap::Models::CrystalBallDS, 1) ;
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
        RooAbsReal&          alphaL    ,  // 
        RooAbsReal&          nL        ,  //   
        RooAbsReal&          alphaR    ,  // 
        RooAbsReal&          nR        ); //  
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& m0     () const { return m_m0     .arg() ; }
      const RooAbsReal& sigma  () const { return m_sigma  .arg() ; }
      const RooAbsReal& alphaL () const { return m_alphaL .arg() ; }
      const RooAbsReal& alphaR () const { return m_alphaR .arg() ; }
      const RooAbsReal& nL     () const { return m_nL     .arg() ; }
      const RooAbsReal& nR     () const { return m_nR     .arg() ; }
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
     *
     *  @see Ostap::Math::Needham
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-05-13
     */
    class  Needham : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Needham, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Needham
      ( const char*          name        ,
        const char*          title       ,
        RooAbsReal&          x           ,
        RooAbsReal&          m0          ,
        RooAbsReal&          sigma       ,
	const double         c0          ,
	const double         c1          ,
	const double         c2          ,
	const double         n    = 0.00 ,
	const double         amin = 0.01 ) ;
      /// constructor from all parameters
      Needham
      ( const char*          name        ,
        const char*          title       ,
        RooAbsReal&          x           ,
        RooAbsReal&          m0          ,
        RooAbsReal&          sigma       ,
        RooAbsReal&          c0          ,
        RooAbsReal&          c1          ,
        RooAbsReal&          c2          , 
        RooAbsReal&          n           ,
	const double         amin = 0.01 ) ;
      /// constructor from all parameters
      Needham
      ( const char*          name        ,
        const char*          title       ,
        RooAbsReal&          x           ,
        RooAbsReal&          m0          ,
        RooAbsReal&          sigma       ,
        RooAbsReal&          c0          ,
        RooAbsReal&          c1          ,
        RooAbsReal&          c2          , 
	const double         amin = 0.01 ) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg  () ; }
      const RooAbsReal& m0     () const { return m_m0     .arg  () ; }
      const RooAbsReal& sigma  () const { return m_sigma  .arg  () ; }
      const RooAbsReal& c0     () const { return m_c0     .arg  () ; }
      const RooAbsReal& c1     () const { return m_c1     .arg  () ; }
      const RooAbsReal& c2     () const { return m_c2     .arg  () ; }
      const RooAbsReal& n      () const { return m_n      .arg  () ; }
      inline double     amin   () const { return m_needham.amin () ; }      
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_c0     ;
      RooRealProxy m_c1     ;
      RooRealProxy m_c2     ;
      RooRealProxy m_n      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Needham m_needham ;                  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallA
     *  The special parametrization of `Crystal Ball-function' with asymmetric core 
     *  @see Ostap::Math::CrystalBallA
     *  @see Ostap::Math::CrystalBall
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-05-13
     */
    class  CrystalBallA : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::CrystalBallA, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      CrystalBallA
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigmaL    ,
        RooAbsReal&          sigmaR    ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          n         ) ;  // n-1
      /// "copy" constructor
      CrystalBallA ( const CrystalBallA& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~CrystalBallA () ;
      /// clone
      CrystalBallA* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      CrystalBallA () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::CrystalBallA& function() const { return m_cb ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& m0     () const { return m_m0     .arg() ; }
      const RooAbsReal& sigmaL () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR () const { return m_sigmaR .arg() ; }
      const RooAbsReal& alpha  () const { return m_alpha  .arg() ; }
      const RooAbsReal& n      () const { return m_n      .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigmaL ;
      RooRealProxy m_sigmaR ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_n      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CrystalBallA m_cb ;                  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallDSA
     *  double-sided `Crystal Ball-function' with asymmetric core 
     *  for description of gaussian with the tail
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @see Ostap::Math::CrystalBallDoubleSided
     *  @date 2011-05-25
     */
    class  CrystalBallDSA : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::CrystalBallDSA, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      CrystalBallDSA
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigmaL    ,  //
        RooAbsReal&          sigmaR    ,  //
        RooAbsReal&          alphaL    ,  // alpha_L
        RooAbsReal&          nL        ,  //     n_L - 1
        RooAbsReal&          alphaR    ,  // alpha_R - 1
        RooAbsReal&          nR        ); //     n_R
      /// "copy" constructor
      CrystalBallDSA ( const CrystalBallDSA& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~CrystalBallDSA () ;
      /// clone
      CrystalBallDSA* clone ( const char* name ) const override;
      // ======================================================================
    public : // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      CrystalBallDSA() {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::CrystalBallDoubleSidedA& function() const
      { return m_cb2 ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& m0     () const { return m_m0     .arg() ; }
      const RooAbsReal& sigmaL () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR () const { return m_sigmaR .arg() ; }
      const RooAbsReal& alphaL () const { return m_alphaL .arg() ; }
      const RooAbsReal& alphaR () const { return m_alphaR .arg() ; }
      const RooAbsReal& nL     () const { return m_nL     .arg() ; }
      const RooAbsReal& nR     () const { return m_nR     .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigmaL ;
      RooRealProxy m_sigmaR ;
      RooRealProxy m_alphaL ;
      RooRealProxy m_nL     ;
      RooRealProxy m_alphaR ;
      RooRealProxy m_nR     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CrystalBallDoubleSidedA m_cb2 ;       // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CrystalBallDSE
     *  double-sided Crystal Ball-like function: 
     *  - asymmetric core 
     *  - left power-law tail 
     *  - right exponential tail
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @see Ostap::Math::CrystalBallDoubleSided
     *  @date 2011-05-25
     */
    class  CrystalBallDSE : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::CrystalBallDSE, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      CrystalBallDSE
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigmaL    ,  //
        RooAbsReal&          sigmaR    ,  //
        RooAbsReal&          alphaL    ,  //
        RooAbsReal&          nL        ,  //
        RooAbsReal&          alphaR    ) ;
      /// "copy" constructor
      CrystalBallDSE ( const CrystalBallDSE& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~CrystalBallDSE () ;
      /// clone
      CrystalBallDSE* clone ( const char* name ) const override;
      // ======================================================================
    public : // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      CrystalBallDSE () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::CrystalBallDoubleSidedE& function() const
      { return m_cb2 ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& m0     () const { return m_m0     .arg() ; }
      const RooAbsReal& sigmaL () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR () const { return m_sigmaR .arg() ; }
      const RooAbsReal& alphaL () const { return m_alphaL .arg() ; }
      const RooAbsReal& alphaR () const { return m_alphaR .arg() ; }
      const RooAbsReal& nL     () const { return m_nL     .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigmaL ;
      RooRealProxy m_sigmaR ;
      RooRealProxy m_alphaL ;
      RooRealProxy m_nL     ;
      RooRealProxy m_alphaR ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CrystalBallDoubleSidedE m_cb2 ;       // the function
      // ======================================================================
    } ;

    
    // ========================================================================
    /** @class Apollonios
     *  An asymetric Apollonious function 
     *
     *  The function is modiication of the original function proposed by Diego Martinez Santos
     *  @see http://arxiv.org/abs/1312.5000
     *
     *  @see Ostap::Math::Apollonios
     *  @author Vanya BELYAEV Ivane.BElyaev@itep.ru
     *  @date 2013-12-01
     */
    // ========================================================================
    class  Apollonios : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Apollonios, 2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Apollonios
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          m0        ,
        RooAbsReal&          sigmaL    ,
        RooAbsReal&          sigmaR    ,
        RooAbsReal&          beta      ) ; 
      /// "copy" constructor
      Apollonios  ( const Apollonios& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Apollonios () ;
      /// clone
      Apollonios* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Apollonios () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Apollonios& function() const { return m_apo ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& m0     () const { return m_m0     .arg() ; }
      const RooAbsReal& sigmaL () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR () const { return m_sigmaR .arg() ; }
      const RooAbsReal& beta   () const { return m_beta   .arg() ; }
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
      mutable Ostap::Math::Apollonios m_apo ;                // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ApolloniosL
     *  A modified gaussian with exponential
     *  tails on low-side
     *
     *  @see Ostap::Math::Apollonios2
     *  @author Vanya BELYAEV Ivane.BElyaev@itep.ru
     *  @date 2013-12-01
     */
    // ========================================================================
    class  ApolloniosL : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::ApolloniosL, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      ApolloniosL
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mode      ,
        RooAbsReal&          sigmaL    ,
        RooAbsReal&          sigmaR    ,
        RooAbsReal&          beta      ,		
        RooAbsReal&          alpha     ,
        RooAbsReal&          n         ) ;
      /// constructor from all parameters
      ApolloniosL
      ( const char*                name  ,
        const char*                title ,
	Ostap::Models::Apollonios& core  , 
        RooAbsReal&                alpha ,
        RooAbsReal&                n     ) ;
      /// "copy" constructor
      ApolloniosL  ( const ApolloniosL& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~ApolloniosL () ;
      /// clone
      ApolloniosL* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      ApolloniosL () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::ApolloniosL& function() const { return m_apoL ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& m0     () const { return m_m0     .arg() ; }
      const RooAbsReal& sigmaL () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR () const { return m_sigmaR .arg() ; }
      const RooAbsReal& beta   () const { return m_beta   .arg() ; }
      const RooAbsReal& alpha  () const { return m_alpha  .arg() ; }
      const RooAbsReal& n      () const { return m_n      .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_m0     ;
      RooRealProxy m_sigmaL ;
      RooRealProxy m_sigmaR ;
      RooRealProxy m_beta   ;
      RooRealProxy m_alpha  ;
      RooRealProxy m_n      ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ApolloniosL m_apoL ;                // the function
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
      ClassDefOverride(Ostap::Models::BifurcatedGauss, 1) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& peak   () const { return m_peak   .arg() ; }
      const RooAbsReal& sigmaL () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR () const { return m_sigmaR .arg() ; }
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
      ClassDefOverride(Ostap::Models::GenGaussV1, 1) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& mu     () const { return m_mu     .arg() ; }
      const RooAbsReal& alpha  () const { return m_alpha  .arg() ; }
      const RooAbsReal& beta   () const { return m_beta   .arg() ; }
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
      ClassDefOverride(Ostap::Models::GenGaussV2, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& xi     () const { return m_xi     .arg() ; }
      const RooAbsReal& alpha  () const { return m_alpha  .arg() ; }
      const RooAbsReal& kappa  () const { return m_kappa  .arg() ; }
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
      ClassDefOverride(Ostap::Models::SkewGauss, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& xi     () const { return m_xi     .arg() ; }
      const RooAbsReal& omega  () const { return m_omega  .arg() ; }
      const RooAbsReal& alpha  () const { return m_alpha  .arg() ; }
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
    /** @class ExGauss 
     *  Exponentially modified Gaussian function, EMG
     *  @see https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
     *  @see Ostap::Math::ExGauss 
     *  It is a distibution for the varibale that is a 
     *  sum (or difference for negative \f$ k\f$) 
     *  of a Gaussian and exponential variables: \f$ X \sim Y + sign(k) Z \f$,  
     *  where 
     *  - \f$ Y \sim N(\mu,\sigma) \f$
     *  - \f$ Z \sim  \frac{1}{k\sigma}\mathrm{e}^{-\frac{x}{k\sigma}} \f$ 
     *  
     *  For \f$ k=0\f$ one gets a Gaussian distrobution
     *  - \f$ k>0\f$ corresponds to the rigth tail  
     *  - \f$ kM0\f$ corresponds to the left tail  
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
    class ExGauss: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::ExGauss, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      ExGauss
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          varsigma  , 
        RooAbsReal&          k         ) ;
      /// "copy" constructor
      ExGauss ( const ExGauss& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~ExGauss () ;
      /// clone
      ExGauss* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      ExGauss () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::ExGauss& function () const { return m_eg ; }
      const Ostap::Math::ExGauss& exgauss  () const { return m_eg ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& mu       () const { return m_mu       .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma .arg() ; }
      const RooAbsReal& k        () const { return m_k        .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        {} ;
      RooRealProxy m_mu       {} ;
      RooRealProxy m_varsigma {} ;
      RooRealProxy m_k        {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ExGauss m_eg ;                     // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ExGauss2
     *  Variant of 
     *  Exponentially modified Gaussian function
     *  with mu parameter being the mode of th edistribution 
     *  @see Ostap::Math::ExGauss
     *  @see Ostap::Math::ExGauss2
     *  @see Ostap::Math::NormalLaplace 
     *  @see Ostap::Models::ExGauss
     */
    class ExGauss2: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::ExGauss2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      ExGauss2
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          varsigma  , 
        RooAbsReal&          k         ) ;
      /// "copy" constructor
      ExGauss2 ( const ExGauss2& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~ExGauss2 () ;
      /// clone
      ExGauss2* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      ExGauss2 () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::ExGauss2& function () const { return m_eg ; }
      const Ostap::Math::ExGauss2& exgauss2 () const { return m_eg ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& mu       () const { return m_mu       .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma .arg() ; }
      const RooAbsReal& k        () const { return m_k        .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        {} ;
      RooRealProxy m_mu       {} ;
      RooRealProxy m_varsigma {} ;
      RooRealProxy m_k        {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ExGauss2 m_eg ;                     // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bukin2 
     *  Compound PDF - sum of two ExGauss2 pdfs with common mode 
     *  @see Ostap::Math::Bukin2
     *  @see Ostap::Math::ExGauss
     *  @see Ostap::Math::ExGauss2
     *  @see Ostap::Math::NormalLaplace 
     *  @see Ostap::Models::ExGauss
     */
    class Bukin2: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Bukin2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Bukin2
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          varsigmaA , 
        RooAbsReal&          varsigmaB , 
        RooAbsReal&          kA        , 
        RooAbsReal&          kB        , 
        RooAbsReal&          phi       ) ;
      /// "copy" constructor
      Bukin2 ( const Bukin2& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Bukin2 () ;
      /// clone
      Bukin2* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Bukin2 () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Bukin2& function () const { return m_b2 ; }
      const Ostap::Math::Bukin2& bukin2   () const { return m_b2 ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x         () const { return m_x        .arg() ; }
      const RooAbsReal& mu        () const { return m_mu       .arg() ; }
      const RooAbsReal& varsigmaA () const { return m_varsigmaA.arg() ; }
      const RooAbsReal& varsigmaB () const { return m_varsigmaB.arg() ; }
      const RooAbsReal& kA        () const { return m_kA       .arg() ; }
      const RooAbsReal& kB        () const { return m_kB       .arg() ; }
      const RooAbsReal& phi       () const { return m_phi      .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x         {} ;
      RooRealProxy m_mu        {} ;
      RooRealProxy m_varsigmaA {} ;
      RooRealProxy m_varsigmaB {} ;
      RooRealProxy m_kA        {} ;
      RooRealProxy m_kB        {} ;
      RooRealProxy m_phi       {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Bukin2 m_b2 ;                     // the function
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
     *  @see Ostap::Math::mills_normal 
     *  @see Ostap::Math::NormalLaplace
     */    
    class NormalLaplace: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::NormalLaplace, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      NormalLaplace
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          varsigma  , 
        RooAbsReal&          kL        , 
        RooAbsReal&          kR        ) ;
      /// "copy" constructor
      NormalLaplace ( const NormalLaplace& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~NormalLaplace () ;
      /// clone
      NormalLaplace* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      NormalLaplace () {} ;
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
      const Ostap::Math::NormalLaplace& function       () const { return m_nl ; }
      const Ostap::Math::NormalLaplace& normallaplace  () const { return m_nl ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& mu       () const { return m_mu       .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma .arg() ; }
      const RooAbsReal& kL       () const { return m_kL       .arg() ; }
      const RooAbsReal& kR       () const { return m_kR       .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        {} ;
      RooRealProxy m_mu       {} ;
      RooRealProxy m_varsigma {} ;
      RooRealProxy m_kL       {} ;
      RooRealProxy m_kR       {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::NormalLaplace m_nl ;                     // the function
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
     */
    class  Novosibirsk : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Novosibirsk, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Novosibirsk
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          peak      ,   // peak position
        RooAbsReal&          sigma     ,   // "width"
        RooAbsReal&          tau       ) ; // tail/asymmetry pareamneter
      /// "copy" constructor
      Novosibirsk ( const Novosibirsk& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~Novosibirsk () ;
      /// clone
      Novosibirsk * clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Novosibirsk () {} ;
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
      const Ostap::Math::Novosibirsk& function    () const { return m_novosibirsk ; }
      const Ostap::Math::Novosibirsk& novosibirsk () const { return m_novosibirsk ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& peak     () const { return m_peak     .arg() ; }
      const RooAbsReal& sigma    () const { return m_sigma    .arg() ; }
      const RooAbsReal& tau      () const { return m_tau      .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x      {} ;
      RooRealProxy m_peak   {} ;
      RooRealProxy m_sigma  {} ;
      RooRealProxy m_tau    {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Novosibirsk m_novosibirsk ;          // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bukin
     *  "Bukin"-function, aka "Modified Novosibirsk function"
     *  @see http://arxiv.org/abs/1107.5751
     *  @see https://doi.org/10.1007/JHEP06(2012)141
     *  @see Ostap::Math::Bukin
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-12-05
     */
    class  Bukin : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Bukin, 1) ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::Bukin& function () const { return m_bukin ; }
      const Ostap::Math::Bukin& bukin    () const { return m_bukin ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& peak     () const { return m_peak     .arg() ; }
      const RooAbsReal& sigma    () const { return m_sigma    .arg() ; }
      const RooAbsReal& xi       () const { return m_xi       .arg() ; }
      const RooAbsReal& rhoL     () const { return m_rhoL     .arg() ; }
      const RooAbsReal& rhoR     () const { return m_rhoR     .arg() ; }
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
      ClassDefOverride(Ostap::Models::StudentT, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      StudentT
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          sigma     ,
        RooAbsReal&          n         ) ;
      /// "copy constructor"
      StudentT
      ( const StudentT&      right     ,
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x     .arg() ; }
      const RooAbsReal& mu       () const { return m_mu    .arg() ; }
      const RooAbsReal& sigma    () const { return m_sigma .arg() ; }
      const RooAbsReal& n        () const { return m_n     .arg() ; }
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
      ClassDefOverride(Ostap::Models::BifurcatedStudentT, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BifurcatedStudentT
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          sigmaL    ,
        RooAbsReal&          sigmaR    ,
        RooAbsReal&          nL        ,
        RooAbsReal&          nR        ) ;
      /// "copy constructor"
      BifurcatedStudentT
      ( const BifurcatedStudentT& right     ,
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
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
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x      .arg() ; }
      const RooAbsReal& mu       () const { return m_mu     .arg() ; }
      const RooAbsReal& sigmaL   () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR   () const { return m_sigmaR .arg() ; }
      const RooAbsReal& nL       () const { return m_nL     .arg() ; }
      const RooAbsReal& nR       () const { return m_nR     .arg() ; }
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
    /** @class PearsonIV 
     *  Pearson Type IV distribution 
     *   Pearson Type IV distribution  
     *   \f$ f(x;\mu, n, \kappa) = 
     *   C \left( 1 + y^{2}\right)^{-(\frac{1}{2}+n)}
     *   \mathrm{e}^{ -\kappa \atan y }}\f$, where 
     *   - \f$  y = \frac{x-\mu}{\sigma}\f$,
     *   - \f$ 0 < n \f$  
     *  @see https://en.wikipedia.org/wiki/Pearson_distribution
     *  For $\kappa=0\f$ one gets Student's t-distribution
     *  @see Ostap::Math::PearsonIV 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2033-07-10
     */
    class PearsonIV : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::PearsonIV, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PearsonIV 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          varsigma  ,
        RooAbsReal&          n         ,
        RooAbsReal&          kappa     );
      /// "copy" constructor
      PearsonIV ( const PearsonIV& right , const char* name = 0  ) ;
      /// virtual destructor
      virtual ~PearsonIV () ;
      /// clone
      PearsonIV* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PearsonIV () {} ;
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
    public: // max-values (for generation)
      // ======================================================================
      Int_t  getMaxVal ( const RooArgSet& vars ) const override ;
      double maxVal    ( Int_t            code ) const override ;
      // ====================================================================== 
    public:
      // ======================================================================
      /// set all parameters
      void setPars () const ; // set all parameters
      // ======================================================================
    public:
      // ======================================================================
      /// access to underlying function
      const Ostap::Math::PearsonIV& function  () const { return m_p4 ; }
      /// access to underlying function
      const Ostap::Math::PearsonIV& pearsonIV () const { return m_p4 ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& mu       () const { return m_mu       .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma .arg() ; }
      const RooAbsReal& n        () const { return m_n        .arg() ; }
      const RooAbsReal& kappa    () const { return m_kappa    .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        {} ;
      RooRealProxy m_mu       {} ;
      RooRealProxy m_varsigma {} ;
      RooRealProxy m_n        {} ;
      RooRealProxy m_kappa    {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PearsonIV m_p4 ;                // the function
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
      ClassDefOverride(Ostap::Models::GramCharlierA, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& m0       () const { return m_m0       .arg() ; }
      const RooAbsReal& sigma    () const { return m_sigma    .arg() ; }
      const RooAbsReal& kappa3   () const { return m_kappa3   .arg() ; }
      const RooAbsReal& kappa4   () const { return m_kappa4   .arg() ; }
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
      ClassDefOverride(Ostap::Models::PhaseSpace2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpace2 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        const double         m1        ,
        const double         m2        ) ;
      /// "copy constructor"
      PhaseSpace2 
      ( const PhaseSpace2& right     ,
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
    public:
      // ======================================================================
      const RooAbsReal& x  () const { return m_x        .arg() ; }
      double            m1 () const { return m_ps2.m1 () ; }
      double            m2 () const { return m_ps2.m2 () ; }
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
     *  Simple model for left-edge of N-body phase-space
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpaceLeft : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::PhaseSpaceLeft, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpaceLeft 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        RooAbsReal&          scale     ,
        const Ostap::Math::PhaseSpaceLeft& left ) ;
      /// constructor from all parameters
      PhaseSpaceLeft
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        const Ostap::Math::PhaseSpaceLeft& left ) ;      
      /// constructor from all parameters
      PhaseSpaceLeft
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        RooAbsReal&          scale     ,
        const Ostap::Math::PhaseSpace2& left ) ;
      /// constructor from all parameters
      PhaseSpaceLeft 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        const Ostap::Math::PhaseSpace2& left ) ;
      /// constructor from all parameters
      PhaseSpaceLeft 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        RooAbsReal&          scale     ,
        const Ostap::Math::PhaseSpace3& left ) ;
      /// constructor from all parameters
      PhaseSpaceLeft 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        const Ostap::Math::PhaseSpace3& left ) ;
      /// constructor from all parameters
      PhaseSpaceLeft
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        RooAbsReal&          scale     ,
        const Ostap::Math::PhaseSpace3s& left ) ;
      /// constructor from all parameters
      PhaseSpaceLeft 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        const Ostap::Math::PhaseSpace3s& left ) ;
      /// constructor from all parameters
      PhaseSpaceLeft
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        RooAbsReal&          scale     ,
        const unsigned short N         ) ;
      /// constructor from all parameters
      PhaseSpaceLeft 
      ( const char*          name      ,
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
      const Ostap::Math::PhaseSpaceLeft& left    () const { return m_left ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x         () const { return m_x         .arg() ; }
      const RooAbsReal& threshold () const { return m_threshold .arg() ; }
      const RooAbsReal& scale     () const { return m_scale     .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x         ;
      RooRealProxy m_threshold ;
      RooRealProxy m_scale     ;
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
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpaceRight : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::PhaseSpaceRight, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpaceRight 
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          threshold ,
        const unsigned short L         ,
        const unsigned short N         ) ;
      /// "copy constructor"
      PhaseSpaceRight 
      ( const PhaseSpaceRight& right     ,
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
    public:
      // ======================================================================
      const RooAbsReal& x         () const { return m_x         .arg() ; }
      const RooAbsReal& threshold () const { return m_threshold .arg() ; }
      unsigned short    L         () const { return m_right.L()        ; }
      unsigned short    N         () const { return m_right.N()        ; }
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
      ClassDefOverride(Ostap::Models::PhaseSpaceNL, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      PhaseSpaceNL 
      ( const char*          name      ,
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
    public:
      // ======================================================================
      const RooAbsReal& x         () const { return m_x         .arg() ; }
      const RooAbsReal& low       () const { return m_low       .arg() ; }
      const RooAbsReal& high      () const { return m_high      .arg() ; }
      unsigned short    N         () const { return m_ps .N()        ; }
      unsigned short    L         () const { return m_ps .L()        ; }
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
     *  The mass-ditribution of L-particles from N-body phase space decays,
     *  modulate with non-negative polynomial
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-10-05
     */
    class PhaseSpacePol : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::PhaseSpacePol, 1) ;
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
      const Ostap::Math::PhaseSpacePol& pspol   () const { return m_ps ; }
      const Ostap::Math::PhaseSpaceNL&  psNL    () const { return m_ps.phasespace ()  ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x         () const { return m_x   .arg() ; }
      const RooArgList& phis      () const { return m_phis       ; }
      // ======================================================================
   private:
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual phase space function
      mutable Ostap::Math::PhaseSpacePol m_ps ;  // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceLeftExpoPol
     *  The mass-ditribuion of L-particles 
     *  modulate with non-negative polynomial end the  exponent 
     *  \f[ f(x) \propto
     *      \Phi_{l}(x;x_{low}) \mathrm{e}^{-\left|\tau\right| x } P_{N}(x) \f]
     *  where :
     *  -  \f$  \Phi_{l}(x;x_{low}) \f$  is a phase space of 
     *     l-particles near the threshold 
     *  -  \f$ P_{N}(x) \f$ is a positive polynomial of degree N
     *  @see Ostap::Math::PhaseSpaceLeftExpoPol
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2018-10-21
     */
    class PhaseSpaceLeftExpoPol : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::PhaseSpaceLeftExpoPol, 1) ;
      // ======================================================================
      /// constructor from all parameters
      PhaseSpaceLeftExpoPol
      ( const char*                        name   ,
        const char*                        title  ,
        RooRealVar&                        x      ,
        const Ostap::Math::PhaseSpaceLeft& ps     ,
        RooAbsReal&                        tau    ,
        RooAbsReal&                        scale  ,
        RooArgList&                        phis   ) ;
      /// constructor from all parameters
      PhaseSpaceLeftExpoPol
      ( const char*                        name   ,
        const char*                        title  ,
        RooRealVar&                        x      ,
        const Ostap::Math::PhaseSpace2&    ps     ,
        RooAbsReal&                        tau    ,
        RooAbsReal&                        scale  ,
        RooArgList&                        phis   ) ;
      /// constructor from all parameters
      PhaseSpaceLeftExpoPol
      ( const char*                        name   ,
        const char*                        title  ,
        RooRealVar&                        x      ,
        const Ostap::Math::PhaseSpace3&    ps     ,
        RooAbsReal&                        tau    ,
        RooAbsReal&                        scale  ,
        RooArgList&                        phis   ) ;
      /// constructor from all parameters
      PhaseSpaceLeftExpoPol
      ( const char*                        name   ,
        const char*                        title  ,
        RooRealVar&                        x      ,
        const Ostap::Math::PhaseSpace3s&   ps     ,
        RooAbsReal&                        tau    ,
        RooAbsReal&                        scale  ,
        RooArgList&                        phis   ) ;
      /// constructor from all parameters
      PhaseSpaceLeftExpoPol
      ( const char*                        name   ,
        const char*                        title  ,
        RooRealVar&                        x      ,
        const Ostap::Math::PhaseSpaceNL&   ps     ,
        RooAbsReal&                        tau    ,
        RooAbsReal&                        scale  ,
        RooArgList&                        phis   ) ;
      /// constructor from all parameters
      PhaseSpaceLeftExpoPol
      ( const char*                        name   ,
        const char*                        title  ,
        RooRealVar&                        x      ,
        const unsigned short               ps     ,
        RooAbsReal&                        tau    ,
        RooAbsReal&                        scale  ,
        RooArgList&                        phis   ) ;
      // ======================================================================
      // "copy" constructor
      // ======================================================================
      PhaseSpaceLeftExpoPol ( const PhaseSpaceLeftExpoPol& right    ,
                              const char*                  name = 0 ) ;
      /// destructor
      virtual ~PhaseSpaceLeftExpoPol () ;
      /// clone
      PhaseSpaceLeftExpoPol* clone( const char* name )  const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PhaseSpaceLeftExpoPol () {} ;
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
      const Ostap::Math::PhaseSpaceLeftExpoPol& function () const { return m_ps              ; }
      const Ostap::Math::PhaseSpaceLeft&        psleft   () const { return m_ps.phasespace() ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x         () const { return m_x    .arg() ; }
      const RooArgList& phis      () const { return m_phis        ; }
      const RooAbsReal& tau       () const { return m_tau  .arg() ; }
      const RooAbsReal& scale     () const { return m_scale.arg() ; }
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy m_x     ;
      RooListProxy m_phis  ;
      RooRealProxy m_tau   ;
      RooRealProxy m_scale ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual phase space function
      mutable Ostap::Math::PhaseSpaceLeftExpoPol m_ps ;  // the actual function
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
     *   \f$\pi^+\pi^-\f$-mass from \f$B^0\rightarrow J/\psi\pi^+\pi^-\f$ decay.
     *
     *  @see Ostap::Math::PhaseSpace23L
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-04-01
     */
    class PhaseSpace23L : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::PhaseSpace23L, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param dalitz Dalizt plot configuration 
       *  @param L  the angular momentum between the first pair and the third particle
       *  @param l  the angular momentum between the first and the second particle
       */
      PhaseSpace23L 
      ( const char*                       name      ,
        const char*                       title     ,
        RooAbsReal&                       x         ,
        const Ostap::Kinematics::Dalitz&  dalitz    ,  
        const unsigned short              L     = 0 ,
        const unsigned short              l     = 0 ) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x         () const { return m_x    .arg() ; }
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
      ClassDefOverride(Ostap::Models::PolyPositive, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x    .arg()      ; }
      const RooArgList& phis () const { return m_phis             ; }
      double            xmin () const { return m_positive.xmin () ; }
      double            xmax () const { return m_positive.xmax () ; }    
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
      ClassDefOverride(Ostap::Models::PolyPositiveEven, 1) ;
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
      const Ostap::Math::PositiveEven& function () const { return m_even ; }
      const Ostap::Math::PositiveEven& polynom  () const { return m_even ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x    .arg() ; }
      const RooArgList& phis () const { return m_phis        ; }
      double            xmin () const { return m_even.xmin () ; }
      double            xmax () const { return m_even.xmax () ; }    
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
    /** @class PolyMonotonic
     *  positive monotonic  polynomial
     *  @see Ostap::Math::Monotonic
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  PolyMonotonic: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::PolyMonotonic, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      PolyMonotonic
      ( const char*          name       ,
        const char*          title      ,
        RooAbsReal&          x          ,
        const RooArgList&    coeffs     ,
        const double         xmin       ,
        const double         xmax       ,
        const bool           increasing ) ;
      /// copy
      PolyMonotonic
        ( const PolyMonotonic&     right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~PolyMonotonic() ;
      /// clone
      PolyMonotonic* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PolyMonotonic () {} ;
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
      const Ostap::Math::Monotonic& function() const { return m_monotonic ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x    .arg() ; }
      const RooArgList& phis () const { return m_phis        ; }
      double            xmin () const { return m_monotonic.xmin () ; }
      double            xmax () const { return m_monotonic.xmax () ; }    
      bool        increasing () const { return m_monotonic.increasing () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Monotonic m_monotonic ;            // the function
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
      ClassDefOverride(Ostap::Models::PolyConvex, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x    .arg() ; }
      const RooArgList& phis () const { return m_phis        ; }
      double            xmin () const { return m_convex.xmin () ; }
      double            xmax () const { return m_convex.xmax () ; }    
      bool        increasing () const { return m_convex.increasing () ; }
      bool        convex     () const { return m_convex.convex     () ; }
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
      ClassDefOverride(Ostap::Models::PolyConvexOnly, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x    .arg() ; }
      const RooArgList& phis () const { return m_phis        ; }
      double            xmin () const { return m_convex.xmin () ; }
      double            xmax () const { return m_convex.xmax () ; }    
      bool        convex     () const { return m_convex.convex     () ; }
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
      ClassDefOverride(Ostap::Models::ExpoPositive, 1) ;
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
      const Ostap::Math::ExpoPositive& expopol () const { return m_positive ; }
      double xmin () const { return m_positive.xmin () ; }
      double xmax () const { return m_positive.xmax () ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x    .arg() ; }
      const RooAbsReal& tau  () const { return m_tau  .arg() ; }
      const RooArgList& phis () const { return m_phis        ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_tau  ;
      RooListProxy m_phis ;
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
      ClassDefOverride(Ostap::Models::PolySigmoid, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooArgList& phis  () const { return m_phis          ; }
      const RooAbsReal& alpha () const { return m_alpha  .arg() ; }
      const RooAbsReal& x0    () const { return m_x0     .arg() ; }
      double            xmin  () const { return m_sigmoid.xmin () ; }
      double            xmax  () const { return m_sigmoid.xmax () ; }
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
      ClassDefOverride(Ostap::Models::TwoExpoPositive, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha .arg() ; }
      const RooAbsReal& delta () const { return m_delta .arg() ; }
      const RooAbsReal& x0    () const { return m_x0    .arg() ; }
      const RooArgList& phis  () const { return m_phis         ; }
      double            xmin  () const { return m_2expopos.xmin () ; }
      double            xmax  () const { return m_2expopos.xmax () ; }
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

    // ========================================================================
    // Generic Math distributions 
    // ========================================================================
    
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
      ClassDefOverride(Ostap::Models::GammaDist, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GammaDist 
      ( const char*          name      ,
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& k     () const { return m_k     .arg() ; }
      const RooAbsReal& theta () const { return m_theta .arg() ; }
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
      ClassDefOverride(Ostap::Models::GenGammaDist, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& k     () const { return m_k     .arg() ; }
      const RooAbsReal& theta () const { return m_theta .arg() ; }
      const RooAbsReal& p     () const { return m_p     .arg() ; }
      const RooAbsReal& low   () const { return m_low   .arg() ; }
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
      ClassDefOverride(Ostap::Models::Amoroso, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& theta () const { return m_theta .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha .arg() ; }
      const RooAbsReal& beta  () const { return m_beta  .arg() ; }
      const RooAbsReal& a     () const { return m_a     .arg() ; }
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
      ClassDefOverride(Ostap::Models::LogGammaDist, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      LogGammaDist 
      ( const char*          name      ,
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& k     () const { return m_k     .arg() ; }
      const RooAbsReal& theta () const { return m_theta .arg() ; }
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
      ClassDefOverride(Ostap::Models::Log10GammaDist, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Log10GammaDist 
      ( const char*           name      ,
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& k     () const { return m_k     .arg() ; }
      const RooAbsReal& theta () const { return m_theta .arg() ; }
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
      ClassDefOverride(Ostap::Models::LogGamma, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      LogGamma
      ( const char*          name      ,
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooAbsReal& nu    () const { return m_nu     .arg() ; }
      const RooAbsReal& lambd () const { return m_lambda .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha  .arg() ; }
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
    /** @class Beta
     *  http://en.wikipedia.org/wiki/Beta_prime_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2025-08-10
     *  @see Ostap::Math::Beta
     */
    class  Beta : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Beta, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Beta
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           alpha     ,
        RooAbsReal&           beta      ,
        RooAbsReal&           scale     ,
        RooAbsReal&           shift     ) ;
      /// "copy constructor"
      Beta
      ( const Beta&           right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~Beta () ;
      /// clone
      Beta* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Beta () {} ;
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
      const Ostap::Math::Beta& function() const { return m_bfun ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha  .arg() ; }
      const RooAbsReal& beta  () const { return m_beta   .arg() ; }
      const RooAbsReal& scale () const { return m_scale  .arg() ; }
      const RooAbsReal& shift () const { return m_shift  .arg() ; }
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
      mutable Ostap::Math::Beta m_bfun ; // the actual function
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
      ClassDefOverride(Ostap::Models::BetaPrime, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha  .arg() ; }
      const RooAbsReal& beta  () const { return m_beta   .arg() ; }
      const RooAbsReal& scale () const { return m_scale  .arg() ; }
      const RooAbsReal& shift () const { return m_shift  .arg() ; }
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
    /** @class GenBetaPrime
     *  Generalized beta-prime distribution 
     *  http://en.wikipedia.org/wiki/Beta_prime_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2024-11-04
     *  @see Ostap::Math::GenBetaPrime
     */
    class  GenBetaPrime : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::GenBetaPrime, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GenBetaPrime
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           alpha     ,
        RooAbsReal&           beta      ,
        RooAbsReal&           p         ,
        RooAbsReal&           q         ,
        RooAbsReal&           scale     ,
        RooAbsReal&           shift     ) ;
      /// "copy constructor"
      GenBetaPrime
      ( const GenBetaPrime&   right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~GenBetaPrime () ;
      /// clone
      GenBetaPrime* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GenBetaPrime () {} ;
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
      const Ostap::Math::GenBetaPrime& function() const { return m_betap ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha  .arg() ; }
      const RooAbsReal& beta  () const { return m_beta   .arg() ; }
      const RooAbsReal& p     () const { return m_p      .arg() ; }
      const RooAbsReal& q     () const { return m_q      .arg() ; }
      const RooAbsReal& scale () const { return m_scale  .arg() ; }
      const RooAbsReal& shift () const { return m_shift  .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_alpha    ;
      RooRealProxy m_beta     ;
      RooRealProxy m_p        ;
      RooRealProxy m_q        ;      
      RooRealProxy m_scale    ;
      RooRealProxy m_shift    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenBetaPrime m_betap ; // the actual function
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
      ClassDefOverride(Ostap::Models::Landau, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooAbsReal& scale () const { return m_scale  .arg() ; }
      const RooAbsReal& shift () const { return m_shift  .arg() ; }
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
      ClassDefOverride(Ostap::Models::SinhAsinh, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& sigma   () const { return m_sigma   .arg() ; }
      const RooAbsReal& epsilon () const { return m_epsilon .arg() ; }
      const RooAbsReal& delta   () const { return m_delta   .arg() ; }
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
      ClassDefOverride(Ostap::Models::JohnsonSU, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& xi      () const { return m_xi      .arg() ; }
      const RooAbsReal& lambd   () const { return m_lambda  .arg() ; }
      const RooAbsReal& delta   () const { return m_delta   .arg() ; }
      const RooAbsReal& gamma   () const { return m_gamma   .arg() ; }
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
     *  \f$  f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\delta x/2}}}{2})\f$,
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
      ClassDefOverride(Ostap::Models::Atlas, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x     .arg() ; }
      const RooAbsReal& mu      () const { return m_mu    .arg() ; }
      const RooAbsReal& sigma   () const { return m_sigma .arg() ; }
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
     *  \f[ f(x,\mu,\sigma) \propto 
     *   \frac{1}{2} {\mathrm{sech}}
     *   \left( \frac{\pi}{2}\frac{x-\mu}{\sigma} \right)
     *  \f]
     *
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
      ClassDefOverride(Ostap::Models::Sech, 1) ;
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
          const char*          name  = 0 )  ;
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
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& sigma   () const { return m_sigma   .arg() ; }
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
    /** @class Losev
     *  Asymmetric variant of hyperbolic secant distribution
     *  \f[ f(x;\mu,\alpha,\beta) \equiv 
     *   \frac{A}{\mathrm{e}^{-\left|\alpha\right| (x-\mu)} + 
     *                         \mathrm{e}^{\left|\beta\right|(x-mu)}}, \f]
     *  where \f$ A = \frac{\left|\alpha\right|+\left|\beta\right|}{\pi}
     *  \sin \frac{\pi\left| \beta\right| }{\left|\alpha\right|+\left|\beta\right|}\f$ 
     *  - Leptokurtic distribution with exponential tails 
     *  @see Losev, A., "A new lineshape for fitting x‐ray photoelectron peaks", 
     *           Surf. Interface Anal., 14: 845-849. doi:10.1002/sia.740141207
     *  @see  https://doi.org/10.1002/sia.740141207
     *  @see  https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-11-25
     *  @see Ostap::Math::Losev
     */
    class  Losev : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Losev, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Losev
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           mu        ,
        RooAbsReal&           alpha     ,
        RooAbsReal&           beta      ) ;
      /// "copy constructor"
      Losev
      ( const Losev&          right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~Losev () ;
      /// clone
      Losev* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Losev () {} ;
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
      const Ostap::Math::Losev& function() const { return m_losev ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& alpha   () const { return m_alpha   .arg() ; }
      const RooAbsReal& beta    () const { return m_beta    .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_alpha    ;
      RooRealProxy m_beta     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Losev m_losev ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Logistic
     *  aka "Sech-square"
     *  \f[ f(x;\mu;s) = \frac{1}{4s}sech^2\left(\frac{x-\mu}{2s}\right)\f]
     *  where
     *  - \f$  s = \sigma \frac{\sqrt{3}}{\pi}\f$
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
      ClassDefOverride(Ostap::Models::Logistic, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& sigma   () const { return m_sigma   .arg() ; }
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
    /** @class GenLogisticIV
     *  - Type I   : beta  = 1 
     *  - Type II  : alpha = 1 
     *  - Type III : alpha = beta         
     *  @see Ostap::Math::GenLogisticIV
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-14
     */
    class  GenLogisticIV: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::GenLogisticIV, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GenLogisticIV
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           mu        ,
        RooAbsReal&           sigma     , 
        RooAbsReal&           alpha     ,
        RooAbsReal&           beta      ) ;
      /// "copy constructor"
      GenLogisticIV
      ( const GenLogisticIV&  right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~GenLogisticIV () ;
      /// clone
      GenLogisticIV* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GenLogisticIV () {} ;
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
      const Ostap::Math::GenLogisticIV& function() const { return m_gl4 ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& sigma   () const { return m_sigma   .arg() ; }
      const RooAbsReal& alpha   () const { return m_alpha   .arg() ; }
      const RooAbsReal& beta    () const { return m_beta    .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_mu       ;
      RooRealProxy m_sigma    ;
      RooRealProxy m_alpha    ;
      RooRealProxy m_beta     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenLogisticIV m_gl4 ; // the actual function
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
      ClassDefOverride(Ostap::Models::Argus, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Argus
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           mu        ,
        RooAbsReal&           c         ,
        RooAbsReal&           chi       ) ;
      /// constructor from all parameters
      Argus
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           c         ,
        RooAbsReal&           chi       ) ;
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
      const Ostap::Math::Argus& function () const { return m_argus ; }
      const Ostap::Math::Argus& argus    () const { return m_argus ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& c       () const { return m_c       .arg() ; }
      const RooAbsReal& chi     () const { return m_chi     .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x   {} ;
      RooRealProxy m_mu  {} ;
      RooRealProxy m_c   {} ;
      RooRealProxy m_chi {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Argus m_argus ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenArgus
     *  Generelised Argus distribution
     *  http://en.wikipedia.org/wiki/ARGUS_distribution
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-06-07
     *  @see Ostap::Math::GenArgus
     *  @see Ostap::Models::Argus
     *  @see Ostap::Math::Argus
     */
    class GenArgus : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::GenArgus, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      GenArgus
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           mu        ,
        RooAbsReal&           c         ,
        RooAbsReal&           chi       ,
        RooAbsReal&           dp        ) ;
      /// constructor from all parameters
      GenArgus
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           c         ,
        RooAbsReal&           chi       , 
        RooAbsReal&           dp      ) ;
      /// "copy constructor"
      GenArgus
      ( const GenArgus&       right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~GenArgus  () ;
      /// clone
      GenArgus* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GenArgus () {} ;
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
      const Ostap::Math::GenArgus& function () const { return m_argus ; }
      const Ostap::Math::GenArgus& argus    () const { return m_argus ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& c       () const { return m_c       .arg() ; }
      const RooAbsReal& chi     () const { return m_chi     .arg() ; }
      const RooAbsReal& dp      () const { return m_dp      .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x   {} ;
      RooRealProxy m_mu  {} ;
      RooRealProxy m_c   {} ;
      RooRealProxy m_chi {} ;
      RooRealProxy m_dp  {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenArgus m_argus ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  Slash 
     *  ``Slash''-distribution -  symmetric peak with very heavy tails
     *  @see https://en.wikipedia.org/wiki/Slash_distribution
     *  Tails arew so heavy that moments (e.g. variance) do not exist 
     *  @see Ostap::Math::Slash
     */
    class Slash : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Slash, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& scale   () const { return m_scale   .arg() ; }
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
      ClassDefOverride(Ostap::Models::AsymmetricLaplace, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x       .arg() ; }
      const RooAbsReal& mu      () const { return m_mu      .arg() ; }
      const RooAbsReal& lambdaL () const { return m_lambdaL .arg() ; }
      const RooAbsReal& lambdaR () const { return m_lambdaR .arg() ; }
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
    /** @class BatesShape
     *  Modified Bates distribution such that it has mean of \f$ \mu \f$
     *  and rms of \f$ \sigma \f$, \f$ n>0\f$ is just a shape parameters 
     *  @see https://en.wikipedia.org/wiki/Bates_distribution
     *  Essentially it is a scaled version of Irwin-Hall distribution 
     *  @see https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution
     *  @see Ostap::Math::BatesShape  
     *  @see Ostap::Math::Bates  
     *  @see Ostap::Math::IrwinHall 
     */
    class BatesShape : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::BatesShape, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      BatesShape
      ( const char*          name  ,
        const char*          title ,
        RooAbsReal&          x     ,
        RooAbsReal&          mu    ,
        RooAbsReal&          sigma , 
        const unsigned short n     ) ;
      /// "copy constructor"
      BatesShape       
      ( const BatesShape&    right     ,
        const char*          name  = 0 )  ;
      /// destructor
      virtual ~BatesShape() ;
      /// clone
      BatesShape* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      BatesShape () {} ;
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
      const Ostap::Math::BatesShape& function   () const { return m_bs ; }
      const Ostap::Math::BatesShape& batesshape () const { return m_bs ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x       .arg () ; }
      const RooAbsReal& mu       () const { return m_mu      .arg () ; }
      const RooAbsReal& sigma    () const { return m_sigma   .arg () ; }
      unsigned short    n        () const { return m_bs      .n   () ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_mu    {} ;
      RooRealProxy m_sigma {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::BatesShape m_bs ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Hat
     *  Finite smooth functon
     *  \f$ f(x;\mu\sigma) = \frac{C}{\sigma} 
     *   \mathrm{e}^{  - frac{1}{1-y^2} }\f$, where 
     *  \F$ y = \frac{m-\mu}{\sigma}\f$ 
     *  @see Ostap::Math::hat 
     *  @see Ostap::Models::Hat 
     */
    class Hat : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Hat, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Hat
      ( const char* name      ,
        const char* title     ,
        RooAbsReal& x         ,
        RooAbsReal& mu        ,
        RooAbsReal& varsigma  ) ;
      /// "copy constructor"
      Hat 
      ( const Hat&  right     ,
        const char* name  = 0 )  ;
      /// destructor
      virtual ~Hat() ;
      /// clone
      Hat* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Hat () {} ;
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
      const Ostap::Math::Hat& function () const { return m_hat ; }
      const Ostap::Math::Hat& hat      () const { return m_hat ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x       .arg() ; }
      const RooAbsReal& mu       () const { return m_mu      .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma.arg() ; }
      const RooAbsReal& sigma    () const { return m_varsigma .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        {} ;
      RooRealProxy m_mu       {} ;
      RooRealProxy m_varsigma {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Hat m_hat ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Up 
     *  Finite stomic functon <code>up</code>,  a finite soltuion  
     *  of the equation 
     *  \f[ f^{\prime(x) = 2 \left( f( 2x+1) - f(2x-1)\right) }\f] with
     *  \f$ f(0) = 1 \f$ 
     *  @see Ostap::Math::up_F 
     *  @see Ostap::Models::Up
     */
    class Up : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Up, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Up
      ( const char* name      ,
        const char* title     ,
        RooAbsReal& x         ,
        RooAbsReal& mu        ,
        RooAbsReal& varsigma  ) ;
      /// "copy constructor"
      Up 
      ( const Up&   right     ,
        const char* name  = 0 )  ;
      /// destructor
      virtual ~Up() ;
      /// clone
      Up* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Up () {} ;
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
      const Ostap::Math::Up& function () const { return m_up ; }
      const Ostap::Math::Up& up       () const { return m_up ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& mu       () const { return m_mu       .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma .arg() ; }
      const RooAbsReal& sigma    () const { return m_varsigma .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        {} ;
      RooRealProxy m_mu       {} ;
      RooRealProxy m_varsigma {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Up m_up ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FupN 
     *  Finite stomic functon <code>fup_N</code>
     *  @see Ostap::Math::fupN_F 
     *  @see Ostap::Models::FupN
     */
    class FupN : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::FupN, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      FupN
      ( const char*          name     ,
        const char*          title    ,
        RooAbsReal&          x        ,
        const unsigned short N        , 
        RooAbsReal&          mu       ,
        RooAbsReal&          varsigma ) ;
      /// "copy constructor"
      FupN 
      ( const FupN& right     ,
        const char* name  = 0 )  ;
      /// destructor
      virtual ~FupN() ;
      /// clone
      FupN* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      FupN () {} ;
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
      const Ostap::Math::FupN& function () const { return m_fupN ; }
      const Ostap::Math::FupN& fupN     () const { return m_fupN ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& mu       () const { return m_mu       .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma .arg() ; }
      unsigned short    N        () const { return m_fupN.N ()       ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        {} ;
      RooRealProxy m_mu       {} ;
      RooRealProxy m_varsigma {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::FupN m_fupN ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Tsallis
     *  Useful function to describe pT-spectra of particles
     *
     *  @see C. Tsallis, ``Possible generalization of Boltzmann-Gibbs statistics'',
     *                   J. Statist. Phys. 52 (1988) 479.
     *
     *  @see C. Tsallis, ``Nonextensive statistics: theoretical, experimental and computational
     *  evidences and connections'', Braz. J. Phys. 29 (1999) 1.
     *
     *  \f[ \frac{d\sigma}{dp_T} \propto
     *    \frac { p_T}{\left( 1 + \frac{E_{kin}}{Tn}\right)^{n}}\f],
     *  where 
     *  - \f$E_{kin} = \sqrt{p_T^2-M^2}-M\f$ is transverse kinetic energy
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
      ClassDefOverride(Ostap::Models::Tsallis, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x        .arg() ; }
      const RooAbsReal& pt   () const { return m_x        .arg() ; }
      const RooAbsReal& n    () const { return m_n        .arg() ; }
      const RooAbsReal& T    () const { return m_T        .arg() ; }
      const RooAbsReal& mass () const { return m_mass     .arg() ; }
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
      ClassDefOverride(Ostap::Models::QGSM, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x        .arg() ; }
      const RooAbsReal& b    () const { return m_b        .arg() ; }
      const RooAbsReal& mass () const { return m_mass     .arg() ; }
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
    /** @class Hagedorn 
     *  Useful function to describe pT-spectra of particles
     *
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
     *  @see Ostap::Math::Hagedorn 
     *  @date 2022-12-06
     */
    class Hagedorn : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Hagedorn, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      Hagedorn
      ( const char*           name      ,
        const char*           title     ,
        RooAbsReal&           x         ,
        RooAbsReal&           beta      ,   // parameter beta 
        RooAbsReal&           mass      ) ; // particle mass (fixed!)
      /// "copy constructor"
      Hagedorn 
      ( const Hagedorn&       right     ,
        const char*           name  = 0 )  ;
      /// destructor
      virtual ~Hagedorn () ;
      /// clone
      Hagedorn* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Hagedorn () {} ;
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
      const Ostap::Math::Hagedorn& function () const { return m_hagedorn ; }
      const Ostap::Math::Hagedorn& hagedorn () const { return m_hagedorn ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x        .arg() ; }
      const RooAbsReal& beta () const { return m_beta     .arg() ; }
      const RooAbsReal& mass () const { return m_mass     .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x        ;
      RooRealProxy m_beta     ;
      RooRealProxy m_mass     ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Hagedorn m_hagedorn ; // the actual function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class TwoExpos
     *  simple difference of two exponents
     *  \f[ f \propto
     *        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} =
     *        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f]
     *  @author Vanya BElyaev Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  TwoExpos : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::TwoExpos, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& alpha () const { return m_alpha .arg() ; }
      const RooAbsReal& delta () const { return m_delta .arg() ; }
      const RooAbsReal& x0    () const { return m_x0    .arg() ; }
      // ======================================================================
    protected:
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_alpha {} ;
      RooRealProxy m_delta {} ;
      RooRealProxy m_x0    {} ;
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
      ClassDefOverride(Ostap::Models::DoubleGauss, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      DoubleGauss 
      ( const char*          name      , 
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
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& sigma    () const { return m_sigma    .arg() ; }
      const RooAbsReal& fraction () const { return m_fraction .arg() ; }
      const RooAbsReal& scale    () const { return m_scale    .arg() ; }
      const RooAbsReal& mean     () const { return m_mean     .arg() ; }
      // ======================================================================
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
     *  \f[  G(x;\mu,\beta) = \frac{1}{\left|\beta\right|} e^{-(z+e^{-z}}\f], 
     *  where 
     *  - \f$ z = \frac{x-\mu}{\beta}\f$.
     *  Important  cases if \f$ E(x) = e^{-\tau x}\f$, and:
     *  - \f$ z \equiv  \log(x)\f$, then \f$ F(z) = E(x) = G(z, -log(\tau) , 1 ) \f$, 
     *  - \f$ z \equiv -\log(x)\f$, then \f$ F(z) = E(x) = G(z, -log(\tau) , 1 ) \f$.
     *
     *  As a direct sequence,  a sum of exponential componets is transformed to 
     *  a sum of ``peak-like'' Gumbel  stuctures
     *  @see Ostap::Math::Gumbel
     */
    class  Gumbel : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Gumbel, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  z    the  variable 
       *  @param  eta  the shape eta-parameter 
       *  @param  b    the scale b-parameter 
       *  @param  xmin the bias  parameter
       */
      Gumbel
      ( const char*          name      , 
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x    .arg() ; }
      const RooAbsReal& mu    () const { return m_mu   .arg() ; }
      const RooAbsReal& beta  () const { return m_beta .arg() ; }
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
     *  \f[ f(x,\lambda,k,x_0) = \frac{k}{\lambda}  y^{k-1} e^{-y^k}\f]
     *  where 
     *  - \f$ y \equiv \frac{x-x_0}{\lambda}\f$
     *  @see https://en.wikipedia.org/wiki/Weibull_distribution
     *  @see Ostap::Math::Weibull
     */
    class Weibull : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Weibull, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  x    the  variable 
       *  @param  scale the scale parameter
       *  @param  shape the shape parameter 
       *  @param  shift the shift parameter 
       */
      Weibull 
      ( const char*          name      , 
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& scale () const { return m_scale .arg() ; }
      const RooAbsReal& shape () const { return m_shape .arg() ; }
      const RooAbsReal& shift () const { return m_shift .arg() ; }
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
     *  \f[ f(x,\mu,s) = \frac{1}{2s}   \left( 1   +\cos \pi y \right)  \f], 
     *  where 
     *  - \f$  y  \equiv = \frac{x-\mu}{s}\f$ 
     *  @see https://en.wikipedia.org/wiki/Raised_cosine_distribution
     *  @see Ostap::Math::RaisngCosine 
     */
    class RaisingCosine: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::RaisingCosine, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  x      the variable 
       *  @param  mean   the mean/mode/median/location 
       *  @param  scale  the scale parameter 
       */
      RaisingCosine 
      ( const char*          name      , 
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& mean  () const { return m_mean  .arg() ; }
      const RooAbsReal& scale () const { return m_scale .arg() ; }
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
     *  \f[ f(x) = \frac{ \sqrt{\beta}}{C_q} e_q (-\beta (x-\mu)^2) \f], 
     *  where  
     *  - \f$ e_q (x) = \left( 1 + (1-q)x\right)^{\frac{1}{1-q}}\f$ 
     *  @see https://en.wikipedia.org/wiki/Q-Gaussian_distribution
     *  If is equal to 
     *  - scaled version of Student' t-distribution for 1<q<3
     *  - Gaussian distribution for q = 1 
     *  - Cauchy   distribution for q = 2 
     *  - has finite  support for q<1 
     *  @see Ostap::Math::QGaussian
     *  Here we use \f$ \beta = \frac{1}{2\sigma^2}\f$
     */
    class QGaussian: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::QGaussian, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  x      the variable 
       *  @param  mean   the mean/mode/median/location 
       *  @param  scale  the scale parameter 
       *  @param  q      the q-value 
       */
      QGaussian
      ( const char*          name      , 
        const char*          title     ,
        RooAbsReal&          x         ,   // observable 
        RooAbsReal&          mean      ,   // mean/mode/location
        RooAbsReal&          scale     ,   // scale parameter/sigma  
        RooAbsReal&          q         ) ; // q-parameter/shape 
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
      const Ostap::Math::QGaussian& qgauss   () const { return m_qgauss ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& mean  () const { return m_mean  .arg() ; }
      const RooAbsReal& scale () const { return m_scale .arg() ; }
      const RooAbsReal& q     () const { return m_q     .arg() ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mean  ;
      RooRealProxy m_scale ;
      RooRealProxy m_q     ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::QGaussian m_qgauss ;
      // =====================================================================
    } ;
    // ========================================================================



    // ========================================================================
    /** @class KGaussian
     *  k-Gaussian distribution:
     *  @see Ostap::Math::KGaussian
     *  @see https://en.wikipedia.org/wiki/Kaniadakis_Gaussian_distribution
     *  Here we use \f$ k = \tanh {\kappa} \f$
     */
    class KGaussian: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::KGaussian, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param  x      the variable 
       *  @param  mean   the mean/mode/median/location 
       *  @param  scale  the scale parameter 
       *  @param  kappa  kappa parameter \f$ k = \tanh \kappa \f$
       */
      KGaussian
      ( const char*          name      , 
        const char*          title     ,
        RooAbsReal&          x         ,   // observable 
        RooAbsReal&          mean      ,   // mean/mode/location
        RooAbsReal&          scale     ,   // scale parameter/sigma  
        RooAbsReal&          kappa     ) ; // kappa-parameter/shape 
      /// "copy" constructor 
      KGaussian ( const KGaussian& , const char* name = 0 ) ;
      /// clone 
      KGaussian* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      KGaussian () {} ;
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
      const Ostap::Math::KGaussian& function () const { return m_kgauss ; }
      const Ostap::Math::KGaussian& kgauss   () const { return m_kgauss ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& mean  () const { return m_mean  .arg() ; }
      const RooAbsReal& scale () const { return m_scale .arg() ; }
      const RooAbsReal& kappa () const { return m_kappa .arg() ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mean  ;
      RooRealProxy m_scale ;
      RooRealProxy m_kappa ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::KGaussian m_kgauss ;
      // =====================================================================
    } ;
    // ========================================================================

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
     *  - location parameter \f$\mu\f$
     *  - parameter               \f$\sigma \gt  0 \f$, related to the width;
     *  - dimensionless parameter \f$\kappa\f$,         related to the asymmetry;
     *  - dimensionless parameter \f$\zeta   \ge 0 \f$, related to the kurtosis 
     *
     * The parameters are defined as:
     * \f[\begin{array}{lcl}
     *     \sigma^2 & \equiv & \gamma^{-2} \zeta \frac{K_2(\zeta)}{\zetaK_1(zeta)} \\
     *     \kappa   & \equiv & \frac{\beta}{\sigma} \                   \
     *     \zeta\equiv\delta \gamma \end{array} \f]
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
     *     \frac{ A^*(\zeta) } { 2 \sigma \sqrt{\kappa^2+A^*(\zeta)} \zeta K^*_1(\zeta) } 
     *     \mathrm{e}^{\zeta - \sqrt{ (\kappa^2+A(\zeta))  \left( \frac{\zeta^2}{A(\zeta)}  +  
     *      \left( \frac{x-\mu}{\sigma}\right)^2  \right) } } 
     *  \f]
     *  where \f$ K^*_n(x)\f$ is a scaled modified Bessel functon to th eseodn kind 
     *   \f$ K^*_n(x) = \mathrm{e}^{x}K_1(x) \f$ 
     *
     *  In all expressions \f$ \left| \sigma \right|\f$ and 
     *  \f$ \left| \zeta \right|\f$ are used instead of \f$\sigma\f$ and \f$\zeta\f$ 
     * 
     *  @see Ostap::Math::Hyperbolic
     */
    class Hyperbolic: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Hyperbolic, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param name  name of PDF
       *  @param title name of PDF
       *  @param x     observable 
       *  @param mu    related to location 
       *  @param sigma related to width
       *  @param zeta  related to shape 
       *  @param kappa related to asymmetry 
       */
      Hyperbolic 
      ( const char*          name  , 
        const char*          title ,
        RooAbsReal&          x     ,   // observable 
        RooAbsReal&          mu    ,   // location
        RooAbsReal&          sigma ,   // wodth 
        RooAbsReal&          zeta  ,   // zeta 
        RooAbsReal&          kappa ) ; // kappa 
      /// "copy" constructor 
      Hyperbolic ( const Hyperbolic& , const char* name = 0 ) ;
      /// clone 
      Hyperbolic* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      Hyperbolic () {} ;
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
      const Ostap::Math::Hyperbolic& function () const { return m_hyperbolic ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& mu    () const { return m_mu    .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma .arg() ; }
      const RooAbsReal& zeta  () const { return m_zeta  .arg() ; }
      const RooAbsReal& kappa () const { return m_kappa .arg() ; }
      // ======================================================================
    public: // canonical parameters 
      // ======================================================================
      double alpha () const ;
      double beta  () const ;
      double gamma () const ;      
      double delta () const ;      
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_mu    ;
      RooRealProxy m_sigma ;
      RooRealProxy m_zeta  ;
      RooRealProxy m_kappa ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::Hyperbolic m_hyperbolic ;
      // =====================================================================
    } ;
    // ========================================================================

    // ========================================================================
    /** @class GenHyperbolic 
     *  Generalized Hyperboilic distribution
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
     * @see Ostap::Math::Hyperbolic
     *  Useful subclasses 
     *  - \f$ \lambda=1\f$ : Hyperbolic distribution
     *  - \f$ \lambda=-\frac{1}{2}\f$ : Normal Inverse Gaussian distributtion
     *  - \f$ \lambda=-\frac{n}{2}, \zeta\rightarrow+0\f$ : Student's t-distibution 
     *  - \f$ \lambda \rightarrow \pm\infty, \kappa=0\f$ : Gaussian distribution 
     *  - \f$ \zeta \rightarrow +\infty, \kappa=0\f$ : Gaussian distribution 
     * 
     *  @see Ostap::Math::Hyperbolic
     *  @see Ostap::Math::GenHyperbolic
     *  @see Ostap::Models::Hyperbolic
     */
    class GenHyperbolic: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::GenHyperbolic, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param name  name of PDF
       *  @param title name of PDF
       *  @param x      observable 
       *  @param mu     related to location 
       *  @param sigma  related to width
       *  @param zeta   related to shape 
       *  @param kappa  related to asymmetry 
       *  @param lambda related to shape 
       */
      GenHyperbolic 
      ( const char*          name   , 
        const char*          title  ,
        RooAbsReal&          x      ,   // observable 
        RooAbsReal&          mu     ,   // location
        RooAbsReal&          sigma  ,   // width 
        RooAbsReal&          zeta   ,   // zeta 
        RooAbsReal&          kappa  ,   // kappa 
        RooAbsReal&          lambda ) ; // lambda 
      /// "copy" constructor 
      GenHyperbolic ( const GenHyperbolic& , const char* name = 0 ) ;
      /// clone 
      GenHyperbolic* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      GenHyperbolic () {} ;
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
      const Ostap::Math::GenHyperbolic& function () const { return m_hyperbolic ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooAbsReal& mu    () const { return m_mu     .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma  .arg() ; }
      const RooAbsReal& zeta  () const { return m_zeta   .arg() ; }
      const RooAbsReal& kappa () const { return m_kappa  .arg() ; }
      const RooAbsReal& lambd () const { return m_lambda .arg() ; }
      // ======================================================================
    public: // canonical parameters 
      // ======================================================================
      double alpha () const ;
      double beta  () const ;
      double gamma () const ;      
      double delta () const ;      
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_mu     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_zeta   ;
      RooRealProxy m_kappa  ;
      RooRealProxy m_lambda ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::GenHyperbolic m_hyperbolic ;
      // =====================================================================
    } ;
    // ========================================================================
    
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
     * @see Ostap::Math::Das 
     */
    class Das: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Das, 2) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param name  name of PDF
       *  @param title name of PDF
       *  @param x      observable 
       *  @param mu     location 
       *  @param sigma  width for Gaussian core 
       *  @param alphaL left tail 
       *  @param alphaR right tail 
       */
      Das 
      ( const char*          name   , 
        const char*          title  ,
        RooAbsReal&          x      ,   // observable 
        RooAbsReal&          mu     ,   // location
        RooAbsReal&          sigma  ,   // width 
        RooAbsReal&          alphaL ,   // left tail 
        RooAbsReal&          alphaR ) ; // right tail 
      /** constructor from all parameters  (symmetric)
       *  @param name  name of PDF
       *  @param title name of PDF
       *  @param x      observable 
       *  @param mu     location 
       *  @param sigma  width for Gaussian core 
       *  @param alpha  left tail == rigth tail
       */
      Das 
      ( const char*          name   , 
        const char*          title  ,
        RooAbsReal&          x      ,   // observable 
        RooAbsReal&          mu     ,   // location
        RooAbsReal&          sigma  ,   // width  
        RooAbsReal&          alpha  ) ; // left & right tails 
      /// "copy" constructor 
      Das ( const Das& , const char* name = 0 ) ;
      /// clone 
      Das* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      Das () {} ;
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
      const Ostap::Math::Das& function () const { return m_das ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& mu     () const { return m_mu     .arg() ; }
      const RooAbsReal& sigma  () const { return m_sigma  .arg() ; }
      const RooAbsReal& alphaL () const { return m_alphaL .arg() ; }
      const RooAbsReal& alphaR () const { return m_alphaR .arg() ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_mu     ;
      RooRealProxy m_sigma  ;
      RooRealProxy m_alphaL ;
      RooRealProxy m_alphaR ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::Das m_das ;
      // =====================================================================
    };
   // ========================================================================
    /** @class ADas
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
     * @see Ostap::Math::Das 
     */
    class ADas: public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::ADas, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param name  name of PDF
       *  @param title name of PDF
       *  @param x      observable 
       *  @param mu     location 
       *  @param sigmaL left width for Gaussian core 
       *  @param sigmaR right width for Gaussian core 
       *  @param alphaL left tail 
       *  @param alphaR  right tail 
       */
      ADas 
      ( const char*          name   , 
        const char*          title  ,
        RooAbsReal&          x      ,   // observable 
        RooAbsReal&          mu     ,   // location
        RooAbsReal&          sigmaL ,   // width 
        RooAbsReal&          sigmaR ,   // width 
        RooAbsReal&          alphaL ,   // left tail 
        RooAbsReal&          alphaR ) ; // right tail 
      /** constructor from all parameters  (symmetric)
       *  @param name  name of PDF
       *  @param title name of PDF
       *  @param x      observable 
       *  @param mu     location 
       *  @param sigma  width for Gaussian core 
       *  @param kL     left tail == rigth tail
       */
      ADas 
      ( const char*          name   , 
        const char*          title  ,
        RooAbsReal&          x      ,   // observable 
        RooAbsReal&          mu     ,   // location
        RooAbsReal&          sigma  ,   // width  
        RooAbsReal&          alpha  ) ; // left & right tails 
      /// "copy" constructor 
      ADas ( const ADas& , const char* name = 0 ) ;
      /// clone 
      ADas* clone ( const char* name ) const override ; 
      // =====================================================================
    public: // some fake functionality
      // =====================================================================
      // fake default contructor, needed just for proper (de)serialization 
      ADas () {} ;
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
      const Ostap::Math::ADas& function () const { return m_adas ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x      .arg() ; }
      const RooAbsReal& mu     () const { return m_mu     .arg() ; }
      const RooAbsReal& sigmaL () const { return m_sigmaL .arg() ; }
      const RooAbsReal& sigmaR () const { return m_sigmaR .arg() ; }
      const RooAbsReal& alphaL () const { return m_alphaL .arg() ; }
      const RooAbsReal& alphaR () const { return m_alphaR .arg() ; }
      // ======================================================================
    protected:
      // =====================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_mu     ;
      RooRealProxy m_sigmaL ;
      RooRealProxy m_sigmaR ;
      RooRealProxy m_alphaL ;
      RooRealProxy m_alphaR ;
      // =====================================================================
    protected : // the function itself 
      // =====================================================================
      mutable Ostap::Math::ADas m_adas ;
      // =====================================================================
    };



    // ========================================================================
    // 1D-splines
    // ========================================================================

    // ========================================================================
    /** @class PositiveSpline
     *  The special spline for non-negative function
     *  Actually it is a sum of M-splines with non-negative coefficients
     *  \f[ f(x) = \sum_i \alpha_i * M_i^k(x) \f],
     *  with constraints \f$  \sum_i \alpha_i=1\f$
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
      ClassDefOverride(Ostap::Models::PositiveSpline, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooArgList& phis  () const { return m_phis          ; }
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
    /** @class MonotonicSpline
     *  The special spline for non-negative monotonic function
     *  @see http://en.wikipedia.org/wiki/I-spline
     *  @see http://en.wikipedia.org/wiki/M-spline
     *  @see http://en.wikipedia.org/wiki/B-spline
     *  @see Ostap::Math::PositiveSpline
     */
    class   MonotonicSpline : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::MonotonicSpline, 1) ;
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
      MonotonicSpline
        ( const char*                          name,
          const char*                          title     ,
          RooAbsReal&                          x         ,
          const Ostap::Math::MonotonicSpline&  spline    ,   // the spline
          RooArgList&                          phis      ) ; // parameters
      /// copy
      MonotonicSpline
        ( const MonotonicSpline& right     ,
          const char*             name = 0  ) ;
      /// destructor
      virtual ~MonotonicSpline() ;
      /// clone
      MonotonicSpline* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      MonotonicSpline  () {} ;
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
      const Ostap::Math::MonotonicSpline& function() const { return m_spline ; }
      const Ostap::Math::MonotonicSpline& spline  () const { return m_spline ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooArgList& phis  () const { return m_phis          ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::MonotonicSpline m_spline ;          // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class ConvexOnlySpline
     *  The special spline for non-negative
     *  convex or concave function
     *  @see Ostap::Math::ConvexSpline
     */
    class   ConvexOnlySpline : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::ConvexOnlySpline, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooArgList& phis  () const { return m_phis          ; }
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
     *  The special spline for non-negative monotonic
     *  convex or concave function
     *  @see Ostap::Math::ConvexSpline
     */
    class   ConvexSpline : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::ConvexSpline, 1) ;
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
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x      .arg() ; }
      const RooArgList& phis  () const { return m_phis          ; }
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
    /** @class CutOffGauss 
     *  Useful function for smooth Gaussian cut-off:
     *  \f[ f(x;x_0;\sigma) = \left\{ 
     *    \begin{array}{ll}
     *    1  & \mathrm{for~} x \le x_0  \                               \
     *    \mathrm{e}^{-\frac{1}{2}\left( \frac{ (x-x_0)^2}{\sigma^2} \right)}
     *       & \mathrm{for~} x >   x_0 
     *    \end{array}\right. \f] 
     *  @see Ostap::Math::CutOffGauss
     */
    class CutOffGauss : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::CutOffGauss, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      CutOffGauss 
      ( const char* name  , 
        const char* title ,
        RooAbsReal& x     , // observable 
        const bool  right , 
        RooAbsReal& x0    , 
        RooAbsReal& sigma ) ;
      // ======================================================================
      CutOffGauss
      ( const char* name  , 
        const char* title ,
        RooAbsReal& x     , // observable 
        RooAbsReal& x0    , 
        RooAbsReal& sigma , 
        const Ostap::Math::CutOffGauss& cutoff ) ;
      // Copy
      CutOffGauss ( const CutOffGauss& right    , 
                    const char*        name = 0 ) ;
      // destructor 
      virtual ~CutOffGauss() ;      
      /// clone
      CutOffGauss* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for the proper (de)serialization
      CutOffGauss () {} ;
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
      const Ostap::Math::CutOffGauss& function() const { return m_cutoff ; }
      const Ostap::Math::CutOffGauss& cutoff  () const { return m_cutoff ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& x0    () const { return m_x0    .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_x0    ;
      RooRealProxy m_sigma ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CutOffGauss m_cutoff ;          // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class CutOffStudent
     *  Useful function for smooth Student's t=-like (power-law) cut-off:
     *  \f[ f(x;x_0;\sigma) = \left\{ 
     *    \begin{array}{ll}
     *    1  & \mathrm{for~} x \le x_0  \                               \
     *    \left( \frac{1}{\nu} \left( \frac{(x-x_0)}{\sigma^2} \right)^{ - \frac{\nu+1}{2}} \right) 
     *       & \mathrm{for~} x >   x_0 
     *    \end{array}\right. \f] 
     *  @see Ostap::Math::CutOffStudent
     */
    class CutOffStudent : public RooAbsPdf 
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::CutOffStudent, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      CutOffStudent 
      ( const char* name  , 
        const char* title ,
        RooAbsReal& x     , // observable 
        const bool  right , 
        RooAbsReal& x0    , 
        RooAbsReal& nu    , 
        RooAbsReal& sigma ) ;
      // ======================================================================
      CutOffStudent
      ( const char* name  , 
        const char* title ,
        RooAbsReal& x     , // observable 
        RooAbsReal& x0    , 
        RooAbsReal& nu    , 
        RooAbsReal& sigma , 
        const Ostap::Math::CutOffStudent& cutoff ) ;
      // Copy
      CutOffStudent ( const CutOffStudent& right    , 
                      const char*        name = 0 ) ;
      // destructor 
      virtual ~CutOffStudent() ;      
      /// clone
      CutOffStudent* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for the proper (de)serialization
      CutOffStudent () {} ;
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
      const Ostap::Math::CutOffStudent& function() const { return m_cutoff ; }
      const Ostap::Math::CutOffStudent& cutoff  () const { return m_cutoff ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x     .arg() ; }
      const RooAbsReal& x0    () const { return m_x0    .arg() ; }
      const RooAbsReal& nu    () const { return m_nu    .arg() ; }
      const RooAbsReal& sigma () const { return m_sigma .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     ;
      RooRealProxy m_x0    ;
      RooRealProxy m_nu    ;
      RooRealProxy m_sigma ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::CutOffStudent m_cutoff ;          // the function
      // ======================================================================
    };
    // ========================================================================
    /** @class Uniform
     *  The trivial model: flat/uniform dsitribution in 1,2,3-dimensions 
     */
    class Uniform : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Uniform, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// flat in 1D 
      Uniform
      ( const char*       name  ,
        const char*       title ,
        RooAbsReal&       x     ) ;
      /// flat in 2D 
      Uniform
      ( const char*       name  ,
        const char*       title ,
        RooAbsReal&       x     , 
        RooAbsReal&       y     ) ;
      /// flat in 3D 
      Uniform
      ( const char*       name  ,
        const char*       title ,
        RooAbsReal&       x     , 
        RooAbsReal&       y     ,
        RooAbsReal&       z     ) ;
      /// copy
      Uniform
        ( const Uniform& right     ,
          const char*    name = 0  ) ;
      /// destructor
      virtual ~Uniform () ;
      /// clone
      Uniform* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Uniform () {} ;
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
      //  dimensionality of the PDF
      unsigned short    dim  () const { return m_dim      ; }      
      const RooAbsReal& x    () const { return m_x .arg() ; }
      const RooAbsReal& y    () const { return m_y .arg() ; }
      const RooAbsReal& z    () const { return m_z .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// variables/observables  
      unsigned short m_dim { 0 } ;
      RooRealProxy   m_x   {   } ;
      RooRealProxy   m_y   {   } ;
      RooRealProxy   m_z   {   } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Rice
     *  @see Ostap::Math::Rice
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
    class Rice : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Rice, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      Rice
      ( const char*  name      , 
        const char*  title     , 
        RooAbsReal&  x         ,
        RooAbsReal&  nu        ,
        RooAbsReal&  varsigma  ,
        RooAbsReal&  shift     ) ;
      Rice
      ( const char*  name      , 
        const char*  title     , 
        RooAbsReal&  x         ,
        RooAbsReal&  nu        ,
        RooAbsReal&  varsigma  ) ;
      /// copy constructor 
      Rice ( const Rice& right , const char* name = nullptr ) ;
      /// clone method
      Rice* clone ( const char* name ) const override ;
      // virtual destructor  
      virtual ~Rice() ;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Rice () {} ;
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
      const Ostap::Math::Rice& function() const { return m_rice ; }
      const Ostap::Math::Rice& rice    () const { return m_rice ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x        () const { return m_x        .arg() ; }
      const RooAbsReal& nu       () const { return m_nu       .arg() ; }
      const RooAbsReal& varsigma () const { return m_varsigma .arg() ; }
      const RooAbsReal& shift    () const { return m_shift    .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy   m_x         { } ;
      RooRealProxy   m_nu        { } ;
      RooRealProxy   m_varsigma  { } ;
      RooRealProxy   m_shift     { } ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Rice m_rice ;          // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenInvGauss
     *  @see Ostap::Math::GenInvGauss
     *  Generalised Inverse Gaussian distribution using
     *  \f$ (\theta,\eta) \f$ parameterisation  
     *  - |f$ \theta = \sqrt{ab}\$ 
     *  - |f$ \eta   = \sqrt{\frac{b}{a}}\$ 
     *  @see https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution
     */
    class GenInvGauss : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::GenInvGauss, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      GenInvGauss
      ( const char*  name      , 
        const char*  title     , 
        RooAbsReal&  x         ,
        RooAbsReal&  theta     ,
        RooAbsReal&  eta       ,
        RooAbsReal&  p         ,
        RooAbsReal&  shift     ) ;
      GenInvGauss 
      ( const char*  name      , 
        const char*  title     , 
        RooAbsReal&  x         ,
        RooAbsReal&  theta     ,
        RooAbsReal&  eta       ,
        RooAbsReal&  p         ) ;
      /// copy constructor 
      GenInvGauss ( const GenInvGauss& right , const char* name = nullptr ) ;
      /// clone method
      GenInvGauss* clone ( const char* name ) const override ;
      // virtual destructor  
      virtual ~GenInvGauss() ;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      GenInvGauss () {} ;
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
      const Ostap::Math::GenInvGauss& function() const { return m_gig ; }
      const Ostap::Math::GenInvGauss& gig     () const { return m_gig ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x     .arg() ; }
      const RooAbsReal& theta   () const { return m_theta .arg() ; }
      const RooAbsReal& eta     () const { return m_eta   .arg() ; }
      const RooAbsReal& p       () const { return m_p     .arg() ; }
      const RooAbsReal& shift   () const { return m_shift .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy   m_x        { } ;
      RooRealProxy   m_theta    { } ;
      RooRealProxy   m_eta      { } ;
      RooRealProxy   m_p        { } ;
      RooRealProxy   m_shift    { } ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenInvGauss m_gig ;          // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class SkewGenT 
     *  Skewwed Generalised t-distribution
     *  @see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution
     *  Original function is parameterised in terms of parameters 
     *  - \f$ \mu \$ related to locartion 
     *  - \f$ \sigma \$ related to width/scale 
     *  - \f$ -1 < \lambda < 1 \f$ related to asymmetry/skewness  
     *  - \f$ 0<p, 0<q \f$ related to kutsosis
     *
     *  Mean value is defined if \f$ 1 < pq \f$ 
     *  RMS si defined for \f$ 2 < pq \f$
     * 
     *  In this view here we adopt sligth reparameterisation in terms of 
     *  - \f$ 0 < r \f$, such as  \f$  r = \frac{1}{p} 
     *  - \f$ 0< \zeta \f$, such as \f$ pq = \zeta + 4 \f$
     *  - \f$ -\infty < \xi < +\infty \f$, such as \f$ \lambda  = \tanh \xi \f$   
     *
     *  Usage of \f$ \zeta\f$ ensures the existance of the  mean, RMS, sewness & kurtosis
     * 
     *  Special limitnig cases:
     *  - \f$ q\rigtharrow +\infty (\zeta \rightarrow +\infty) \f$ 
     *     Generalized Error Distribution 
     *  - \f$ \lambda=0 (\xi = 0)  \f$ Generalized t-distribution 
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
     * @see Ostap::Math::SkewGenT 
     * @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class SkewGenT : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::SkewGenT, 1) ;
      // ======================================================================
    public:
      // ======================================================================      
     SkewGenT 
     ( const char*  name  , 
       const char*  title , 
       RooAbsReal&  x     ,
       RooAbsReal&  mu    ,   // location/mean  
       RooAbsReal&  sigma ,   // scale/rms 
       RooAbsReal&  psi   ,   // related to asymmetry 
       RooAbsReal&  r     ,   // shape parameter 
       RooAbsReal&  zeta  ) ; // shape parameter 
      /// copy constructor 
      SkewGenT  ( const SkewGenT& right , const char* name = nullptr ) ;
      /// clone method
      SkewGenT* clone ( const char* name ) const override ;
      /// virtual destructor
      virtual ~SkewGenT() ;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      SkewGenT () {} ;
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
      const Ostap::Math::SkewGenT& function () const { return m_sgt ; }
      const Ostap::Math::SkewGenT& sgt      () const { return m_sgt ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x     .arg() ; }
      const RooAbsReal& mu      () const { return m_mu    .arg() ; }
      const RooAbsReal& sigma   () const { return m_sigma .arg() ; }
      const RooAbsReal& psi     () const { return m_psi   .arg() ; }
      const RooAbsReal& r       () const { return m_r     .arg() ; }
      const RooAbsReal& zeta    () const { return m_zeta  .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy   m_x      {} ;
      RooRealProxy   m_mu     {} ;
      RooRealProxy   m_sigma  {} ;
      RooRealProxy   m_psi    {} ;
      RooRealProxy   m_r      {} ;
      RooRealProxy   m_zeta   {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::SkewGenT m_sgt ;          // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class SkewGenError
     *  Skewed gheneralised error districbution 
     *  @see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution#Skewed_generalized_error_distribution
     *
     *  The Special  case of Skewwed Generaliaed T-distribution 
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
     *  - \f$ -\infty < \xi < +\infty \f$, such as \f$ \lambda  = \tanh \xi \f$   
     * 
     *  special cases: 
     *  - \f$ \xi=0 (\lambda=0), p=2\$ corresponds to Gaussian function 
     *  - \f$ \xi=0 (\lambda=0), p=1\$ corresponds to Laplace case 
     *
     *  @see Ostap::Math::SkewGenError
     *  @see Ostap::Math::SkewGenT 
     *  @see Ostap::Models::SkewGenT 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class SkewGenError : public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::SkewGenError, 1) ;
      // ======================================================================
    public:
      // ======================================================================      
     SkewGenError 
     ( const char*  name  , 
       const char*  title , 
       RooAbsReal&  x     ,
       RooAbsReal&  mu    ,   // location/mean  
       RooAbsReal&  sigma ,   // scale/rms 
       RooAbsReal&  xi    ,   // related to asymmetry 
       RooAbsReal&  p     ) ; // shape parameter 
      /// copy constructor 
      SkewGenError  ( const SkewGenError& right , const char* name = nullptr ) ;
      /// clone method
      SkewGenError* clone ( const char* name ) const override ;
      /// virtual destructor
      virtual ~SkewGenError() ;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      SkewGenError () {} ;
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
      const Ostap::Math::SkewGenError& function () const { return m_sge ; }
      const Ostap::Math::SkewGenError& sge      () const { return m_sge ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x     .arg() ; }
      const RooAbsReal& mu      () const { return m_mu    .arg() ; }
      const RooAbsReal& sigma   () const { return m_sigma .arg() ; }
      const RooAbsReal& xi      () const { return m_xi    .arg() ; }
      const RooAbsReal& p       () const { return m_p     .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy   m_x      {} ;
      RooRealProxy   m_mu     {} ;
      RooRealProxy   m_sigma  {} ;
      RooRealProxy   m_xi     {} ;
      RooRealProxy   m_p      {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::SkewGenError m_sge ;          // the function
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
     *  @see Ostap::Math::HORNSdini 
     */
    class HORNSdini : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::HORNSdini, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      HORNSdini 
      ( const char*  name  , 
        const char*  title , 
        RooAbsReal&  x     ,
        RooAbsReal&  a     , 
        RooAbsReal&  delta ,
        RooAbsReal&  phi   ) ;      
      /// copy constructor 
      HORNSdini ( const HORNSdini & right , const char* name = nullptr ) ;
      /// clone method 
      HORNSdini* clone ( const char* name ) const override ;
      // virtual destructor  
      virtual ~HORNSdini() ;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      HORNSdini  () {} ;
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
      const Ostap::Math::HORNSdini& function () const { return m_horns ; }
      const Ostap::Math::HORNSdini& horns    () const { return m_horns ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x     .arg() ; }
      const RooAbsReal& a       () const { return m_a     .arg() ; }
      const RooAbsReal& delta   () const { return m_delta .arg() ; }
      const RooAbsReal& phi     () const { return m_phi   .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy   m_x     {} ;
      RooRealProxy   m_a     {} ;
      RooRealProxy   m_delta {} ;
      RooRealProxy   m_phi   {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::HORNSdini m_horns ;          // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class HILLdini 
     *  \f[ f(x;a,\delta, \phi) = 
     *  \frac{3}{2\delta}\left( 1 - z^2 \right)
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
     *  @see Ostap::Math::HORNSdini 
     */
    class HILLdini : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::HILLdini, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      HILLdini 
      ( const char*  name  , 
        const char*  title , 
        RooAbsReal&  x     ,
        RooAbsReal&  a     , 
        RooAbsReal&  delta ,
        RooAbsReal&  phi   ) ;
      /// copy constructor 
      HILLdini ( const HILLdini & right , const char* name = nullptr ) ;
      /// clone method 
      HILLdini* clone ( const char* name ) const override ;
      // virtual destructor  
      virtual ~HILLdini() ;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      HILLdini  () {} ;
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
      const Ostap::Math::HILLdini& function () const { return m_hill ; }
      const Ostap::Math::HILLdini& hill     () const { return m_hill ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x     .arg() ; }
      const RooAbsReal& a       () const { return m_a     .arg() ; }
      const RooAbsReal& delta   () const { return m_delta .arg() ; }
      const RooAbsReal& phi     () const { return m_phi   .arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy   m_x     {} ;
      RooRealProxy   m_a     {} ;
      RooRealProxy   m_delta {} ;
      RooRealProxy   m_phi   {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::HILLdini m_hill ;          // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class KarlinShapley
     *  Positive polinomial on the interval 
     *  @see Ostap::Math::KarlinStudden
     *  @see Ostap::Math::KarlinShapley
     *  @see Ostap::Math::Posititive 
     *  @see Ostap::Models::PplyPosititive 
     *  @see Ostap::Models::KarlinShapley
     *  - note that Ostap::Models::PolyPositive shodu be better here 
     */
    class KarlinShapley : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::KarlinShapley, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      KarlinShapley
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        const RooArgList&    coeffs    ,
        const double         xmin      ,
        const double         xmax      ) ;
      /// general
      KarlinShapley
      ( const char*          name      ,
        const char*          title     ,
        RooRealVar&          x         ,
        const RooArgList&    coeffs    ) ;
      /// copy
      KarlinShapley
      ( const KarlinShapley& right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~KarlinShapley() ;
      /// clone
      KarlinShapley* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      KarlinShapley () {} ;
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
      const Ostap::Math::KarlinShapley& function       () const { return m_positive ; }
      const Ostap::Math::KarlinShapley& positive       () const { return m_positive ; }
      const Ostap::Math::KarlinShapley& karlin_shapley () const { return m_positive ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x    .arg()      ; }
      const RooArgList& phis () const { return m_phis             ; }
      double            xmin () const { return m_positive.xmin () ; }
      double            xmax () const { return m_positive.xmax () ; }    
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    {} ;
      RooListProxy m_phis {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::KarlinShapley  m_positive {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class KarlinStudden
     *  Positive polinomial on the interval
     *  @see Ostap::Math::KarlinStudden
     *  @see Ostap::Math::KarlinShapley
     *  @see Ostap::Math::Posititive 
     *  @see Ostap::Models::PplyPosititive 
     *  @see Ostap::Models::KarlinShapley
     *  - note that Ostap::Models::PolyPositive shodul be better choice  
     */
    class KarlinStudden : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::KarlinStudden, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      KarlinStudden
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        const RooArgList&    coeffs    ,
        const double         xmin      ,
        const double         scale = 1 ) ;
      /// general
      KarlinStudden 
      ( const char*          name      ,
        const char*          title     ,
        RooRealVar&          x         ,
        const RooArgList&    coeffs    ,
        const double         scale = 1 ) ;
      /// copy
      KarlinStudden
      ( const KarlinStudden& right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~KarlinStudden() ;
      /// clone
      KarlinStudden* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      KarlinStudden () {} ;
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
      const Ostap::Math::KarlinStudden& function       () const { return m_positive ; }
      const Ostap::Math::KarlinStudden& positive       () const { return m_positive ; }
      const Ostap::Math::KarlinStudden& karlin_studden () const { return m_positive ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x     () const { return m_x    .arg()       ; }
      const RooArgList& phis  () const { return m_phis              ; }
      double            xmin  () const { return m_positive.xmin  () ; }
      double            scale () const { return m_positive.scale () ; }    
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    {} ;
      RooListProxy m_phis {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::KarlinStudden  m_positive {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenPareto
     *  Generalized Pareto Distribution
     *  @see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
     *  @see Ostap::Math::GenPareto
     */
    class GenPareto : public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::GenPareto, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      GenPareto
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          scale     ,
        RooAbsReal&          shape     ) ;
      /// copy
      GenPareto
      ( const GenPareto&     right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~GenPareto () ;
      /// clone
      GenPareto* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      GenPareto  () {} ;
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
      const Ostap::Math::GenPareto& function () const { return m_gpd ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg () ; }
      const RooAbsReal& mu     () const { return m_mu    .arg () ; }
      const RooAbsReal& scale  () const { return m_scale .arg () ; }
      const RooAbsReal& shape  () const { return m_shape .arg () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_mu    {} ;
      RooRealProxy m_scale {} ;
      RooRealProxy m_shape {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GenPareto  m_gpd {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ExGenPareto
     *  Reparameterised Exponentiated Generalized Pareto Distribution
     *  @see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
     *  @see Ostap::Math::ExGenPareto
     */
    class ExGenPareto : public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::ExGenPareto, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      ExGenPareto
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          scale     ,
        RooAbsReal&          shape     ) ;
      /// copy
      ExGenPareto
      ( const ExGenPareto&   right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~ExGenPareto () ;
      /// clone
      ExGenPareto* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      ExGenPareto  () {} ;
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
      const Ostap::Math::ExGenPareto& function () const { return m_egpd ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg () ; }
      const RooAbsReal& mu     () const { return m_mu    .arg () ; }
      const RooAbsReal& scale  () const { return m_scale .arg () ; }
      const RooAbsReal& shape  () const { return m_shape .arg () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_mu    {} ;
      RooRealProxy m_scale {} ;
      RooRealProxy m_shape {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::ExGenPareto  m_egpd {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Benini
     *  Benini Distribution
     *  @see https://en.wikipedia.org/wiki/Benini_distribution
     *  @see Ostap::Math::Benini
     */
    class Benini : public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Benini, 2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// Modified Benini
      Benini
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooArgList&          shape     , 
        RooAbsReal&          scale     , 
        RooAbsReal&          shift     ) ;
      /// Modified Benini
      Benini
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          beta      ,
        RooAbsReal&          gamma     ,
        RooAbsReal&          delta     ,
        RooAbsReal&          scale     , 
        RooAbsReal&          shift     ) ;
      /// standard Benini
      Benini
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          alpha     ,
        RooAbsReal&          beta      ,
        RooAbsReal&          scale     ) ;
      /// copy
      Benini
      ( const Benini&        right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Benini () ;
      /// clone
      Benini* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      Benini  () {} ;
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
      const Ostap::Math::Benini& function () const { return m_benini ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg () ; }
      const RooArgList& shape  () const { return m_shape         ; }
      const RooAbsReal& scale  () const { return m_scale .arg () ; }
      const RooAbsReal& shift  () const { return m_shift .arg () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooListProxy m_shape {} ;
      RooRealProxy m_scale {} ;
      RooRealProxy m_shift {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Benini  m_benini {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GEV
     *  Generalized extreme value distribution
     *  @see https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
     *  @see Ostap::Math::GEV
     */
    class GEV : public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::GEV, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      GEV
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          scale     ,
        RooAbsReal&          shape     ) ;
      /// copy
      GEV
      ( const GEV&           right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~GEV () ;
      /// clone
      GEV* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      GEV  () {} ;
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
      const Ostap::Math::GEV& function () const { return m_gev ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg () ; }
      const RooAbsReal& mu     () const { return m_mu    .arg () ; }
      const RooAbsReal& scale  () const { return m_scale .arg () ; }
      const RooAbsReal& shape  () const { return m_shape .arg () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_mu    {} ;
      RooRealProxy m_scale {} ;
      RooRealProxy m_shape {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::GEV  m_gev {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FisherZ
     *  Fisher's Z-distirbution with additional location-scale parameters 
     *  @see https://en.wikipedia.org/wiki/Fisher%27s_z-distribution
     *  @see Ostap::Math::FisherZ
     */
    class FisherZ: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::FisherZ, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      FisherZ
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,    // mode 
        RooAbsReal&          d1        ,    // d1-shape 
        RooAbsReal&          d2        ,    // d2-shape 
        RooAbsReal&          scale     ) ;  // scale 
      /// general
      FisherZ
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,    // mode 
        RooAbsReal&          d1        ,    // d1-shape 
        RooAbsReal&          d2        ,    // d2-shape 
        const double         scale = 1 ) ;  // scale 
      /// copy
      FisherZ
      ( const FisherZ&       right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~FisherZ() ;
      /// clone
      FisherZ* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      FisherZ () {} ;
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
      const Ostap::Math::FisherZ& function () const { return m_fz ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg  () ; }
      const RooAbsReal& mu     () const { return m_mu    .arg  () ; }
      const RooAbsReal& scale  () const { return m_scale .arg  () ; }
      const RooAbsReal& d1     () const { return m_d1    .arg  () ; }
      const RooAbsReal& d2     () const { return m_d2    .arg  () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_mu    {} ;
      RooRealProxy m_scale {} ;
      RooRealProxy m_d1    {} ;
      RooRealProxy m_d2    {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::FisherZ  m_fz {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BirnbaumSaunders
     *  Birnbaum-Saunders distribution 
     *  @see https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution
     *  \f[ f(x;\mu, \beta,\gamma) = 
     *   \frac{ z + z^{-1}}{2\gamma(x-\mu)}\phi( \frac{1}{\gamma}(z-z^{-1}) \f]
     *  where
     *   - \f$ z=\frac{x-\mu}{\beta}\f$
     *   - \f$ \phi\f$ is Gaussian PDF 
     *  @see Ostap::Math::BirnbaumSaunders
     */
    class BirnbaumSaunders: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::BirnbaumSaunders, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      BirnbaumSaunders
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,    // location 
        RooAbsReal&          beta      ,    // scale 
        RooAbsReal&          gamma     ) ;  // shape 
      /// copy
      BirnbaumSaunders
      ( const BirnbaumSaunders& right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~BirnbaumSaunders() ;
      /// clone
      BirnbaumSaunders* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      BirnbaumSaunders() {} ;
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
      const Ostap::Math::BirnbaumSaunders& function () const { return m_bs ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg  () ; }
      const RooAbsReal& mu     () const { return m_mu    .arg  () ; }
      const RooAbsReal& beta   () const { return m_beta  .arg  () ; }
      const RooAbsReal& gamma  () const { return m_gamma .arg  () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_mu    {} ;
      RooRealProxy m_beta  {} ;
      RooRealProxy m_gamma {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::BirnbaumSaunders m_bs {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class MPERT
     *  Modified PERT distribution 
     *  @see https://en.wikipedia.org/wiki/PERT_distribution
     *  @see https://www.vosesoftware.com/riskwiki/ModifiedPERTdistribution.php
     *  @see Ostap::Math::MPERT
     */
    class MPERT: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::MPERT, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general
      MPERT
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          xi        ,
        RooAbsReal&          gamma     , 
        const double         xmin      , 
        const double         xmax      ) ;
      /// get xmin/xmax from variable limits 
      MPERT
      ( const char*          name      ,
        const char*          title     ,
        RooAbsRealLValue&    x         ,
        RooAbsReal&          xi        ,
        RooAbsReal&          gamma     ) ;
      /// copy
      MPERT
      ( const MPERT&         right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~MPERT () ;
      /// clone
      MPERT* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      MPERT  () {} ;
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
      const Ostap::Math::MPERT& function () const { return m_mpert ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg  () ; }
      const RooAbsReal& xi     () const { return m_xi    .arg  () ; }
      const RooAbsReal& gamma  () const { return m_gamma .arg  () ; }
      double            xmin   () const { return m_mpert .xmin () ; }
      double            xmax   () const { return m_mpert .xmax () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_xi    {} ;
      RooRealProxy m_gamma {} ;
      RooRealProxy m_shape {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::MPERT  m_mpert {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Rational
     *  Ratio of two positive Bernstein polynomials 
     *  @see Ostap::Math::RationalPositive 
     */
    class Rational : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Rational, 2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      Rational 
      (  const char*          name  , 
         const char*          title ,
         RooAbsReal&          x     , // observable 
         const RooArgList&    p     , // parameters for numerator 
         const RooArgList&    q     , // parameters for denumerator 
         const double         xmin  , 
         const double         xmax  ) ;
      /// constructor
      Rational 
      (  const char*          name  , 
         const char*          title ,
         RooAbsReal&          x     , // observable 
         const unsigned short p     , // degree of numerator 
         const RooArgList&    a     , // all parameters
         const double         xmin  , 
         const double         xmax  ) ;
      /// copy 
      Rational
      ( const Rational& right          , 
	const char*     name = nullptr ) ;
      /// destructor 
      virtual ~Rational() ;
      /// clone method
      Rational* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// fake default constructor 
      Rational () {} ;
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
      /// access to underlying rational function
      inline const Ostap::Math::RationalPositive& rational    () const
      { setPars () ; return m_rational ; } 
      /// access to numerator of underlying   rational function
      inline const Ostap::Math::Positive&         numerator   () const
      { setPars () ; return m_rational.numerator() ; } 	
      /// access to denomerator of underlying rational function
      inline const Ostap::Math::Positive&         denominator () const 
      { setPars () ; return m_rational.denominator() ; } 	
      /// access to underlying rational function
      inline const Ostap::Math::RationalPositive& function    () const
      { setPars () ; return m_rational ; } 	
      // ======================================================================      
    public:
      // ======================================================================
      /// parameters 
      const RooArgList& pars () const { return m_pars               ; }
      /// observable 
      const RooAbsReal& x    () const { return m_x.arg()            ; }
      /// degree of nuemrator 
      unsigned short    p    () const { return m_rational.pnpars () ; }
      /// degree of denominator  
      unsigned short    q    () const { return m_rational.qnpars () ; }
      // ======================================================================      
    private:
      // ======================================================================
      /// observable 
      RooRealProxy                          m_x        {} ; // observable 
      /// parameters  
      RooListProxy                          m_pars     {} ; // parameters
      // the function 
      mutable Ostap::Math::RationalPositive m_rational {} ; // fnuction 
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class Shape1D
     *  simple generic PDF
     */
    class Shape1D final : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Shape1D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
      Shape1D
      ( const std::string& name   , 
        const std::string& title  , 
        RooAbsReal&        x      ,
        FUNCTION           f      , 
        const std::size_t tag = 0 )
        : RooAbsPdf   (  name.c_str() ,  title.c_str() ) 
        , m_x         ( "!x"   , "Variable" , this , x ) 
        , m_function  ( f   ) 
        , m_tag       ( tag ) 
      {}
      // ======================================================================
      Shape1D
      ( const std::string&            name    , 
        const std::string&            title   , 
        RooAbsReal&                   x       ,
        std::function<double(double)> f       , 
        const std::size_t             tag = 0 ) ;
      /// copy constructor 
      Shape1D ( const Shape1D& right , const char* name = nullptr ) ;
      /// clone method
      Shape1D* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
      static inline Shape1D 
      create
      ( const std::string& name     , 
        const std::string& title    , 
        RooAbsReal&        x        ,
        FUNCTION           f        , 
        const std::size_t  tag = 0  ) 
      { return Shape1D  ( name , title , x , f , tag ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the PDF 
      Double_t evaluate () const override
      { const double x = m_x ; return std::max ( m_function ( x ) , 0.0 ) ; ; }
      // ======================================================================        
    public:
      // ======================================================================
      /// evaluate the function
      double func  ( const double x ) const 
      { return std::max ( m_function ( x ) , 0.0 ) ; }
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
    private :
      // ======================================================================
      /// variable 
      RooRealProxy                   m_x            ; // variable 
      /// the function itself 
      std::function<double(double)>  m_function     ; // function 
      /// helper (hopefully unique) tag 
      std::size_t                    m_tag          ; // tag 
      // ======================================================================
    } ;    
    // ========================================================================
    /** @class Shape2D
     *  simple generic PDF
     */
    class Shape2D final : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Shape2D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
      Shape2D
      ( const std::string& name    ,  
        const std::string& title   ,  
        RooAbsReal&       x       ,
        RooAbsReal&       y       ,
        FUNCTION          f       , 
        const std::size_t tag = 0 ) 
        : RooAbsPdf  (  name.c_str () ,  title.c_str()  ) 
        , m_x        ( "!x"   , "x-variable" , this , x ) 
        , m_y        ( "!y"   , "y-variable" , this , y ) 
        , m_function ( f   ) 
        , m_tag      ( tag ) 
      {}
      // ======================================================================
      Shape2D
      ( const std::string&                   name    ,  
        const std::string&                   title   ,  
        RooAbsReal&                          x       ,
        RooAbsReal&                          y       ,
        std::function<double(double,double)> f       , 
        const std::size_t                    tag = 0 ) ;
      /// copy constructor 
      Shape2D ( const Shape2D& right , const char* name = nullptr ) ;
      /// clone method
      Shape2D* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
        static inline Shape2D
        create
        ( const std::string& name    ,  
          const std::string& title   ,  
          RooAbsReal&        x       ,
          RooAbsReal&        y       ,
          FUNCTION           f       ,
          const std::size_t  tag = 0 ) 
      { return Shape2D ( name , title , x , y , f , tag ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the PDF 
      Double_t evaluate () const override
      { 
        const double x = m_x ; 
        const double y = m_y ; 
        return std::max ( m_function ( x , y ) , 0.0 ) ; 
      }
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
      /// evaluate the function
      double func  ( const double x , 
                     const double y ) const 
      { return std::max ( m_function ( x , y ) , 0.0 ) ; }
      // ======================================================================        
    private :
      // ======================================================================
      /// x-variable 
      RooRealProxy                          m_x        ; // x-variable 
      /// y-variable 
      RooRealProxy                          m_y        ; // y-variable 
      /// the function itself 
      std::function<double(double,double)>  m_function ; // function 
      /// unique tag 
      std::size_t                           m_tag      ; // the tag      
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Shape3D
     *  simple generic PDF
     */
    class Shape3D final : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Shape3D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
      Shape3D
      ( const std::string& name    , 
        const std::string& title   , 
        RooAbsReal&        x       ,
        RooAbsReal&        y       ,
        RooAbsReal&        z       ,
        FUNCTION           f       , 
        const std::size_t  tag = 0 )
        : RooAbsPdf  (  name.c_str()  ,  title.c_str() ) 
        , m_x        ( "!x"   , "x-variable" , this , x ) 
        , m_y        ( "!y"   , "y-variable" , this , y ) 
        , m_z        ( "!z"   , "z-variable" , this , z ) 
        , m_function ( f ) 
        , m_tag      ( tag ) 
      {}
      // ======================================================================
      Shape3D
      ( const std::string&                          name    , 
        const std::string&                          title   , 
        RooAbsReal&                                 x       ,
        RooAbsReal&                                 y       ,
        RooAbsReal&                                 z       ,
        std::function<double(double,double,double)> f       , 
        const std::size_t                           tag = 0 ) ;
      /// copy constructor 
      Shape3D ( const Shape3D& right , const char* name = nullptr ) ;
      /// clone method
      Shape3D* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
      static inline Shape3D 
      create
      ( const std::string& name    , 
        const std::string& title   , 
        RooAbsReal&        x       ,
        RooAbsReal&        y       ,
        RooAbsReal&        z       ,
        FUNCTION           f       , 
        const std::size_t  tag = 0 )
      { return Shape3D ( name , title , x , y , z , f , tag ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the PDF 
      Double_t evaluate () const override
      { 
        const double x = m_x ; 
        const double y = m_y ; 
        const double z = m_z ; 
        return std::max ( m_function ( x , y , z ) , 0.0 ) ; 
      }
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
      /// evaluate the function
      double func  ( const double x , 
                     const double y , 
                     const double z ) const 
      { return std::max ( m_function ( x , y , z ) , 0.0 ) ; }
      // ======================================================================        
    private :
      // ======================================================================
      /// x-variable 
      RooRealProxy                          m_x        ; // x-variable 
      /// y-variable 
      RooRealProxy                          m_y        ; // y-variable 
      /// z-variable 
      RooRealProxy                          m_z        ; // z-variable 
      /// the function itself 
      std::function<double(double,double,double)>  m_function ; // function 
      /// unuqie  tag 
      std::size_t                                  m_tag      ; // unuqire  tag 
      // ======================================================================      
    } ;        
    // ========================================================================
    /** @class Histo1D
     *  simple generic PDF from the histogram 
     */
    class Histo1D final : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Histo1D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      Histo1D 
        ( const char*                 name  , 
          const char*                 title , 
          RooAbsReal&                 x     ,
          const Ostap::Math::Histo1D& histo ) ;
      /// copy constructor 
      Histo1D ( const Histo1D& right , const char* name = nullptr ) ;
      /// clone method
      Histo1D* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // fake default contructor, needed just for the proper (de)serialization
      Histo1D () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the PDF 
      Double_t evaluate () const override { return func ( m_x ) ; }
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
      /// the function itself 
      const Ostap::Math::Histo1D& histo () const { return   m_histo ; }
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function
      double func  
        ( const double x ) const 
      { return std::max ( m_histo ( x ) , 0.0 ) ; }
      // ======================================================================        
    public:
      // ======================================================================
      const RooAbsReal& x () const { return m_x     .arg() ; }
      // ======================================================================
    private :
      // ======================================================================
      /// variable 
      RooRealProxy                   m_x            ; // variable 
      /// the function itself 
      Ostap::Math::Histo1D           m_histo        ; // function 
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Histo2D
     *  simple generic PDF from the histogram 
     */
    class Histo2D final : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Histo2D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      Histo2D 
        ( const char*                 name  , 
          const char*                 title , 
          RooAbsReal&                 x     ,
          RooAbsReal&                 y     ,
          const Ostap::Math::Histo2D& histo ) ;
      /// copy constructor 
      Histo2D ( const Histo2D& right , const char* name = nullptr ) ;
      /// clone method
      Histo2D* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // fake default contructor, needed just for the proper (de)serialization
      Histo2D () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the PDF 
      Double_t evaluate () const override { return func ( m_x , m_y ) ; }
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
      /// the function itself 
      const Ostap::Math::Histo2D& histo () const { return  m_histo ; }
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function
      double func
        ( const double x , 
          const double y ) const 
      { return std::max ( m_histo ( x , y ) , 0.0 ) ; }
      // ======================================================================        
    public:
      // ======================================================================
      const RooAbsReal& x () const { return m_x     .arg() ; }
      const RooAbsReal& y () const { return m_y     .arg() ; }
      // ======================================================================
    private :
      // ======================================================================
      /// x-variable 
      RooRealProxy                   m_x     ; // x-variable 
      /// y-variable 
      RooRealProxy                   m_y     ; // y-variable 
      /// the function itself 
      Ostap::Math::Histo2D           m_histo ; // function 
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Histo3D
     *  simple generic PDF from the histogram 
     */
    class Histo3D final : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::Histo3D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      Histo3D
        ( const char*                 name  , 
          const char*                 title , 
          RooAbsReal&                 x     ,
          RooAbsReal&                 y     ,
          RooAbsReal&                 z     ,
          const Ostap::Math::Histo3D& histo ) ;
      /// copy constructor 
      Histo3D ( const Histo3D& right , const char* name = nullptr ) ;
      /// clone method
      Histo3D* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // fake default contructor, needed just for the proper (de)serialization
      Histo3D () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the PDF 
      Double_t evaluate () const override { return func ( m_x , m_y , m_z ) ; }
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
      /// evaluate the function
      double func 
        ( const double x , 
          const double y , 
          const double z ) const 
      { return std::max ( m_histo ( x , y , z ) , 0.0 ) ; }
      // ======================================================================        
    public:
      // ======================================================================        
      /// the function itself 
      const Ostap::Math::Histo3D& histo () const { return  m_histo ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x () const { return m_x     .arg() ; }
      const RooAbsReal& y () const { return m_y     .arg() ; }
      const RooAbsReal& z () const { return m_z     .arg() ; }
      // ======================================================================
    private :
      // ======================================================================
      /// x-variable 
      RooRealProxy                   m_x     ; // x-variable 
      /// y-variable 
      RooRealProxy                   m_y     ; // y-variable 
      /// z-variable 
      RooRealProxy                   m_z     ; // z-variable 
      /// the function itself 
      Ostap::Math::Histo3D           m_histo ; // function 
      // ======================================================================      
    } ;
    // ========================================================================


   // ========================================================================
    /** @class Meixner
     *  Modified PERT distribution 
     *  @see https://en.wikipedia.org/wiki/PERT_distribution
     *  @see https://www.vosesoftware.com/riskwiki/ModifiedPERTdistribution.php
     *  @see Ostap::Math::MExner
     * 
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
    class Meixner: public RooAbsPdf 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Meixner, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// general case
      Meixner
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          sigma     , 
        RooAbsReal&          psi       , 
        RooAbsReal&          shape     ) ;
      /// symmetric case: psi = 0 
      Meixner
      ( const char*          name      ,
        const char*          title     ,
        RooAbsReal&          x         ,
        RooAbsReal&          mu        ,
        RooAbsReal&          sigma     ,  
        RooAbsReal&          shape     ) ;
      /// copy
      Meixner
      ( const Meixner&       right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Meixner () ;
      /// clone
      Meixner* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default constructor, needed just for proper (de)serialization
      Meixner   () {} ;
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
      const Ostap::Math::Meixner& function () const { return m_meixner; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x      () const { return m_x     .arg () ; }
      const RooAbsReal& mu     () const { return m_mu    .arg () ; }
      const RooAbsReal& sigma  () const { return m_sigma .arg () ; }
      const RooAbsReal& psi    () const { return m_psi   .arg () ; }
      const RooAbsReal& shape  () const { return m_shape .arg () ; }
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x     {} ;
      RooRealProxy m_mu    {} ;
      RooRealProxy m_sigma {} ;
      RooRealProxy m_psi   {} ;
      RooRealProxy m_shape {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Meixner m_meixner {} ;  // the function
      // ======================================================================
    } ;
    // ========================================================================
  } //                                           end of namespace Ostap::Models
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // ANALYSIS_MODELS_H
// ============================================================================
