// ============================================================================
#ifndef OSTAP_MOREVARS_H 
#define OSTAP_NOREVARS_H 1
// ============================================================================
// Include  files 
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
// ============================================================================
// Ostap
// ============================================================================
#include  "Ostap/Bernstein.h"
#include  "Ostap/Bernstein1D.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace MoreRooFit 
  {
    // ========================================================================
    /** @class Bernstein
     *  Simple Bernstein polynomial 
     *  \f[ p(x) = \sum_{k=0}^{n} a_i B_n^k(x)  \f],
     *  where \f$ B_n^k(x)\f$ is basic Bernstein polynomial   
     *  @see Ostap::Math::Bernstein
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-06     
     */
    class Bernstein : public RooAbsReal 
    {
      // ======================================================================
      ClassDef ( Ostap::MoreRooFit::Bernstein , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and the list of coefficients
      Bernstein ( const std::string& name  ,
                  const std::string& title ,
                  RooAbsReal&        xvar  ,
                  const double       xmin  , 
                  const double       xmax  ,
                  const RooArgList&  pars  ) ;
      // ======================================================================
      /// copy constructor 
      Bernstein ( const Bernstein& right , 
                  const char*      name = nullptr ) ;
      // ======================================================================
      Bernstein () ;
      virtual ~Bernstein() ;
      // ======================================================================
      Bernstein* clone ( const char* newname ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&  allVars             , 
        RooArgSet&  analVars            , 
        const char* rangeName = nullptr ) const ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Bernstein& bernstein () const { return m_bernstein ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars       () const ;
      // ======================================================================
    public: 
      // ======================================================================
      Double_t evaluate  () const override ;
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy                   m_xvar      ;
      RooListProxy                   m_pars      ;
      mutable Ostap::Math::Bernstein m_bernstein ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Monotonic 
     *  Simple monotonical polynomial 
     *  \f[ p(x) = a + b P(x) \f], 
     *  where \f$ P(x)\f$ is normalized positive monotonic polynomial:
     *  - \f$ \int_{x_{\mathrm{min}}}^{x_{\mathrm{max}}} P(x) dx = 1 \f$ 
     *  - \f$ P (x) \ge              0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
     *  - \f$ P^{\prime}(x) \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
     *  @see Ostap::Math::Monotonic 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-06     
     */
    class Monotonic : public RooAbsReal 
    {
      // ======================================================================
      ClassDef ( Ostap::MoreRooFit::Monotonic , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      Monotonic ( const std::string& name       ,
                  const std::string& title      ,
                  RooAbsReal&        xvar       ,
                  const bool         increasing , 
                  const double       xmin       , 
                  const double       xmax       ,
                  RooAbsReal&        a          ,
                  RooAbsReal&        b          ,
                  const RooArgList&  pars       ) ;
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      Monotonic ( const std::string& name       ,
                  const std::string& title      ,
                  RooAbsReal&        xvar       ,
                  const bool         increasing , 
                  const double       xmin       , 
                  const double       xmax       ,
                  const RooArgList&  pars       ) ;
      // ======================================================================
      /// copy constructor 
      Monotonic ( const Monotonic& right , 
                  const char*      name = nullptr ) ;
      // ======================================================================
      Monotonic  () ;
      virtual ~Monotonic () ;
      // ======================================================================
      Monotonic* clone ( const char* newname ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&  allVars             , 
        RooArgSet&  analVars            , 
        const char* rangeName = nullptr ) const ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Monotonic& monotonic () const { return m_monotonic ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars       () const ;
      // ======================================================================
    public: 
      // ======================================================================
      Double_t evaluate  () const override ;
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy                   m_xvar      ;
      RooRealProxy                   m_a         ;
      RooRealProxy                   m_b         ;
      RooListProxy                   m_pars      ;
      mutable Ostap::Math::Monotonic m_monotonic ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Convex
     *  Simple convex/concave polynomial 
     *  \f[ p(x) = a + b P(x) \f], 
     *  where \f$ P(x)\f$ is normalized positive monotonic polynomial:
     *  - \f$ \int_{x_{\mathrm{min}}}^{x_{\mathrm{max}}} P(x) dx = 1 \f$ 
     *  - \f$ P (x) \ge                    0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
     *  - \f$ P^{\prime}(x)       \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
     *  - \f$ P^{\prime\prime}(x) \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
     *  @see Ostap::Math::Convex
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-06     
     */
    class Convex : public RooAbsReal 
    {
      // ======================================================================
      ClassDef ( Ostap::MoreRooFit::Convex  , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      Convex    ( const std::string& name       ,
                  const std::string& title      ,
                  RooAbsReal&        xvar       ,
                  const bool         increasing , 
                  const bool         convex     , 
                  const double       xmin       , 
                  const double       xmax       ,
                  RooAbsReal&        a          ,
                  RooAbsReal&        b          ,
                  const RooArgList&  pars       ) ;
      // ======================================================================
      /// copy constructor 
      Convex    ( const Convex&    right , 
                  const char*      name = nullptr ) ;
      // ======================================================================
      Convex    () ;
      virtual ~Convex () ;
      // ======================================================================
      Convex* clone ( const char* newname ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&  allVars             , 
        RooArgSet&  analVars            , 
        const char* rangeName = nullptr ) const ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Convex& convex () const { return m_convex ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars       () const ;
      // ======================================================================
    public: 
      // ======================================================================
      Double_t evaluate  () const override ;
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy                   m_xvar      ;
      RooRealProxy                   m_a         ;
      RooRealProxy                   m_b         ;
      RooListProxy                   m_pars      ;
      mutable Ostap::Math::Convex    m_convex    ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ConvexOnly
     *  Simple convex/concave polynomial 
     *  \f[ p(x) = a + b P(x) \f], 
     *  where \f$ P(x)\f$ is normalized positive monotonic polynomial:
     *  - \f$ \int_{x_{\mathrm{min}}}^{x_{\mathrm{max}}} P(x) dx = 1 \f$ 
     *  - \f$ P (x) \ge                    0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
     *  - \f$ P^{\prime\prime}(x) \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
     *  @see Ostap::Math::ConvexOnly
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-06     
     */
    class ConvexOnly : public RooAbsReal 
    {
      // ======================================================================
      ClassDef ( Ostap::MoreRooFit::ConvexOnly  , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      ConvexOnly ( const std::string& name       ,
                   const std::string& title      ,
                   RooAbsReal&        xvar       ,
                   const bool         convex     , 
                   const double       xmin       , 
                   const double       xmax       ,
                   RooAbsReal&        a          ,
                   RooAbsReal&        b          ,
                   const RooArgList&  pars       ) ;
      // ======================================================================
      /// copy constructor 
      ConvexOnly ( const ConvexOnly&    right , 
                   const char*      name = nullptr ) ;
      // ======================================================================
      ConvexOnly () ;
      virtual ~ConvexOnly () ;
      // ======================================================================
      ConvexOnly* clone ( const char* newname ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&  allVars             , 
        RooArgSet&  analVars            , 
        const char* rangeName = nullptr ) const ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::ConvexOnly& convex () const { return m_convex ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars       () const ;
      // ======================================================================
    public: 
      // ======================================================================
      Double_t evaluate  () const override ;
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy                    m_xvar      ;
      RooRealProxy                    m_a         ;
      RooRealProxy                    m_b         ;
      RooListProxy                    m_pars      ;
      mutable Ostap::Math::ConvexOnly m_convex    ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                   The end of namespace Ostap::MoreRooFit
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MOREVARS_H
// ============================================================================
