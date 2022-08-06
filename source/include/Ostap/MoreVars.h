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
#include  "Ostap/BSpline.h"
// ============================================================================
/// froward declarations 
// ============================================================================
class RooAddPdf   ; // ROOT,RooFit 
class RooProdPdf  ; // ROOT,RooFit 
class RooGaussian ; // ROOT,RooFit 
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
      ClassDefOverride ( Ostap::MoreRooFit::Bernstein , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and the list of coefficients
      Bernstein
      ( const std::string& name  ,
        const std::string& title ,
        RooAbsReal&        xvar  ,
        const RooArgList&  pars  ,
        const double       xmin  , 
        const double       xmax  ) ;
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
        const char* rangeName = nullptr ) const override ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the varibale 
      const RooAbsReal& xvar   () const { return m_xvar.arg() ; }
      /// get parameters 
      const RooArgList& pars   () const { return m_pars       ; }
      /// xmin for bernstein 
      double            xmin   () const { return m_bernstein.xmin() ; }
      /// xmax for bernstein 
      double            xmax   () const { return m_bernstein.xmax() ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars       () const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Bernstein& function  () const { return m_bernstein ; }
      const Ostap::Math::Bernstein& bernstein () const { return m_bernstein ; }
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
      ClassDefOverride ( Ostap::MoreRooFit::Monotonic , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      Monotonic 
      ( const std::string& name       ,
        const std::string& title      ,
        RooAbsReal&        xvar       ,
        const RooArgList&  pars       ,
        const bool         increasing , 
        const double       xmin       , 
        const double       xmax       ,
        RooAbsReal&        a          ,
        RooAbsReal&        b          ) ;
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      Monotonic 
      ( const std::string& name       ,
        const std::string& title      ,
        RooAbsReal&        xvar       ,
        const RooArgList&  pars       ,
        const bool         increasing , 
        const double       xmin       , 
        const double       xmax       ,
        const double       a  = 0     , 
        const double       b  = 1     ) ;
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
        const char* rangeName = nullptr ) const override ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Monotonic& monotonic () const { return m_monotonic ; }
      const Ostap::Math::Monotonic& function  () const { return m_monotonic ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the varibale 
      const RooAbsReal& xvar   () const { return m_xvar.arg() ; }
      /// get parameters 
      const RooArgList& pars   () const { return m_pars       ; }
      /// get the shift  
      const RooAbsReal& a      () const { return m_a.arg()    ; }
      /// get the scale   
      const RooAbsReal& b      () const { return m_b.arg()    ; }
      /// xmin for bernstein 
      double xmin () const { return m_monotonic.xmin() ; }
      /// xmax for bernstein 
      double xmax () const { return m_monotonic.xmax() ; }
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
      ClassDefOverride ( Ostap::MoreRooFit::Convex  , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      Convex
      ( const std::string& name       ,
        const std::string& title      ,
        RooAbsReal&        xvar       ,
        const RooArgList&  pars       ,
        const bool         increasing , 
        const bool         convex     , 
        const double       xmin       , 
        const double       xmax       ,
        RooAbsReal&        a          ,
        RooAbsReal&        b          ) ;
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      Convex
      ( const std::string& name       ,
        const std::string& title      ,
        RooAbsReal&        xvar       ,
        const RooArgList&  pars       ,
        const bool         increasing , 
        const bool         convex     , 
        const double       xmin       , 
        const double       xmax       ,
        const double       a = 0      , 
        const double       b = 1      ) ;
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
        const char* rangeName = nullptr ) const override ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the varibale 
      const RooAbsReal& xvar   () const { return m_xvar.arg() ; }
      /// get parameters 
      const RooArgList& pars   () const { return m_pars       ; }
      /// get the shift  
      const RooAbsReal& a      () const { return m_a.arg()    ; }
      /// get the scale   
      const RooAbsReal& b      () const { return m_b.arg()    ; }
      /// xmin for bernstein 
      double xmin () const { return m_convex.xmin() ; }
      /// xmax for bernstein 
      double xmax () const { return m_convex.xmax() ; }
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::Convex& function () const { return m_convex ; }
      const Ostap::Math::Convex& convex   () const { return m_convex ; }
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
      ClassDefOverride ( Ostap::MoreRooFit::ConvexOnly  , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and list of coefficients
      ConvexOnly
      ( const std::string& name       ,
        const std::string& title      ,
        RooAbsReal&        xvar       ,
        const RooArgList&  pars       ,
        const bool         convex     , 
        const double       xmin       , 
        const double       xmax       ,
        RooAbsReal&        a          ,
        RooAbsReal&        b          ) ;
      /// constructor from the variable, range and list of coefficients
      ConvexOnly
      ( const std::string& name       ,
        const std::string& title      ,
        RooAbsReal&        xvar       ,
        const RooArgList&  pars       ,
        const bool         convex     , 
        const double       xmin       , 
        const double       xmax       ,
        const double       a = 0      , 
        const double       b = 1      ) ;
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
        const char* rangeName = nullptr ) const override ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::ConvexOnly& function () const { return m_convex ; }
      const Ostap::Math::ConvexOnly& convex   () const { return m_convex ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the variable
      const RooAbsReal& xvar   () const { return m_xvar.arg() ; }
      /// get parameters 
      const RooArgList& pars   () const { return m_pars       ; }
      /// get the shift  
      const RooAbsReal& a      () const { return m_a.arg()    ; }
      /// get the scale   
      const RooAbsReal& b      () const { return m_b.arg()    ; }
      /// xmin for bernstein 
      double xmin () const { return m_convex.xmin() ; }
      /// xmax for bernstein 
      double xmax () const { return m_convex.xmax() ; }
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
    /** @class BSpline
     *  The basic spline   ("B-spline")
     *  @see http://en.wikipedia.org/wiki/B-spline
     *  @see http://link.springer.com/chapter/10.1007%2F978-3-0348-7692-6_6
     *  @see Ostap::Math::Bspline 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2022-07-27     
     */
    class BSpline : public RooAbsReal 
    {
      // ======================================================================
      ClassDefOverride ( Ostap::MoreRooFit::BSpline , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the variable, range and the list of coefficients
      BSpline
      ( const std::string&         name  ,
        const std::string&         title ,
        RooAbsReal&                xvar  ,
        const std::vector<double>& knots ,
        const RooArgList&          pars  ) ;
      // ======================================================================
      /// copy constructor 
      BSpline
      ( const BSpline& right , 
        const char*    name = nullptr ) ;
      // ======================================================================
      BSpline () ;
      virtual ~BSpline () ;
      // ======================================================================
      BSpline* clone ( const char* newname ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&  allVars             , 
        RooArgSet&  analVars            , 
        const char* rangeName = nullptr ) const override ;
      Double_t    analyticalIntegral
      ( Int_t code , 
        const char* rangeName = nullptr ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the varibale 
      const RooAbsReal& xvar   () const { return m_xvar.arg() ; }
      /// get parameters 
      const RooArgList& pars   () const { return m_pars       ; }
      /// vector of knots 
      const std::vector<double>& knots () const { return m_bspline.knots() ; }
      /// xmin
      double            xmin   () const { return m_bspline.xmin() ; }
      /// xmax 
      double            xmax   () const { return m_bspline.xmax() ; }
      /// degree 
      unsigned short    degree () const { return m_bspline.degree () ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars       () const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::BSpline& function () const { return m_bspline ; }
      const Ostap::Math::BSpline& bspline  () const { return m_bspline ; }
      // ======================================================================
    public: 
      // ======================================================================
      Double_t evaluate  () const override ;
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy                  m_xvar    {} ;
      RooListProxy                  m_pars    {} ;
      mutable Ostap::Math::BSpline  m_bspline {} ; 
      // ======================================================================
    } ; //                          The end of class Ostap::MoreRooFit::BSpline 
    // ========================================================================
    /** Helper method to check if recursive fractions were 
     *  used for creation of RooAddPdf object
     *  @see RooAddPdf 
     */
    bool       recursive 
    ( const RooAddPdf& pdf ) ;
    // ========================================================================
    /** get the original fractions from the <code>RooAddPdf</code>
     *  @see RooAddPdf
     */
    RooArgList fractions
    ( const RooAddPdf& pdf       , 
      bool&            recursive ) ;  
    // ========================================================================
    /** get the original fractions from the <code>RooAddPdf</code>
     *  @see RooAddPdf
     */
    RooArgList fractions
    ( const RooAddPdf& pdf       ) ;
    // ========================================================================
    /** get x-observable
     *  @see RooGauissian
     */
    const RooAbsReal& getX    ( const RooGaussian& pdf ) ;
    // ========================================================================
    /** get mean value 
     *  @see RooGauissian
     */
    const RooAbsReal& getMean ( const RooGaussian& pdf ) ;
    // ========================================================================
    /** get sigma 
     *  @see RooGauissian
     */
    const RooAbsReal& getSigma ( const RooGaussian& pdf ) ;
    // ========================================================================
  } //                                   The end of namespace Ostap::MoreRooFit
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MOREVARS_H
// ============================================================================
