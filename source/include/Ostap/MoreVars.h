// ============================================================================
#ifndef OSTAP_MOREVARS_H 
#define OSTAP_MOREVARS_H 1
// ============================================================================
// Include  files 
// ============================================================================
// ROOT
// ============================================================================
#include "TVectorDfwd.h"
#include "RVersion.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooProfileLL.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein1D.h"
#include "Ostap/BSpline.h"
#include "Ostap/Rational.h"
#include "Ostap/HistoInterpolators.h"
// ============================================================================
/// forward declarations 
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
      Bernstein 
      ( const Bernstein& right , 
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
      const RooAbsReal& x      () const { return m_xvar.arg() ; }
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
      const RooAbsReal& x      () const { return m_xvar.arg() ; }
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
      const RooAbsReal& x      () const { return m_xvar.arg() ; }
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
      const RooAbsReal& x      () const { return m_xvar.arg() ; }
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
      const RooAbsReal& x      () const { return m_xvar.arg() ; }
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
    /** @class Rational 
     *  A simple pole-free rational function at interval \f$ x_{min} \le x \le x_{max}\f$
     *  \f[ F(x) = \frac{p(x)}{q(x)} \f]
     *  Actually internally it uses 
     *  the Floater-Hormann rational barycentric interpolant 
     *  and parameters are the function valeus at Chebyshev's nodes  
     *   
     *  @see Ostap::Math::Rational
     *  @see Ostap::Math::FloaterHormann
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2023-09-21
     */
    class Rational final : public RooAbsReal 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::Rational, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      Rational 
        ( const std::string&   name  , 
          const std::string&   title , 
          RooAbsReal&          xvar  ,
          const RooArgList&    pars  ,
          const unsigned short d     ,
          const double         xmin  , 
          const double         xmax  ) ;
      /// copy constructor 
      Rational ( const Rational& right , const char* name = nullptr ) ;
      /// default 
      Rational () ;
      /// clone method
      Rational* clone ( const char* name ) const override ;
      /// virtual destructor 
      virtual ~Rational() ;
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
      /// get the variable/observable 
      const RooAbsReal& x      () const { return m_xvar.arg() ; }
      /// get the variable/observable 
      const RooAbsReal& xvar   () const { return m_xvar.arg() ; }
      /// get parameters 
      const RooArgList& pars   () const { return m_pars       ; }
      /// get n
      unsigned short    n      () const { return m_rational.n    () ; }
      /// get d 
      unsigned short    d      () const { return m_rational.d    () ; }
      /// get p (==d)  
      unsigned short    p      () const { return m_rational.d    () ; }
      /// xmin for Rational
      double            xmin   () const { return m_rational.xmin () ; }
      /// xmax for Rational 
      double            xmax   () const { return m_rational.xmax () ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars () const ;
      // ======================================================================
   public:
      // ======================================================================
      const Ostap::Math::Rational& function () const { return m_rational ; }
      const Ostap::Math::Rational& rational () const { return m_rational ; }
      // ======================================================================
    public: 
      // ======================================================================
      Double_t evaluate  () const override ;
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy                   m_xvar     {} ;
      RooListProxy                   m_pars     {} ;
      mutable Ostap::Math::Rational  m_rational {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class RationalBernstein 
     *  Rational fnuction as ratio of Bernstein polyhomial and 
     *  positive Bernstein polynomial
     *  \f[ R ( x ) = \frac{B(x)}{P(x)\frac{1} {x_{max} - x_{min} } \f]
     *  @see Ostap::Math::RationalBernstein
     *  @see Ostap::Math::Bernstein
     *  @see Ostap::Math::Positive 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2023-09-21
     */
    class RationalBernstein final : public RooAbsReal 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::RationalBernstein, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      RationalBernstein 
        ( const std::string&   name  , 
          const std::string&   title , 
          RooAbsReal&          xvar  ,
          const RooArgList&    p     , // numerator 
          const RooArgList&    q     , // denominator 
          const double         xmin  , 
          const double         xmax  ) ;
      RationalBernstein 
        ( const std::string&   name  , 
          const std::string&   title , 
          RooAbsReal&          xvar  ,
          const RooArgList&    pars  , // all pars 
          const unsigned short p     , // degree of numerator 
          const double         xmin  , 
          const double         xmax  ) ;
      /// copy constructor 
      RationalBernstein ( const RationalBernstein& right , const char* name = nullptr ) ;
      /// default 
      RationalBernstein () ;
      /// clone method
      RationalBernstein* clone ( const char* name ) const override ;
      /// virtual destructor 
      virtual ~RationalBernstein() ;
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
      /// get the variable/observable 
      const RooAbsReal& x      () const { return m_xvar.arg() ; }
      /// get the variable/observable 
      const RooAbsReal& xvar   () const { return m_xvar.arg() ; }
      /// get parameters 
      const RooArgList& pars   () const { return m_pars       ; }
      /// get p 
      unsigned short    p      () const { return m_rational.pdegree () ; }
      /// get q 
      unsigned short    q      () const { return m_rational.qdegree () ; }
      /// xmin for Rational
      double            xmin   () const { return m_rational.xmin    () ; }
      /// xmax for Rational 
      double            xmax   () const { return m_rational.xmax    () ; }
      // ======================================================================
    public:
      // ======================================================================
      void setPars () const ;
      // ======================================================================
   public:
      // ======================================================================
      const Ostap::Math::RationalBernstein& function () const { return m_rational ; }
      const Ostap::Math::RationalBernstein& rational () const { return m_rational ; }
      // ======================================================================
    public: 
      // ======================================================================
      Double_t evaluate  () const override ;
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy                            m_xvar     {} ;
      RooListProxy                            m_pars     {} ;
      mutable Ostap::Math::RationalBernstein  m_rational {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Shape1D
     *  The generic "fixed shape" function 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2023-01-30     
     */
    class Shape1D final : public RooAbsReal 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::Shape1D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
        Shape1D ( const char*  name  , 
                  const char*  title , 
                  RooAbsReal&  x     ,
                  FUNCTION     f     )
        : RooAbsReal (  name ,  title ) 
        , m_x        ( "x"   , "Variable" , this , x ) 
        , m_function ( f ) 
      {}
      /// copy constructor 
      Shape1D ( const Shape1D& right , const char* name = nullptr ) ;
      /// clone method
      Shape1D* clone ( const char* name ) const override ;
      /// virtual destructor 
      virtual ~Shape1D() ;
      // ======================================================================
    public:
      // ======================================================================
      /// templated constructor 
      template <class FUNCTION> 
        static inline Shape1D 
        create
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        x     ,
          FUNCTION           f     ) 
      { return Shape1D ( name.c_str () , title.c_str () , x , f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the PDF 
      Double_t evaluate () const override
      { const double x = m_x ; return m_function ( x ) ; }
      // ======================================================================        
    public:
      // ======================================================================
      /// evaluate the function
      double func  ( const double x ) const { return m_function ( x ) ; }
      // ======================================================================        
    private :
      // ======================================================================
      /// variable 
      RooRealProxy                   m_x        ; // variable 
      /// the function itself 
      std::function<double(double)>  m_function ; // function 
      // ======================================================================      
    } ; //                          The end of class Ostap::MoreRooFit::Shape1D
    // ========================================================================
    /** @class Histo1D
     *  simple generic function from the histogram 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2023-01-30     
     */
    class Histo1D final : public RooAbsReal 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::Histo1D, 1) ;
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
      /// virtual destructor 
      virtual ~Histo1D() ;
      // ======================================================================
    public:
      // ======================================================================
      // fake default contructor, needed just for the proper (de)serialization
      Histo1D () {} ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      Double_t evaluate () const override { return func ( m_x ) ; }
      // ======================================================================        
    public:
      // ======================================================================
      /// the function itself 
      const Ostap::Math::Histo1D& histo () const { return   m_histo ; }
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function
      double func  ( const double x ) const { return m_histo ( x )  ; }
      // ======================================================================        
    public:
      // ======================================================================
      const RooAbsReal& x () const { return m_x     .arg() ; }
      // ======================================================================
    private :
      // ======================================================================
      /// variable 
      RooRealProxy                   m_x     ; // variable 
      /// the function itself 
      Ostap::Math::Histo1D           m_histo ; // function 
      // ======================================================================      
    } ; //                          The end of class Ostap::MoreRooFit::Histo1D
    // ========================================================================
    /** @class ProfileLL 
     *  Slight extension for class RooProfileLL
     *  - Do no subtract the minimum...
     *  @see RooProfileLL
     */
    class ProfileLL : public RooProfileLL
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::ProfileLL,0) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      ProfileLL
      ( const char*      name        , 
	const char*      title       , 
	RooAbsReal&      nll         ,  
	const RooArgSet& observables ) ;
      /// "copy" constructor
      ProfileLL
      ( const RooProfileLL& right    , 
	const char*         name = nullptr ) ;
      /// virual destructor 
      virtual ~ProfileLL() ;
      /// clone method 
      ProfileLL* clone ( const char* newname ) const override ;
      // ======================================================================
#if ROOT_VERSION_CODE<ROOT_VERSION(6,35,1)
      // ============================================================================
    public:
      // ======================================================================
      /// default constructor 
      ProfileLL () {} ;
      // ============================================================================
#endif 
      // ======================================================================
    public:
      // ======================================================================
      using RooProfileLL::nll ;
      const RooAbsReal& nll () const { return _nll.arg()  ; }
      const RooArgSet&  obs () const { return _obs        ; }
      const RooArgSet&  par () const { return _par        ; }
      // ======================================================================
      /// min-value 
      double abs_min        () const { return _absMin       ; }
      /// Is min-valeu valid ? 
      bool   abs_min_valid  () const { return _absMinValid  ; }
      // ======================================================================
      /// main method: do not subtract min-value! 
      double evaluate  () const override ;
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
