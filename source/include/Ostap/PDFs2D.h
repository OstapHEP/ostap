// ============================================================================
#ifndef OSTAP_PDFS2D_H 
#define OSTAP_PDFS2D_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Models2D.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/BSpline.h"
#include "Ostap/Peaks.h"
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
/** @file Ostap/PDFs2D.h
 *  Collection of non-facrorizeable 2D-models 
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Models
  {
    // ========================================================================
    /** @class Poly2DPositive 
     *  Poly2DPositive polynomial
     *  @see Ostap::Math::Positive2D
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  Poly2DPositive: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Poly2DPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// linear
      Poly2DPositive
      ( const char*          name      ,
        const char*          title     ,
        RooRealVar&          x         ,
        RooRealVar&          y         ,
        const unsigned short nX        ,
        const unsigned short nY        ,
        RooArgList&          phis      ) ; // at least (n+1)*(n+2)-1 elements
      /// copy
      Poly2DPositive
      ( const Poly2DPositive&     right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Poly2DPositive() ;
      /// clone
      Poly2DPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Poly2DPositive  () {} ;
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
      const Ostap::Math::Positive2D&  function  () const { return m_positive ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Positive2D m_positive ;              // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Poly2DSymPositive
     *  Poly2DSymPositive polynomial
     *  @see Ostap::Math::Positive2D
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  Poly2DSymPositive: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Poly2DSymPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// linear
      Poly2DSymPositive
      ( const char*          name      ,
        const char*          title     ,
        RooRealVar&          x         ,
        RooRealVar&          y         ,
        const unsigned short n         ,
        RooArgList&          phis      ) ; // (n+1)*(n+2)/2-1 elements
      /// copy
      Poly2DSymPositive
      ( const Poly2DSymPositive&     right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Poly2DSymPositive() ;
      /// clone
      Poly2DSymPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Poly2DSymPositive  () {} ;
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
      const Ostap::Math::Positive2DSym& function() const { return m_positive ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Positive2DSym m_positive ;           // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PS2DPol
     *  The 2D-function, that represent a cross-product of two phase-space factors,
     *  \f$ Ps_x(x)\f$ and \f$ Ps_y(y)\f$,  modulated by the 2D-positive polynomial
     *  The function is:
     *  \f[ f(x,y) = Ps_{x}(x) Ps_{y}(y) P_{pos}(x,y)\f], where 
     *  - \f$ Ps_x(x)\f$ is 1D phase-space function 
     *  - \f$ Ps_y(y)\f$ is 1D phase-space function 
     *  - \f$ P_{pos}(x,y) \f$ is 2D positive Bernstein polynomial 
     *  @see Ostap::Math::PS2DPol
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2D
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  PS2DPol: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PS2DPol,2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      PS2DPol
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PhaseSpaceNL& psx        ,
        const Ostap::Math::PhaseSpaceNL& psy        ,
        const unsigned short             nX         ,
        const unsigned short             nY         ,
        RooArgList&                      phis       ) ; // at least (nX+1)*(nY+1)-1 elements
      /// constructor
      PS2DPol
      ( const char*                      name      ,
        const char*                      title     ,
        RooRealVar&                      x         ,
        RooRealVar&                      y         ,
        const Ostap::Math::PS2DPol&      ps        ,
        RooArgList&                      phis      ) ; // at least (nX+1)*(nY+1)-1 elements
      /// copy
      PS2DPol
      ( const PS2DPol&       right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~PS2DPol() ;
      /// clone
      PS2DPol* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PS2DPol  () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::PS2DPol&      function    () const { return m_function ; }
      const Ostap::Math::Positive2D&   positive    () const { return m_function.positive   () ; }
      const Ostap::Math::Positive2D&   polynom     () const { return m_function.positive   () ; }
      const Ostap::Math::PhaseSpaceNL& psX         () const { return m_function.phasespaceX() ; }
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_function.phasespaceY() ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::PS2DPol m_function ;              // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PS2DPolSym
     *  The symmetric 2D-function, that represent a cross-product of two identical 
     *  phase-space factors,
     *  \f$ Ps(x)\f$ and \f$ Ps(y)\f$,  modulated by the symmetric 2D-positive 
     *  polynomial.
     *
     *  The function is:
     *  \f[ f(x,y) = Ps(x) Ps(y) P_{pos}(x,y)\f], where 
     *  - \f$ Ps(x)\f$ is 1D phase-space function
     *  - \f$ P_{pos}(x,y) \f$ is symmetric 2D positive Bernstein polynomial 
     *
     * Clearly the function is symmetric under 
     * \f$ x\leftrightarrow y \f$ transformmation: \f$f(x,y) = f(y,x) \f$ 
     *  @see Ostap::Math::PS2DPolSym
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2DSym
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  PS2DPolSym: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PS2DPolSym, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      PS2DPolSym
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PhaseSpaceNL& ps         ,
        const unsigned short N                      , 
        RooArgList&                      phis       ) ; // at least (nX+1)*(nY+1)-1 elements
      /// constructor
      PS2DPolSym
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PS2DPolSym&   ps         ,
        RooArgList&                      phis       ) ; // at least (nX+1)*(nY+1)-1 elements
      /// copy
      PS2DPolSym
        ( const PS2DPolSym&    right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~PS2DPolSym() ;
      /// clone
      PS2DPolSym* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PS2DPolSym () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::PS2DPolSym&    function    () const { return m_function ; }
      const Ostap::Math::Positive2DSym& positive    () const { return m_function.positive   () ; }
      const Ostap::Math::Positive2DSym& polynom     () const { return m_function.positive   () ; }
      const Ostap::Math::PhaseSpaceNL&  psX         () const { return m_function.phasespaceX() ; }
      const Ostap::Math::PhaseSpaceNL&  psY         () const { return m_function.phasespaceY() ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceY () const { return psY () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::PS2DPolSym m_function ;              // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PS2DPol2
     *  The 2D-function, that represent a cross-product of two phase-space factors,
     *  \f$ Ps_x(x)\f$ and \f$ Ps_y(y)\f$,  modulated by the 2D-positive polynomial
     *  The function is:
     *  \f[ f(x,y) = Ps_{x}(x) Ps_{y}(y) P_{pos}(x,y)\f], where 
     *  - \f$ Ps_x(x)\f$ is 1D phase-space function 
     *  - \f$ Ps_y(y)\f$ is 1D phase-space function 
     *  - \f$ P_{pos}(x,y) \f$ is 2D positive Bernstein polynomial 
     *  @see Ostap::Math::PS2DPol
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2D
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  PS2DPol2: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PS2DPol2,2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      PS2DPol2
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PhaseSpaceNL& psx        ,
        const Ostap::Math::PhaseSpaceNL& psy        ,
        const double                     mmax       ,
        const unsigned short             nX         ,
        const unsigned short             nY         ,
        RooArgList&                      phis       ) ; // at least (nX+1)*(nY+1)-1 elements
      /// constructor
      PS2DPol2
      ( const char*                      name      ,
        const char*                      title     ,
        RooRealVar&                      x         ,
        RooRealVar&                      y         ,
        const Ostap::Math::PS2DPol2&     ps        ,
        RooArgList&                      phis      ) ; // at least (nX+1)*(nY+1)-1 elements
      /// copy
      PS2DPol2
      ( const PS2DPol2&      right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~PS2DPol2() ;
      /// clone
      PS2DPol2* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PS2DPol2 () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::PS2DPol2&     function    () const { return m_function ; }
      const Ostap::Math::Positive2D&   positive    () const { return m_function.positive   () ; }
      const Ostap::Math::Positive2D&   polynom     () const { return m_function.positive   () ; }
      const Ostap::Math::PhaseSpaceNL& psX         () const { return m_function.phasespaceX() ; }
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_function.phasespaceY() ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::PS2DPol2 m_function ;                // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PS2DPol2Sym
     *  The symmetric 2D-function, that represent a cross-product of two identical 
     *  phase-space factors,
     *  \f$ Ps(x)\f$ and \f$ Ps(y)\f$,  modulated by the symmetric 2D-positive 
     *  polynomial.
     *
     *  The function is:
     *  \f[ f(x,y) = Ps(x) Ps(y) P_{pos}(x,y)\f], where 
     *  - \f$ Ps(x)\f$ is 1D phase-space function
     *  - \f$ P_{pos}(x,y) \f$ is symmetric 2D positive Bernstein polynomial 
     *
     * Clearly the function is symmetric under 
     * \f$ x\leftrightarrow y \f$ transformmation: \f$f(x,y) = f(y,x) \f$ 
     *  @see Ostap::Math::PS2DPolSym
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2DSym
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  PS2DPol2Sym: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PS2DPol2Sym, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      PS2DPol2Sym
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PhaseSpaceNL& ps         ,
        const double                     mmax       ,
        const unsigned short N                      , 
        RooArgList&                      phis       ) ; // at least (nX+1)*(nY+1)-1 elements
      /// constructor
      PS2DPol2Sym
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PS2DPol2Sym&  ps         ,
        RooArgList&                      phis       ) ; // at least (nX+1)*(nY+1)-1 elements
      /// copy
      PS2DPol2Sym
        ( const PS2DPol2Sym&   right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~PS2DPol2Sym() ;
      /// clone
      PS2DPol2Sym* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PS2DPol2Sym () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::PS2DPol2Sym&   function    () const { return m_function ; }
      const Ostap::Math::Positive2DSym& positive    () const { return m_function.positive   () ; }
      const Ostap::Math::Positive2DSym& polynom     () const { return m_function.positive   () ; }
      const Ostap::Math::PhaseSpaceNL&  psX         () const { return m_function.phasespaceX() ; }
      const Ostap::Math::PhaseSpaceNL&  psY         () const { return m_function.phasespaceY() ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceY () const { return psY () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::PS2DPol2Sym m_function ;             // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PS2DPol3
     *  The 2D-function, that represent a cross-product of two phase-space factors,
     *  \f$ Ps_x(x)\f$ and \f$ Ps_y(y)\f$,  modulated by the 2D-positive polynomial
     *  The function is:
     *  \f[ f(x,y) = Ps_{x}(x) Ps_{y}(y) P_{pos}(x,y)\f], where 
     *  - \f$ Ps_x(x)\f$ is 1D phase-space function 
     *  - \f$ Ps_y(y)\f$ is 1D phase-space function 
     *  - \f$ P_{pos}(x,y) \f$ is 2D positive Bernstein polynomial 
     *  @see Ostap::Math::PS2DPol
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2D
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  PS2DPol3: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PS2DPol3,2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      PS2DPol3
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PhaseSpaceNL& psx        ,
        const Ostap::Math::PhaseSpaceNL& psy        ,
        const double                     mmax       ,
        const unsigned short             nX         ,
        const unsigned short             nY         ,
        RooArgList&                      phis       ) ; // at least nX+nY elements
      /// constructor
      PS2DPol3
      ( const char*                       name      ,
        const char*                       title     ,
        RooRealVar&                       x         ,
        RooRealVar&                       y         ,
        const Ostap::Math::PS2DPol3&      ps        ,
        RooArgList&                       phis      ) ; // at least nX+nY elements
      /// copy
      PS2DPol3
      ( const PS2DPol3&      right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~PS2DPol3() ;
      /// clone
      PS2DPol3* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PS2DPol3 () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::PS2DPol3&      function    () const { return m_function ; }
      const Ostap::Math::PhaseSpacePol& psX         () const { return m_function.phasespaceX() ; }
      const Ostap::Math::PhaseSpacePol& psY         () const { return m_function.phasespaceY() ; }
      const Ostap::Math::PhaseSpacePol& phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpacePol& phasespaceY () const { return psY () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::PS2DPol3 m_function ;                // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PS2DPol3Sym
     *  The symmetric 2D-function, that represent a cross-product of two identical 
     *  phase-space factors,
     *  \f$ Ps(x)\f$ and \f$ Ps(y)\f$,  modulated by the symmetric 2D-positive 
     *  polynomial.
     *
     *  The function is:
     *  \f[ f(x,y) = Ps(x) Ps(y) P_{pos}(x,y)\f], where 
     *  - \f$ Ps(x)\f$ is 1D phase-space function
     *  - \f$ P_{pos}(x,y) \f$ is symmetric 2D positive Bernstein polynomial 
     *
     * Clearly the function is symmetric under 
     * \f$ x\leftrightarrow y \f$ transformmation: \f$f(x,y) = f(y,x) \f$ 
     *  @see Ostap::Math::PS2DPolSym
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2DSym
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  PS2DPol3Sym: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::PS2DPol3Sym, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      PS2DPol3Sym
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PhaseSpaceNL& ps         ,
        const double                     mmax       ,
        const unsigned short             N          , 
        RooArgList&                      phis       ) ; // at least N elements
      /// constructor
      PS2DPol3Sym
      ( const char*                      name       ,
        const char*                      title      ,
        RooRealVar&                      x          ,
        RooRealVar&                      y          ,
        const Ostap::Math::PS2DPol3Sym&  ps         ,
        RooArgList&                      phis       ) ; // at least N elements
      /// copy
      PS2DPol3Sym
        ( const PS2DPol3Sym&   right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~PS2DPol3Sym() ;
      /// clone
      PS2DPol3Sym* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      PS2DPol3Sym () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::PS2DPol3Sym&   function    () const { return m_function ; }
      const Ostap::Math::PhaseSpacePol& psX         () const { return m_function.phasespaceX() ; }
      const Ostap::Math::PhaseSpacePol& psY         () const { return m_function.phasespaceY() ; }
      const Ostap::Math::PhaseSpacePol& phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpacePol& phasespaceY () const { return psY () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::PS2DPol3Sym m_function ;             // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ExpoPS2DPol
     *
     *  F(x,y) = exp(x)*PS(y)*PPOL(x,y)
     *
     *  @see Ostap::Math::ExpoPS2DPol
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  ExpoPS2DPol: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::ExpoPS2DPol, 2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      ExpoPS2DPol
        ( const char*                      name      ,
          const char*                      title     ,
          RooRealVar&                      x         ,
          RooRealVar&                      y         ,
          RooAbsReal&                      tau       ,
          const Ostap::Math::PhaseSpaceNL& psy       ,
          const unsigned short nX                    ,
          const unsigned short nY                    ,
          RooArgList&          phis      ) ; // at least (nX+1)*(nY+1)-1 elements
      /// constructor
      ExpoPS2DPol
        ( const char*                      name      ,
          const char*                      title     ,
          RooRealVar&                      x         ,
          RooRealVar&                      y         ,
          RooAbsReal&                      tau       ,
          const Ostap::Math::ExpoPS2DPol&  ps        ,
          RooArgList&          phis      ) ; // at least (nX+1)*(nY+1)-1 elements
      /// copy
      ExpoPS2DPol
        ( const ExpoPS2DPol&   right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~ExpoPS2DPol() ;
      /// clone
      ExpoPS2DPol* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      ExpoPS2DPol () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::ExpoPS2DPol&  function    () const { return m_function ; }
      const Ostap::Math::Positive2D&   positive    () const { return m_function.positive   () ; }
      const Ostap::Math::Positive2D&   polynom     () const { return m_function.positive   () ; }
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_function.phasespaceY() ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooRealProxy m_tau  ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::ExpoPS2DPol m_function ;             // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Expo2DPol
     *
     *  F(x,y) = exp(x)*exp(y)*PPOL(x,y)
     *
     *  @see Ostap::Math::ExpoPS2DPol
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  Expo2DPol: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Expo2DPol, 2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      Expo2DPol
        ( const char*                      name      ,
          const char*                      title     ,
          RooRealVar&                      x         ,
          RooRealVar&                      y         ,
          RooAbsReal&                      taux      ,
          RooAbsReal&                      tauy      ,
          const unsigned short nX                    ,
          const unsigned short nY                    ,
          RooArgList&          phis      ) ; // at least (nX+1)*(nY+1)-1 elements
      /// copy
      Expo2DPol
        ( const Expo2DPol&     right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~Expo2DPol() ;
      /// clone
      Expo2DPol* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Expo2DPol () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::Expo2DPol&    function    () const { return m_function ; }
      const Ostap::Math::Positive2D&   positive    () const { return m_function.positive   () ; }
      const Ostap::Math::Positive2D&   polynom     () const { return m_function.positive   () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooRealProxy m_taux ;
      RooRealProxy m_tauy ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::Expo2DPol m_function ;               // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Expo2DPolSym
     *
     *  F(x,y) = exp(x)*exp(y)*SPOL(x,y)
     *
     *  @see Ostap::Math::Expo2DPolSym
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class  Expo2DPolSym: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Expo2DPolSym, 2) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor
      Expo2DPolSym
        ( const char*                      name      ,
          const char*                      title     ,
          RooRealVar&                      x         ,
          RooRealVar&                      y         ,
          RooAbsReal&                      tau       ,
          const unsigned short n                     ,
          RooArgList&          phis      ) ; // at least (nX+1)*(nY+1)-1 elements
      /// copy
      Expo2DPolSym
        ( const Expo2DPolSym&  right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~Expo2DPolSym() ;
      /// clone
      Expo2DPolSym* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Expo2DPolSym () {} ;
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
      /// access to underlying function(s)
      const Ostap::Math::Expo2DPolSym&  function    () const { return m_function ; }
      const Ostap::Math::Positive2DSym& positive    () const { return m_function.positive   () ; }
      const Ostap::Math::Positive2DSym& polynom     () const { return m_function.positive   () ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooRealProxy m_tau  ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual functions()
      mutable Ostap::Math::Expo2DPolSym m_function ;           // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Spline2D
     *  Positive 2D-spline
     *  @see Ostap::Math::Spline2D
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class Spline2D: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Spline2D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// linear
      Spline2D
      ( const char*                          name      ,
        const char*                          title     ,
        RooRealVar&                          x         ,
        RooRealVar&                          y         ,
        const Ostap::Math::PositiveSpline2D& spline    ,
        RooArgList&                          phis      ) ;
      /// copy
      Spline2D
      ( const Spline2D&      right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Spline2D() ;
      /// clone
      Spline2D* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Spline2D () {} ;
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
      const Ostap::Math::PositiveSpline2D& function() const { return m_spline ; }
      const Ostap::Math::PositiveSpline2D& spline  () const { return m_spline ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PositiveSpline2D m_spline ;         // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Spline2DSym
     *  Positive symmetric 2D-spline
     *  @see Ostap::Math::Spline2DSym
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2011-05-25
     */
    class Spline2DSym: public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDef(Ostap::Models::Spline2DSym, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// linear
      Spline2DSym
      ( const char*                             name      ,
        const char*                             title     ,
        RooRealVar&                             x         ,
        RooRealVar&                             y         ,
        const Ostap::Math::PositiveSpline2DSym& spline    ,
        RooArgList&                             phis      ) ;
      /// copy
      Spline2DSym
      ( const Spline2DSym&   right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Spline2DSym() ;
      /// clone
      Spline2DSym* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Spline2DSym () {} ;
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
      const Ostap::Math::PositiveSpline2DSym& function() const { return m_spline ; }
      const Ostap::Math::PositiveSpline2DSym& spline  () const { return m_spline ; }
      // ======================================================================
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const ;
      // ======================================================================
#endif
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::PositiveSpline2DSym m_spline ;       // the function
      // ======================================================================
    } ;
    // ========================================================================
  } //                                      The end of  namespace Ostap::Models
  // ==========================================================================
} //                                              The end of namespace Analysis 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PDFS2D_H
// ============================================================================
