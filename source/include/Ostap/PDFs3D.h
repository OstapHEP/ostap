// ============================================================================
#ifndef OSTAP_PDFS3D_H 
#define OSTAP_PDFS3D_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Bernstein3D.h"
#include "Ostap/Models3D.h"
// ============================================================================
// ROOT
// ============================================================================
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsReal.h"
// ============================================================================
/** @file Ostap/PDFs3D.h
 *  Collection of non-facrorizeable 3D-models 
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Models 
  {
    // ========================================================================
    /** @class Poly3DPositive
     *  The 3D-polynomial of order Nx*Ny*Nz, that is constrained 
     *  to be non-negative over the  defined range      
     *  \f[  P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n_x}_i(x) B^{n_y}_j(y) B^{n_z}_k(z)\f] 
     *  where all coefficients \f$a_{ijk}\f$ are non-negative and 
     *  \f$ \sum_{i,j,k} a_{ijk}=1 \f$ 
     *  @author Vanya BELYAEV Ivan.Belayev@itep.ru
     *  @date 2017-11-14
     *  @see Ostap::Math::Positive3D
     *  @see Ostap::Math::Bernstein3D
     */
    class  Poly3DPositive : public RooAbsPdf
    {
      // ======================================================================
      ClassDefOverride(Ostap::Models::Poly3DPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// linear
      Poly3DPositive
      ( const char*          name      ,
        const char*          title     ,
        RooRealVar&          x         ,
        RooRealVar&          y         ,
        RooRealVar&          z         ,
        const unsigned short nX        ,
        const unsigned short nY        ,
        const unsigned short nZ        ,
        RooArgList&          phis      ) ; // at least (nx+1)*(ny+1)*(nz+1)-1 elements
      /// copy
      Poly3DPositive
        ( const Poly3DPositive&     right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~Poly3DPositive() ;
      /// clone
      Poly3DPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Poly3DPositive  () {} ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public:  // integrals
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
      const Ostap::Math::Positive3D&  function  () const { return m_positive ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x  .arg() ; }
      const RooAbsReal& y    () const { return m_y  .arg() ; }
      const RooAbsReal& z    () const { return m_z  .arg() ; }
      const RooArgList& phis () const { return m_phis      ; }      
      //
      unsigned short    nX   () const { return m_positive.nX() ; }
      unsigned short    nY   () const { return m_positive.nY() ; }
      unsigned short    nZ   () const { return m_positive.nZ() ; }      
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooRealProxy m_z    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Positive3D m_positive ;              // the function
      // ======================================================================
    } ;
    // ========================================================================   
    /** @class Poly3DSymPositive
     *  The 3D-polynomial of order N*N*N, that is constrained 
     *  to be non-negative ans symmetric over the  defined range      
     *  \f[  P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n}_k(z)\f] 
     *  where all coefficients \f$a_{ijk}\f$ are:
     * - non-negative: \f$ a_{ijk}\ge0 \f$
     * - symmetric:    \f$ a_{ijk}=a_{jik}=a_{ikj}\f$
     * - constrainted: \f$ \sum_{i,j,k} a_{ijk}=1 \f$ 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2017-11-14
     *  @see Ostap::Math::Positive3DSym
     *  @see Ostap::Math::Bernstein3DSym
     */
    class  Poly3DSymPositive : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Poly3DSymPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// linear
      Poly3DSymPositive
        ( const char*          name      ,
          const char*          title     ,
          RooRealVar&          x         ,
          RooRealVar&          y         ,
          RooRealVar&          z         ,
          const unsigned short n         ,
          RooArgList&          phis      ) ; // at least (n+1)*(n+2)*(n+3)/6-1 elements
      /// copy
      Poly3DSymPositive
        ( const Poly3DSymPositive&     right     ,
          const char*          name = 0  ) ;
      /// destructor
      virtual ~Poly3DSymPositive() ;
      /// clone
      Poly3DSymPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Poly3DSymPositive  () {} ;
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
      const Ostap::Math::Positive3DSym& function() const { return m_positive ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x  .arg() ; }
      const RooAbsReal& y    () const { return m_y  .arg() ; }
      const RooAbsReal& z    () const { return m_z  .arg() ; }
      const RooArgList& phis () const { return m_phis      ; }      
      // ======================================================================
      unsigned short    nX   () const { return m_positive.nX() ; }
      unsigned short    nY   () const { return m_positive.nY() ; }
      unsigned short    nZ   () const { return m_positive.nZ() ; }      
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooRealProxy m_z    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Positive3DSym m_positive ;           // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Poly3DMixPositive
     *  The 3D-polynomial of order N*N*Nz, that is constrained 
     *  to be non-negative and symmetric for \f$ x \leftrightarrow y\f$ interchange 
     *  over the  defined range      
     *  \f[  P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n_z}_k(z)\f] 
     *  where all coefficients \f$a_{ijk}\f$ are:
     * - non-negative: \f$ a_{ijk}\ge0 \f$
     * - symmetric for \f$ x \leftrightarrow y \f$ interchange: \f$a_{ijk}=a_{jik}\f$
     * - constrainted: \f$ \sum_{i,j,k} a_{ijk}=1 \f$ 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2017-11-14
     *  @see Ostap::Math::Positive3DMix
     *  @see Ostap::Math::Bernstein3DMix
     */
    // ========================================================================
    class  Poly3DMixPositive : public RooAbsPdf
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Poly3DMixPositive, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// linear
      Poly3DMixPositive
      ( const char*          name      ,
        const char*          title     ,
        RooRealVar&          x         ,
        RooRealVar&          y         ,
        RooRealVar&          z         ,
        const unsigned short n         ,
        const unsigned short nz        ,
        RooArgList&          phis      ) ; // at least (n+1)*(n+2)*(nz+1)/2-1 elements
      /// copy
      Poly3DMixPositive
      ( const Poly3DMixPositive&     right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Poly3DMixPositive() ;
      /// clone
      Poly3DMixPositive* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Poly3DMixPositive  () {} ;
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
      const Ostap::Math::Positive3DMix& function() const { return m_positive ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x    () const { return m_x  .arg() ; }
      const RooAbsReal& y    () const { return m_y  .arg() ; }
      const RooAbsReal& z    () const { return m_z  .arg() ; }
      const RooArgList& phis () const { return m_phis      ; }      
      // ======================================================================
      unsigned short    nX   () const { return m_positive.nX() ; }
      unsigned short    nY   () const { return m_positive.nY() ; }
      unsigned short    nZ   () const { return m_positive.nZ() ; }      
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x    ;
      RooRealProxy m_y    ;
      RooRealProxy m_z    ;
      RooListProxy m_phis ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Positive3DMix m_positive ;           // the function
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Gauss3D 
     *  Simple 3D rotated Gaussian function 
     *  @see Ostap::Math::Gauss3D
     *  @date 2022-06-22
     */
    class Gauss3D: public RooAbsPdf
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Models::Gauss3D, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      ///
      Gauss3D
      ( const char* name      ,
        const char* title     ,
        RooAbsReal& x         ,
        RooAbsReal& y         ,
        RooAbsReal& z         ,
        RooAbsReal& muX       , 
        RooAbsReal& muY       , 
        RooAbsReal& muZ       , 
        RooAbsReal& sigmaX    , 
        RooAbsReal& sigmaY    , 
        RooAbsReal& sigmaZ    , 
        RooAbsReal& phi       ,
        RooAbsReal& theta     , 
        RooAbsReal& psi       ) ;
      /// copy
      Gauss3D
      ( const Gauss3D&      right     ,
        const char*          name = 0  ) ;
      /// destructor
      virtual ~Gauss3D() ;
      /// clone
      Gauss3D* clone ( const char* name ) const override;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      // fake default contructor, needed just for proper (de)serialization
      Gauss3D () {} ;
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
      const Ostap::Math::Gauss3D& function() const { return m_gauss3D ; }
      /// access to underlying function
      const Ostap::Math::Gauss3D& gauss3d () const { return m_gauss3D ; }
      // ======================================================================
    public:
      // ======================================================================
      const RooAbsReal& x       () const { return m_x     .arg() ; }
      const RooAbsReal& y       () const { return m_y     .arg() ; }
      const RooAbsReal& z       () const { return m_z     .arg() ; }
      const RooAbsReal& muX     () const { return m_muX   .arg() ; }
      const RooAbsReal& muY     () const { return m_muY   .arg() ; }
      const RooAbsReal& muZ     () const { return m_muZ   .arg() ; }
      const RooAbsReal& sigmaX  () const { return m_sigmaX.arg() ; }      
      const RooAbsReal& sigmaY  () const { return m_sigmaY.arg() ; }      
      const RooAbsReal& sigmaZ  () const { return m_sigmaZ.arg() ; }      
      const RooAbsReal& phi     () const { return m_phi   .arg() ; }      
      const RooAbsReal& theta   () const { return m_theta .arg() ; }      
      const RooAbsReal& psi     () const { return m_psi   .arg() ; }      
      // ======================================================================
    protected :
      // ======================================================================
      RooRealProxy m_x      ;
      RooRealProxy m_y      ;
      RooRealProxy m_z      ;
      RooRealProxy m_muX    ;
      RooRealProxy m_muY    ;
      RooRealProxy m_muZ    ;
      RooRealProxy m_sigmaX ;
      RooRealProxy m_sigmaY ;
      RooRealProxy m_sigmaZ ;
      RooRealProxy m_phi    ;
      RooRealProxy m_theta  ;
      RooRealProxy m_psi    ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual function
      mutable Ostap::Math::Gauss3D m_gauss3D ;           // the function
      // ======================================================================
    } ;
    // ========================================================================
  } //                                    The end of namespace Ostap::Models
  // ==========================================================================
} //                                              The end of namespace Analysis
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // ANALYSIS_MODELS3D_H
// ============================================================================
