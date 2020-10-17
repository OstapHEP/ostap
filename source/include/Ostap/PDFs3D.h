// ============================================================================
#ifndef OSTAP_PDFS3D_H 
#define OSTAP_PDFS3D_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Bernstein3D.h"
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
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const override ;
      // ======================================================================
#endif
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
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const override ;
      // ======================================================================
#endif
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
    private:
      // ======================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
      // ======================================================================
      RooSpan<double> evaluateBatch 
      ( std::size_t begin     , 
        std::size_t batchSize ) const override ;
      // ======================================================================
#endif
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
  } //                                    The end of namespace Ostap::Models
  // ==========================================================================
} //                                              The end of namespace Analysis
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // ANALYSIS_MODELS3D_H
// ============================================================================
