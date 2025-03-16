// ============================================================================
#ifndef OSTAP_SPLOT_H 
#define OSTAP_SPLOT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/COWs.h"
#include "Ostap/FitResult.h"
// ============================================================================
// forward declaratios 
// ============================================================================
// ROOT/RooFit 
// ============================================================================
class RooAddPdf  ; // ROOT/RooFit 
class RooAbsData ; // ROOT/Roofit 
class RooArgSet  ; // ROOT/Roofit 
class RooArgList ; // ROOT/Roofit 
class RooAbsData ; // ROOT/Roofit 
class RooAbsData ; // ROOT/Roofit 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Utils
  {
    // ========================================================================
   /** @class SPLOT
     *  Hellper class to add sPlot results to the TTree 
     */
    class SPLOT : public COWs
    {
      // ======================================================================
    public :
      // ======================================================================
      /** @param addpdf input extended RooAddPdf 
       *  @param observabels list of observables 
       *  @param normalzation normalisation set 
       */
      SPLOT
      ( const RooAddPdf&        addpdf                  ,
        const RooAbsCollection& observables             ,
        const RooFitResult&     fitresult               , 
        const RooAbsCollection* normalization = nullptr ) ;
      // ======================================================================
      // copy constructor 
      SPLOT ( const SPLOT&  right ) ;
      // move constructor 
      SPLOT (       SPLOT&& right ) = default ;
      // ======================================================================
      // destructore 
      virtual ~SPLOT () ;
      // ======================================================================
      SPLOT* clone () const override ;
      // ======================================================================       
    public:
      // ======================================================================
      /// get list of coefficiencts  
      const RooArgList&  coefficients () const { return *m_coefs  ; }
      /// fit result
      const RooFitResult& fitresult   () const { return *m_result ; }
      // ======================================================================      
    private:
      // ======================================================================
      /// coefficients
      std::unique_ptr<RooArgList>               m_coefs  {} ; // coefficients
      /// fir result
      std::unique_ptr<Ostap::Utils::FitResults> m_result {} ; // 
      // ======================================================================    
    } ; // ====================================================================
    // ========================================================================
  } //                                   The END of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The END of namespace Ostap
// ============================================================================
#endif // OSTAP_SPLOT_H 
// ============================================================================
//                                                                      The END 
// ============================================================================
