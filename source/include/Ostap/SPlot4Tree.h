// ============================================================================
#ifndef OSTAP_SPLOT4TREE_H 
#define OSTAP_SPLOT4TREE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/RooFun.h"
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
  /** @namespace Ostap::MoreRooFit   Ostap/MoreRooFit.h
   *  Collection of small additions to RooFit 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
   *  @date 2019-11-21
   */
  namespace MoreRooFit 
  {
    // ========================================================================
   /** @class SPlot4Tree 
     *  Hellper class to add sPlot results to the TTree 
     */
    class SPlot4Tree : public RooFun 
    {
      // ======================================================================
    public :
      // ======================================================================
      /** @param addpdf input extended RooAddPdf 
       *  @param observabels list of observables 
       *  @param normalzation normalisation set 
       */
      SPlot4Tree
      ( const RooAddPdf&        addpdf                  ,
        const RooAbsCollection& observables             ,
        const RooFitResult&     fitresult               , 
        const RooAbsCollection* normalization = nullptr ) ;
      // ======================================================================
      // copy constructor 
      SPlot4Tree ( const SPlot4Tree&  right ) ;
      // move constructor 
      SPlot4Tree (       SPlot4Tree&& right ) = default ;
      // ======================================================================
      // destructore 
      virtual ~SPlot4Tree() ;
      // ======================================================================      
    public:
      // ======================================================================
      /// get the pdf 
      const RooAddPdf&    pdf          () const ; 
      /// g`et list of components
      const RooArgList&   components   () const { return *m_cmps   ; }
      /// get list of coefficiencts  
      const RooArgList&   coefficients () const { return *m_coefs  ; }
      /// fit result
      const RooFitResult& fitresult    () const { return *m_result ; }
      /// size of object: number of componnents
      std::size_t         size         () const ; 
      // ======================================================================      
    private:
      // ======================================================================
      /// components 
      std::unique_ptr<RooArgList>               m_cmps   {} ; // components 
      /// coefficients
      std::unique_ptr<RooArgList>               m_coefs  {} ; // coefficients 
      /// fir result
      std::unique_ptr<Ostap::Utils::FitResults> m_result {} ; // 
      // ======================================================================    
    } ; // ====================================================================
    // ========================================================================
  } //                                   The END of namespace Ostap::MoreRooFit
  // ==========================================================================
} //                                                 The END of namespace Ostap
// ============================================================================
#endif // OSTAP_SPLOT4TREE_H 
// ============================================================================
//                                                                      The END 
// ============================================================================
