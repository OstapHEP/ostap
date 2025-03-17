// ============================================================================
#ifndef OSTAP_COWS_H 
#define OSTAP_COWS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
// ROOT 
// ============================================================================
#include "TMatrixDSym.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/RooFun.h"
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
  namespace Utils
  {
    // ========================================================================
    // forward declaration 
    class ProgressConf ;
    // ========================================================================
    /** @class COWs 
     * (over) simplified version of COWS
     *
     * @see H.Dembinski, M.Kenzie, C.Langenbruch and M.Schmelling, 
     *      "Custom Orthogonal Weight functions (COWs) for Event Classification",
     *      NIMA 1040 (2022) 167270, arXiv:2112.04574
     * @see https://doi.org/10.48550/arXiv.2112.04574
     * @see https://doi.org/10.1016/j.nima.2022.167270
     */
    class COWs : public RooFun 
    {
      // ======================================================================
    public :
      // ======================================================================
      /** Constuctor
       *  @param addpdf    input extended RooAddPdf 
       *  @param data      input data
       *  @param normalzation normalisation set
       *  @param progress  show progress bar? 
       */
      COWs
      ( const RooAddPdf&                  addpdf                  ,
	const RooAbsData&                 data                    ,
	const RooAbsCollection*           normalization           , 
        const Ostap::Utils::ProgressConf& progress      = false   ) ;
      // ======================================================================= 
      /** Constuctor
       *  @param addpdf    input extended RooAddPdf 
       *  @param data      input data
       *  @param normalzation normalisation set
       *  @param progress  show progress bar? 
       */
      COWs
      ( const RooAddPdf&                  addpdf                  ,
	const RooAbsData&                 data                    ,
        const Ostap::Utils::ProgressConf& progress      = false   ) ; 
      //   ======================================================================
      /** "Recovery constructor
       *  @param addpdf input extended RooAddPdf 
       *  @param observables  observables set  
       *  @param normalzation normalisation set 
       *  @param cows         the symmetric matrix of weights 
       */
      COWs
      ( const RooAddPdf&        addpdf        ,
	const RooAbsCollection& observables   , 
	const RooAbsCollection* normalization , 
	const TMatrixDSym&      cows          ) ;
      // ======================================================================
      /** "Recovery constructor
       *  @param addpdf input extended RooAddPdf 
       *  @param observables  observables set  
       *  @param normalzation normalisation set 
       *  @param cows         the symmetric matrix of weights 
       */
      COWs
      ( const RooAddPdf&        addpdf                  ,
	const RooAbsCollection& observables             , 
	const TMatrixDSym&      cows                    , 
	const RooAbsCollection* normalization = nullptr ) ; 	 
      // ======================================================================
      // copy constructor 
      COWs ( const COWs& right ) ;
      // move constructor 
      COWs (      COWs&& right ) = default ;
      // ======================================================================
      // destructore 
      virtual ~COWs() ;
      // ======================================================================
      /// clone funcntion/virtual constructire 
      COWs* clone() const override ;  
      // ======================================================================
    protected:    
      // ======================================================================      
      /** protected constructor fro SPLOT 
       *  @param addpdf input extended RooAddPdf 
       *  @param observables  observables set  
       *  @param normalzation normalisation set 
       *  @param cows         the symmetric matrix of weights 
       */
      COWs
      ( const RooAddPdf&        addpdf        ,
	const RooAbsCollection& observables   , 
	const RooAbsCollection* normalization ) ; 
      // ======================================================================
     public:
      // ======================================================================
      /// get the pdf 
      const RooAddPdf&   pdf          () const ; 
      /// get list of components
      const RooArgList&  components   () const { return *m_cmps   ; }
      // ======================================================================  
      /// get he matrix
      const TMatrixDSym& A            () const { return m_A       ; }
      /// size of object: number of componnents
      std::size_t        size         () const ;
      // ======================================================================      
    protected :
       // ======================================================================
      /// components 
      std::unique_ptr<RooArgList> m_cmps   {} ; // components 
      /// the matrix A
      TMatrixDSym                 m_A      {} ; // the matix A 
      // ======================================================================    
    } ; //                              The end of the class Ostap::Utils::COWs
    // ========================================================================
  } //                                        The END of namespace Ostap::Utils 
  // ==========================================================================
} //                                                 The END of namespace Ostap
// ============================================================================
#endif // OSTAP_COWS_H 
// ============================================================================
//                                                                      The END 
// ============================================================================
