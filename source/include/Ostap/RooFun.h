// ============================================================================
#ifndef OSTAP_ROOFUN_H 
#define OSTAP_ROOFUN_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
// ============================================================================
// ROOT/RooFit 
// ============================================================================
// #include "RooAbsReal.h"
// #include "RooArgSet.h"
// ============================================================================
// forward declaratios 
// ============================================================================
// ROOT/RooFit 
// ============================================================================
class RooAbsData       ; // ROOT/Roofit 
class RooAbsReal       ; // ROOT/Roofit 
class RooAbsCollection ; // ROOT/Roofit 
class RooArgSet        ; // ROOT/Roofit 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** @class RooFun
     *  Helper class to keep RooAbsPdf&observables&normalzation 
     */
    class RooFun 
    {
      // ======================================================================
    public :
      // ======================================================================
      /** @param fun the function 
       *  @param observabels list of observables 
       *  @param normalzation normalisation set 
       */
      RooFun 
      ( const RooAbsReal&       fun                     ,
        const RooAbsCollection& observables             ,
        const RooAbsCollection* normalization = nullptr ) ;
      // ======================================================================
      /** @param fun the function 
       *  @param observabels list of observables 
       *  @param normalzation normalisation set 
       */
      RooFun 
      ( const RooAbsReal&       fun                     ,
        const RooAbsData&       data                    ,
        const RooAbsCollection* normalization = nullptr ) ;
      // ======================================================================
      RooFun ( const RooFun&  right ) ;
      RooFun (       RooFun&& right ) = default ;
      RooFun (                      ) = default ;
      // ======================================================================
      virtual ~RooFun() ;  
      /// clone funcntion/virtual constructire 
      virtual RooFun* clone() const ;  
      // ======================================================================      
    public: // constant getters 
      // ======================================================================
      /// get the pdf 
      const  RooAbsReal& fun           () const { return *m_fun            ; }
      /// get normalzation 
      const  RooArgSet*  normalization () const { return  m_normset.get () ; }
      // ======================================================================      
    public: // non-constant getters 
      // ======================================================================
      /// get list of observables
      inline RooArgSet&  observables   () const { return *m_observables    ; }
      /// get list of parameters
      inline RooArgSet&  parameters    () const { return *m_parameters     ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// assign observables 
      void set_observables ( const RooAbsCollection& obs  ) const ;
      /// assign parameters 
      void set_parameters  ( const RooAbsCollection& pars ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate it! 
      double evaluate () const ; 
      // ======================================================================      
    private :
      // ======================================================================
      /// perform initialization 	  
      void Init
      ( const RooAbsCollection& observables   ,
        const RooAbsCollection* normalization ) ;
      // ======================================================================      
    protected : 
      // ======================================================================
      std::unique_ptr<RooAbsReal> m_fun         {} ; // function 
      /// observables
      std::unique_ptr<RooArgSet>  m_observables {} ; // observables 
      /// parameters
      std::unique_ptr<RooArgSet>  m_parameters  {} ; // parameters 
      /// normalization
      std::unique_ptr<RooArgSet>  m_normset     {} ; // normalization
      // ======================================================================
    };
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namesapce Ostap
// ============================================================================
#endif // OSTAP_ROOFUN_H
// ============================================================================
//                                                                      The END
// ============================================================================
