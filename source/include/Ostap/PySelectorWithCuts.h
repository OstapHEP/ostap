// ============================================================================
#ifndef ANALYSIS_PYSELECTORWITHCUTS_H 
#define ANALYSIS_PYSELECTORWITHCUTS_H 1
// ============================================================================
// Include files
// ============================================================================
//  STD & STL 
// ============================================================================
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PySelector.h"
// ============================================================================
// forward decalrations 
class TCut ; // ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  class Formula ;
  // ==========================================================================
  /** @class PySelectorWithCuts Ostap/PySelectorWithCuts.h
   *  @author Vanya Belyaev
   *  @date   2013-05-06
   */
  class SelectorWithCuts : public Selector 
  {
    // ========================================================================
  public:
    // ========================================================================
    ClassDef (Ostap::SelectorWithCuts , 1 ) ;
    // ========================================================================
  public:
    // ========================================================================
    /// constructor 
    SelectorWithCuts
    ( const std::string& cuts = "" , 
      TTree*             tree = 0  , 
      PyObject*          self = 0  ) ;
    SelectorWithCuts
    ( const TCut&        cuts      , 
      TTree*             tree = 0  , 
      PyObject*          self = 0  ) ;
    /// virtual destructor 
    virtual ~SelectorWithCuts () ; // virtual destructor 
    // ========================================================================
  public:
    // ========================================================================
    Bool_t         Notify       () override ;
    // ========================================================================
    virtual void   Init         ( TTree*   tree  ) ;
    virtual void   Begin        ( TTree*   tree  ) ;
    virtual void   SlaveBegin   ( TTree*   tree  ) ;
    virtual Bool_t Process      ( Long64_t entry ) ;
    // ========================================================================
  public:
    // ========================================================================    
    /// is formula OK ? 
    bool ok () const  ; // is formula OK ? 
    /// get the formula 
    Ostap::Formula*    formula () const { return fMyformula.get() ; }
    /// get the formula
    const std::string& cuts    () const { return fMycuts    ; }
    /// event counter (useless for PROOF, useful for interactive python) 
    unsigned long long event   () const { return m_event    ; }
    // ========================================================================
  private:
    // ========================================================================
    /// the selection formula 
    std::string                        fMycuts    ; 
    std::unique_ptr<Ostap::Formula>    fMyformula ;
    /// event counter 
    unsigned long long  m_event    ; // event counter: useless for PROOF
    // ========================================================================    
  };
  // ==========================================================================
} //                                                  end of namespace Analysis 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // ANALYSIS_PYSELECTORWITHCUTS_H
// ============================================================================
