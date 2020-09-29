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
// ROOT 
// ============================================================================
#include "RVersion.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PySelector.h"
// ============================================================================
// forward declarations 
// ============================================================================
class TCut  ; // ROOT 
class TTree ; // ROOT 
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
  class SelectorWithCuts : public Ostap::Selector 
  {
    // ========================================================================
  public:
    // ========================================================================
    ClassDef (Ostap::SelectorWithCuts , 2) ;
    // ========================================================================
  public:
    // ========================================================================
    /// constructor 
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
    // ========================================================================
    SelectorWithCuts
    ( PyObject*          self           , 
      const std::string& cuts = ""      , 
      TTree*             tree = nullptr ) ;
    SelectorWithCuts
    ( PyObject*          self           , 
      const TCut&        cuts           ,      
      TTree*             tree = nullptr ) ;
    // ========================================================================
#else 
    // ========================================================================
    SelectorWithCuts 
    ( const std::string& cuts = ""      , 
      TTree*             tree = nullptr ) ;
    SelectorWithCuts 
    ( const TCut&        cuts           ,
      TTree*             tree = nullptr ) ;
    // ========================================================================
#endif 
    // ========================================================================
    /// virtual destructor 
    virtual ~SelectorWithCuts () ; // virtual destructor 
    // ========================================================================
  public:
    // ========================================================================
    Bool_t Notify       () override ;
    // ========================================================================
    void   Init         ( TTree*   tree  ) override ;
    void   Begin        ( TTree*   tree  ) override ;
    void   SlaveBegin   ( TTree*   tree  ) override ;
    Bool_t Process      ( Long64_t entry ) override ;
    // ========================================================================
  public:
    // ========================================================================    
    /// make formula cuts 
    bool make_formula  ( TTree*   tree  ) ;
    /// check the entry 
    bool good_entry    ( Long64_t entry ) ;
    // ========================================================================
    /// reset formula unisy new tree 
    void reset_formula ( TTree*   tree  ) ;
    // ========================================================================
    /// process a good entry
    bool process_entry () override ;
    // ========================================================================
  public:
    // ========================================================================    
    /// set new cuts 
    void set_cuts ( const std::string& cuts ) ;
    /// set new cuts 
    void set_cuts ( const TCut&        cuts ) ;
    // ========================================================================
  public:
    // ========================================================================    
    /// is formula OK ? 
    bool ok () const  ; // is formula OK ? 
    /// get the formula 
    Ostap::Formula*    formula () const { return fMyFormula.get() ; }
    /// get the formula
    const std::string& cuts    () const { return fMyCuts ; }
    /** good event counter 
     *  - useless for PROOF
     *  - useful for interactive python
     *  incremented in Ostap::SelectorWithCuts::Process 
     */
    unsigned long long good    () const { return m_good  ; }
    // ========================================================================
  private:
    // ========================================================================
    /// the selection formula 
    std::string                        fMyCuts    ; 
    std::unique_ptr<Ostap::Formula>    fMyFormula ;
    /// event counter 
    unsigned long long  m_good  { 0 } ; // event counter: useless for PROOF
    // ========================================================================    
  };
  // ==========================================================================
} //                                                  end of namespace Analysis 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // ANALYSIS_PYSELECTORWITHCUTS_H
// ============================================================================
