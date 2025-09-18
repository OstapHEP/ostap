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
    ClassDefOverride (Ostap::SelectorWithCuts , 3) ;
    // ========================================================================
  public:
    // ========================================================================
    /// constructor 
    // ========================================================================
    /** full contructor 
     *  @param cuts cuts 
     *  @param tree input TTree/TChain
     *  @param progress configuration of progress bar
     */
    SelectorWithCuts 
    ( const std::string&                cuts     , 
      TTree*                            tree     , 
      const Ostap::Utils::ProgressConf& progress ) ;
    // ========================================================================
    /** full contructor 
     *  @param cuts cuts 
     *  @param tree input TTree/TChain
     *  @param progress configuration of progress bar
     */
    SelectorWithCuts 
    ( const TCut&                       cuts     , 
      TTree*                            tree     , 
      const Ostap::Utils::ProgressConf& progress ) ;
    // ========================================================================      
    SelectorWithCuts 
    ( const std::string&                cuts = ""      , 
      TTree*                            tree = nullptr ) ;
    // ========================================================================      
    SelectorWithCuts 
    ( const TCut&                       cuts           ,
      TTree*                            tree = nullptr ) ;
    // ========================================================================      
    SelectorWithCuts 
    ( const std::string&                cuts           , 
      const Ostap::Utils::ProgressConf& progress       , 
      TTree*                            tree = nullptr ) ;    
    // ========================================================================
    SelectorWithCuts 
    ( const TCut&                       cuts           ,
      const Ostap::Utils::ProgressConf& progress       , 
      TTree*                            tree = nullptr ) ;
    // ========================================================================      
    /// virtual destructor 
    virtual ~SelectorWithCuts () ; // virtual destructor 
    // ========================================================================
  public:
    // ========================================================================
    Bool_t Notify       () override ;
    // ========================================================================\
    /// Init 
    void   Init           ( TTree*   tree  )      override ;
    /// begin 
    void   Begin          ( TTree*   tree  )      override ;
    /// initialize the slave 
    void   SlaveBegin     ( TTree*   tree  )      override ;
    /// process 
    Bool_t Process        ( Long64_t entry )      override ;
    /// terminate the slave 
    void   SlaveTerminate ()                      override ;
    /// terminate
    void   Terminate      ()                      override ;
    /// get entry 
    Int_t  GetEntry       ( Long64_t entry      , 
                            Int_t    getall = 0 ) override ;
    /// Version 
    Int_t  Version        () const                override ;
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
    Ostap::Formula*    formula () const { return m_formula.get() ; }
    /// get the formula
    const std::string& cuts    () const { return m_cuts ; }
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
    std::string                     m_cuts    {} ; 
    std::unique_ptr<Ostap::Formula> m_formula {} ;
    /// event counter 
    unsigned long long              m_good    { 0 } ;
    // ========================================================================    
  };
  // ==========================================================================
} //                                                  end of namespace Analysis 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // ANALYSIS_PYSELECTORWITHCUTS_H
// ============================================================================
