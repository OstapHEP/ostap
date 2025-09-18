// ============================================================================
#ifndef OSTAP_PYSELECTOR_H 
#define OSTAP_PYSELECTOR_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "TSelector.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ProgressConf.h" 
#include "Ostap/ProgressBar.h" 
// ============================================================================
// Forward declaratios
// ============================================================================
class TTree  ;            // ROOT 
class TChain ;            // ROOT
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class Selector PySelector.h Ostap/PySelector.h
   *  Helper class for implementation of "python TSelector".
   *  The fix has been kindly provided by Wim Lavrijsen
   *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
   *  @date   2011-01-21
   */
  // ==========================================================================
  class Selector : public TSelector
  {
    // ========================================================================
  public:
    // ========================================================================
    ClassDefOverride(Ostap::Selector, 4) ;
    // ========================================================================
  public:
    // ========================================================================
    /// constructor 
    // ========================================================================
    /// constructor from the TTree & progress bar configuration 
    Selector ( TTree* tree , const Ostap::Utils::ProgressConf& progress ) ;
    /// constructor from the TTree 
    Selector ( TTree* tree ) ;
    /// constructor progress bar configuration 
    Selector (               const Ostap::Utils::ProgressConf& progress ) ;
    /// default constructor 
    Selector () ;
    // ========================================================================
    /// destructor
    virtual ~Selector() ;
    // ========================================================================
  public: // the basic innterface 
    // ========================================================================
    /// init 
    void   Init           ( TTree*   tree       ) override ;
    /// begin 
    void   Begin          ( TTree*   tree       ) override ;
    /// initialize the slave 
    void   SlaveBegin     ( TTree*   tree       ) override ;
    // =======================================================================
    /** the main method: process entry 
     *  Note:
     *  - internally  calls  <code>GetEntry</code>
     *  - increment the event counter 
     *  - advance the progress bar 
     *  @see Ostap::Selector::GetEntry 
     *  @see Ostap::Selector::GetEntry 
     *  @see Ostap::Selector::process_entry 
     */  
    Bool_t Process        ( Long64_t entry      ) override ;
    /// notify 
    Bool_t Notify         ()                      override ;
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
  public: // Ostap-specific 
    // ========================================================================
    /// process an entry 
    virtual bool process_entry () ; 
    // ========================================================================
  public:
    // ========================================================================
    /// get the tree 
    TTree* get_tree () const { return m_tree ; } 
    ///  set the tree 
    void   set_tree  ( TTree* tree ) ;
    // ========================================================================
    /** event counter 
     *  - useless for PROOF
     *  - useful for interactive python
     *  incremented in Ostap::Selector::Process 
     */
    inline unsigned long long event                () const { return   m_event ; }
    // ========================================================================
  public: 
    // ========================================================================
    /** reset the selector:
     *  - set event count to zero 
     *  - set the maximunm count for progress bar to new  max-count
     */
    void reset ( const unsigned long long maxevents = 0 ) ;
    // ========================================================================
  protected:
    // ========================================================================
    /// increment the event counter & advance the progress bar  
    inline unsigned long long increment_event ()
    {
      if ( m_progress.enabled() ) { ++m_progress ; }
      return ++m_event ;
    }
    // ========================================================================
  protected:
    // ========================================================================
    /// progress bar
    Ostap::Utils::ProgressBar m_progress {} ;
    // ========================================================================
  private:
    // ========================================================================
    /// number of processed events 
    unsigned long long  m_event { 0 } ; // number of processed events
    // ========================================================================
  private : 
    // ========================================================================
    /// the tree 
    TTree* m_tree { nullptr } ; // the tree 
    // ========================================================================
  } ;
  // ==========================================================================
} //                                              the end of namespace Analysis 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class Process Ostap/PySelector.h
   *  Helper class to fix strange "feature" of (Py)ROOT:
   *     - ROOT.TTree  does have method <c>TTree.Process</c>  with TSelector as argument 
   *     - ROOT.TChain has *NO*  method <c>TChain.Process</c> with TSelector as argument 
   *  This trick allows to access these methods indirectly.
   *  It is due to "python-unfriendly" signature of TTree::Process method  
   *  @author Vanya BELYAEV Ivan.Belyaev Ivan.Belyaev@cern.ch
   *  @date 2010-11-21
   */
  namespace Utils 
  {
    // ========================================================================
    /** helper function to use TTree::Process in python 
     * 
     *  @param tree      root-tree 
     *  @param selector  the selector 
     *  
     *  @see TTree 
     *  @see TTree::Process 
     *  @see TSelector 
     *
     *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
     *  @date   2011-01-21
     */
    long process
    ( TTree*             tree      ,
      TSelector*         selector  ) ;
    // ========================================================================
    /** helper function to use TTree::Process in python 
     * 
     *  @param tree      root-tree 
     *  @param selector  the selector 
     *  @param events    events to be processed 
     *  
     *  @see TTree 
     *  @see TTree::Process 
     *  @see TSelector 
     *
     *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
     *  @date   2013-02-10
     */
    long process
    ( TTree*              tree         ,
      TSelector*          selector     , 
      const unsigned long events       , 
      const unsigned long first    = 0 ) ;
    // ========================================================================
    /** helper function to use TChain::Process in python 
     * 
     *  @param chain     root-tree/chain
     *  @param selector  the selector 
     *  
     *  @see TTree 
     *  @see TTree::Process 
     *  @see TSelector 
     *
     *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
     *  @date   2011-01-21
     */
    long process 
    ( TChain*    chain     ,
      TSelector* selector  ) ;
    // ========================================================================
    /** helper function to use TChain::Process in python 
     * 
     *  @param chain     root-tree/chain
     *  @param selector  the selector 
     *  @param events    events to be processed 
     *  
     *  @see TTree 
     *  @see TTree::Process 
     *  @see TSelector 
     *
     *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
     *  @date   2013-02-10
     */
    long process 
    ( TChain*             chain        ,
      TSelector*          selector     ,
      const unsigned long events       , 
      const unsigned long first    = 0 ) ;
    // ========================================================================
  };  
  // ==========================================================================
} //                                                 the end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYSELECTOR_H
// ============================================================================
