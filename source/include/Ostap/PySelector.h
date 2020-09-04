// ============================================================================
#ifndef OSTAP_PYSELECTOR_H 
#define OSTAP_PYSELECTOR_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
// ============================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,22,0) 
// ============================================================================
#include "TPySelector.h"
// ============================================================================
#else 
// ============================================================================
#include "TPyArg.h"
#include "TSelector.h"
// ============================================================================
#endif 
// ============================================================================
// Forward declaratios
// ============================================================================
class TTree  ;            // ROOT 
class TChain ;            // ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
  typedef TPySelector ROOT_Selector ;
#else 
  typedef   TSelector ROOT_Selector ;
#endif
  // ==========================================================================
  /** @class Selector PySelector.h Ostap/PySelector.h
   *  Helper class for implementation of "python TSelector".
   *  The fix has been kindly provided by Wim Lavrijsen
   *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
   *  @date   2011-01-21
   */
  // ==========================================================================
  class Selector : public ROOT_Selector
  {
    // ========================================================================
  public:
    // ========================================================================
    ClassDef(Ostap::Selector, 2) ;
    // ========================================================================
  public:
    // ========================================================================
    /// constructor 
    // ========================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
    // ========================================================================
    Selector ( PyObject* self           , 
               TTree*    tree = nullptr ) ;
    // ========================================================================
#else 
    // ========================================================================
    Selector ( TTree*    tree = nullptr ) ;
    // ========================================================================
#endif 
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
    /// process 
    Bool_t Process        ( Long64_t entry      ) override ;
    /// notify 
    Bool_t Notify         ()                      override ;
    /// teminnate the slave 
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
    /// get the tree 
    TTree* tree() const ;
    ///  set the tree 
    void   set_tree  ( TTree* tree ) ;
    // ========================================================================
  private: 
    // ========================================================================
    /// the tree 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
    TTree* m_tree { nullptr } ; // the tree 
#endif
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
   *
   *  It is due to "python-unfriendly" signature of TTree::Process method  
   *
   *  @author Vanya BELYAEV Ivan.Belyaev Ivan.Belyaev@cern.ch
   *  @date 2010-11-21
   */
  class Process 
  {
  public:
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
    static
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
    static
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
    static 
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
    static 
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
