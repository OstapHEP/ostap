// $Id$
// ============================================================================
#ifndef OSTAP_PYSELECTOR_H 
#define OSTAP_PYSELECTOR_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT
// ============================================================================
#include "TPySelector.h"
// ============================================================================
// Forward declaratios
// ============================================================================
class TTree  ;            // ROOT 
class TChain ;            // ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class Selector PySelector.h Analysis/PySelector.h
   *
   *  Helper class for implementation of "python TSelector"
   *
   * The fix has been kindly provided by Wim Lavrijsen
   *
   *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
   *  @date   2011-01-21
   * 
   *                    $Revision$
   *  Last modification $Date$
   *                 by $Author$
   */
  class Selector : public  TPySelector
  {
    // ========================================================================
  public:
    // ========================================================================
    ClassDef(Ostap::Selector, 1) ;
    // ========================================================================
  public:
    // ========================================================================
    /// constructor 
    Selector
    ( TTree*    tree = 0 , 
      PyObject* self = 0 ) ;
    // ========================================================================
    /// destructor
    virtual ~Selector() ;
    // ========================================================================
  } ;
  // ==========================================================================
} //                                              the end of namespace Analysis 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class Process
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
    ( TTree*              tree      ,
      TSelector*          selector  , 
      const unsigned long events    ) ;
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
    ( TChain*             chain     ,
      TSelector*          selector  ,
      const unsigned long events    ) ;
    // ========================================================================
  };  
  // ==========================================================================
} //                                                 the end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYSELECTOR_H
// ============================================================================
