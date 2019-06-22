// ============================================================================
#ifndef OSTAP_ADDBRANCH_H 
#define OSTAP_ADDBRANCH_H 1
// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
// ============================================================================
// Forward declarations 
// ============================================================================
class TTree   ; // from ROOT 
class TBranch ; // from ROOT 
class TH1     ;
class TH2     ;
class TH3     ;
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Trees
  {
    // ========================================================================
    /**  add new branch with name <code>name</code> to the tree
     *   the value of the branch is taken from  function <code>func</code>
     *   @param tree    input tree 
     *   @param name    the name for new branch 
     *   @param func    the function to be used to fill the branch 
     *   @return new  branch 
     *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *   @date 2019-05-14
     */
    TBranch* add_branch 
    ( TTree*                  tree ,  
      const std::string&      name , 
      const Ostap::IFuncTree& func ) ;
    // ========================================================================
    /**  add new branch with name <code>name</code> to the tree
     *   the value of the branch is taken from  function <code>func</code>
     *   @param tree    input tree 
     *   @param name    the name for new branch 
     *   @param formula the fomula use to calculate new  branch
     *   @return new  branch 
     *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *   @date 2019-05-14
     */
    TBranch* add_branch 
    ( TTree*             tree    ,  
      const std::string& name    , 
      const std::string& formula ) ;
    // ========================================================================
    /** add new branch to TTree, sampling it from   the 1D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param name   name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return new  branch 
     *  @see TH1::GetRandom 
     */
    TBranch* add_branch 
    ( TTree*               tree  , 
      const std::string&   name  , 
      const TH1&           histo ) ;
    // =========================================================================
    /** add new branch to TTree, sampling it from   the 1D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param namex  name of the new branch 
     *  @param namey  name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return new  brances 
     *  @see TH2::GetRandom2 
     */
    TBranch* add_branch 
    ( TTree*               tree  , 
      const std::string&   namex , 
      const std::string&   namey , 
      const TH2&           histo ) ;
    // ========================================================================
    /** add new branch to TTree, sampling it from   the 1D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param namex  name of the new branch 
     *  @param namey  name of the new branch 
     *  @param namez  name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return new  brances 
     *  @see TH3::GetRandom3 
     */
    TBranch* add_branch 
    ( TTree*               tree  , 
      const std::string&   namex , 
      const std::string&   namey , 
      const std::string&   namez , 
      const TH3&           histo ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADDBRANCH_H
// ============================================================================
