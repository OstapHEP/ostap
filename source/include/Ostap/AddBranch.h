// ============================================================================
#ifndef OSTAP_ADDBRANCH_H 
#define OSTAP_ADDBRANCH_H 1
// ============================================================================
// Include files 
// ============================================================================
//  STD&STL
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
// ============================================================================
// Forward declarations 
// ============================================================================
class TTree   ; // from ROOT 
class TBranch ; // from ROOT 
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
    TBranch* add_branch ( TTree*                  tree ,  
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
    TBranch* add_branch ( TTree*             tree    ,  
                          const std::string& name    , 
                          const std::string& formula ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADDBRANCH_H
// ============================================================================
