// ============================================================================
#ifndef OSTAP_ADDBRANCH_H 
#define OSTAP_ADDBRANCH_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================ 
#include <string>
#include <map>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
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
    /** @typedef FUNCTREEEMAP 
     *  helper type to deal with map of functions 
     */
    typedef std::map<std::string,const Ostap::IFuncTree*>        FUNCTREEMAP  ;
    // ========================================================================
    /**  add new branch with name <code>name</code> to the tree
     *   the value of the branch is taken from  function <code>func</code>
     *   @param tree    input tree 
     *   @param name    the name for new branch 
     *   @param func    the function to be used to fill the branch 
     *   @return status code 
     *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *   @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                  tree ,  
      const std::string&      name , 
      const Ostap::IFuncTree& func ) ;
    // ========================================================================
    /** add new branch with name <code>name</code> to the tree
     *  the value of the branch is taken from  function <code>func</code>
     *  @param tree    input tree 
     *  @param name    the name for new branch 
     *  @param formula the fomula use to calculate new  branch
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*             tree    ,  
      const std::string& name    , 
      const std::string& formula ) ;
    // ========================================================================
    /** add new branches to the tree
     *  the value of the branch each  is taken from <code>branches</code>
     *  @param tree     input tree 
     *  @param name     the name for new branch 
     *  @param branches the map name->formula use to calculate newbranch
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                                   tree     ,  
      const std::map<std::string,std::string>& branches ) ;
    // ========================================================================
    /** add new branches to the tree
     *  the value of the branch each  is taken from <code>branches</code>
     *  @param tree     input tree 
     *  @param name     the name for new branch 
     *  @param branches the map name->function use to calculate new branch
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*             tree     ,  
      const FUNCTREEMAP& branches ) ;
    // ========================================================================
    /** add new branch to TTree, sampling it from   the 1D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param name   name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return status code 
     *  @see TH1::GetRandom 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*               tree  , 
      const std::string&   name  , 
      const TH1&           histo ) ;
    // =========================================================================
    /** add new branch to TTree, sampling it from   the 1D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param namex  name of the new branch 
     *  @param namey  name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return status code 
     *  @see TH2::GetRandom2 
     */
    Ostap::StatusCode 
    add_branch 
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
     *  @return status code 
     *  @see TH3::GetRandom3 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*               tree  , 
      const std::string&   namex , 
      const std::string&   namey , 
      const std::string&   namez , 
      const TH3&           histo ) ;
    // ========================================================================
#if ROOT_VERSION(6,24,0)<=ROOT_VERSION_CODE
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data fuffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree        , 
      const std::string&   vname       ,  
      const double*        data        , 
      const unsigned long  size        , 
      const double         value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data fuffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const std::string&   vname     ,  
      const float*         data      , 
      const unsigned long  size      , 
      const float          value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data fuffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const std::string&   vname     ,  
      const int*           data      , 
      const unsigned long  size      , 
      const int            value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data fuffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const std::string&   vname     ,  
      const long*          data      , 
      const unsigned long  size      , 
      const long           value = 0 ) ;
    // ========================================================================
#endif
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree    The tree 
     *  @param namex   name of the new branch 
     *  @param value   the value 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree        , 
      const std::string&   vname       ,  
      const double         value       ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree    The tree 
     *  @param namex   name of the new branch 
     *  @param value   the value 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree        , 
      const std::string&   vname       ,  
      const int            value       ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADDBRANCH_H
// ============================================================================
