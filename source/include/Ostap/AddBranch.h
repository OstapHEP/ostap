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
#include <functional>
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
  namespace Utils 
  {
    // ========================================================================
    /// progress bar configrutaion 
    class ProgressConf ; // progress bar configuration 
    // ========================================================================
  }
  // ==========================================================================
  namespace Trees
  { 
    // ========================================================================
    /** @typedef FUNCTREEEMAP 
     *  helper type to deal with map of functions 
     */
    typedef std::map<std::string,const Ostap::IFuncTree*>        FUNCTREEMAP  ;
    // =======================================================================
    /** @typedef FUNCTREEEMAP 
     *  helper type to deal with map of functions 
     */
    typedef FUNCTREEMAP::value_type                              FUNCTREEPAIR ;
    // ========================================================================
    class Branches 
    {
      // ===================================================================== 
    public:
      // ======================================================================
      bool add
      ( const std::string&      name ,
	const Ostap::IFuncTree& func ) ;      
      bool add
      ( const Ostap::IFuncTree& func ,		 
	const std::string&      name ) ;
      // ======================================================================
    public:
      // ======================================================================
      FUNCTREEMAP::const_iterator begin () const { return m_map.begin () ; }
      FUNCTREEMAP::const_iterator end   () const { return m_map.end   () ; }
      std::size_t                 size  () const { return m_map.size  () ; }
      bool                        empty () const { return m_map.empty () ; }
      const FUNCTREEMAP&          map   () const { return m_map          ; }
      // ======================================================================
    private :
      // ======================================================================
      /// the actual storage 
      FUNCTREEMAP m_map {} ; // the actual storage 
      // ======================================================================
    };
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
    /**  add new branch with name <code>name</code> to the tree
     *   the value of the branch is taken from  function <code>func</code>
     *   @param tree     input tree 
     *   @param progress configuration of the progress bar
     *   @param name     the name for new branch 
     *   @param func     the function to be used to fill the branch 
     *   @return status code 
     *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *   @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree     ,  
      const Ostap::Utils::ProgressConf& progress , 
      const std::string&                name     , 
      const Ostap::IFuncTree&           func     ) ;
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
    /** add new branch with name <code>name</code> to the tree
     *  the value of the branch is taken from  function <code>func</code>
     *  @param tree    input tree 
     *  @param progress configuration of the progress bar
     *  @param name    the name for new branch 
     *  @param formula the fomula use to calculate new  branch
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree     ,  
      const Ostap::Utils::ProgressConf& progress , 
      const std::string&                name     , 
      const std::string&                formula  ) ;
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
     *  @param progress configuration of the progress bar
     *  @param name     the name for new branch 
     *  @param branches the map name->formula use to calculate newbranch
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                                   tree     ,  
      const Ostap::Utils::ProgressConf&        progress , 
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
    /** add new branches to the tree
     *  the value of the branch each  is taken from <code>branches</code>
     *  @param tree     input tree 
     *  @param progress configuration of the progress bar
     *  @param name     the name for new branch 
     *  @param branches the map name->function use to calculate new branch
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree     ,  
      const Ostap::Utils::ProgressConf& progress , 
      const FUNCTREEMAP&                branches ) ;
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
      const Branches&    branches ) ;
    // ========================================================================
    /** add new branches to the tree
     *  the value of the branch each  is taken from <code>branches</code>
     *  @param tree     input tree 
     *  @param progress configuration of the progress bar
     *  @param name     the name for new branch 
     *  @param branches the map name->function use to calculate new branch
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree     ,  
      const Ostap::Utils::ProgressConf& progress , 
      const Branches&                   branches ) ;
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
    // ========================================================================
    /** add new branch to TTree, sampling it from   the 1D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param progress configuration of the progress bar
     *  @param name   name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return status code 
     *  @see TH1::GetRandom 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree     , 
      const Ostap::Utils::ProgressConf& progress , 
      const std::string&                name     , 
      const TH1&                        histo    ) ;
    // =========================================================================
    /** add new branch to TTree, sampling it from   the 2D-histogram
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
    /** add new branch to TTree, sampling it from   the 2D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param progress configuration of the progress bar
     *  @param namex  name of the new branch 
     *  @param namey  name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return status code 
     *  @see TH2::GetRandom2 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree     , 
      const Ostap::Utils::ProgressConf& progress , 
      const std::string&                namex    , 
      const std::string&                namey    , 
      const TH2&                        histo    ) ;
    // ========================================================================
    /** add new branch to TTree, sampling it from   the 3D-histogram
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
    /** add new branch to TTree, sampling it from   the 3D-histogram
     *  @param tree (UPFATE) input tree 
     *  @param progress configuration of the progress bar
     *  @param namex  name of the new branch 
     *  @param namey  name of the new branch 
     *  @param namez  name of the new branch 
     *  @param histo  the historgam to be  sampled
     *  @return status code 
     *  @see TH3::GetRandom3 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree     , 
      const Ostap::Utils::ProgressConf& progress , 
      const std::string&                namex    , 
      const std::string&                namey    , 
      const std::string&                namez    , 
      const TH3&                        histo    ) ;
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
     *  @param progress configuration of the progress bar
     *  @param data   input data fuffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*                            tree        , 
      const Ostap::Utils::ProgressConf& progress    , 
      const std::string&                vname       ,  
      const double*                     data        , 
      const unsigned long               size        , 
      const double                      value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
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
     *  @param progress configuration of the progress bar
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*                            tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&                vname     ,  
      const float*                      data      , 
      const unsigned long               size      , 
      const float                       value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const std::string&   vname     ,  
      const char*          data      , 
      const unsigned long  size      , 
      const char           value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname     ,  
      const char*          data      , 
      const unsigned long  size      , 
      const char           value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*                tree      , 
      const std::string&    vname     ,  
      const unsigned char*  data      , 
      const unsigned long   size      , 
      const unsigned char   value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*                tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&    vname     ,  
      const unsigned char*  data      , 
      const unsigned long   size      , 
      const unsigned char   value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const std::string&   vname     ,  
      const short*         data      , 
      const unsigned long  size      , 
      const short          value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname     ,  
      const short*         data      , 
      const unsigned long  size      , 
      const short          value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*                tree      , 
      const std::string&    vname     ,  
      const unsigned short* data      , 
      const unsigned long   size      , 
      const unsigned short  value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*                tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&    vname     ,  
      const unsigned short* data      , 
      const unsigned long   size      , 
      const unsigned short  value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
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
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname     ,  
      const int*           data      , 
      const unsigned long  size      , 
      const int            value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const std::string&   vname     ,  
      const unsigned int*  data      , 
      const unsigned long  size      , 
      const unsigned int   value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname     ,  
      const unsigned int*  data      , 
      const unsigned long  size      , 
      const unsigned int   value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
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
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname     ,  
      const long*          data      , 
      const unsigned long  size      , 
      const long           value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const std::string&   vname     ,  
      const unsigned long* data      , 
      const unsigned long  size      , 
      const unsigned long  value = 0 ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree   The tree 
     *  @param data   input data buffer 
     *  @param size   length of the buffer
     *  @param value  default value (used for short buffers) 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname     ,  
      const unsigned long* data      , 
      const unsigned long  size      , 
      const unsigned long  value = 0 ) ;
    // ========================================================================
#endif
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree    The tree 
     *  @param vname   name of the new branch 
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
     *  @param vname   name of the new branch 
     *  @param value   the value 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree        , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname       ,  
      const double         value       ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree    The tree 
     *  @param vname   name of the new branch 
     *  @param value   the value 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree        , 
      const std::string&   vname       ,  
      const int            value       ) ;
    // ========================================================================
    /** copy data from buffer into new branch 
     *  @param tree    The tree 
     *  @param vname   name of the new branch 
     *  @param value   the value 
     *  @return status code 
     */
    Ostap::StatusCode
    add_branch 
    ( TTree*               tree        , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&   vname       ,  
      const int            value       ) ;
    // ========================================================================
    // Add branch from generic 1D function 
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param xname (INPUT) name of input variable 
     *  @param fun   (INPUT) the function 
     *  @return status code 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                        tree  , 
      const std::string&            bname , 
      const std::string&            xname , 
      std::function<double(double)> fun   ) ;
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree     (INPUT) The tree 
     *  @param progress (INPUT) configuration of the progress bar
     *  @param bname    (INPUT) branch name 
     *  @param xname    (INPUT) name of input variable 
     *  @param fun      (INPUT) the function 
     *  @return status code 
     */
    Ostap::StatusCode 
    add_branch
    ( TTree*                            tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&                bname     , 
      const std::string&                xname     , 
      std::function<double(double)>     fun       ) ;
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param fun   (INPUT) the function 
     *  @param xname (INPUT) name of input variable 
     *  @return status code 
     */
    template <class FUNCTION>
    inline Ostap::StatusCode 
    add_branch_ 
    ( TTree*                        tree  , 
      const std::string&            bname , 
      FUNCTION                      fun   ,
      const std::string&            xname ) 
    { return add_branch ( tree , bname , xname , std::cref ( fun ) ) ; }        
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param fun   (INPUT) the function 
     *  @param xname (INPUT) name of input variable 
     *  @return status code 
     */
    template <class FUNCTION>
    inline Ostap::StatusCode 
    add_branch_ 
    ( TTree*                            tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&                bname     , 
      FUNCTION                          fun       ,
      const std::string&                xname     ) 
    { return add_branch ( tree , progress , bname , xname , std::cref ( fun ) ) ; }        
    // ========================================================================
    // Add branch from generic 2D function 
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param xname (INPUT) name of input variable 
     *  @param yname (INPUT) name of input variable 
     *  @param fun   (INPUT) the function 
     *  @return status code 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                               tree  , 
      const std::string&                   bname , 
      const std::string&                   xname , 
      const std::string&                   yname , 
      std::function<double(double,double)> fun   ) ;
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree     (INPUT) The tree 
     *  @param progress (INPUT) configuration of the progress bar
     *  @param bname    (INPUT) branch name 
     *  @param xname    (INPUT) name of input variable 
     *  @param yname    (INPUT) name of input variable 
     *  @param fun      (INPUT) the function 
     *  @return status code 
     */
    Ostap::StatusCode 
    add_branch
    ( TTree*                               tree      , 
      const Ostap::Utils::ProgressConf&    progress  , 
      const std::string&                   bname     , 
      const std::string&                   xname     , 
      const std::string&                   yname     , 
      std::function<double(double,double)> fun       ) ;
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param fun   (INPUT) the function 
     *  @param xname (INPUT) name of input variable 
     *  @param yname (INPUT) name of input variable 
     *  @return status code 
     */
    template <class FUNCTION>
    inline Ostap::StatusCode 
    add_branch_ 
    ( TTree*                        tree  , 
      const std::string&            bname , 
      FUNCTION                      fun   ,
      const std::string&            xname , 
      const std::string&            yname ) 
    { return add_branch ( tree , bname , xname , yname , std::cref ( fun ) ) ; }        
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param fun   (INPUT) the function 
     *  @param xname (INPUT) name of input variable 
     *  @param yname (INPUT) name of input variable 
     *  @return status code 
     */
    template <class FUNCTION>
    inline Ostap::StatusCode 
    add_branch_ 
    ( TTree*                            tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&                bname     , 
      FUNCTION                          fun       ,
      const std::string&                xname     , 
      const std::string&                yname     ) 
    { return add_branch ( tree , progress , bname , xname , yname , std::cref ( fun ) ) ; }
    // ========================================================================
    // Add branch from generic 3D function 
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param xname (INPUT) name of input variable 
     *  @param yname (INPUT) name of input variable 
     *  @param zname (INPUT) name of input variable 
     *  @param fun   (INPUT) the function 
     *  @return status code 
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                                      tree  , 
      const std::string&                          bname , 
      const std::string&                          xname , 
      const std::string&                          yname , 
      const std::string&                          zname , 
      std::function<double(double,double,double)> fun   ) ;
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree     (INPUT) The tree 
     *  @param progress (INPUT) configuration of the progress bar
     *  @param bname    (INPUT) branch name 
     *  @param xname    (INPUT) name of input variable 
     *  @param yname    (INPUT) name of input variable 
     *  @param zname    (INPUT) name of input variable 
     *  @param fun      (INPUT) the function 
     *  @return status code 
     */
    Ostap::StatusCode 
    add_branch
    ( TTree*                                      tree      , 
      const Ostap::Utils::ProgressConf&           progress  , 
      const std::string&                          bname     , 
      const std::string&                          xname     , 
      const std::string&                          yname     , 
      const std::string&                          zname     , 
      std::function<double(double,double,double)> fun       ) ;
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param fun   (INPUT) the function 
     *  @param xname (INPUT) name of input variable 
     *  @param yname (INPUT) name of input variable 
     *  @param zname (INPUT) name of input variable 
     *  @return status code 
     */
    template <class FUNCTION>
    inline Ostap::StatusCode 
    add_branch_ 
    ( TTree*                        tree  , 
      const std::string&            bname , 
      FUNCTION                      fun   ,
      const std::string&            xname , 
      const std::string&            yname , 
      const std::string&            zname ) 
    { return add_branch ( tree , bname , xname , yname , zname , std::cref ( fun ) ) ; }        
    // ========================================================================
    /** add new branch to the tree, calculated from the function  
     *  @param tree  (INPUT) The tree 
     *  @param bname (INPUT) branch name 
     *  @param fun   (INPUT) the function 
     *  @param xname (INPUT) name of input variable 
     *  @param yname (INPUT) name of input variable 
     *  @param zname (INPUT) name of input variable 
     *  @return status code 
     */
    template <class FUNCTION>
    inline Ostap::StatusCode 
    add_branch_ 
    ( TTree*                            tree      , 
      const Ostap::Utils::ProgressConf& progress  , 
      const std::string&                bname     , 
      FUNCTION                          fun       ,
      const std::string&                xname     , 
      const std::string&                yname     ,
      const std::string&                zname     ) 
    { return add_branch ( tree , progress , bname , xname , yname , zname , std::cref ( fun ) ) ; }
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADDBRANCH_H
// ============================================================================
