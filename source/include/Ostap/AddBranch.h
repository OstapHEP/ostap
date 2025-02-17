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
#include <list>
#include <functional>
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/StatusCode.h"
#include "Ostap/ProgressConf.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
// ============================================================================
// Forward declarations from ROOT 
// ============================================================================
class TTree      ; 
class TBranch    ; 
class TH1        ;
class TH2        ;
class TH3        ;
// ============================================================================
// Forward declarations from ROOT/RooFit 
// ============================================================================
class RooAbsReal       ;
class RooAbsCollection ;
// ============================================================================
// Forward declarations from Ostap 
// ============================================================================
namespace  Ostap { namespace MoreRooFit { class SPlot4Tree ; } }  // Ostap 
namespace  Ostap { namespace MoreRooFit { class RooFun     ; } }  // Ostap 
namespace  Ostap { namespace Math       { class Histo1D    ; } }  // Ostap 
namespace  Ostap { namespace Math       { class Histo2D    ; } }  // Ostap 
namespace  Ostap { namespace Math       { class Histo3D    ; } }  // Ostap 
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
    /** valid name for branch ?
     *  - not empty 
     *  - no blanks 
     *  - no special symbols 
     */
    bool valid_name_for_branch ( const std::string& name ) ;
    // ========================================================================
    /** @typdef DCT
     *  string->string dictionary 
     */
    typedef  std::map<std::string,std::string>                    DCT ;
    // ========================================================================
    /** @class Branches 
     *  helper map-like class to keep informatio abotu the new branches 
     */
    class Branches 
    {
      // =====================================================================
    public :
      // =====================================================================
      Branches () ; 
      /// copy constructor 
      Branches ( const Branches&  right ) ;
      /// move constructor 
      Branches (       Branches&& right ) ;
      /// desctructor 
      ~Branches () ; // destructor 
      // ===================================================================== 
    public:
      // ===================================================================== 
      ///  add a fnuction to calculate new branch 
      bool add
      ( const std::string&      name                  ,
        const Ostap::IFuncTree& func                  ,
        const TTree*            tree        = nullptr ) ;
      // ======================================================================
      ///  add an expression to calculate new branch 
      bool add
      ( const std::string&      name                  ,
        const std::string&      expression  = ""      ,
        const TTree*            tree        = nullptr ) ;
      // ======================================================================
    public:
      // ======================================================================
      bool                      empty  () const { return m_map.empty () ; }
      std::size_t               size   () const { return m_map.size  () ; }
      // ======================================================================     
    public: 
      // ======================================================================
      typedef const Ostap::IFuncTree*                  BRANCH          ;
      typedef std::vector<std::string>                 NAMES           ; 
      typedef NAMES::const_iterator                    const_iterator; ; 
      typedef const_iterator                           iterator        ; 
      // ======================================================================
      typedef std::map<std::string,BRANCH>             BRANCHES        ;  
      // ======================================================================
    public:
      // =====================================================================      
      const_iterator begin  () const { return m_names.begin () ; }
      const_iterator end    () const { return m_names.end   () ; }
      // ======================================================================
      bool has_key ( const std::string& name ) const ; 
      // ======================================================================
      /// get new branch by name 
      const Ostap::IFuncTree*   branch ( const std::string& key ) const ;
      // ======================================================================
    private :
      // ======================================================================
      /// the actual storage
      BRANCHES m_map   {} ; // the actual storage
      /// branch  
      NAMES    m_names {} ;
      // ======================================================================
    }; //                               The end of class Ostap::Trees::Branches
    // ========================================================================
    /**  add new branch with name <code>name</code> to the tree
     *   the value of the branch is taken from  function <code>func</code>
     *   @param tree     input tree 
     *   @param name     the name for new branch 
     *   @param func     the function to be used to fill the branch 
     *   @param progress configuration of the progress bar
     *   @return status code 
     *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *   @date 2019-05-14
     */
    Ostap::StatusCode 
    add_branch 
    ( TTree*                            tree             ,  
      const std::string&                name             , 
      const Ostap::IFuncTree&           func             , 
      const Ostap::Utils::ProgressConf& progress = false ) ;
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
    ( TTree*                            tree             ,  
      const std::string&                name             , 
      const std::string&                formula          , 
      const Ostap::Utils::ProgressConf& progress = false ) ;
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
    ( TTree*                                   tree         ,  
      const std::map<std::string,std::string>& branches     , 
      const Ostap::Utils::ProgressConf&        progress = 0 ) ;
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
    ( TTree*                                   tree             ,  
      const Branches&                          branches         , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;
    // ========================================================================
    // sample variabales from 1D,2D&3D histograms 
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
    ( TTree*                                   tree             , 
      const std::string&                       name             , 
      const TH1&                               histo            , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;
    // =========================================================================
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
    ( TTree*                                   tree             , 
      const std::string&                       namex            , 
      const std::string&                       namey            , 
      const TH2&                               histo            , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;
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
    ( TTree*                                   tree             , 
      const std::string&                       namex            , 
      const std::string&                       namey            , 
      const std::string&                       namez            , 
      const TH3&                               histo            , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;
    // ========================================================================
    
    // ========================================================================
    // sample varibiales from (smoothed) 1D,2D&3D smoothed historgams 
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
    ( TTree*                                   tree             , 
      const std::string&                       name             , 
      const Ostap::Math::Histo1D&              histo            , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;
    // =========================================================================
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
    ( TTree*                                   tree             , 
      const std::string&                       namex            , 
      const std::string&                       namey            , 
      const Ostap::Math::Histo2D&              histo            , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;
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
    ( TTree*                                   tree             , 
      const std::string&                       namex            , 
      const std::string&                       namey            , 
      const std::string&                       namez            , 
      const Ostap::Math::Histo3D&              histo            , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;
    // ========================================================================
    
    // ========================================================================
    // Add branch from generic 1D function 
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
    ( TTree*                            tree             , 
      const std::string&                bname            , 
      const std::string&                xname            , 
      std::function<double(double)>     fun              , 
      const Ostap::Utils::ProgressConf& progress = false ) ;
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
    ( TTree*                            tree             , 
      const std::string&                bname            , 
      FUNCTION                          fun              ,
      const std::string&                xname            , 
      const Ostap::Utils::ProgressConf& progress = false )       
    { return add_branch ( tree  ,
                          bname ,
                          xname ,
                          std::function<double(double)> ( fun ) ,                          
                          progress ) ; }
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
    ( TTree*                               tree             , 
      const std::string&                   bname            , 
      const std::string&                   xname            , 
      const std::string&                   yname            , 
      std::function<double(double,double)> fun              ,
      const Ostap::Utils::ProgressConf&    progress = false ) ;       
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
    ( TTree*                            tree             , 
      const std::string&                bname            , 
      FUNCTION                          fun              ,
      const std::string&                xname            ,  
      const std::string&                yname            ,
      const Ostap::Utils::ProgressConf& progress = false )        
    { return add_branch ( tree  ,
                          bname ,
                          xname ,
                          yname ,
                          std::function<double(double,double)>  ( fun ) ,
                          progress ) ; }        
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
    ( TTree*                                      tree             , 
      const std::string&                          bname            , 
      const std::string&                          xname            , 
      const std::string&                          yname            , 
      const std::string&                          zname            , 
      std::function<double(double,double,double)> fun              , 
      const Ostap::Utils::ProgressConf&           progress = false ) ; 
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
    ( TTree*                            tree             , 
      const std::string&                bname            , 
      FUNCTION                          fun              ,
      const std::string&                xname            , 
      const std::string&                yname            , 
      const std::string&                zname            , 
      const Ostap::Utils::ProgressConf& progress = false )  
    { return add_branch ( tree  ,
                          bname ,
                          xname ,
                          yname ,
                          zname ,
                          std::function<double(double,double,double)>  ( fun ) ,
                          progress ) ; }        
    // ========================================================================
    /** add new branch to the tree from RooFit function
     *  @param tree input tree 
     *  @param bname branch name 
     *  @param fun   the function 
     *  @param mapping : observables -> brnaches 
     */
    Ostap::StatusCode
    add_branch
    ( TTree*                            tree              ,
      const std::string&                bname             , 
      const Ostap::MoreRooFit::RooFun&  fun               ,
      const DCT&                        mapping  = DCT () , 
      const Ostap::Utils::ProgressConf& progress = false  ) ;  
    // ============================================================
    /** add new branch to the tree from RooFir function
     *  @param tree input tree 
     *  @param bname branch name 
     *  @param fun   the function 
     *  @param observables   function observables 
     *  @oaram normalization normalization set 
     *  @param mapping : observables -> brnaches 
     */
    Ostap::StatusCode
    add_branch
    ( TTree*                            tree                    ,
      const std::string&                bname                   , 
      const RooAbsReal&                 fun                     ,
      const RooAbsCollection&           observables             ,
      const RooAbsCollection*           normalization = nullptr , 
      const DCT&                        mapping       = DCT ()  , 
      const Ostap::Utils::ProgressConf& progress      = false   ) ;
    // ========================================================================
    /** Add sPlot information to the tree 
     *  @param tree  input tree 
     *  @param splot sPlot object 
     *  @param prefix prefix for the branch names 
     *  @param suffix suffix for the branch names 
     *  @return StatusCode
     */
    StatusCode
    add_branch
    ( TTree*                               tree             ,
      const Ostap::MoreRooFit::SPlot4Tree& splot            ,
      const std::string&                   prefix  = ""     ,
      const std::string&                   suffix  = "_sw"  , 
      const DCT&                           mapping = DCT () , 
      const Ostap::Utils::ProgressConf&   progress = false  ) ; 
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADDBRANCH_H
// ============================================================================
