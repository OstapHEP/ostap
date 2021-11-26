// ============================================================================
#ifndef OSTAP_TREEGETTER_H 
#define OSTAP_TREEGETTER_H 1
// ============================================================================
// Include files 
// ============================================================================
/// STD&STL
// ============================================================================
#include <memory>
#include <vector>
#include <string>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Notifier.h"
#include "Ostap/Formula.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Forward declarations
// ============================================================================
class TTree ; // ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Trees
  { 
    // ========================================================================
    /** @class Getter 
     *  Helper class to get value for the certains expressions 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2021-09-22
     */
    class Getter
    {
    public:
      // ======================================================================
      /// constuctor form the tree and list of expressions
      Getter ( TTree*                           tree        , 
               const std::vector<std::string>& expressions ) ;
      /// constructor from the tree and expression 
      Getter ( TTree*             tree        ,
               const std::string& expression  ) ;
      /// constructor from the tree and expressions 
      Getter ( TTree*             tree        ,
               const std::string& expression1 ,
               const std::string& expression2 ) ;
      /// constructor from the tree and expressions 
      Getter ( TTree*             tree        ,
               const std::string& expression1 ,
               const std::string& expression2 ,
               const std::string& expression3 ) ;
      /// constructor from the tree and expressions 
      Getter ( TTree*             tree        ,
               const std::string& expression1 ,
               const std::string& expression2 ,
               const std::string& expression3 ,
               const std::string& expression4 ) ;
      /// constructor from the tree and expressions 
      Getter ( TTree*             tree        ,
               const std::string& expression1 ,
               const std::string& expression2 ,
               const std::string& expression3 ,
               const std::string& expression4 ,
               const std::string& expression5 ) ;
      // ======================================================================
      /// copy constructor 
      Getter ( const Getter&  ) = delete ;
      /// move constructor 
      Getter (       Getter&& ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tree 
      TTree* tree() const { return m_tree ; }
      // ======================================================================
      /// get the results 
      Ostap::StatusCode   eval ( std::vector<double>& result ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// The tree 
      TTree*                                        m_tree     { nullptr } ;
      /// notifier for the formulas and the tree  
      std::unique_ptr<Ostap::Utils::Notifier>       m_notifier { nullptr } ; 
      // list of formulas 
      std::vector<std::unique_ptr<Ostap::Formula> > m_formulas {}          ;         
      // ======================================================================
    } ;  
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_TREEGETTER_H
// ============================================================================

