// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/TreeGetter.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Trees::Getter
 *  @date 2021-09-22 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  enum {
    INVALID_TREE          = 750 , 
    CANNOT_CREATE_BRANCH  = 751 , 
    CANNOT_CREATE_FORMULA = 752 , 
    CANNOT_READ_TREE      = 753 , 
  };
  // ==========================================================================
}
// ============================================================================
//  constructor from the tree and expression 
// ============================================================================
Ostap::Trees::Getter::Getter
( TTree*             tree        ,
  const std::string& expression  ) 
  : Getter ( tree , std::vector<std::string>( 1 , expression ) ) 
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( TTree*             tree        ,
  const std::string& expression1 ,
  const std::string& expression2 ) 
  : Getter ( tree , std::vector<std::string>( { expression1 , expression2 } ) ) 
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( TTree*             tree        ,
  const std::string& expression1 ,
  const std::string& expression2 , 
  const std::string& expression3 )  
  : Getter ( tree , std::vector<std::string>( { expression1 , 
          expression2 , 
          expression3 } ) ) 
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( TTree*             tree        ,
  const std::string& expression1 ,
  const std::string& expression2 , 
  const std::string& expression3 ,  
  const std::string& expression4 )  
  : Getter ( tree , std::vector<std::string>( { expression1 , 
          expression2 , 
          expression3 , 
          expression4 } ) ) 
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( TTree*             tree        ,
  const std::string& expression1 ,
  const std::string& expression2 , 
  const std::string& expression3 ,  
  const std::string& expression4 , 
  const std::string& expression5 )  
  : Getter ( tree , std::vector<std::string>( { expression1 , 
          expression2 , 
          expression3 , 
          expression4 ,
          expression5 } ) ) 
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( TTree*                          tree        ,
  const std::vector<std::string>& expressions )
  : m_tree ( tree ) 
  , m_notifier ( nullptr ) 
  , m_formulas ( )
{
  // ==========================================================================
  Ostap::Assert  ( tree                   ,
                   "Invalid tree "        , 
                   "Ostap::Tree::Getter"  ,
                   INVALID_TREE           ) ;
  // ==========================================================================
  m_notifier.reset ( new Ostap::Utils::Notifier ( tree ) ) ;
  // ==========================================================================
  for ( const auto& e : expressions  ) 
  {
    auto p = std::make_unique<Ostap::Formula>( e , tree ) ;
    Ostap::Assert  ( p && p->ok() ,
                     "Invalid expression:" + e , 
                     "Ostap::Tree::Getter"      , 
                     CANNOT_CREATE_FORMULA     ) ;
    m_formulas.push_back ( std::move ( p ) ) ;  
    m_notifier->add ( m_formulas.back().get() ) ;
  }
  // ==========================================================================
}
// ============================================================================
// get the results 
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::Getter::eval ( std::vector<double>& result ) const 
{
  if ( m_tree == nullptr ) { return INVALID_TREE ; }
  //
  result.clear() ;
  result.reserve ( m_formulas.size() ) ;
  //
  /// read the tree entry if needed 
  if ( m_tree->GetReadEntry() < 0 ) 
  { 
    if ( m_tree->GetEntry ( 0 ) < 0 ) { return CANNOT_READ_TREE ; } 
  }
  // ==========================================================================
  std::vector<double> res ;
  for ( const auto& f : m_formulas ) 
  {
    f->evaluate ( res ) ;
    for ( const double r : res ) { result.push_back ( r ) ; }
  }
  // 
  return Ostap::StatusCode::SUCCESS ;
}




// ============================================================================
//                                                                      The END 
// ============================================================================
