// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "TChain.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/TreeGetter.h"
#include "Ostap/ToStream.h"
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
( const std::string& expression , 
  TTree*             tree       ) 
  : Getter ( std::vector<std::string>( 1 , expression ) , tree ) 
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::string& expression1 ,
  const std::string& expression2 , 
  TTree*             tree        ) 
  : Getter ( std::vector<std::string>( { expression1 , expression2 } ) , tree )
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::string& expression1 ,
  const std::string& expression2 , 
  const std::string& expression3 ,
  TTree*             tree        ) 
  : Getter ( std::vector<std::string>( { expression1 , 
        expression2 , 
        expression3 } ) , tree )
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::string& expression1 ,
  const std::string& expression2 , 
  const std::string& expression3 ,  
  const std::string& expression4 ,
  TTree*             tree        ) 
  : Getter ( std::vector<std::string>( { expression1 , 
        expression2 , 
        expression3 , 
        expression4 } ) , tree ) 
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::string& expression1 ,
  const std::string& expression2 , 
  const std::string& expression3 ,  
  const std::string& expression4 , 
  const std::string& expression5 ,
  TTree*             tree        ) 
  : Getter ( std::vector<std::string>( { expression1 , 
        expression2 , 
        expression3 , 
        expression4 ,
        expression5 } ) , tree )
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::vector<std::string>& expressions ,
  TTree*                          tree        )
  : m_tree     ( tree    ) 
  , m_map      ()
  , m_formulae ()   
{
  unsigned short index = 0 ;
  for ( auto item : expressions ) { m_map [ "var" + Ostap::Utils::toString ( ++index ) ] = item ; }
  if ( nullptr != tree ) { make_formulae () ; }
}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::map<std::string,std::string>& expressions , 
  TTree*                                   tree        )
  : m_tree     ( tree    ) 
  , m_map      ( expressions )
  , m_formulae ()   
{
  if ( nullptr != tree ) { make_formulae () ; } 
}
// ===========================================================================
// notify 
// ============================================================================
Bool_t Ostap::Trees::Getter::Notify () 
{
  /// attention! here  we delete the variables instead of notify/reset 
  m_formulae.clear() ;
  return true ;
}
// ============================================================================
/// make formulae 
// ============================================================================
Ostap::StatusCode
Ostap::Trees::Getter::make_formulae () const // make formulae
{
  m_formulae.clear() ;
  //
  Ostap::Assert ( nullptr != m_tree      ,
                  "Invalid Tree"         ,
                  "Ostap::Trees::Getter" , 
                  Ostap::StatusCode(701) ) ;
  
  TTree* t  = const_cast<TTree*> ( m_tree ) ;
  for ( SIT it = m_map.begin() ; m_map.end() != it ; ++it  )
    {
      auto formula = std::make_unique<Ostap::Formula> ( it->first , it->second , t ) ; 
      Ostap::Assert ( formula && formula->ok ()       ,
                      "Invalid Formula:" + it->second ,
                      "Ostap::Trees::Getter"          , 
                      Ostap::StatusCode(701)          ) ;
      formula->Notify() ;
      m_formulae [ it->first ].reset ( formula.release() ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
// is everything OK ? 
// ============================================================================
Ostap::StatusCode
Ostap::Trees::Getter::ok ( TTree* tree  ) const
{
  //
  if ( nullptr != tree ) 
    {
      const TChain* chain = dynamic_cast<const TChain*> ( tree ) ;
      if ( chain ) { tree = chain->GetTree() ; }
    }
  //
  if ( nullptr != tree && m_tree != tree )
    {
      m_tree = tree ;
      Ostap::Assert ( make_formulae ().isSuccess() ,
                      "Invalid Formulae"           ,
                      "Ostap::Trees::Getter"       , 
                      Ostap::StatusCode(701)       );                      
    }
  //
  if ( m_formulae.empty() && nullptr != m_tree )
    {
      Ostap::Assert ( make_formulae ().isSuccess() ,
                      "Invalid Formulae"           ,
                      "Ostap::Trees::Getter"       , 
                      Ostap::StatusCode(701)       );                      
    }
  //
  Ostap::Assert ( nullptr != m_tree      ,
                  "Invalid Tree"         ,
                  "Ostap::Trees::Getter" , 
                  Ostap::StatusCode(701) ) ;
  Ostap::Assert ( !m_formulae.empty()    ,
                  "Invalid Formulae"     ,
                  "Ostap::Trees::Getter" , 
                  Ostap::StatusCode(701) ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
// get the results 
// ============================================================================
Ostap::StatusCode
Ostap::Trees::Getter::eval
( std::vector<double>& result ,
  TTree*               tree   ) const
{
  // trivial case 
  if ( m_map.empty() )
    {
      result.clear() ;
      return Ostap::StatusCode::SUCCESS ;
    }
  //
  Ostap::Assert ( ok ( tree ).isSuccess () ,
                  "Invalid Getter"         ,
                  "Ostap::Trees::Getter"   ,  
                  Ostap::StatusCode(701)   ) ;
  //              
  result.resize ( m_formulae.size() ) ;
  unsigned short index = 0 ; 
  for ( FIT it = m_formulae.begin() ; m_formulae.end() != it ; ++it ) 
    { result [ index++ ] = it->second->evaluate() ; }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
// get the results 
// ============================================================================
Ostap::StatusCode
Ostap::Trees::Getter::eval
( std::map<std::string,double>& result ,
  TTree*                        tree  ) const
{
  // 
  result.clear() ;
  // trivial case 
  if ( m_map.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Assert ( ok ( tree ).isSuccess () ,
                  "Invalid Getter"         ,
                  "Ostap::Trees::Getter"   ,  
                  Ostap::StatusCode(701)   ) ;
  //
  for ( FIT it = m_formulae.begin() ; m_formulae.end() != it ; ++it ) 
    { result [ it->first ] = it->second->evaluate() ; }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
/// ROOT technicalities 
// ============================================================================
ClassImp(Ostap::Trees::Getter) // ROOT technicalities 
// ============================================================================
//                                                                      The END 
// ============================================================================
