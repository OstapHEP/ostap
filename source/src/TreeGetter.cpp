// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
#include "TTree.h"
#include "TChain.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsCollection.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/ToStream.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_roofit.h"
#include"status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Trees::Getter
 *  @date 2021-09-22 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::vector<std::string>& expressions ,
  const TTree*                    tree        )
  : m_tree     ( tree    ) 
  , m_map      ()
  , m_formulae ()   
{
  unsigned short index = 0 ;
  for ( auto item : expressions ) { m_map [ "var" + Ostap::Utils::toString ( ++index ) ] = item ; }
  //
  if ( nullptr != tree )
    {
      const Ostap::StatusCode sc = make_formulae () ;
      Ostap::Assert ( sc.isSuccess() , "Invalid Getter" ,  "Ostap::Trees::Getter" , sc , __FILE__ , __LINE__) ;
    }
}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::Getter::Getter
( const std::map<std::string,std::string>& expressions , 
  const TTree*                             tree        )
  : m_tree     ( tree    ) 
  , m_map      ( expressions )
  , m_formulae ()   
{
  if ( nullptr != tree )
    {
      const Ostap::StatusCode sc = make_formulae () ;
      Ostap::Assert ( sc.isSuccess() , "Invalid Getter" ,  "Ostap::Trees::Getter" , sc , __FILE__ , __LINE__ ) ;
    }
}
// ===========================================================================
// copy constructor 
// ===========================================================================
Ostap::Trees::Getter::Getter
( const Ostap::Trees::Getter&  right )
  : TObject ( right )
  , m_tree     ( right.m_tree )
  , m_map      ( right.m_map  )
  , m_formulae () 
{
  if ( nullptr != m_tree )
    {
      const Ostap::StatusCode sc = make_formulae () ;
      Ostap::Assert ( sc.isSuccess() , "Invalid Getter" ,  "Ostap::Trees::Getter" , sc , __FILE__ , __LINE__ ) ;
    }
}
// ===========================================================================
// clone method 
// ===========================================================================
Ostap::Trees::Getter*
Ostap::Trees::Getter::Clone ( const char* /* newname */ ) const
{ return new Getter ( *this ) ; }
// ===========================================================================
Ostap::Trees::Getter::~Getter() {}
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
                  "Ostap::Trees::Getter::make_formulae" , 
                  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  for ( SIT it = m_map.begin() ; m_map.end() != it ; ++it  )
    {
      auto formula = std::make_unique<Ostap::Formula> ( it->first , it->second , m_tree  ) ; 
      Ostap::Assert ( formula && formula->ok ()       ,
                      "Invalid Formula:" + it->second ,
                      "Ostap::Trees::Getter::make_formulae" , 
                      INVALID_FORMULA , __FILE__ , __LINE__ ) ;
      
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
Ostap::Trees::Getter::ok ( const TTree* tree  ) const
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
      const Ostap::StatusCode sc =  make_formulae () ;
      Ostap::Assert ( sc.isSuccess() , "Invalid Formulae" , "Ostap::Trees::Getter::ok" , sc , __FILE__ , __LINE__ ) ;
    }
  //
  if ( m_formulae.empty() && nullptr != m_tree )
    {      
      const Ostap::StatusCode sc = make_formulae () ;
      Ostap::Assert ( sc.isSuccess () , "Invalid Formulae" , "Ostap::Trees::Getter::ok"   ,  sc , __FILE__ , __LINE__ ) ;
    }
  //
  Ostap::Assert ( nullptr != m_tree  ,
                  "Invalid Tree"     ,
                  "Ostap::Trees::Getter::ok"  , 
                  INVALID_TREE     , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( !m_formulae.empty()         ,
                  "Invalid Formulae/3"        ,
                  "Ostap::Trees::Getter::ok"  , 
                  INVALID_FORMULAE  , __FILE__ , __LINE__ ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
// get the results 
// ============================================================================
Ostap::StatusCode
Ostap::Trees::Getter::eval
( Ostap::Trees::Getter::RVCT& result ,
  const TTree*                tree   ) const
{
  // trivial case 
  if ( m_map.empty() )
    {
      result.clear() ;
      return Ostap::StatusCode::SUCCESS ;
    }
  //
  const Ostap::StatusCode sc = ok ( tree ) ;
  Ostap::Assert ( sc.isSuccess () , "Invalid Getter" , "Ostap::Trees::Getter::eval" ,  sc , __FILE__ , __LINE__ ) ;
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
( Ostap::Trees::Getter::RMAP& result ,
  const TTree*                tree  ) const
{
  // 
  result.clear() ;
  // trivial case 
  if ( m_map.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::StatusCode sc = ok ( tree ) ;
  Ostap::Assert ( sc.isSuccess () , "Invalid Getter" , "Ostap::Trees::Getter::eval" ,  sc , __FILE__ , __LINE__ ) ;
  //
  for ( FIT it = m_formulae.begin() ; m_formulae.end() != it ; ++it ) 
    { result [ it->first ] = it->second->evaluate() ; }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
/*  add entry mapping 
 *  @attention no replacement!!! 
 */
// ============================================================================
bool
Ostap::Trees::Getter::add
( const std::string& item       ,
  const std::string& expression )
{
  if ( expression.empty() ) { return add ( item , item ) ; }
  //
  if ( m_map.end() != m_map.find ( item ) ) { return false ; } // RETURN 
  //
  m_map.insert ( DCT::value_type ( item , expression ) ) ;
  return true ;                                                // RETURN
}
// ============================================================================
/// copy constructor 
// ============================================================================
Ostap::Trees::RooGetter::RooGetter
( const Ostap::Trees::RooGetter& right )
  : Getter ( right )
{}
// ===========================================================================
// clone method 
// ===========================================================================
Ostap::Trees::RooGetter*
Ostap::Trees::RooGetter::Clone ( const char* /* newname */ ) const
{ return new RooGetter ( *this ) ; }
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::RooGetter::RooGetter
( const std::vector<std::string>& expressions ,
  const TTree*                    tree        )
  : Getter ( expressions , tree )
{}
// ============================================================================
//  constructor from the tree and expressions 
// ============================================================================
Ostap::Trees::RooGetter::RooGetter
( const std::map<std::string,std::string>& expressions , 
  const TTree*                             tree        )
  : Getter ( expressions , tree )
{}
// ============================================================================
// get the results
// ============================================================================
Ostap::StatusCode
Ostap::Trees::RooGetter::assign
( const TTree&      tree   , 
  RooAbsCollection& result ) const 
{ return assign ( result , &tree ) ; }
// ===========================================================================
Ostap::Trees::RooGetter::~RooGetter() {}
// ============================================================================
// get the results
// ============================================================================
Ostap::StatusCode
Ostap::Trees::RooGetter::assign
( RooAbsCollection& result ,
  const TTree*      tree   ) const
{
  // 
  if ( 0 == ::size ( result ) ) { return Ostap::StatusCode::SUCCESS ; }
  if ( m_map.empty()          ) { return Ostap::StatusCode::SUCCESS ; }
  //
  typedef RMAP::const_iterator RIT ;
  RMAP rmap {} ;
  Ostap::StatusCode sc = eval ( rmap , tree ) ;
  Ostap::Assert ( sc.isSuccess () , "Invalid Getter" , "Ostap::Trees::RooGetter::assign" , sc , __FILE__ , __LINE__ ) ;
  // ==========================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0) // ===============================
  // ==========================================================================
  //
  Ostap::Utils::Iterator iter( result  ) ; // only for ROOT < 6.18
  RooAbsArg* o = 0 ;
  while ( o = (RooAbsArg*) iter .next() )
    {
      // =====================================================================
#else  // ====================================================================
      // =====================================================================
  for ( RooAbsArg* o : result )
    {
      // =====================================================================
#endif // ====================================================================
      // =====================================================================
      if ( nullptr == o ) { continue ; }
      RooAbsRealLValue*       rlv = dynamic_cast<RooAbsRealLValue*>     ( o ) ;
      RooAbsCategoryLValue*   clv = nullptr ;
      if ( nullptr == rlv ) { clv = dynamic_cast<RooAbsCategoryLValue*> ( o ) ; }
      //
      if ( !rlv && !clv        ) { continue ; }
      RIT found = rmap.find ( o->GetName() ) ;
      if ( rmap.end() == found ) { continue ; }
      //
      if      ( rlv ) { rlv->setVal    (                      found->second   ) ; }
      else if ( clv ) { clv->setIndex  ( Ostap::Math::round ( found->second ) ) ; }
      // ======================================================================
    } //                                       The end of the loop over results
  // ==========================================================================
  //
  return Ostap::StatusCode::SUCCESS ;
  //
}    
// ============================================================================
/// ROOT technicalities 
// ============================================================================
ClassImp(Ostap::Trees::Getter)   // ROOT technicalities 
ClassImp(Ostap::Trees::RooGetter) // ROOT technicalities 
// ============================================================================
//                                                                      The END 
// ============================================================================
