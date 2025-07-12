// ============================================================================
// Include files 
// ============================================================================
#include <string>
#include <tuple>
#include <functional>
#include <map>
#include <typeinfo> 
// ============================================================================
// ROOT
// ============================================================================
#include "RVersion.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Names.h"
#include "Ostap/StatusCode.h"
#include "Ostap/AddBranch.h"
#include "Ostap/Funcs.h"
#include "Ostap/Math.h"
#include "Ostap/Notifier.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/HistoInterpolators.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// Local stuff 
// ============================================================================
#include "local_roofit.h" 
#include "status_codes.h" 
// ============================================================================
/** @file
 *  Implementation file for function Ostap::AddBranch
 *  @see Ostap::AddBranch::add_branch 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
// valid name for branch ?
// ============================================================================
bool Ostap::Trees::valid_name_for_branch ( const std::string& name )
{ return Ostap::primitive ( name ) ; }
// ============================================================================
// default constructor 
// ============================================================================
Ostap::Trees::Branches::Branches () {}  
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Trees::Branches::Branches
( const Ostap::Trees::Branches&  right )
  : m_names () 
  , m_map   () 
{ 
  for ( const std::string& name : right.names() )
    {
      const Ostap::IFuncTree* fun = right.branch ( name ) ;
      Ostap::Assert ( nullptr != fun  ,
                      "Inconsistent names/map structure: " + name ,
                      "Ostap::Trees::Branches"                    ,
                      INVALID_TREEFUNCTION , __FILE__ , __LINE__  ) ;                  
      add ( name , *fun ) ;
    }
  Ostap::Assert ( m_names.size() == m_map.size() ,
                  "Inconsistent names/map structure"      ,
                  "Ostap::Trees::Branches"                ,
                  INVALID_BRANCH   , __FILE__ , __LINE__  ) ;                  
  Ostap::Assert ( size() == right.size() ,
                  "Inconsistent copy"                     ,
                  "Ostap::Trees::Branches"                 ,
                  INVALID_BRANCH   , __FILE__ , __LINE__  ) ;                  
  
}
// ============================================================================
// move consructor
// ============================================================================
Ostap::Trees::Branches::Branches
( Ostap::Trees::Branches&& right )
  : m_map   () 
  , m_names () 
{
  std::swap ( m_map   , right.m_map   ) ;
  std::swap ( m_names , right.m_names ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Trees::Branches::~Branches()
{
  // ==========================================================================
  for ( BRANCHES::const_iterator it = m_map.begin() ; m_map.end() != it ; ++it )
    {
      const Ostap::IFuncTree* f = it->second ;
      if ( nullptr != f ) { delete f ; } 
    }  
  // ==========================================================================
}
// ============================================================================
bool
Ostap::Trees::Branches::has_key
( const std::string& name ) const 
{
  if ( m_names.end() == std::find ( m_names.begin() , m_names.end() , name ) ) { return false ; }
  const Ostap::IFuncTree* br = branch ( name ) ;
  Ostap::Assert ( nullptr != br                               ,
                  "Inconsistent names/map structure: " + name ,
                  "Ostap::Trees::Branches"                    ,
                  INVALID_TREEFUNCTION , __FILE__ , __LINE__      ) ;
  return nullptr != br ;
}   
// ============================================================================
bool Ostap::Trees::Branches::add
( const std::string&         name    ,
  const Ostap::IFuncTree&    func    ,
  const TTree*            /* tree */ ) 
{
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name ) ,
                  "Invalid name for branch:\'" + name + "\'" ,
                  "Ostap::Trees::Branches"                   ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  Ostap::Assert ( !has_key ( name )                          , 
                  "Branch already defined :\'" + name + "\'" ,
                  "Ostap::Trees::Branches"                   ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  // ==========================================================================
  m_names.push_back ( name ) ;
  m_map [ name ] = func.clone() ; //
  Ostap::Assert ( m_names.size() == m_map.size() ,
                  "Inconsistent names/map structure: " + name ,
                  "Ostap::Trees::Branches"                    ,
                  INVALID_BRANCH   , __FILE__ , __LINE__      ) ;                  
  // ==========================================================================
  return true ;
}
// ============================================================================
bool Ostap::Trees::Branches::add
( const std::string& name       ,
  const std::string& expression ,
  const TTree*       tree       ) 
{
  if ( expression.empty() ) { return add ( name , name , tree ) ; }   
  // create formula and add it into the map
  const Ostap::Functions::FuncFormula formula { expression , tree } ;
  Ostap::Assert ( !tree || formula.ok()                 ,
                  "Invalid formula: " + expression      , 
                  "Ostap::Trees::Branches::add"         ,                  
                  INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  return add ( name , formula ) ;
}
// =============================================================================
// get the branch 
// ===========================================================================
const Ostap::IFuncTree*
Ostap::Trees::Branches::branch
( const std::string& name ) const 
{
  if ( m_names.end() == std::find ( m_names.begin() , m_names.end() , name ) ) { return nullptr ; }
  BRANCHES::const_iterator found = m_map.find ( name ) ;
  if ( m_map.end  () == found )                                                { return nullptr ; }
  return found->second ; // FIX ME! 
}
// ============================================================================



// ============================================================================
// constructor with ProgressBar configuration
// =============================================================================
Ostap::AddBranch::AddBranch
( const Ostap::Utils::ProgressConf& progress)
  : m_progress ( progress )
{}
// ============================================================================
/* add new branch with name <code>name</code> to the tree
 * the value of the branch is taken from  function <code>func</code>
 * @param tree    input tree 
 * @param progress configuration of the progress bar
 * @param name    the name for new branch 
 * @param func    the function to be used to fill the branch 
 * @return status code 
 * @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 * @date 2019-05-14
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     ,  
  const std::string&                name     , 
  const Ostap::IFuncTree&           func     ) const 
{
  if ( !tree   ) { return INVALID_TREE ; }
  Ostap::Trees::Branches branches {} ;
  branches.add ( name , func , tree ) ;
  return add_branch ( tree , branches ) ;
}
// =============================================================================
/*   add new branch with name <code>name</code> to the tree
 *   the value of the branch is taken from the function <code>func</code>
 *   @param tree    input tree 
 *   @param progress configurtaion of the progress bar
 *   @param name    the name for new branch 
 *   @param formula the fomula use to calculate new  branch
 *   @return status code 
 *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *   @date 2019-05-14
 */
// =============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     ,  
  const std::string&                name     , 
  const std::string&                formula  ) const 
{
  if ( !tree ) { return INVALID_TREE  ; }
  Ostap::Trees::Branches branches {} ;
  branches.add ( name , formula , tree ) ;
  return add_branch ( tree , branches ) ;
}
// =============================================================================
/*  add new branches to the tree
 *  the value of the branch each  is taken from <code>branches</code>
 *  @param tree     input tree 
 *  @param progress configuration of the progress bar
 *  @param name     the name for new branch 
 *  @param branches the map name->formula use to calculate newbranch
 *  @return status code 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                          tree     ,  
  const Ostap::Dict<std::string>& branches ) const 
{
  //
  if      ( !tree            ) { return INVALID_TREE               ; }
  else if ( branches.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Trees::Branches brs {} ;
  for ( const auto& entry : branches ) { brs.add ( entry.first , entry.second , tree ) ; }
  return add_branch ( tree , brs ) ;
}
// ============================================================================
/* add new branches to the tree
 *  the value of the branch each  is taken from <code>branches</code>
 *  @param tree     input tree 
 *  @param progress configuration of the progress bar
 *  @param name     the name for new branch 
 *  @param branches the map name->function use to calculate new branch
 *  @return status code 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     ,  
  const Ostap::Trees::Branches&     branches ) const 
{
  //
  if      ( !tree            ) { return INVALID_TREE               ; }
  else if ( branches.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //

  /// (1) keep branches locally 
  const Ostap::Trees::Branches lbranches { branches } ;
  
  /// #of brached to be added 
  const std::size_t N = lbranches.size() ;  
  
  // ==========================================================================
  //                 branch   function          store 
  typedef std::tuple<TBranch*,const Ostap::IFuncTree*,double> ITEM  ;
  typedef std::vector<ITEM>                                   ITEMS ;
  ITEMS items {} ; items.reserve ( N ) ;
  // ==========================================================================
  
  // ==========================================================================
  // Notifier (needed for chain processing 
  Ostap::Utils::Notifier notifier { tree } ;
  /// create branches# 
  for ( const std::string& name : lbranches.names() ) 
    {
      // ======================================================================
      const Ostap::IFuncTree* func = lbranches.branch ( name ) ;;
      // ======================================================================
      Ostap::Assert ( func                                       ,
                      "Invalid IFuncTree"                        ,
                      "Ostap::AddBranch::add_branch"             ,
                      INVALID_TREEFUNCTION , __FILE__ , __LINE__ ) ;
      //
      const TObject* o = dynamic_cast<const TObject*> ( func ) ;
      if ( nullptr != o ) { notifier.add ( const_cast<TObject*> ( o ) ) ; }
      //
      items.emplace_back ( nullptr , func , 0.0 ) ;
      //
      TBranch* branch   = tree->Branch
        ( name.c_str()                   ,   // name 
          &std::get<2> ( items.back() )  ,   // storage std::get<2> - store 
          ( name + "/D" ).c_str()        ) ; // spoecification 
      Ostap::Assert ( branch ,
                      "Cannot create branch: " + name            ,
                      "Ostap::AddBranch::add_branch"             ,
                      CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
      //
      std::get<0> ( items.back() ) = branch ;  // TBranch 
      std::get<1> ( items.back() ) = func   ;  // IFuncTree
    }
  //
  // due to some very strange reasons we need to invoke the Notifier explicitely.
  // - otherwise crash could happen  
  notifier.Notify() ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , m_progress ) ;
  for ( Long64_t entry = 0 ; entry < nentries ; ++entry , ++bar  )
    {
      if ( tree->GetEntry ( entry ) < 0 ) { break ; };
      // (A) evaluate the functions. and ut  results into store 
      for ( auto& item : items ) { std::get<2> ( item ) = ( *std::get<1> ( item ) ) ( tree ) ; }
      // (B) fill the branches
      for ( auto& item : items ) { std::get<0> ( item ) -> Fill () ; }
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================
/*  add new branch to TTree, sampling it from   the 1D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param progress configuration of the progress bar
 *  @param name   name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH1::GetRandom 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     , 
  const std::string&                name     , 
  const TH1&                        histo    ) const 
{
  if ( !tree ) { return INVALID_TREE ; }
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name ) ,
                  "Invalid name for branch:\"" + name + "\"" ,
                  "Ostap::AddBracnh::add_branch"             ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  //
  Ostap::Assert  ( 1 == histo.GetDimension()                     ,
                   "Invalid TH1 type:"  + std::string ( typeid( histo ).name() ) ,
                   "Ostap::AddBranch::add_branch"                ,
                   INVALID_TH1  , __FILE__ , __LINE__            ) ; 
  //
  TH1D hh {} ;
  histo.Copy ( hh ) ;
  hh.SetDirectory ( nullptr ) ;
  //
  Double_t value  = 0  ;
  TBranch* branch = tree->Branch( name.c_str() , &value , (name + "/D").c_str() );
  Ostap::Assert ( branch ,
                  "Cannot create branch: " + name            ,
                  "Ostap::AddBranch::add_branch"             ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , m_progress ) ;
  for ( Long64_t entry = 0 ; entry < nentries ; ++entry , ++bar )
    {
      if ( tree->GetEntry ( entry ) < 0 ) { break ; };
      //
      value = hh.GetRandom() ;
      //
      branch -> Fill (       ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/** add new branch to TTree, sampling it from   the 2D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param namex  name of the new branch 
 *  @param namey  name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH2::GetRandom2 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     , 
  const std::string&                namex    , 
  const std::string&                namey    , 
  const TH2&                        histo    ) const 
{
  if ( !tree ) { return INVALID_TREE ; }
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namex ) ,
                  "Invalid name for branch:\"" + namex + "\""   ,
                  "Ostap::AddBranch::add_branch"                ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__     ) ;
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namey ) ,
                  "Invalid name for branch:\"" + namey + "\""   ,
                  "Ostap::AddBranch::add_branch"                ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__     ) ;
  //
  Ostap::Assert  ( 2 == histo.GetDimension()                     ,
                   "Invalid TH2 type:"  + std::string ( typeid( histo ).name() ) ,
                   "Ostap::AddBranch::add_branch"                    ,
                   INVALID_TH2  , __FILE__ , __LINE__            ) ; 
  //
  TH2D hh {} ;
  histo.Copy ( hh ) ;
  hh.SetDirectory ( nullptr ) ;
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  Ostap::Assert ( branch_x                                   ,
                  "Cannot create branch: " + namex           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() ) ;
  Ostap::Assert ( branch_y                                   ,
                  "Cannot create branch: " + namey           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , m_progress ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i, ++bar )
    {
      if ( tree->GetEntry ( i ) < 0 ) { break ; };
      //
      hh.GetRandom2 ( value_x , value_y ) ;
      //
      branch_x -> Fill (       ) ;
      branch_y -> Fill (       ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/** add new branch to TTree, sampling it from   the 3D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param progress configuration of the progress bar
 *  @param namex  name of the new branch 
 *  @param namey  name of the new branch 
 *  @param namez  name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH2::GetRandom2 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     , 
  const std::string&                namex    , 
  const std::string&                namey    , 
  const std::string&                namez    , 
  const TH3&                        histo    ) const 
{
  if ( !tree ) { return INVALID_TREE ; }
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namex ) ,
                  "Invalid name for branch:\"" + namex + "\""   ,
                  "Ostap::AddBranch::add_branch"                ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__     ) ;
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namey ) ,
                  "Invalid name for branch:\"" + namey + "\""   ,
                  "Ostap::AddBranch::add_branch"                ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__     ) ;
  //
  Ostap::Assert  ( 3 == histo.GetDimension()                     ,
                   "Invalid TH3 type:"  + std::string ( typeid( histo ).name() ) ,
                   "Ostap::AddBranch::add_branch"                    ,
                   INVALID_TH3  , __FILE__ , __LINE__            ) ;   
  //
  TH3D hh {} ;
  histo.Copy ( hh ) ;
  hh.SetDirectory ( nullptr ) ;
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  Ostap::Assert ( branch_x                                   ,
                  "Cannot create branch: " + namex           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;    
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() );
  Ostap::Assert ( branch_y                                   ,
                  "Cannot create branch: " + namey           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  Double_t value_z   = 0  ;
  TBranch* branch_z  = tree->Branch( namez.c_str() , &value_z , (namez + "/D").c_str() );
  Ostap::Assert ( branch_z                                   ,
                  "Cannot create branch: " + namez           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , m_progress ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar )
    {
      if ( tree->GetEntry ( i ) < 0 ) { break ; };
      //
      hh.GetRandom3 ( value_x , value_y , value_z ) ;
      //
      branch_x -> Fill (       ) ;
      branch_y -> Fill (       ) ;
      branch_z -> Fill (       ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================

// ============================================================================
/*  add new branch to TTree, sampling it from   the 1D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param progress configuration of the progress bar
 *  @param name   name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH1::GetRandom 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     , 
  const std::string&                name     , 
  const Ostap::Math::Histo1D&       histo    ) const 
{
  if ( !tree ) { return INVALID_TREE ; }
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name )             ,
                  "Invalid name for branch:\"" + name + "\"" ,
                  "Ostap::AddBranch::add_branch"                 ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  //
  const Ostap::Math::Histo1D hh { histo } ;
  //
  Double_t value  = 0  ;
  TBranch* branch = tree->Branch( name.c_str() , &value , (name + "/D").c_str() );
  Ostap::Assert ( branch ,
                  "Cannot create branch: " + name            ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , m_progress ) ;
  for ( Long64_t entry = 0 ; entry < nentries ; ++entry , ++bar )
    {
      if ( tree->GetEntry ( entry ) < 0 ) { break ; };
      //
      hh.random ( value ) ;
      //
      branch -> Fill ( ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/** add new branch to TTree, sampling it from   the 2D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param namex  name of the new branch 
 *  @param namey  name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH2::GetRandom2 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     , 
  const std::string&                namex    , 
  const std::string&                namey    , 
  const Ostap::Math::Histo2D&       histo    ) const 
{
  if ( !tree ) { return INVALID_TREE ; }
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namex )             ,
                  "Invalid name for branch:\"" + namex + "\"" ,
                  "Ostap::AddBranch::add_branch"                  ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namey )             ,
                  "Invalid name for branch:\"" + namey + "\"" ,
                  "Ostap::AddBranch::add_branch"                  ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  //
  const Ostap::Math::Histo2D hh { histo } ;
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  Ostap::Assert ( branch_x                                   ,
                  "Cannot create branch: " + namex           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() ) ;
  Ostap::Assert ( branch_y                                   ,
                  "Cannot create branch: " + namey           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  std::array<double,2> result {} ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , m_progress ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i, ++bar )
    {
      if ( tree->GetEntry ( i ) < 0 ) { break ; };
      //
      hh.random ( value_x , value_y ) ;
      //
      branch_x -> Fill (       ) ;
      branch_y -> Fill (       ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/** add new branch to TTree, sampling it from   the 3D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param progress configuration of the progress bar
 *  @param namex  name of the new branch 
 *  @param namey  name of the new branch 
 *  @param namez  name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH2::GetRandom2 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch 
( TTree*                            tree     , 
  const std::string&                namex    , 
  const std::string&                namey    , 
  const std::string&                namez    , 
  const Ostap::Math::Histo3D&       histo    ) const 
{
  if ( !tree ) { return INVALID_TREE ; }
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namex )             ,
                  "Invalid name for branch:\"" + namex + "\"" ,
                  "Ostap::AddBranch::add_branch"                  ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  //
  Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( namey )             ,
                  "Invalid name for branch:\"" + namey + "\"" ,
                  "Ostap::AddBranch::add_branch"                  ,
                  INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
  //
  const Ostap::Math::Histo3D hh { histo } ;
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  Ostap::Assert ( branch_x                                   ,
                  "Cannot create branch: " + namex           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;    
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() );
  Ostap::Assert ( branch_y                                   ,
                  "Cannot create branch: " + namey           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  Double_t value_z   = 0  ;
  TBranch* branch_z  = tree->Branch( namez.c_str() , &value_z , (namez + "/D").c_str() );
  Ostap::Assert ( branch_z                                   ,
                  "Cannot create branch: " + namez           ,
                  "Ostap::AddBranch::add_branch"                 ,
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;  
  //
  std::array<double,3> result {} ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , m_progress ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar )
    {
      if ( tree->GetEntry ( i ) < 0 ) { break ; };
      //
      hh.random ( value_x , value_y , value_z ) ;
      //
      branch_x -> Fill (       ) ;
      branch_y -> Fill (       ) ;
      branch_z -> Fill (       ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================

// ============================================================================
// Generic 1D-functions 
// ============================================================================
/*  add new branch to the tree, calculated from the function  
 *  @param tree     (INPUT) The tree 
 *  @param progress (INPUT) configuration of the progress bar
 *  @param bname    (INPUT) branch name 
 *  @param xname    (INPUT) name of input variable 
 *  @param fun      (INPUT) the function 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch
( TTree*                            tree      , 
  const std::string&                bname     , 
  const std::string&                xname     , 
  std::function<double(double)>     fun       ) const 
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  // create the function 
  const Ostap::Functions::Func1D fun1d { std::cref ( fun ) , xname  , tree } ;
  return add_branch ( tree , bname , fun1d ) ;
}
// ============================================================================
// Generic 2D-functions 
// ============================================================================
/*  add new branch to the tree, calculated from the function  
 *  @param tree     (INPUT) The tree 
 *  @param progress (INPUT) configuration of the progress bar
 *  @param bname    (INPUT) branch name 
 *  @param xname    (INPUT) name of input variable 
 *  @param yname    (INPUT) name of input variable 
 *  @param fun      (INPUT) the function 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch
( TTree*                               tree      , 
  const std::string&                   bname     , 
  const std::string&                   xname     , 
  const std::string&                   yname     , 
  std::function<double(double,double)> fun       ) const 
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  // create the function 
  const Ostap::Functions::Func2D fun2d { std::cref ( fun ) , xname  , yname , tree } ;
  return add_branch ( tree , bname , fun2d ) ;
}
// ============================================================================
// Generic 3D-functions 
// ============================================================================
/*  add new branch to the tree, calculated from the function  
 *  @param tree     (INPUT) The tree 
 *  @param progress (INPUT) configuration of the progress bar
 *  @param bname    (INPUT) branch name 
 *  @param xname    (INPUT) name of input variable 
 *  @param yname    (INPUT) name of input variable 
 *  @param zname    (INPUT) name of input variable 
 *  @param fun      (INPUT) the function 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::AddBranch::add_branch
( TTree*                                      tree      , 
  const std::string&                          bname     , 
  const std::string&                          xname     , 
  const std::string&                          yname     , 
  const std::string&                          zname     , 
  std::function<double(double,double,double)> fun       ) const 
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  // create the function 
  const Ostap::Functions::Func3D fun3d { std::cref ( fun ) , xname  , yname , zname , tree } ;
  return add_branch ( tree , bname , fun3d ) ;
}

  
// ============================================================================
//                                                                      The END 
// ============================================================================
