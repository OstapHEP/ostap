// ============================================================================
// Include files 
// ============================================================================
#include <string>
#include <tuple>
#include <functional>
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
#include "Ostap/StatusCode.h"
#include "Ostap/AddBranch.h"
#include "Ostap/Funcs.h"
#include "Ostap/Math.h"
#include "Ostap/Notifier.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// Local stuff 
// ============================================================================
#include "local_roofit.h" 
#include "status_codes.h" 
// ============================================================================
/** @file
 *  Implementation file for function Ostap::Trees::add_branch 
 *  @see Ostap::Trees::add_branch 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
// Destructor 
// ============================================================================
Ostap::Trees::Branches::~Branches(){} 
// ============================================================================
bool Ostap::Trees::Branches::add
( const std::string&      name ,
  const Ostap::IFuncTree& func )
{
  const TObject* o = dynamic_cast<const TObject*> ( &func ) ;
  m_map [ name ] = &func ; 
  return true ;
}
// =============================================================================
bool Ostap::Trees::Branches::add
( const Ostap::IFuncTree& func ,
  const std::string&      name )
{ return add ( name , func ) ; }
// ===========================================================================
// get the branch 
// ===========================================================================
const Ostap::IFuncTree*
Ostap::Trees::Branches::branch
( const std::string&key ) const 
{
  // explicit loop 
  FUNCTREEMAP::const_iterator it = m_map.find ( key ) ;
  return ( m_map.end() != it ) ? it->second : nullptr ; 
}
// ===========================================================================
/*  add new branches to the tree
 *  the value of the branch each  is taken from <code>branches</code>
 *  @param tree     input tree 
 *  @param name     the name for new branch 
 *  @param branches the map name->function use to calculate new branch
 *  @return status code 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ===========================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*                        tree     ,  
  const Ostap::Trees::Branches& branches )
{ return add_branch ( tree , branches.map() ) ; }
// ========================================================================
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
Ostap::Trees::add_branch 
( TTree*                            tree     ,  
  const Ostap::Utils::ProgressConf& progress , 
  const Ostap::Trees::Branches&     branches )
{ return add_branch ( tree , progress , branches.map() ) ; }
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
Ostap::Trees::add_branch 
( TTree*                            tree     ,  
  const Ostap::Utils::ProgressConf& progress , 
  const std::string&                name     , 
  const Ostap::IFuncTree&           func     ) 
{
  if ( !tree   ) { return INVALID_TREE ; }
  //
  Double_t value    = 0  ;
  TBranch* branch   = tree->Branch( name.c_str() , &value , (name + "/D").c_str() );
  if ( !branch ) { return CANNOT_CREATE_BRANCH ; }
  //
  const TObject* o = dynamic_cast<const TObject*>( &func ) ;
  Ostap::Utils::Notifier notifier { tree , o ? const_cast<TObject*>(o) : nullptr } ;
  //
  // due to  some strange reasons we need to invoke the Notifier explicitely.
  // - otherwise FuncTH1::Notify causes crash..
  notifier.Notify() ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress  ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar )
    {
      if ( tree->GetEntry ( i ) < 0 ) { break ; };
      //
      value  =  func ( tree  ) ;
      //
      branch -> Fill (       ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/* add new branch with name <code>name</code> to the tree
 * the value of the branch is taken from  function <code>func</code>
 * @param tree    input tree 
 * @param name    the name for new branch 
 * @param func    the function to be used to fill the branch 
 * @return status code 
 * @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 * @date 2019-05-14
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*                  tree ,  
  const std::string&      name , 
  const Ostap::IFuncTree& func ) 
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another mehtod 
  return add_branch ( tree , progress , name , func ) ;
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
Ostap::Trees::add_branch 
( TTree*                            tree    ,  
  const Ostap::Utils::ProgressConf& progress , 
  const std::string&                name    , 
  const std::string&                formula ) 
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  //
  auto func = std::make_unique<Ostap::Functions::FuncFormula>( formula ,  tree ) ;
  if ( !func ) { return Ostap::StatusCode ( CANNOT_CREATE_FORMULA ) ; }
  //
  Ostap::Utils::Notifier notifier ( tree , func.get () ) ;
  //
  return add_branch ( tree , progress , name , *func ) ;
}
// =============================================================================
/*   add new branch with name <code>name</code> to the tree
 *   the value of the branch is taken from the function <code>func</code>
 *   @param tree    input tree 
 *   @param name    the name for new branch 
 *   @param formula the fomula use to calculate new  branch
 *   @return status code 
 *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *   @date 2019-05-14
 */
// =============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*             tree    ,  
  const std::string& name    , 
  const std::string& formula ) 
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another mehtod 
  return add_branch ( tree , progress , name , formula ) ;
}
// ============================================================================
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
Ostap::Trees::add_branch 
( TTree*                                   tree     ,  
  const Ostap::Utils::ProgressConf&        progress , 
  const std::map<std::string,std::string>& branches ) 
{
  //
  if      ( !tree            ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  else if ( branches.empty() ) { return Ostap::StatusCode::SUCCESS         ; }
  //
  typedef FUNCTREEMAP                                                  MAP   ;
  typedef std::vector<std::unique_ptr<Ostap::Functions::FuncFormula> > STORE ;
  //
  STORE store ;
  MAP   map ;
  //
  store.reserve( branches.size () ) ;
  //
  Ostap::Utils::Notifier notifier ( tree ) ;
  //
  for ( const auto& entry : branches )
  {
    store.emplace_back ( std::make_unique<Ostap::Functions::FuncFormula>( entry.second ,  tree ) ) ;
    auto const& func = store.back() ;
    if ( !func ) { return Ostap::StatusCode ( CANNOT_CREATE_FORMULA ) ; }
    //
    notifier.add ( func.get() )   ;
    map [ entry.first ] = func.get () ;    
    //
  }
  //
  // due to some very strange reasons we need to invoke the Notifier explicitely.
  // - otherwise crash could happen causes crash..
  notifier.Notify() ;
  //
  return add_branch ( tree , progress , map ) ;
}
// ============================================================================
/*  add new branches to the tree
 *  the value of the branch each  is taken from <code>branches</code>
 *  @param tree     input tree 
 *  @param name     the name for new branch 
 *  @param branches the map name->formula use to calculate newbranch
 *  @return status code 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*                                   tree     ,  
  const std::map<std::string,std::string>& branches ) 
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another mehtod 
  return add_branch ( tree , progress , branches ) ;
}
// ============================================================================
/*  add new branches to the tree
 *  the value of the branch each  is taken from <code>branches</code>
 *  @param tree     input tree 
 *  @param progress configuration of the progress bar
 *  @param branches the map name->function use to calculate new branch
 *  @return status code 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*                            tree     ,  
  const Ostap::Utils::ProgressConf& progress , 
  const Ostap::Trees::FUNCTREEMAP&  branches ) 
{
  //
  if      ( !tree            ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  else if ( branches.empty() ) { return Ostap::StatusCode::SUCCESS         ; }
  //
  typedef std::vector<const Ostap::IFuncTree*> FUNCTIONS ;
  typedef std::vector<double>                  VALUES    ;
  typedef std::vector<TBranch*>                BRANCHES  ;
  //
  const std::size_t N = branches.size() ;
  VALUES    values    ( N ) ;
  BRANCHES  tbranches ( N ) ;
  FUNCTIONS functions ( N ) ;
  //
  // Notifier
  Ostap::Utils::Notifier notifier { tree } ;
  //
  unsigned int index = 0 ;
  for ( const auto& entry : branches ) 
  {
    const auto&             name = entry.first  ;
    const Ostap::IFuncTree* func = entry.second ;
    //
    if ( !func   ) { return Ostap::StatusCode ( INVALID_TREEFUNCTION ) ; }
    functions [ index ] = func ;
    //
    const TObject* o = dynamic_cast<const TObject*> ( func ) ;
    if ( nullptr != o ) { notifier.add ( const_cast<TObject*> ( o ) ) ; }
    //
    TBranch* branch   = tree->Branch
      ( name.c_str() ,
        &values[index] , (name + "/D").c_str() );
    if ( !branch ) { return Ostap::StatusCode ( CANNOT_CREATE_BRANCH ) ; }
    tbranches [ index ] = branch ;
    //
    ++index ;
  }
  //
  // due to some very strange reasons we need to invoke the Notifier explicitely.
  // - otherwise crash could happen  
  //
  notifier.Notify() ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress ) ;
  //
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar  )
  {
    //
    if ( tree->GetEntry ( i ) < 0 ) { break ; };
    //
    // evaluate the functions
    for ( unsigned int k = 0 ; k < N ; ++k ) { values    [ k ] = (*functions[ k ]) ( tree ) ; }
    //
    // fill the branches
    for ( unsigned int j = 0 ; j < N ; ++j ) { tbranches [ j ] -> Fill ()                   ; }
  }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================
/*  add new branches to the tree
 *  the value of the branch each  is taken from <code>branches</code>
 *  @param tree     input tree 
 *  @param name     the name for new branch 
 *  @param branches the map name->function use to calculate new branch
 *  @return status code 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*                           tree     ,  
  const Ostap::Trees::FUNCTREEMAP& branches ) 
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another mehtod 
  return add_branch ( tree , progress , branches ) ;
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
Ostap::Trees::add_branch 
( TTree*                            tree     , 
  const Ostap::Utils::ProgressConf& progress , 
  const std::string&                name     , 
  const TH1&                        histo    )
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  //
  const TH1* h1 = &histo ;
  if ( nullptr != dynamic_cast<const TH2*>( h1 ) ) 
  { return Ostap::StatusCode ( INVALID_TH1 ) ; }
  //
  Double_t value    = 0  ;
  TBranch* branch   = tree->Branch( name.c_str() , &value , (name + "/D").c_str() );
  if ( !branch ) { return Ostap::StatusCode ( CANNOT_CREATE_BRANCH ) ; }
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar )
  {
    if ( tree->GetEntry ( i ) < 0 ) { break ; };
    //
    value = histo.GetRandom() ;
    //
    branch -> Fill (       ) ;
  }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/*  add new branch to TTree, sampling it from   the 1D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param name   name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH1::GetRandom 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   name  , 
  const TH1&           histo )
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another mehtod 
  return add_branch ( tree , progress , name , histo ) ;
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
Ostap::Trees::add_branch 
( TTree*                            tree     , 
  const Ostap::Utils::ProgressConf& progress , 
  const std::string&                namex    , 
  const std::string&                namey    , 
  const TH2&                        histo    )
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  //
  const TH2* h2 = &histo ;
  if ( nullptr != dynamic_cast<const TH3*>( h2 ) ) 
  { return Ostap::StatusCode ( INVALID_TH2 ) ; }
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  if ( !branch_x ) { return Ostap::StatusCode ( CANNOT_CREATE_BRANCH ) ; }
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() );
  if ( !branch_y ) { return Ostap::StatusCode ( CANNOT_CREATE_BRANCH ) ; }
  //
  TH2& h = const_cast<TH2&> ( histo ) ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i, ++bar )
  {
    if ( tree->GetEntry ( i ) < 0 ) { break ; };
    //
    h.GetRandom2 ( value_x , value_y ) ;
    //
    branch_x -> Fill (       ) ;
    branch_y -> Fill (       ) ;
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
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   namex , 
  const std::string&   namey , 
  const TH2&           histo )
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another mehtod 
  return add_branch ( tree , progress , namex , namey  , histo ) ;
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
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const Ostap::Utils::ProgressConf& progress , 
  const std::string&   namex , 
  const std::string&   namey , 
  const std::string&   namez , 
  const TH3&           histo )
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  if ( !branch_x ) { return Ostap::StatusCode ( CANNOT_CREATE_BRANCH ) ; }
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() );
  if ( !branch_y ) { return Ostap::StatusCode ( CANNOT_CREATE_BRANCH ) ; }
  //
  Double_t value_z   = 0  ;
  TBranch* branch_z  = tree->Branch( namez.c_str() , &value_z , (namez + "/D").c_str() );
  if ( !branch_z ) { return Ostap::StatusCode ( CANNOT_CREATE_BRANCH ) ; }
  //
  TH3& h = const_cast<TH3&> ( histo ) ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar )
  {
    if ( tree->GetEntry ( i ) < 0 ) { break ; };
    //
    h.GetRandom3 ( value_x , value_y , value_z ) ;
    //
    branch_x -> Fill (       ) ;
    branch_y -> Fill (       ) ;
    branch_z -> Fill (       ) ;
  }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/** add new branch to TTree, sampling it from   the 3D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param namex  name of the new branch 
 *  @param namey  name of the new branch 
 *  @param namez  name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return status code 
 *  @see TH2::GetRandom2 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   namex , 
  const std::string&   namey , 
  const std::string&   namez , 
  const TH3&           histo )
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another mehtod 
  return add_branch ( tree , progress , namex , namey , namez , histo ) ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  template <class DATA>
  inline Ostap::StatusCode 
  _add_branch_ 
  ( TTree*                            tree     , 
    const Ostap::Utils::ProgressConf& progress , 
    const std::string&                vname    ,                
    const std::string&                vtype    , 
    const DATA*                       data     ,
    const unsigned long               size     , 
    const DATA                        value    ) 
  {
    //
    if ( !tree ) { return INVALID_TREE   ; }
    if ( !data ) { return INVALID_BUFFER ; }
    //    
    DATA  bvalue  { value } ;
    TBranch* branch = tree->Branch( vname.c_str() , &bvalue , ( vname + vtype ).c_str() );
    //
    const Long64_t total    = tree->GetEntries() ;
    const Long64_t nentries = total < size ? total : size ;
    //
    for ( Long64_t i = 0 ; i < nentries ; ++i ) 
    {
      bvalue =  *(data + i) ;
      branch->Fill() ;
    }
    //
    Ostap::Utils::ProgressBar bar ( nentries , progress ) ;
    for ( Long64_t i = nentries ; i < total  ; ++i , ++bar ) 
    {
      bvalue = value ;
      branch->Fill() ;
    }
    //
    return Ostap::StatusCode::SUCCESS ;  
  }
  // ==========================================================================
  template <class DATA>
  inline Ostap::StatusCode 
  _add_branch_ 
  ( TTree*                            tree     , 
    const std::string&                vname    ,                
    const std::string&                vtype    , 
    const DATA*                       data     ,
    const unsigned long               size     , 
    const DATA                        value    ) 
  {
    /// create fake progress bar
    Ostap::Utils::ProgressConf progress { 0 } ;
    /// delegate to another mehtod 
    return _add_branch_ ( tree , progress , vname , vtype , data , size , value ) ;
  }
  // ==========================================================================
}
// ============================================================================
#if ROOT_VERSION(6,24,0)<=ROOT_VERSION_CODE
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   vname ,  
  const double*        data  , 
  const unsigned long  size  , 
  const double         value ) 
{ return _add_branch_ ( tree , vname , "/D" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*                            tree      , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&                vname     ,  
  const double*                     data      , 
  const unsigned long               size      , 
  const double                      value     ) 
{ return _add_branch_ ( tree , progress , vname , "/D" , data , size , value ) ; }
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   vname ,  
  const float*         data  , 
  const unsigned long  size  , 
  const float          value ) 
{ return _add_branch_ ( tree , vname , "/F" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname ,  
  const float*         data  , 
  const unsigned long  size  , 
  const float          value ) 
{ return _add_branch_ ( tree , progress , vname , "/F" , data , size , value ) ; }
// ========================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   vname ,  
  const char*          data  , 
  const unsigned long  size  , 
  const char           value ) 
{ return _add_branch_ ( tree , vname , "/B" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname ,  
  const char*          data  , 
  const unsigned long  size  , 
  const char           value ) 
{ return _add_branch_ ( tree , progress , vname , "/B" , data , size , value ) ; }
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*                tree  , 
  const std::string&    vname ,  
  const unsigned char*  data  , 
  const unsigned long   size  , 
  const unsigned char   value ) 
{ return _add_branch_ ( tree , vname , "/b" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*                tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&    vname ,  
  const unsigned char*  data  , 
  const unsigned long   size  , 
  const unsigned char   value ) 
{ return _add_branch_ ( tree , progress , vname , "/b" , data , size , value ) ; }
// ========================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   vname ,  
  const short*         data  , 
  const unsigned long  size  , 
  const short          value ) 
{ return _add_branch_ ( tree , vname , "/S" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname ,  
  const short*         data  , 
  const unsigned long  size  , 
  const short          value ) 
{ return _add_branch_ ( tree , progress , vname , "/S" , data , size , value ) ; }
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*                tree  , 
  const std::string&    vname ,  
  const unsigned short* data  , 
  const unsigned long   size  , 
  const unsigned short  value ) 
{ return _add_branch_ ( tree , vname , "/s" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*                tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&    vname ,  
  const unsigned short* data  , 
  const unsigned long   size  , 
  const unsigned short  value ) 
{ return _add_branch_ ( tree , progress , vname , "/s" , data , size , value ) ; }
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   vname ,  
  const int*           data  , 
  const unsigned long  size  , 
  const int            value ) 
{ return _add_branch_ ( tree , vname , "/I" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname ,  
  const int*           data  , 
  const unsigned long  size  , 
  const int            value ) 
{ return _add_branch_ ( tree , progress , vname , "/I" , data , size , value ) ; }
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   vname ,  
  const unsigned int*  data  , 
  const unsigned long  size  , 
  const unsigned int   value ) 
{ return _add_branch_ ( tree , vname , "/i" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname ,  
  const unsigned int*  data  , 
  const unsigned long  size  , 
  const unsigned int   value ) 
{ return _add_branch_ ( tree , progress , vname , "/i" , data , size , value ) ; }
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  ,
  const std::string&   vname ,
  const long*          data  , 
  const unsigned long  size  , 
  const long           value ) 
{ return _add_branch_ ( tree , vname , "/L" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  ,
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname ,
  const long*          data  , 
  const unsigned long  size  , 
  const long           value ) 
{ return _add_branch_ ( tree , progress , vname , "/L" , data , size , value ) ; }
// ============================================================================
/*  copy data from buffer into new branch 
 *  @param tree   The tree 
 *  @param data   input data fuffer 
 *  @param size   length of the buffer
 *  @param value  default value (used for short buffers) 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   vname ,  
  const unsigned long* data  , 
  const unsigned long  size  , 
  const unsigned long value ) 
{ return _add_branch_ ( tree , vname , "/l" , data , size , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname ,  
  const unsigned long* data  , 
  const unsigned long  size  , 
  const unsigned long value ) 
{ return _add_branch_ ( tree , progress , vname , "/l" , data , size , value ) ; }
// ============================================================================
#endif
// ============================================================================
/** copy data from buffer into new branch 
 *  @param tree    The tree 
 *  @param namex   name of the new branch 
 *  @param value   the value 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree        , 
  const std::string&   vname       ,  
  const double         value       ) 
{ return _add_branch_ ( tree , vname , "/D" , &value , 1 , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree        , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname       ,  
  const double         value       ) 
{ return _add_branch_ ( tree , progress , vname , "/D" , &value , 1 , value ) ; }
// ============================================================================
/** copy data from buffer into new branch 
 *  @param tree    The tree 
 *  @param namex   name of the new branch 
 *  @param value   the value 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree        , 
  const std::string&   vname       ,  
  const int            value       ) 
{ return _add_branch_ ( tree , vname , "/I" , &value , 1 , value ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch 
( TTree*               tree        , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&   vname       ,  
  const int            value       ) 
{ return _add_branch_ ( tree , progress , vname , "/I" , &value , 1 , value ) ; }
// ============================================================================


// ============================================================================
// Generic 1D-functions 
// ============================================================================
/*  add new branch to the tree, calculated from the function  
 *  @param tree  (INPUT) The tree 
 *  @param bname (INPUT) branch name 
 *  @param xname (INPUT) name of input variable 
 *  @param fun   (INPUT) the function 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch
( TTree*                        tree  , 
  const std::string&            bname , 
  const std::string&            xname , 
  std::function<double(double)> fun   ) 
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another method 
  return add_branch  ( tree              ,
                       progress          , 
                       bname             , 
                       xname             , 
                       std::cref ( fun ) ) ;                     
}
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
Ostap::Trees::add_branch
( TTree*                            tree      , 
  const Ostap::Utils::ProgressConf& progress  , 
  const std::string&                bname     , 
  const std::string&                xname     , 
  std::function<double(double)>     fun       ) 
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  // create the function 
  const Ostap::Functions::Func1D fun1d { std::cref ( fun ) , xname  , tree } ;
  //
  return add_branch ( tree , progress , bname , fun1d ) ;
}
// ============================================================================
// Generic 2D-functions 
// ============================================================================
/*  add new branch to the tree, calculated from the function  
 *  @param tree  (INPUT) The tree 
 *  @param bname (INPUT) branch name 
 *  @param xname (INPUT) name of input variable 
 *  @param yname (INPUT) name of input variable 
 *  @param fun   (INPUT) the function 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch
( TTree*                               tree  , 
  const std::string&                   bname , 
  const std::string&                   xname , 
  const std::string&                   yname , 
  std::function<double(double,double)> fun   ) 
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another method 
  return add_branch  ( tree              ,
                       progress          , 
                       bname             , 
                       xname             , 
                       yname             , 
                       std::cref ( fun ) ) ;                     
}
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
Ostap::Trees::add_branch
( TTree*                               tree      , 
  const Ostap::Utils::ProgressConf&    progress  , 
  const std::string&                   bname     , 
  const std::string&                   xname     , 
  const std::string&                   yname     , 
  std::function<double(double,double)> fun       ) 
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  // create the function 
  const Ostap::Functions::Func2D fun2d { std::cref ( fun ) , xname  , yname , tree } ;
  //
  return add_branch ( tree , progress , bname , fun2d ) ;
}
// ============================================================================
// Generic 3D-functions 
// ============================================================================
/*  add new branch to the tree, calculated from the function  
 *  @param tree  (INPUT) The tree 
 *  @param bname (INPUT) branch name 
 *  @param xname (INPUT) name of input variable 
 *  @param yname (INPUT) name of input variable 
 *  @param zname (INPUT) name of input variable 
 *  @param fun   (INPUT) the function 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Trees::add_branch
( TTree*                                      tree  , 
  const std::string&                          bname , 
  const std::string&                          xname , 
  const std::string&                          yname , 
  const std::string&                          zname , 
  std::function<double(double,double,double)> fun   ) 
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another method 
  return add_branch  ( tree              ,
                       progress          , 
                       bname             , 
                       xname             , 
                       yname             , 
                       zname             , 
                       std::cref ( fun ) ) ;                     
}
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
Ostap::Trees::add_branch
( TTree*                                      tree      , 
  const Ostap::Utils::ProgressConf&           progress  , 
  const std::string&                          bname     , 
  const std::string&                          xname     , 
  const std::string&                          yname     , 
  const std::string&                          zname     , 
  std::function<double(double,double,double)> fun       ) 
{
  if ( !tree ) { return Ostap::StatusCode ( INVALID_TREE ) ; }
  // create the function 
  const Ostap::Functions::Func3D fun3d { std::cref ( fun ) , xname  , yname , zname , tree } ;
  //
  return add_branch ( tree , progress , bname , fun3d ) ;
}



// ============================================================================
/*  add new branch to the tree from RooFir function
 *  @param tree input tree 
 *  @param bname branch name 
 *  @param fun   the function 
 *  @param observables   function observables 
     *  @oaram normalization normalization set 
 *  @param mapping : observables -> brnaches 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                   tree          ,
  const std::string&       bname         , 
  const RooAbsReal&        fun           ,
  const RooAbsCollection&  observables   ,
  const RooAbsCollection*  normalization ,       
  const Ostap::Trees::DCT& mapping       )
{
  /// create fake progress bar
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// delegate to another method 
  return add_branch  ( tree          ,
                       progress      ,
                       bname         ,
                       fun           ,
                       observables   ,
                       normalization , 
                       mapping       ) ;
}
// ============================================================================
/*  add new branch to the tree from RooFir function
 *  @param tree input tree 
 *  @param bname branch name 
 *  @param fun   the function 
 *  @param observables   function observables 
 *  @oaram normalization normalization 
 *  @param mapping : observables -> brnaches 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                            tree          ,
  const Ostap::Utils::ProgressConf& progress      , 
  const std::string&                bname         , 
  const RooAbsReal&                 fun           ,
  const RooAbsCollection&           observables   ,
  const RooAbsCollection*           normalization ,         
  const Ostap::Trees::DCT&          mapping       )
{
  //
  if ( !tree )                       { return INVALID_TREE ; }
  // 
  /// get/copy the observables
  if ( 0 == ::size ( observables ) ) { return INVALID_OBSERVABLES ; }
  // ==========================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,26,0)
  // ==========================================================================
  RooArgSet obsset { observables } ; ::copy ( observables , obsset ) ;
  // ==========================================================================
#else // ======================================================================
  // ==========================================================================
  RooArgSet obsset { observables } ;
  //===========================================================================
#endif // =====================================================================
  // ==========================================================================
  // actual observables 
  const std::unique_ptr<RooArgSet> obsvars { fun.getObservables ( obsset ) } ;
  if ( !obsvars || ::size ( *obsvars )  != ::size ( obsset ) ) { return INVALID_OBSERVABLES ; }
  
  /// rooFit observables <-> Tree branches  mapping 
  DCT dct  { mapping } ;
  
  // check the proper type of all obsevables
  typedef std::vector<RooAbsRealLValue*>     RSET ;
  typedef std::vector<RooAbsCategoryLValue*> CSET ;
  //
  RSET rvars {} ; rvars.reserve ( ::size ( *obsvars ) ) ;
  CSET cvars {} ; cvars.reserve ( ::size ( *obsvars ) ) ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
  //
  Ostap::Utils::Iterator iter (  *obsvars ) ; // only for ROOT < 6.18
  RooAbsArg* o = 0 ;
  while ( o = (RooAbsArg*) iter .next() )
    {
      //
#else
      //
  for ( auto* o : *obsvars )
    {
      //
#endif 
      //
      Ostap::Assert ( nullptr != o ,
                      "Invalid/nullptr observable" ,
                      "Ostap::Trees:add_branch"    ,
                      INVALID_OBSERVABLE           , __FILE__ , __LINE__ ) ;
      //
      DCT::const_iterator found = dct.find ( o->GetName() ) ;
      if ( dct.end() == found ) { dct [ o->GetName() ] = o->GetName()  ; }
      
      RooAbsRealLValue*       rlv = dynamic_cast<RooAbsRealLValue*>     ( o ) ;
      RooAbsCategoryLValue*   clv = nullptr ;
      if ( nullptr == rlv ) { clv = dynamic_cast<RooAbsCategoryLValue*> ( o ) ; } 
      Ostap::Assert (  ( nullptr != rlv ) ||  ( nullptr != clv ) ,
                       "Invalid/nullptr observable" ,
                       "Ostap::Trees:add_branch"    ,
                       INVALID_OBSERVABLE           , __FILE__ , __LINE__ ) ;
      //
      if ( rlv ) { rvars.push_back ( rlv ) ; }
      if ( clv ) { cvars.push_back ( clv ) ; }
    }
    
  /// get/copy normalization set 
  RooArgSet normset {} ;
  if ( nullptr != normalization ) { ::copy_real ( *normalization , normset ) ; }
  
  
  // Create the getter object 
  Ostap::Trees::Getter getter ( dct , tree ) ;
  
  // create the branch 
  Double_t bvalue = 0  ;
  TBranch* branch = tree->Branch ( bname.c_str() , &bvalue , ( bname + "/D" ).c_str() );
  if ( !branch ) { return CANNOT_CREATE_BRANCH ; }

  /// create the notifies 
  Ostap::Utils::Notifier notifier { tree , &getter } ; 
  // due to  some strange reasons we need to invoke the Notifier explicitely.
  notifier.Notify() ;
  
  /// result 
  typedef Ostap::Trees::Getter::RMAP RMAP ;
  typedef RMAP::const_iterator       RIT   ;
  RMAP   result {} ; 

  // normalization 
  const RooArgSet* norm = ( 0 !=! normalization ) ? &normset : nullptr ;
  
  // loop over the TTree entries
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress  ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar )
    {
      if ( tree->GetEntry ( i ) < 0 ) { break ; };
      //
      Ostap::StatusCode sc = getter.eval ( result , tree ) ;
      if ( sc.isFailure() ) { return sc ; } 
      //
      // real variables 
      for ( RSET::const_iterator rit = rvars.begin() ; rvars.end() != rit ; ++rit )
        {
          RooAbsRealLValue* rlv = *rit ;
          RIT found = result.find ( rlv->GetName() ) ;           
          if ( result.end() == found ) { return INVALID_OBSERVABLE ; }
          rlv->setVal ( found->second ) ;          
        }
      // category variables 
      for ( CSET::const_iterator cit = cvars.begin() ; cvars.end() != cit ; ++cit )
        {
          RooAbsCategoryLValue* clv = *cit ;
          RIT found = result.find ( clv->GetName() ) ;           
          if ( result.end() == found ) { return INVALID_OBSERVABLE ; }          
          clv->setIndex ( Ostap::Math::round ( found->second ) ) ;          
        }
      //
      // Evaluate the RooFit function
      bvalue = fun.getVal ( norm ) ;
      //
      branch->Fill () ; 
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}

// ============================================================================
//                                                                      The END 
// ============================================================================
