// ============================================================================
// Include files 
// ============================================================================
#include <string>
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/AddBranch.h"
#include "Ostap/Funcs.h"
#include "Ostap/Notifier.h"
// ============================================================================
/** @file
 *  Implementation file for function Ostap::Trees::add_branch 
 *  @see Ostap::Trees::add_branch 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
/* add new branch with name <code>name</code> to the tree
 * the value of the branch is taken from  function <code>func</code>
 * @param tree    input tree 
 * @param name    the name for new branch 
 * @param func    the function to be used to fill the branch 
 * @return new  branch 
 * @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 * @date 2019-05-14
 */
// ============================================================================
TBranch* Ostap::Trees::add_branch 
( TTree*                  tree ,  
  const std::string&      name , 
  const Ostap::IFuncTree& func ) 
{
  if ( !tree   ) { return nullptr ; }
  //
  Double_t value    = 0  ;
  TBranch* branch   = tree->Branch( name.c_str() , &value , (name + "/D").c_str() );
  if ( !branch ) { return nullptr ; }
  //
  const TObject* o = dynamic_cast<const TObject*>( &func ) ;
  Ostap::Utils::Notifier notifier { tree , o ? const_cast<TObject*>(o) : nullptr } ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  for ( Long64_t i = 0 ; i < nentries ; ++i )
  {
    if ( tree->GetEntry ( i ) < 0 ) { break ; };
    //
    value = func   ( tree  ) ;
    //
    branch -> Fill (       ) ;
  }
  //
  return branch ;
}
// =============================================================================
/*   add new branch with name <code>name</code> to the tree
 *   the value of the branch is taken from the function <code>func</code>
 *   @param tree    input tree 
 *   @param name    the name for new branch 
 *   @param formula the fomula use to calculate new  branch
 *   @return new  branch 
 *   @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *   @date 2019-05-14
 */
// =============================================================================
TBranch* Ostap::Trees::add_branch 
( TTree*             tree    ,  
  const std::string& name    , 
  const std::string& formula ) 
{
  if ( !tree ){ return  nullptr ; }
  //
  auto func = std::make_unique<Ostap::Functions::FuncFormula>( formula ,  tree ) ;
  //
  Ostap::Utils::Notifier notifier ( tree , func.get() ) ;
  //
  return add_branch ( tree , name , *func ) ;
}
// ============================================================================
/*  add new branch to TTree, sampling it from   the 1D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param name   name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return new  branch 
 *  @see TH1::GetRandom 
 */
// ============================================================================
TBranch* Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   name  , 
  const TH1&           histo )
{
  if ( !tree   ) { return nullptr ; }
  //
  const TH1* h1 = &histo ;
  if ( nullptr != dynamic_cast<const TH2*>( h1 ) ) { return nullptr ; }
  //
  Double_t value    = 0  ;
  TBranch* branch   = tree->Branch( name.c_str() , &value , (name + "/D").c_str() );
  if ( !branch ) { return nullptr ; }
  //
  const Long64_t nentries = tree->GetEntries(); 
  for ( Long64_t i = 0 ; i < nentries ; ++i )
  {
    if ( tree->GetEntry ( i ) < 0 ) { break ; };
    //
    value = histo.GetRandom() ;
    //
    branch -> Fill (       ) ;
  }
  //
  return branch ; 
}
// ============================================================================
/** add new branch to TTree, sampling it from   the 1D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param namex  name of the new branch 
 *  @param namey  name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return new  brances 
 *  @see TH2::GetRandom2 
 */
// ============================================================================
TBranch* 
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   namex , 
  const std::string&   namey , 
  const TH2&           histo )
{
  if ( !tree   ) { return nullptr ; }
  //
  const TH2* h2 = &histo ;
  if ( nullptr != dynamic_cast<const TH3*>( h2 ) ) { return nullptr ; }
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  if ( !branch_x ) { return nullptr ; }
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() );
  if ( !branch_y ) { return nullptr ; }
  //
  TH2& h = const_cast<TH2&> ( histo ) ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  for ( Long64_t i = 0 ; i < nentries ; ++i )
  {
    if ( tree->GetEntry ( i ) < 0 ) { break ; };
    //
    h.GetRandom2 ( value_x , value_y ) ;
    //
    branch_x -> Fill (       ) ;
    branch_y -> Fill (       ) ;
  }
  //
  return branch_y ; 
}
// ============================================================================
/** add new branch to TTree, sampling it from   the 1D-histogram
 *  @param tree (UPFATE) input tree 
 *  @param namex  name of the new branch 
 *  @param namey  name of the new branch 
 *  @param namez  name of the new branch 
 *  @param histo  the historgam to be  sampled
 *  @return new  brances 
 *  @see TH2::GetRandom2 
 */
// ============================================================================
TBranch*
Ostap::Trees::add_branch 
( TTree*               tree  , 
  const std::string&   namex , 
  const std::string&   namey , 
  const std::string&   namez , 
  const TH3&           histo )
{
  if ( !tree   ) { return nullptr ; }
  //
  Double_t value_x   = 0  ;
  TBranch* branch_x  = tree->Branch( namex.c_str() , &value_x , (namex + "/D").c_str() );
  if ( !branch_x ) { return nullptr ; }
  //
  Double_t value_y   = 0  ;
  TBranch* branch_y  = tree->Branch( namey.c_str() , &value_y , (namey + "/D").c_str() );
  if ( !branch_y ) { return nullptr ; }
  //
  Double_t value_z   = 0  ;
  TBranch* branch_z  = tree->Branch( namez.c_str() , &value_z , (namez + "/D").c_str() );
  if ( !branch_z ) { return nullptr ; }
  //
  TH3& h = const_cast<TH3&> ( histo ) ;
  //
  const Long64_t nentries = tree->GetEntries(); 
  for ( Long64_t i = 0 ; i < nentries ; ++i )
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
  return branch_z ; 
}
// ============================================================================
//                                                                      The END 
// ============================================================================
