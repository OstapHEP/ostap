// ============================================================================
// Include files 
// ============================================================================
#include <string>
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
#include "TBranch.h"
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
 * the value of the branch is taked from  function <code>func</code>
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
//                                                                      The END 
// ============================================================================
