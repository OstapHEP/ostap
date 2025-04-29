// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
#include "TTree.h"
#include "TChain.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/PySelector.h"
// ============================================================================
// local
// ============================================================================
#include "status_codes.h"
// ============================================================================
/** @file 
 * 
 *  Implementation file for class Ostap::PySelector
 *
 *  @see Ostap::Selector
 *  @see TPySelector 
 * 
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,36,0)
// ============================================================================
ClassImp(Ostap::Selector) ;
// ============================================================================
#endif 
// ============================================================================
// constructor 
// ============================================================================
Ostap::Selector::Selector
( TTree* tree ) 
  : TSelector ()
  , m_event { 0    }
  , m_tree  ( tree )
{
  set_tree ( tree ) ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Selector::~Selector(){}
// ============================================================================
// init 
// ============================================================================
void   Ostap::Selector::Init
( TTree*   tree       )
{
  set_tree ( tree ) ;
  TSelector::Init ( m_tree ) ;
}
// ============================================================================
// beginn
// ============================================================================
void   Ostap::Selector::Begin 
( TTree*   tree       )
{
  set_tree ( tree ) ;
  TSelector::Begin ( m_tree ) ;
}
// ============================================================================
// initialize the slave 
// ============================================================================
void   Ostap::Selector::SlaveBegin   
( TTree*   tree       ) 
{
  set_tree ( tree ) ;
  TSelector::SlaveBegin ( m_tree ) ;
} 
// ============================================================================
// process 
// ============================================================================
Bool_t Ostap::Selector::Process ( Long64_t entry ) 
{ 
  // increment number of processed events  
  increment_event() ;
  //
  if ( Ostap::Selector::GetEntry ( entry ) <= 0 ) 
  {
    Abort ( "" , TSelector::kAbortFile ) ;
    return false ; 
  }
  //
  return process_entry () ;
}
// ============================================================================
// notify 
// ============================================================================
Bool_t Ostap::Selector::Notify         () { return TSelector::Notify ()  ; }
// ============================================================================
// teminate the slave 
// ============================================================================
void   Ostap::Selector::SlaveTerminate () { TSelector::SlaveTerminate () ; }
// ============================================================================
/// terminate
// ============================================================================
void   Ostap::Selector::Terminate      () { TSelector::Terminate () ; }
// ============================================================================
// get entry 
// ============================================================================
Int_t  Ostap::Selector::GetEntry       
( Long64_t entry  , 
  Int_t    getall ) 
{ return  m_tree ? m_tree->GetTree()->GetEntry ( entry , getall ) : 0 ; }
// ============================================================================
// version
// ============================================================================
Int_t Ostap::Selector::Version()const 
{
  //
  return  2 ; // NB! note 2 here!!!
  //
}
// ============================================================================
// get the tree 
// ============================================================================
// TTree* Ostap::Selector::get_tree () const { return  m_tree ; }
// ============================================================================
//  set the tree 
// ============================================================================
void Ostap::Selector::set_tree  ( TTree* tree )
{
  // if ( tree )
  // {
  // TChain* chain = dynamic_cast<TChain*> ( m_tree ) ;
  // if ( chain ) { tree = chain->GetTree() ;}
  // }
  //
  m_tree = tree ;
}
// ============================================================================
/// process an entry 
bool Ostap::Selector::process_entry ()
{
  Ostap::Assert ( false ,
                  "`process_entry` method must be overrided!" , 
                  "Ostap::Selector"                           ,
                  UNDEFINED_METHOD , __FILE__ , __LINE__      ) ;
  return true ;
} 
// ============================================================================
/*  helper function to use TTree::Process in python 
 * 
 *  @param tree      root-tree 
 *  @param selector  the selector 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
long Ostap::Utils::process
( TTree*             tree     ,
  TSelector*         selector )
{
  if ( 0 == tree || 0 == selector ) { return -1 ; }
  return tree->Process ( selector ) ;
}
// ============================================================================
/*  helper function to use TTree::Process in python 
 * 
 *  @param tree      root-tree 
 *  @param selector  the selector 
 *  @param events    events to be processed 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2013-02-10
 */
// ============================================================================
long Ostap::Utils::process
( TTree*              tree      ,
  TSelector*          selector  , 
  const unsigned long events    , 
  const unsigned long first     ) 
{
  if ( 0 == tree || 0 == selector ) { return -1 ; }
  return tree->Process ( selector , "" , events , first ) ;
} 
// ============================================================================
/* helper function to use TChain::Process in python 
 * 
 *  @param chain     root-chain
 *  @param selector  the selector 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
long Ostap::Utils::process
( TChain*            chain    ,
  TSelector*         selector )
{
  if ( 0 == chain || 0 == selector ) { return -1 ; }
  return chain -> Process ( selector ) ;
}
// ============================================================================
/* helper function to use TChain::Process in python 
 * 
 *  @param chain     root-chain
 *  @param selector  the selector 
 *  @param events    events to be processed 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
long Ostap::Utils::process
( TChain*             chain    ,
  TSelector*          selector ,
  const unsigned long events   ,
  const unsigned long first    ) 
{
  if ( 0 == chain || 0 == selector ) { return -1 ; }
  return chain -> Process ( selector , "" , events , first ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
