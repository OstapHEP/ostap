// ============================================================================
// Include files 
// ============================================================================
// ROOT
// ============================================================================
#if   __has_include("ROOT/RVersion.hxx")
#include "ROOT/RVersion.hxx"
#elif __has_include("RVersion.h")
#include "RVersion.h"
#endif
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "TChain.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/PySelector.h"
#include "Ostap/ProgressConf.h"
#include "Ostap/ProgressBar.h"
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
( TTree*                            tree     , 
  const Ostap::Utils::ProgressConf& progress )
  : TSelector  ()
  , m_event    { 0        }
  , m_tree     ( tree     )
  , m_progress ( progress )
{
  set_tree ( tree ) ;
}
// ============================================================================
// constructor 
// ============================================================================
Ostap::Selector::Selector
( TTree*  tree  )
  : Selector ( tree , false )
{}
// ============================================================================
// constructor 
// ============================================================================
Ostap::Selector::Selector
( const Ostap::Utils::ProgressConf& progress )
  : Selector ( nullptr  , progress )
{}
// ============================================================================
// constructor 
// ============================================================================
Ostap::Selector::Selector ()
  : Selector ( nullptr , false )
{}
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
  if ( tree ) { m_progress.reset ( tree->GetEntries () ) ; }
  TSelector::Init ( m_tree ) ;
}
// ============================================================================
// beginn
// ============================================================================
void   Ostap::Selector::Begin 
( TTree*   tree       )
{
  set_tree ( tree ) ;
  if ( tree ) { m_progress.reset ( tree->GetEntries () ) ; } 
  TSelector::Begin ( m_tree ) ;
}
// ============================================================================
// initialize the slave 
// ============================================================================
void   Ostap::Selector::SlaveBegin   
( TTree*   tree       ) 
{
  set_tree ( tree ) ;
  if ( tree ) { m_progress.reset ( tree->GetEntries () ) ; } 
  TSelector::SlaveBegin ( m_tree ) ;
} 
// ============================================================================
// process 
// ============================================================================
Bool_t Ostap::Selector::Process ( Long64_t entry ) 
{ 
  // increment number of processed events  and advance the progress bar 
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
// terminate
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
//  set the tree 
// ============================================================================
void Ostap::Selector::set_tree  ( TTree* tree )
{
  m_tree = tree ;
  if ( m_tree ) { m_progress.reset ( tree->GetEntries () ) ; }
}
// ============================================================================
/// reset the progress bar (and use new  max-count) 
// ============================================================================
void Ostap::Selector::reset
( const unsigned long long maxevents )
{
  m_event = 0 ;
  if ( m_tree ) { m_progress.reset ( maxevents ? maxevents : m_tree -> GetEntries () ) ; } 
}
// ============================================================================
// process an entry 
// ============================================================================
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
( TTree*                  tree     ,
  TSelector*              selector ,
  const Ostap::EventIndex first    ,  
  const Ostap::EventIndex last     ) 
{  
  if ( !tree || !selector ) { return -1 ; }
  if ( last <= first      ) { return  0 ; }
  //
  const Long64_t size = tree -> GetEntries () ;
  if ( size <= first      ) { return  0 ; }
  //
  const Long64_t LAST  = last <= size ? static_cast<Long64_t> ( last  ) : size ;
  const Long64_t FIRST =                static_cast<Long64_t> ( first )        ; 
  //
  return tree -> Process ( selector , "" , LAST - FIRST , FIRST  ) ;
} 
// ============================================================================
//                                                                      The END 
// ============================================================================
