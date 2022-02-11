// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "TChain.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PySelector.h"
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
ClassImp(Ostap::Selector) ;
// ============================================================================
// constructor 
// ============================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
// ============================================================================
Ostap::Selector::Selector ( PyObject* self , 
                            TTree*    tree ) 
  : ROOT_Selector ( tree , self )
  , m_event { 0 }
{}
// ============================================================================
#else
// ============================================================================
Ostap::Selector::Selector ( TTree* tree ) 
  : ROOT_Selector ()
  , m_event { 0 }
  , m_tree  ( tree ) 
{}
// ============================================================================
#endif
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
  ROOT_Selector::Init ( tree ) ;
}
// ============================================================================
// beginn
// ============================================================================
void   Ostap::Selector::Begin 
( TTree*   tree       )
{
  set_tree ( tree ) ;
  ROOT_Selector::Begin ( tree ) ;
}
// ============================================================================
// initialize the slave 
// ============================================================================
void   Ostap::Selector::SlaveBegin   
( TTree*   tree       ) 
{
  set_tree ( tree ) ;
  ROOT_Selector::SlaveBegin ( tree ) ;
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
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
  //
  return ROOT_Selector::Process ( entry ) ; 
  //
#else 
  //
  return process_entry () ;
  //
#endif 
  //
}
// ============================================================================
// notify 
// ============================================================================
Bool_t Ostap::Selector::Notify         () { return ROOT_Selector::Notify ()  ; }
// ============================================================================
// teminate the slave 
// ============================================================================
void   Ostap::Selector::SlaveTerminate () { ROOT_Selector::SlaveTerminate () ; }
// ============================================================================
/// terminate
// ============================================================================
void   Ostap::Selector::Terminate      () { ROOT_Selector::Terminate () ; }
// ============================================================================
// get entry 
// ============================================================================
Int_t  Ostap::Selector::GetEntry       
( Long64_t entry  , 
  Int_t    getall ) 
{
  //
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
  //
  return ROOT_Selector::GetEntry ( entry , getall ) ;
  //
#else 
  //
  return  m_tree ? m_tree->GetTree()->GetEntry ( entry , getall ) : 0 ;
  //
#endif 
  //
}
// ============================================================================
// version
// ============================================================================
Int_t Ostap::Selector::Version()const 
{
  //
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
  //
  return ROOT_Selector::Version ()  ;
  //
#else 
  //
  return  2 ; // NB! note 2 here!!!
  //
#endif 
  //
}
// ============================================================================
// get the tree 
// ============================================================================
TTree* Ostap::Selector::get_tree () const
{
  //
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
  //
  return ROOT_Selector::fChain ;
  //
#else 
  //
  return  m_tree ;
  //
#endif 
}
// ============================================================================
//  set the tree 
// ============================================================================
void Ostap::Selector::set_tree  ( TTree* tree ) 
{
  //
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
  //
  ROOT_Selector::fChain = tree ;
  //
#else 
  //
  m_tree = tree ;
  //
#endif 
  //
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
  if ( 0 == tree || 0 == selector ) { return 0 ; }
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
  if ( 0 == tree || 0 == selector ) { return 0 ; }
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
  if ( 0 == chain || 0 == selector ) { return 0 ; }
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
  if ( 0 == chain || 0 == selector ) { return 0 ; }
  return chain -> Process ( selector , "" , events , first ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
