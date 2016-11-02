// $Id$ 
// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TCut.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/PySelectorWithCuts.h"
// ============================================================================
/** @file 
 *  Implementation file for class Analysis::SelectorWithCuts
 * 
 *  @date 2013-05-06 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *
 *                    $Revision$
 *  Last modification $Date$
 *  by                $Author$
 */
// ============================================================================
ClassImp(Ostap::SelectorWithCuts) ;
// ============================================================================
// constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
( const std::string& cuts , 
  TTree*             tree , 
  PyObject*          self ) 
  : Ostap::Selector ( tree , self ) 
  , fMycuts            ( cuts        ) 
  , fMyformula         () 
  , m_event            ( 0           )            
{
  if ( tree ) { fMyformula.reset ( new Ostap::Formula ( "" , fMycuts , tree ) ) ; }
}
// ============================================================================
// constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
( const TCut&        cuts , 
  TTree*             tree , 
  PyObject*          self ) 
  : Ostap::Selector ( tree , self      ) 
  , fMycuts            ( cuts.GetTitle () ) 
  , fMyformula         () 
  , m_event            ( 0           )            
{
  if ( tree ) { fMyformula.reset ( new Ostap::Formula ( "" , fMycuts , tree ) ) ; }
}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::SelectorWithCuts::~SelectorWithCuts () {}
// ============================================================================
// notify 
// ============================================================================
Bool_t Ostap::SelectorWithCuts::Notify() 
{
  if ( fMyformula ) { fMyformula->Notify() ; }
  return TPySelector::Notify () ;
}
// ============================================================================
// init 
// ============================================================================
void Ostap::SelectorWithCuts::Init ( TTree* tree ) 
{
  /// reset the event counter 
  m_event = 0 ;
  //
  fMyformula.reset ( new Ostap::Formula ( "" , fMycuts , tree ) ) ;
  //
  TPySelector::Init ( tree ) ;
}
// ============================================================================
//  begin 
// ============================================================================
void Ostap::SelectorWithCuts::Begin ( TTree* tree ) 
{ TPySelector::Begin ( tree ) ; }
// ============================================================================
// slave begin 
// ============================================================================
void Ostap::SelectorWithCuts::SlaveBegin ( TTree* tree ) 
{
  //
  fMyformula.reset( new Ostap::Formula ( "" , fMycuts , tree ) ) ;
  //
  TPySelector::SlaveBegin ( tree ) ;
}
// ============================================================================
// process 
// ============================================================================
Bool_t Ostap::SelectorWithCuts::Process      ( Long64_t entry ) 
{
  /// increment the event counter 
  ++m_event  ;
  //
  if ( fMyformula && fMyformula->GetNdim() && !fMyformula ->evaluate() ) 
  { return false ; }
  //
  return TPySelector::Process ( entry ) ;
}
// ============================================================================
// is formula OK?
// ============================================================================
bool Ostap::SelectorWithCuts::ok () const // is formula OK ? 
{ return fMyformula && fMyformula->ok () ; }


// ============================================================================
// The END 
// ============================================================================
