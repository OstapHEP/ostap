// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Notifier.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Utils::Notifier
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2018-04-09 
 */
// ============================================================================
Ostap::Utils::Notifier::Notifier 
( TTree*   tree , 
  TObject* obj0 , 
  TObject* obj1 , 
  TObject* obj2 , 
  TObject* obj3 , 
  TObject* obj4 ,
  TObject* obj5 , 
  TObject* obj6 , 
  TObject* obj7 , 
  TObject* obj8 , 
  TObject* obj9 )
  : TObject   () 
  , m_tree    ( tree    ) 
  , m_old     ( nullptr )
  , m_objects () 
{
  //
  _pre_action () ;
  //
  add ( obj0 ) ;
  add ( obj1 ) ;
  add ( obj2 ) ;
  add ( obj3 ) ;
  add ( obj4 ) ;
  add ( obj5 ) ;
  add ( obj6 ) ;
  add ( obj7 ) ;
  add ( obj8 ) ;
  add ( obj9 ) ;
  //
  _post_action () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Utils::Notifier::~Notifier() { exit() ; }
// ============================================================================
// Notify them 
// ============================================================================
Bool_t Ostap::Utils::Notifier::Notify   () 
{
  for ( TObject* o : m_objects ) { if ( nullptr != o ) { o->Notify() ; } }    
  return kTRUE ;
}
// ============================================================================
void Ostap::Utils::Notifier::_pre_action()
{
  if ( nullptr != m_tree ) { m_old = m_tree->GetNotify ()  ; }
  add ( m_old ) ;
}
// ============================================================================
void Ostap::Utils::Notifier::_post_action()
{ if ( nullptr != m_tree ) { m_tree->SetNotify   ( this ) ; } }
// ============================================================================
bool Ostap::Utils::Notifier::exit()
{
  //
  if ( nullptr == m_tree || m_tree->GetNotify() != this ) { return false ; }
  //
  if ( this != m_old ) { m_tree->SetNotify ( m_old ) ; } // RESTORE OLD NOTIFICATIONS 
  m_tree = nullptr ;
  return true ;
}



// ============================================================================
ClassImp(Ostap::Utils::Notifier) 
// ============================================================================
// The END 
// ============================================================================
