// ============================================================================
// Include files
// ============================================================================
// ROOT&RooFit 
// ============================================================================
#include "TIterator.h"           // ROOT 
#include "TCollection.h"         // ROOT 
#include "RooLinkedList.h"       // RooFit
#include "RooAbsCollection.h"    // RooFit
// ============================================================================
// Local
// ============================================================================
#include "Ostap/Iterator.h"
// ============================================================================
// standard constructor: create and keep the ietrator 
// ============================================================================  
#if ROOT_VERSION_CODE < ROOT_VERSION(6,31,0)
Ostap::Utils::Iterator::Iterator ( const RooAbsCollection& collection ) 
  : m_iterator ( collection.createIterator() ) 
{}
#endif 
// ============================================================================  
// standard constructor: create and keep the ietrator 
// ============================================================================  
Ostap::Utils::Iterator::Iterator ( const TCollection& collection ) 
  : m_iterator ( collection.MakeIterator() ) 
{}
// ============================================================================  
// standard constructor: create and keep the ietrator 
// ============================================================================  
Ostap::Utils::Iterator::Iterator ( const RooLinkedList& collection ) 
  : m_iterator ( collection.MakeIterator() ) 
{}
// ============================================================================  
Ostap::Utils::Iterator::Iterator ( TIterator* iterator ) 
  : m_iterator ( iterator )
{}
// ============================================================================  
// invoke TIterator::Next
// ============================================================================  
TObject* 
Ostap::Utils::Iterator::next  () const   // invoke TIterator::Next
{ return m_iterator ? m_iterator->Next() : 0 ; }
// ============================================================================  
// invoke TIterator::Reset 
// ============================================================================  
bool Ostap::Utils::Iterator::reset () const 
{
  if ( !m_iterator ) { return false ; }
  m_iterator->Reset() ;
  return true ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
