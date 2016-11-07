// $Id$ 
// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
#include <iostream>
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "TCut.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/PyIterator.h"
// ============================================================================
/** @file 
 *  implemetation file for class Analysis::PyIterator
 *  @author Vanya Belyaev
 *  @date   2013-05-06
 *  
 *                    $Revision$
 *  Last modification $Date$
 *                 by $Author$
 */
// ============================================================================
namespace 
{
  static_assert ( std::numeric_limits<unsigned long>::is_specialized , 
                  "std::numeric_limits<unsigned long> is not specialized" ) ;
}
// ============================================================================
// constructor 
// ============================================================================
Ostap::PyIterator::PyIterator 
( TTree*              tree  , 
  const std::string&  cuts  , 
  const unsigned long first , 
  const unsigned long last  )  
  : m_tree     ( tree  ) 
  , m_formula  (       )  
  , m_current  ( first )
  , m_last     ( last  )
{
  // 
  if ( 0 == m_tree ) 
  {
    m_current = 0 ;
    m_last    = 0 ;
  }
  else
  { 
    //
    m_last    = std::min ( m_last , (unsigned long) tree->GetEntries() ) ;
    m_formula.reset ( new Ostap::Formula ( "" ,  cuts , m_tree ) ) ;
    //
    if ( !m_formula->GetNdim() ) { m_formula.reset()                      ; }
    else                         { m_tree->SetNotify ( m_formula.get() )  ; }
    //
  }
  //
  m_tree = next () ;
}
// ============================================================================
// constructor 
// ============================================================================
Ostap::PyIterator::PyIterator 
( TTree*              tree  , 
  const TCut&         cuts  , 
  const unsigned long first , 
  const unsigned long last  )  
  : m_tree     ( tree  )
  , m_formula  (       )  
  , m_current  ( first )
  , m_last     ( last  )
{
  // 
  if ( 0 == m_tree ) 
  {
    m_current = 0 ;
    m_last    = 0 ;
  }
  else
  { 
    //
    m_last    = std::min ( m_last , (unsigned long) tree->GetEntries() ) ;
    m_formula.reset ( new Ostap::Formula ( "" ,  cuts , m_tree ) ) ;
    //
    if ( !m_formula->GetNdim() ) { m_formula.reset()                      ; }
    else                         { m_tree->SetNotify ( m_formula.get() )  ; }
    //
  }
  //
  m_tree = next () ;
}
// ============================================================================
// go to next item 
// ============================================================================
TTree* Ostap::PyIterator::next () const                   // go to next item 
{
  //
  if ( 0 == m_tree ) { return nullptr ; }
  if ( !m_formula  ) { return nullptr ; }
  //
  for ( ; m_current <= m_last ; ++m_current ) 
  {
    //
    const long ievent = m_tree->GetEntryNumber ( m_current ) ;
    if ( 0 > ievent             ) { continue ; }  // CONTINUE
    //
    // ATTENTION! Load here everything!
    const long result = m_tree -> GetEntry ( ievent ) ; 
    if ( 0 >= result ) { return nullptr ; }
    //
    // check the cuts: 
    if ( m_formula && !m_formula->evaluate() ) { continue ; }  // CONTINUE 
    //
    ++m_current   ;   // ADVANCE THE COUNTER   
    return m_tree ;   // return TREE 
  }
  //
  return nullptr ;
}
// ============================================================================
// check if formula is ok 
// ============================================================================
bool Ostap::PyIterator::ok   () const 
{ return m_formula && m_formula->ok() ; }
// ============================================================================
// The END 
// ============================================================================
