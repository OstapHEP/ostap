// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
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
#include "Ostap/ProgressBar.h"
// ============================================================================
/** @file 
 *  implemetation file for class Analysis::PyIterator
 *  @author Vanya Belyaev
 *  @date   2013-05-06
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned long>::is_specialized , 
                  "std::numeric_limits<unsigned long> is not specialized" ) ;
  // ==========================================================================
}
// ============================================================================
// constructor 
// ============================================================================
Ostap::PyIterator::PyIterator 
( TTree*                            the_tree  , 
  const Ostap::Utils::ProgressConf& progress  ,
  const std::string&                cuts      , 
  const unsigned long               first     , 
  const unsigned long               last      ) 
  : m_tree     ( the_tree ) 
  , m_formula  ()  
  , m_last     ( first < last ? last : 0  )
  , m_weight   ( 1       )
  , m_current  ( first   )
  , m_progress ( progress                                            , 
                 nullptr == the_tree || last <= first            ? 0 : 
                 (unsigned long) the_tree->GetEntries() <= first ? 0 :
                 std::min ( last - first , (unsigned long) the_tree->GetEntries() ) - first ) 
{
  init ( cuts ) ;
}
// ============================================================================
// constructor
// ============================================================================
Ostap::PyIterator::PyIterator 
( TTree*                            tree      , 
  const Ostap::Utils::ProgressConf& progress  ,
  const TCut&                       cuts      , 
  const unsigned long               first     , 
  const unsigned long               last      ) 
  : PyIterator ( tree       ,
                 progress   , 
                 std::string ( cuts.GetTitle() ) , 
                 first      , 
                 last       ) 
{}
// ============================================================================
// constructor
// ============================================================================
Ostap::PyIterator::PyIterator 
( TTree*                            tree      , 
  const std::string&                cuts      , 
  const unsigned long               first     , 
  const unsigned long               last      ) 
  : PyIterator ( tree  , Ostap::Utils::ProgressConf ( 0 ) , cuts , first , last ) 
{}
// ============================================================================
// constructor
// ============================================================================
Ostap::PyIterator::PyIterator 
( TTree*                            tree      , 
  const TCut&                       cuts      , 
  const unsigned long               first     , 
  const unsigned long               last      ) 
  : PyIterator ( tree  , Ostap::Utils::ProgressConf ( 0 ) , cuts , first , last ) 
{}
// ============================================================================
// intialize cuts 
// ============================================================================
void Ostap::PyIterator::init  ( const std::string& cuts ) 
{
  if ( nullptr == m_tree || m_tree->GetEntries () <= m_current ) 
    {
      m_current = 0       ;
      m_last    = 0       ;
      m_tree    = nullptr ;
      m_weight  = 0       ;
    }
  else 
    { 
      //
      m_last    = std::min ( m_last , (unsigned long) m_tree->GetEntries() ) ;
      m_formula.reset ( new Ostap::Formula ( cuts , m_tree ) ) ;
      //
      if ( !m_formula->GetNdim() ) { m_formula.reset()                      ; }
      else                         { m_tree->SetNotify ( m_formula.get() )  ; }
      //
    }
  /// advance the Tree to the first good event (if possible) 
  if ( nullptr !=  m_tree ) { m_tree = next () ; } 
}
// ============================================================================
// go to next item 
// ============================================================================
TTree* Ostap::PyIterator::next () const                   // go to next item 
{
  //
  m_weight = 0 ;
  if ( 0 == m_tree ) { return nullptr ; }
  if ( !m_formula  ) { return nullptr ; }
  //
  for ( ; m_current < m_last ; ++m_current , ++m_progress )
    {
      // ========================================================================
      //
      const long ievent = m_tree->GetEntryNumber ( m_current ) ;
      if ( 0 > ievent  ) { return nullptr ; } // RETURN
      // ATTENTION
      const long result = m_tree -> GetEntry ( ievent ) ; 
      if ( 0 >= result ) { return nullptr ; }         // RETURN 
      //
      // check the cuts:
      // if ( m_formula && !m_formula->evaluate() ) { continue ; }  // CONTINUE 
      m_weight = m_formula ? m_formula->evaluate() : 1.0 ;
      if ( !m_weight ) { continue ; }  // CONTINUE 
      //
      ++m_current   ;
      ++m_progress  ;
      return m_tree ;
    }
  //
  m_weight = 0 ;
  return nullptr ;
}
// ============================================================================
// check if formula is ok 
// ============================================================================
bool Ostap::PyIterator::ok   () const 
{ return m_formula && m_formula->ok() ; }
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
