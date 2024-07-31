// $Id$ 
// ============================================================================
// Include files 
// ============================================================================
// STD&STL 
// ============================================================================
#include <limits>
#include <cassert>
#include <sstream>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
#include "Ostap/NStatEntity.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::NStatEntity
 *  @see Ostap::NStatEntity
 *  @see Ostap::WStatEntity
 *  @see Ostap::StatEntity
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 *  @date   2015-04-03
 */
// ===========================================================================
namespace 
{
  // =========================================================================
  static_assert ( std::numeric_limits<unsigned long>::is_specialized   , 
                  "numeric_limits<unsigned long> is not specialized!" ) ;
  // =========================================================================
  /// define the maximum half-width of a window 
  const unsigned long s_max = std::numeric_limits<unsigned long>::max() << 3 ;
  // =========================================================================
}
// ===========================================================================
// constructor with N-parameter  
// ===========================================================================
Ostap::NStatEntity::NStatEntity
( const unsigned long N ) 
  : m_cnt1 (   ) 
  , m_cnt2 (   ) 
  , m_N    ( N ) 
{
  Ostap::Assert ( 2 <= m_N                ,
		  "Invalid window size!" , 
		  "Ostap::NstatEntity"   ) ;
  
}
// ===========================================================================
// constructor from two counters and sliding window
// ===========================================================================
Ostap::NStatEntity::NStatEntity
( const unsigned long N    ,
  const StatEntity&   cnt1 ,
  const StatEntity&   cnt2 )
  : m_cnt1 ( cnt1 )
  , m_cnt2 ( cnt2 )
  , m_N    ( N    ) 
{
  Ostap::Assert ( 2 <= m_N                ,
		  "Invalid window size!" , 
		  "Ostap::NStatEntity"   ) ;
  Ostap::Assert ( m_cnt1.nEntries () <= 2 *  m_N ,
		  "Invalid first counter!"  , 
		  "Ostap::NstatEntity"   ) ;
  Ostap::Assert ( m_cnt2.nEntries () <= 2 *  m_N ,
		  "Invalid second counter!"  , 
		  "Ostap::NStatEntity"   ) ;					  
}
// ======================================================================


// ===========================================================================
/*  printout  to std::ostream
 *  @param s the reference to the output stream
 */
// ===========================================================================
std::ostream& Ostap::NStatEntity::fillStream ( std::ostream& o ) const 
{
  o << "N=" << m_N << " " ;
  // print the longest counter 
  return counter().fillStream ( o ) ;
}
// ===========================================================================
// conversion to string
// ===========================================================================
std::string Ostap::NStatEntity::toString() const 
{ 
  std::ostringstream ost ;
  fillStream ( ost )  ;
  return ost.str () ;
}
// ===========================================================================
// The END 
// ===========================================================================
